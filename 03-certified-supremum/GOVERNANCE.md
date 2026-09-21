# Governance: the production decision record and the operational policy (plan section 9)

tier: measured (a policy document, not a statistical result)

The plan's §9 is explicit that the statistical procedure does **not** by itself solve research-process
multiplicity, and that the two belong to different layers. This file is the layer that wraps the
procedure: what every promotion decision must record, and what the desk is allowed to claim.

---

## 1. Production decision record (plan §9.2) — the template

Every promotion, rejection or override writes one record with these fields. The first block is
produced automatically by `rule_certified.py` (`deposits/rule_certified.json` is one such record per
family); the rest is desk-supplied and must not be omitted.

| field | source | notes |
|---|---|---|
| `strategy_code_hash`, `benchmark_code_hash` | desk | the exact revision tested |
| `candidate_family` and `declaration_date` | desk | the family declared BEFORE the results were seen |
| `data_snapshot` / vintage | desk | the dataset identity and the extraction date |
| `sample_period` | desk | start and end, with any gaps |
| `return_construction`, `costs`, `corporate_actions` | desk | how the return series was built |
| `estimand` | desk, per `ESTIMAND-DECISION.md` | raw / paired / factor alpha / residual alpha / risk-adjusted / utility |
| `factor_model`, `residualisation` | desk | required when the estimand is factor alpha |
| `tier` | automatic | `region_sup_certified` for the theorem-grade rule; `envelope_measured`, `calibrated_measured`, `plugin_measured` for the others |
| `abar`, `gamma`, `alpha0` | automatic | the hard ceiling, the miss budget and the internal level actually used |
| `region`, `coverage_exact` | automatic | the ratio region and its exact quadrature coverage |
| `c_hat`, `c_L`, `c_U`, `gap`, `certified_at_tol` | automatic | the certified threshold and its enclosure |
| `statistics`, `decisions` | automatic | the Welch statistics and the promotions |
| `certificate_id` | automatic | the branch-and-bound certificate identity |
| `error_budget` | automatic | `Delta p` and its c-equivalent |
| `assumption_diagnostics` | desk | normality, independence, serial-correlation and factor-exposure diagnostics |
| `numerical_certificate_id` | automatic | the deposit's own identifier |
| `decision`, `approver`, `override_trail` | desk | who decided, and every override with its reason |

## 2. Operational policy (plan §9.3)

1. **Theorem-grade rule for capital/promotion gates.** A gate that moves capital uses
   `region_sup_certified` only: the threshold is a certified upper bound on the continuous supremum
   over a coverage-certified region, at `alpha0 = abar - gamma`, and the guarantee reads
   `FWER <= alpha0 + gamma = abar`.
2. **Measured and calibrated rules only for triage**, and always labelled. `envelope_measured` (the
   02-bounded-rule vertex envelope), `plugin_measured`, `calibrated_measured` and every simulation-based rule may
   rank research candidates but may not carry the word "guaranteed".
3. **Minimum effective sample size and dependence diagnostics.** Before any theorem-grade decision:
   report the per-group degrees of freedom and a serial-correlation diagnostic. The programme's own
   boundary is explicit: no rule in the tested set survives a lag-one correlation of about 0.1 at
   these sample sizes (02-bounded-rule AR(1) study, `deposits/ar2/`), so a serially dependent series must be
   routed to a dependence-aware route, not to the certified rule.
4. **Mandatory factor-residual or paired route** when the raw-return assumptions fail: if the
   candidate's factor exposure is material, the estimand must be factor alpha (per
   `ESTIMAND-DECISION.md`), and the decision is measured-tier until a certified factor law covers the
   estimated-loading, general-K case.
5. **Escalation when certification is not affordable.** If the certified threshold is too conservative
   for the decision at hand, or the certification cannot be completed within the run-time budget, the
   escalation is *recorded*, not silent: name the region, the achieved gap, the cost, and the
   fallback used. Never widen the tolerance without a written decision.
6. **No silent substitution.** A calibrated or measured threshold must never replace the certified one
   inside a governance gate; if both are computed, the certified one governs and the other is
   advisory.
7. **Tier discipline is enforced mechanically.** `tests_gates.py` gate E4 fails if a deposit lacks a
   valid tier or if measured-tier text uses the reserved words; `run_modes.py quick` prints the gate's
   verdict and returns a non-zero exit status when it fails.

## 3. What this layer does NOT do

It does not model, and must not be represented as modelling: strategy generation and variant
proliferation; feature and hyperparameter search; repeated peeking and rerunning; benchmark selection
after the results are known; data-vintage changes; validation/holdout reuse; discretionary overrides.
Those are the research-process multiplicity the plan assigns to a complementary audit protocol, and
the correct desk statement is that the statistical guarantee is conditional on the family having been
declared in advance.
