# Estimand decision record (plan section 8.3)

tier: measured (a governance choice, not a statistical result)

**Status: decision required from the desk; the statistical machinery supports every option below, and
the choice changes what "outperformance" means.** The plan is explicit that this must be a documented
desk/governance decision rather than a code default, so this file records the options, what each one
costs, and what the current code does.

---

## 1. Why the estimand must be chosen first

The many-to-one test answers one question: *is candidate `i`'s mean above the benchmark's, relative to
the noise in that difference?* Under a factor model the candidates' returns share the factor, so the
answer depends entirely on what has been subtracted:

| estimand | what is tested | what it means for the desk | what the current code does |
|---|---|---|---|
| **raw** candidate mean | `E[r_i - r_0]` on the raw series | "does this strategy beat the sleeve?" - but a strategy that is simply high-beta to the factor will pass whenever the factor did well | supported (the paired statistic) |
| **paired difference** | the same statistic on `d_t = r_{i,t} - r_{0,t}` | the correct *geometry* (the benchmark is removed pairwise), but the differences still share the factor, so the factor exposure remains in the answer | supported |
| **factor alpha** | `alpha_i` in `r_{i,t} = alpha_i + beta_i f_t + eps_{i,t}` | "does this strategy beat the benchmark *after* the factor it is exposed to?" - changes the question from outperformance to alpha relative to the chosen factor model | supported (the residual route in `paired_sim.py`) |
| **residual alpha** | `alpha_i` with a residual-risk adjustment (residual correlation, residual volatility) | the same as factor alpha but with the residual covariance geometry respected in the *test*, not only in the estimate | partially: the residual route is simulated, not certified |
| **risk-adjusted return** | mean return divided by realised volatility (or a Sharpe-like ratio) after costs | what a capital allocator actually compares, but the statistic is no longer the Welch difference and the certified law does not apply to it | not supported; a different law would be required |
| **utility / capacity-adjusted** | the incremental contribution after the cost model and capacity limits | the most honest desk objective, and the furthest from any certified law in this programme | not supported |

## 2. What the programme's evidence says

* **`paired_sim.py`** measures the three *statistical* estimands separately (raw, paired difference,
  factor residual) across 20 scenarios of factor structure, loadings, residual correlation, serial
  correlation, volatility clustering, overlap, unequal lengths and estimated loadings. Its deposit
  (`deposits/paired_sim.json`) reports the realised FWER of the *certified iid rule* under each
  estimand and each scenario, so the desk can see how far each one departs from the theorem's
  assumptions.
* **`factor_general.py`** provides the exact factor-conditioned law for the *known-loading* case
  (validated against Monte-Carlo at K ≤ 3, conservative and stated at K = 5) - the basis on which a
  factor-alpha rule could be certified rather than merely measured.
* **The AR(1) boundary** (02-bounded-rule, `deposits/ar2/`) shows that no rule in the tested set survives a
  lag-one correlation of about 0.1 at these sample sizes, which is why the serial-dependence axis of
  `paired_sim.py` is a *boundary*, not a nuisance parameter.

## 3. Recommendation (to be ratified, not assumed)

1. **For any promotion or capital gate**: state the estimand explicitly in the decision record, and
   use the **factor alpha** estimand when the candidate's factor exposure is material (which it almost
   always is for strategy returns). The certified threshold applies to the raw/paired statistic; if
   the estimand is factor alpha, the threshold must come from the factor-conditioned law, and until
   that law is certified for estimated loadings and general K, the decision is **measured-tier** and
   must be labelled so.
2. **For triage and research ranking**: the raw or paired statistic with the certified threshold is
   acceptable, labelled as the raw-return question it actually answers.
3. **Never** report a risk-adjusted or utility-adjusted comparison as carrying the certified
   guarantee: the certified law is for the mean-difference statistic, and no result in this programme
   covers a ratio statistic.

## 4. What the code will and will not do

* It will not choose for you: `rule_certified.py` takes the statistics as inputs, and every deposit
  records which statistic was used.
* It will refuse to describe a measured-tier route as guaranteed: the tier gate (`tests_gates.py`,
  gate E4) rejects reserved words attached to measured-tier deposits and reports.
* It records, for every decision, the region, the coverage, the threshold, its bracket and its gap -
  so a governance audit can see exactly which claim was made and on what basis.
