# PHASE 6 — factor/paired extension and governance (plan M6: sections 8, 9)

tier: oracle_exact_law for the factor law (with the stated limits); calibrated_measured for the
simulation matrix; measured for the governance documents

**Status: complete.** Every number is read from a strict-JSON deposit under `deposits/`.

---

## 1. What P6 delivered

| file | role |
|---|---|
| `factor_general.py` | the general-K one-factor conditioned exact law, with its validation |
| `paired_sim.py` | the §8.6 paired-return simulation matrix, three estimands measured separately |
| `ESTIMAND-DECISION.md` | the estimand choice, recorded as a desk/governance decision |
| `GOVERNANCE.md` | the production decision-record template (§9.2) and the operational policy (§9.3) |
| deposits | `factor_general.json`, `paired_sim.json` |

## 2. The factor-conditioned law — `deposits/factor_general.json`

**The conditioning that makes it exact.** For `r_{i,t} = alpha_i + beta_i f_t + eps_{i,t}` with one
Gaussian factor, condition on the realised factor path: each group's sample mean is
`N(alpha_i + beta_i fbar, sigma_i^2/n_i)` and its sample variance is an independent scaled chi-square,
so the candidates are conditionally independent and

    F(c) = E_{fbar}[ prod_i E_{S_i^2, S_0^2}[ Phi( (c * SE_i - (beta_i - beta_0) fbar) / sigma_diff_i ) ] ],

with `SE_i = sqrt(S_i^2/n_i + S_0^2/n_0)`.  The factor adds a 1-D outer expectation to the structure
the iid law already has; the inner expectations are over the independent chi-square PAIR, evaluated on
a product Gauss-Legendre grid on the probability scale.

**Validation (three routes):**

| check | result |
|---|---|
| K = 1 against a direct Monte-Carlo of the same model (40,000 draws) | `c* = 2.10981`; simulated FWER `0.024525` against `alpha0 = 0.0245`, **z = +0.03** — PASS |
| `beta = 0` must reproduce the iid law | `2.39736` against the iid oracle `2.36622` (**+1.31%**) — PASS, and the residual is the factor quadrature's grid, documented |
| general K against a direct simulation (40,000 draws) | K = 3: simulated FWER `0.02237`, z = −2.87 → **PASS**; K = 5: `0.02135`, z = −4.36 → **conservative, stated** |

**The language must stay cautious.** The factor-conditioned extension **agrees with the
designated checks at K = 1** (`z = +0.03`) and is **consistent with the simulation reference at K = 3**
(`z = -2.87`, inside four Monte-Carlo standard errors), but it requires **larger independent runs before
being described as fully validated**.  At K = 5 the result is **measured conservatism, not
theorem-grade validation** (the simulated FWER sits ~3e-3 below `alpha0`), and **estimated loadings
remain a plug-in**.  Per the plan's own instruction ("Do not claim a general solution until both
validations pass"), no general-K theorem is claimed.  Two earlier versions of the law were wrong in
instructive ways, both caught by the `beta = 0` limit: the first paired the two chi-square quantile
grids index-by-index (a comonotonic coupling) and produced a critical value **2.5× too large**; the
second used midpoint grids, which are only first-order.

**Estimated loadings are a labelled plug-in**, not a theorem: the rolling-window loading estimates in
`paired_sim.py` show their effect, and no certified claim is attached to them.

## 3. The paired-return simulation matrix — `deposits/paired_sim.json`

17 declared scenarios (the deposited count): factor structure (market / market+sector), loadings (positive / mixed / one
negative), residual correlation, serial correlation `phi` in {0, 0.2, 0.4}, volatility clustering,
signal overlap, unequal lengths, estimated loadings, K in {3, 5}, T in {60, 120}, spread in {1, 2},
plus three combined stresses.  The three estimands are formed **separately in every replication**:
raw, paired difference, and factor residual.

| quantity | raw | paired | residual |
|---|---|---|---|
| max realised FWER (nominal `0.0245`) | 0.0135 | 0.0135 | 0.0045 |
| scenarios above nominal | **0 / 17** | **0 / 17** | **0 / 17** |

**The certified iid rule is conservative in every scenario tested**, and most so under the factor
route — which is the expected direction, because the rule's threshold is the supremum over a wide
ratio region and the factor route removes shared variance that the rule was not calibrated for.
Two honest readings:

1. **Conservative is not correct.**  The rule's FWER being *below* nominal means the desk is
   discarding real candidates, and the size of that conservatism is the price of applying an
   independence-calibrated threshold to a dependent process.  It is reported, not celebrated.
2. **The serial-correlation axis still bites.**  The matrix's `phi` scenarios keep the *realised* FWER
   under nominal here only because the threshold is so conservative; the 02-bounded-rule AR(1) study, which used
   the *tight* plug-in threshold, is the study that shows the boundary
   (`deposits/ar2/`, phi ≥ 0.1 at n = 36).

## 4. Governance — `GOVERNANCE.md`, `ESTIMAND-DECISION.md`

* The **decision-record template** captures every field the plan's §9.2 lists: code hashes, the family
  declaration and its date, the data vintage, the sample period, the return construction, the
  estimand, the factor model, the tier, `abar`/`gamma`/`alpha0`, the region and its coverage, the
  certified threshold with its bracket and gap, the certificate identity, the error budget, the
  assumption diagnostics, and the approval/override trail.  The automatic half is exactly what
  `rule_certified.py` already writes.
* The **operational policy** fixes the theorem tier for capital/promotion gates, measured/calibrated
  rules for triage only, the minimum-diagnostics requirement, the mandatory factor-residual route when
  raw assumptions fail, and the escalation path when certification is not affordable - with the tier
  gate (`tests_gates.py` E4) as the mechanical enforcement.
* The **estimand document** records the six candidate estimands, what each one costs, what the
  programme's evidence says about it, and the recommendation (factor alpha for material exposures,
  labelled measured-tier until a certified factor law covers estimated loadings and general K).

## 5. What is NOT claimed

1. No theorem for the factor route with **estimated loadings**.
2. No theorem for the factor route at **K > 3** (conservative and stated, not validated).
3. No theorem for **serially dependent** returns: the AR(1) boundary stands.
4. No coverage of **risk-adjusted or utility-adjusted** estimands: the certified law is for the
   mean-difference statistic.
5. The simulation matrix is **measured**, not certified: it measures how far each estimand departs
   from the iid theorem's assumptions, and it does not extend the theorem.

## 6. Reproduction

    python factor_general.py     # the factor law + validation   -> factor_general.json
    python paired_sim.py         # the §8.6 matrix               -> paired_sim.json
    # governance and the estimand decision are documents, not runs:
    #   GOVERNANCE.md, ESTIMAND-DECISION.md
