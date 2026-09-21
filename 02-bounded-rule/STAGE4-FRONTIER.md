# Stage 4 — Frontier: region shape, and the factor (D2) boundary

**Question this stage answers.** Two loose ends that bound the guarantee's reach: (i) is the shipped
admission region's *shape* leaving money on the table — could a smarter shape at the same certified
coverage buy a smaller critical value — and (ii) what happens when the observations share a common
factor with **known** loadings, where the conditional law is verifiable but the desk statistic is
no longer a heteroscedastic ratio at all?

Status: **both steps complete** (S4.1 `frontier.py`, S4.2 `factor_d2.py`). Numbers below are read from
`deposits/frontier/shapes.json` and `deposits/w6/factor_d2.json`.

---

## 1. S4.1 — Region shape is at the floor; optimising width is not optimising `c`

Three shapes were compared, each at its own **maximal** level mass such that the exact quadrature
coverage is still `>= 1 - gamma` (all twelve cells land on coverage 0.999500 for every shape, so the
comparison is at equal certified coverage — this is what makes it fair):

- `uniform` — the shipped split (half the budget to the shared benchmark scale, half spread equally
  across candidates), equal-tail chi-square levels;
- `best_eqtail` — the split ratio `rho` chosen to minimise the mean log-width per candidate, equal-tail;
- `best_hpd` — the same, but with **shortest** (equal-density) chi-square intervals, which is the
  narrowest interval shape available at a given mass.

**Headline: the shipped uniform split is at the practical floor.** Across all 12 desk cells the best
shape changes the mean critical value by at most **+0.39%** and typically by about **±0.1%**:

| cell | uniform `c` | best_eqtail | best_hpd | best shape vs uniform |
|---|---|---|---|---|
| A df3 | 3.2707 | 3.2691 | 3.2687 | -0.06% |
| A df5 | 2.7071 | 2.7015 | 2.7020 | -0.21% |
| A df7 | 2.5268 | 2.5215 | 2.5219 | -0.21% |
| A df11 | 2.3910 | 2.3892 | 2.3890 | -0.08% |
| A df19 | 2.3081 | 2.3094 | 2.3094 | +0.06% |
| A df35 | 2.3000 | 2.3030 | 2.3032 | +0.14% |
| B df3 | 3.3158 | 3.3182 | 3.3288 | **+0.39%** |
| B df5 | 2.8236 | 2.8251 | 2.8248 | +0.05% |
| B df7 | 2.6605 | 2.6623 | 2.6621 | +0.07% |
| B df11 | 2.5412 | 2.5413 | 2.5413 | +0.00% |
| B df19 | 2.5262 | 2.5268 | 2.5268 | +0.02% |
| B df35 | 2.5363 | 2.5370 | 2.5370 | +0.03% |

Two things worth more than the table itself:

1. **Width and `c` are not the same objective.** At B df3 the HPD shape has a **larger** width at its own grid-best point
   (9.3662) than the equal-tail shape (8.6831) and than the shipped uniform split (8.8949), yet it *raises*
   the critical value by +0.39% — the largest effect in the table. Width and `c` disagree in sign there,
   so a pure width objective would have picked a shape that costs power; the floor statement must be made
   in `c`, which is what the table does. (The A and B widths are identical by construction — region
   geometry depends on n, K and `gamma`, not on `lambda_true` — so the two families are not independent
   evidence about shape.)
2. **All realised FWERs under the shipped shape are `<= alpha`** in these cells (0.0103-0.0238 at 40
   draws per cell), so the shape comparison is not being confounded by a shape that quietly breaks the
   guarantee.

**Floor claim, stated with its scope.** Among the shapes tried — uniform/optimal split crossed with
equal-tail/HPD levels, all at exactly certified coverage — the gap to the shipped rule is bounded by
**0.4% in `c`**. That is a best-possible statement *within the family*, not a theorem about all
possible regions; the deposit's `note` says so, and this report does not upgrade it.

## 2. S4.2 (D2) — the factor design, and a boundary with an exact reference

Design: `X_i = b_i F + eps_i` with one common factor path `F ~ N(0,1)`, `eps ~ N(0, 1-rho)`, equal n,
so every group has variance 1 and the desk rule is the plug-in law at `lambda = 1`. Three loading
patterns: `same`, `cand_only`, `sign_flip`. Measured at B = 2000 draws per cell.

**What the correction had to be, and why the first two attempts failed (both caught by smoke runs).**

1. *λ-matching (rejected as infeasible).* Matching the misspecification by an effective ratio
   `lambda_eff = 1 + (b_i-b_0)^2 - 2 rho` gives 0 at `rho = 0.5` and **-0.5 at rho = 0.75** for
   same-sign exposure — a negative variance ratio, hence `sqrt(lambda)` is NaN. When the factor
   cancels out of the difference, the desk statistic simply is not a heteroscedastic ratio.
2. *Scale matching (right idea, wrong normaliser).* The implemented correction multiplies the critical
   value by `w = sqrt( v / (sigma_i^2 + sigma_0^2) )` with `v = (b_i-b_0)^2 + 2(1-rho)`. My first
   version divided by the constant 2, which silently assumes both groups have unit variance — false for
   `cand_only` (the benchmark carries no loading, so its variance is `1-rho`). Fixing the normaliser
   makes `w = 1.000` for `cand_only` at every `rho`, and that is *exactly* right: numerator and
   denominator expectations fall together and no correction is warranted. The empirical check agrees
   (naive ≡ corrected to the last digit in all four `cand_only` cells).

**Measured outcome at B = 2000** (`deposits/w6/factor_d2.json`, MCSE ≈ 0.003-0.004):

| mode | rho | naive | corrected | w |
|---|---|---|---|---|
| same | 0.00 | 0.0210 (-16%) | 0.0210 (-16%) | 1.000 |
| same | 0.50 | 0.0020 (-92%) | 0.0260 (+4%) | 0.707 |
| same | 0.75 | 0.0000 (-100%) | 0.0255 (+2%) | 0.500 |
| cand_only | 0.25-0.75 | 0.0190-0.0275 | **identical** | 1.000 |
| sign_flip | 0.25 | 0.0270 (+8%) | 0.0120 (-52%) | 1.118 |
| sign_flip | 0.50 | 0.0355 (+42%) | 0.0125 (-50%) | 1.225 |
| sign_flip | 0.75 | 0.0410 (+64%) | 0.0115 (-54%) | 1.323 |

So the correction is **exact where it should be** (`cand_only`), **fixes drastic over-conservatism**
(`same`: 0.0000 → 0.025), and **overshoots where the naive rule is liberal** (`sign_flip`: +42%/+64%
becomes -50%/-54%). Acceptance therefore resolves to the *failure boundary*, and the boundary is
quantified rather than described, by an independent exact evaluation.

**The exact reference (the sharpest number in this stage).** For K = 1 the null law can be computed
exactly by deterministic quadrature — conditioning on the factor's sufficient statistics
(`Fbar ~ N(0,1/n)`, `Q_F ~ chi2_{n-1}`), with the candidate's variance a non-central chi-square
`nc = (b_1/se)^2 Q_F` and a shifted normal numerator (packaged as a Poisson-mixture over
`chi2_{nu+2k}`, no simulation). That gives, for `sign_flip`:

| rho | exact `c` (K=1) | naive `c` | required factor | my `w` |
| --- | --- | --- | --- | --- |
| 0.25 | 2.5746 | 1.9944 | 1.291 | 1.118 |
| 0.50 | 3.4537 | 1.9944 | 1.732 | 1.225 |
| 0.75 | **5.2733** | 1.9944 | **2.644** | 1.323 |

The desk rule needs a critical value **2.644× larger** at `rho = 0.75` (exact `c` / deposited naive
reference), and the closed-form numerator correction supplies 1.323× — less than half (in log terms) of
what the exact law demands. Both factors are computed from `deposits/cert/k1_reference.json`. The residual
is the correlation the shared factor induces in the *denominators*: given the path, each `S_i^2` is
non-central chi-square whose variability, not just its mean, is inflated. A single multiplicative
factor cannot represent that, which is why the correction overshoots `sign_flip` and why the honest
deliverable is the boundary plus the exact reference, not a repaired rule.

## 3. What Stage 4 adds to the boundary statement

 **What closing D2 for general `K` would take (named, not hand-waved).** The K=1 exact law already conditions on the factor's two sufficient statistics `(Fbar, Q_F)`; conditioned *also* on the benchmark block `(xbar_0, S_0^2)` the candidates are independent, each with a non-central chi-square denominator (`nc_i = (b_i/se)^2 Q_F`) and a shifted normal numerator. So the general-`K` exact law is a 3-D quadrature over `(Fbar, Q_F, S_0^2)` with a product over candidates — the same machinery already implemented for K=1, extended in one dimension. What is missing is (i) the implementation and (ii) its own quadrature-error certification under a refinement test, i.e. exactly the standard this project applies to everything else. The factorisation argument was independently checked by the reviewer: writing `u = (F - F_1)/||F - F_1||`, each candidate's variance is `(b_i sqrt(Q_F) + u.eps_i)^2 + ||eps_i,perp||^2` - a non-central chi-square whose non-centrality is driven by `b_i^2 Q_F` - and although `u` is shared across candidates, the conditional law of `u.eps_i` given `u` does not depend on `u`, so the joint law across candidates factorises exactly. That endorsement removes the mathematical risk from this item: what remains is bounded implementation work with a named refinement certification, not research risk.

- **The shipped region is at the family floor**: no split/shape among those tried lowers `c` by more
  than 0.4%, and the HPD shape can *raise* it while narrowing width. The desk rule is not leaving an
  easy optimisation on the table.
- **Factor structure with known loadings is outside the guarantee's coverage**: at `rho = 0.75` under
  opposite-sign loadings the correct critical value is 2.6× the plug-in's, and the naive rule is
  liberal by +64% of nominal. The D2 deliverable is the *quantified* boundary plus an exact K=1
  reference implementation — the full composite law for K >= 2 remains open, with its own nuisance and
  therefore its own envelope.

## 4. Reproduction

    python frontier.py                # S4.1  (~3.6k s; shapes at equal certified coverage)
    W26_D2_B=2000 python factor_d2.py # S4.2  (~3.2k s; factor grid + K=1 exact quadrature)
