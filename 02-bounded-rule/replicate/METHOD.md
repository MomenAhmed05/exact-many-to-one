# Method specification for independent replication

This document specifies **exactly what to compute**. It is intended for a party that has *not* seen
the original codebase. Implement it from this text alone; do not read `02-bounded-rule/*.py`.
Deliverables and the reporting template are in `README.md`. Expected values are deliberately
**withheld** — report your numbers and they will be compared against the sealed originals.

Notation: `K` candidates compared with a shared benchmark (control). `n_0` benchmark sample size,
`n_i` candidate sample sizes. `lambda_i = sigma_i^2 / sigma_0^2` are the variance ratios;
`lambda_hat_i = s_i^2 / s_0^2` are their estimates. `alpha = 0.025`, `gamma = 0.0005`.

---

## 1. The exact conditional law (the object everything is built on)

For fixed variance ratios `lambda`, the null distribution of the many-to-one max-t statistic

    T_i = (xbar_i - xbar_0) / sqrt( s_i^2/n_i + s_0^2/n_0 ),        i = 1..K

is **exactly**

    F(c, lambda) = E_U [ prod_{i=1..K} Phi( s_i (U + c d_i) ) ],     U ~ N(0, 1/n_0),

    s_i   = sqrt( n_i / lambda_i ),
    d_i   = sqrt( lambda_i / n_i + 1 / n_0 )        (sigma_0^2 = 1 without loss of generality).

Equivalently: condition on the benchmark's sufficient statistics `U = xbar_0 - mu_0` and
`V = s_0^2/n_0`; given them the candidates are independent. Compute the 1-D expectation in `U` by
Gauss-Hermite (we use 40 panels x order 10). The population critical value at `lambda` is the root of
`F(c, lambda) = 1 - alpha`; call it `c*_pop(lambda)`.

**Two objects share this notation, and the blind replication showed that conflating them is the easiest mistake in this kit.** The display above is the law *conditional on the variance ratios* (the `s_i` are known constants). That conditional object is what the certificate maximises over a region of `lambda`, and it is what an ORACLE that knows the true variances obeys; for K=1, `n_0=n_1`, `lambda=1` it returns the NORMAL quantile, not the Student t. The **unconditional** null law - the object referenced by the FWER, by the desk rule, and by checks F1 and F2 - additionally integrates the estimated variances: with `U ~ N(0, 1/n_0)`, `V = S_0^2/n_0` and each `S_i^2` independent chi-square scaled by the true `sigma_i^2 = lambda_i sigma_0^2`,  P(max_i t_i <= c) = E_{U, V, S_1^2..S_K^2}[ prod_i Phi( sqrt(n_i/lambda_i) (U + c sqrt(S_i^2/n_i + V)) ) ]  (with `sigma_0^2 = 1`). Check F1 validates THIS form. Implement the conditional form only where you need the oracle's critical value at a fixed `lambda`.

**Check your implementation** before going further: with `K = 1`, `n_0 = n_1`, `sigma_0 = sigma_1`
(i.e. `lambda = 1`), the statistic is the pooled two-sample t and the root must reproduce the
Student `t_{2n-2}` quantile. Report the worst absolute deviation you see.

## 2. The admission region (data-dependent, certified coverage)

Given estimates `lambda_hat`, admit the box

    Lambda(lambda_hat) = prod_i [ lambda_hat_i * w_lo / y_hi_i , lambda_hat_i * w_hi / y_lo_i ],

with, for `nu_i = n_i - 1` and a **level scale** `s >= 0`,

    w_lo = chi2.ppf( (s*gamma/2) / 2 , nu_0 ) / nu_0
    w_hi = chi2.ppf( 1 - (s*gamma/2)/2 , nu_0 ) / nu_0
    y_lo_i = chi2.ppf( (s*gamma/(2K)) / 2 , nu_i ) / nu_i
    y_hi_i = chi2.ppf( 1 - (s*gamma/(2K))/2 , nu_i ) / nu_i

(these are the equal-tail chi-square levels with the budget split 1/2 to the shared benchmark scale
and 1/2 equally across candidates).

Since `lambda_i = lambda_hat_i * W / Y_i` with `W ~ chi2_{nu_0}/nu_0` and `Y_i ~ chi2_{nu_i}/nu_i`
independent, the **exact** coverage is the 1-D integral

    cov(s) = INT_0^inf  prod_i [ F_chi2( nu_i * w * a_hi_i , nu_i )
                               - F_chi2( nu_i * w * a_lo_i , nu_i ) ]
              * f_chi2( nu_0 * w , nu_0 ) * nu_0  dw ,
    a_lo_i = y_lo_i / w_hi ,   a_hi_i = y_hi_i / w_lo .

Coverage *decreases* in `s`. Use the **largest** scale with `cov(s) >= 1 - gamma + 1e-6` (bisection;
we get scales 12.7-37.8 depending on the design). Report the scale and `cov(s)` per cell — the original values are withheld with the other outputs.

## 3. The two rules

**Envelope rule (measured rule).** Form the candidate set as the Cartesian product over `i` of
`{lambda_lo_i, lambda_hat_i, lambda_hi_i}` when `K <= 3` (3^K vertices; at `K = 5` use the 2^K corners;
at `K = 8` the original used a reduced axis-wise set (`3 + 2K` points); the full `3^K` product or the `2^K` corners are equally valid candidate sets, and a replication should report which one it used). Set

    c_env(lambda_hat) = max over candidate vertices of c*_pop(vertex).

Reject when `max_i T_i > c_env`. This is the rule whose guarantee is
`FWER <= alpha + P(lambda not in Lambda)`.

**Plug-in rule (comparison).** Reject when `max_i T_i > c*_pop(lambda_hat)`.

## 4. What to measure, on which designs

Desk map: two families x six benchmark degrees of freedom.
`df_0 in {3, 5, 7, 11, 19, 35}`, `n = [df_0 + 1, 13, 13]`.
Family A: `lambda = (0.04, 0.04)`; Family B: `lambda = (1.0, 25.0)`.

For each of the 12 cells:

1. Draw 200 independent replication draws from the true model (benchmark variance 1):
   `s_0^2 = chi2_{nu_0}/nu_0`, `s_i^2 = lambda_i * chi2_{nu_i}/nu_i`, `xbar_i - mu_0` given `s_i^2`
   as in §1. **Pin the draw order, or the comparison is not draw-for-draw.** A seed does not by itself fix a replication: the order in which the variances are generated does. Use `rng = numpy.random.default_rng(seed)` and draw, per seed: `s_0^2 = rng.chisquare(nu_0)/nu_0`; then `s_i^2 = lambda_i * rng.chisquare(nu_i)/nu_i` for `i = 1..K` in that order; then any other variates you need (or integrate the means analytically as in the display above). The blind replication that used a different order reproduced every design-level quantity - `c_oracle`, the region scale, the exact coverage - to 4-5 significant figures, but its draw-level quantities (`c_env`, realised FWER, dominance counts) differ in the third digit, which is the signature of unpinned draws rather than of a disagreement about the method. **Use the seeds in `inputs/designs.json`** (listed explicitly per cell, 50 draws per seed) — do **not** reconstruct them from a formula: the original deposits record the *offsets* [101, 202, 303, 404], and a reconstruction that is off by the family index would silently break the draw-for-draw comparison so the comparison is
   draw-for-draw with the original.
2. Per draw compute `T = max_i T_i`, `lambda_hat`, `c_env` and `c*_pop(lambda_hat)`, and the
   **oracle** value `c*_pop(lambda_true)`.
3. Report:
   - realised FWER of the envelope rule and of the plug-in rule (fraction of draws rejected);
   - **paired dominance**: `P( T > c_env ) - P( T > c_oracle )` as a mean over draws of the
     per-draw difference `1{T > c_env} - 1{T > c_oracle}`, with its standard error. This is the
     decisive quantity; the acceptance rule used here is `mean + 2*SE <= 0`;
   - the mean `c_env`, mean `c*_pop(lambda_hat)`, and the population `c*_pop(lambda_true)`;
   - the fraction of draws for which `lambda_true` lies outside the admitted region (a `gamma` event);
   - the cost `c_env / c*_pop(lambda_true) - 1` in percent.

## 5. Robustness checks to report (each must pass or be reported as a failure)

- **F1 (identity):** §1's check against the Student t quantile, `K=1`, `n_0=n_1 in {10, 31, 61}`,
  `alpha in {0.05, 0.01}` — worst absolute deviation.
- **F2 (known-control limit):** `sigma_0 -> 0` (use `sigma_0 = 1e-5`) with `K = 3, 5, 8`: `F(c,.)`
  must approach `prod_i F_{t_{n_i-1}}(c)`; report the worst deviation at `c in {2.0, 2.5, 3.0}`.
- **F3 (grid refinement):** recompute §1's integral with 10, 32 and 48 panels per axis; the root must
  **converge**, not wander. Report the movement between the two finest grids.
- **F4 (Monte Carlo):** simulate the null law of `T` directly (>= 1e6 draws) at one unbalanced design
  and check the empirical size at `c*_pop` against `alpha`. Report the z-score.

## 6. What is *not* required

You do not need to reproduce the region-shape optimisation, the AR(1) extension, the calibration
gate, or the factor (D2) work. Those are separate stages and are reported separately.
