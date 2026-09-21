# SPEC — 03-certified-supremum: theorem-grade many-to-one comparison under estimated variance ratios

Version: **v1 (P0)** — finalised at P5. This document is the standalone specification: it must be
sufficient for a clean-room implementation of the certified rule from scratch.

---

## 1. Model and estimand

One benchmark group (index `0`) and `K` candidate groups (index `1..K`). Independent Gaussian
observations across groups,

    X_{i,1..n_i} iid N(mu_i, sigma_i^2),  i = 0..K,  groups independent,

with `mu_i`, `sigma_i^2` unknown. **This release covers ONE-SIDED promotion tests only.** The
one-sided comparison of each candidate against the benchmark uses the **separate-variance Welch
statistic**

    T_i = (Xbar_i − Xbar_0) / sqrt(S_i^2/n_i + S_0^2/n_0),

with `S_i^2` the usual unbiased sample variance (`n_i − 1` degrees of freedom). The estimand tested is
the mean difference `mu_i − mu_0` (one-sided promotion rule: reject `H_0i: mu_i ≤ mu_0`).

Population variance ratios and their estimates:

    lambda_i = sigma_i^2 / sigma_0^2,      lambda_hat_i = S_i^2 / S_0^2.

Goal: a **data-computable** threshold `c_hat` with

    FWER = P_lambda{ max_i T_i > c_hat } <= abar      for every admissible lambda,

at a user-chosen hard ceiling `abar` (default 0.025), without simulation in the decision path.

The two-sided analogue needs the joint law of `max_i |T_i|`, which is **not implemented and not
approximated** in this release; a two-sided familywise guarantee requires its own derivation and
certified implementation (see section 13, and gate E9 which fails theorem-tier text that implies
otherwise).

## 2. The exact fixed-ratio null law (already established; `02-bounded-rule/core.py`)

Conditioning on the benchmark's **two** sufficient statistics — its mean and its variance estimate —
makes the `K` comparison events **exactly independent**. With

    W = S_0^2/sigma_0^2 ~ chi^2_{n_0−1}/(n_0−1),   Y_i = S_i^2/(lambda_i sigma_0^2) ~ chi^2_{n_i−1}/(n_i−1),

independent of each other, the null distribution of `max_i T_i` reduces to an outer expectation over
`(Z_0, W)` (2 dimensions, independent of `K`) and `K` independent inner expectations. Precisely:

    F(c; lambda) = E_{S}[ F_cond(c; lambda; S) ],
    F_cond(c; lambda; S) = E_{U}[ prod_{i=1}^{K} Phi( sqrt(n_i/lambda_i) (U + c d_i) ) ],
    U ~ N(0, 1/n_0),   d_i = sqrt(S_i^2/n_i + S_0^2/n_0),

where the outer expectation is over the chi-square sampling law of the sample variances
`S = (S_0^2, S_1^2, ..., S_K^2)`, and the inner expectation is the conditional law given the
benchmark's two sufficient statistics (its mean and its variance estimate). `F_cond` is the form
evaluated throughout this programme (Section 3); `core.HeteroMaxT` evaluates the population law
`F` directly by its own composite-panel quadrature, and the draw-average identity
`E_S[F_cond(c*; S)] = 1 - alpha` is checked in `test_enclose.py` (falsifier F8).

Implemented as `core.HeteroMaxT` (production) and cross-validated against `fastlaw.py`, an independent
transformation/panel architecture, in `check_core.py` (kernel regression max abs difference 0.000e+00).
The **oracle** critical value at a known ratio vector is `c_star(lambda) = F^{-1}(1 - alpha; lambda)`.

## 3. Reduced form actually used by the certified path

For a realised draw, let `d_i` be the estimated standard error of contrast `i`,

    d_i = sqrt(S_i^2/n_i + S_0^2/n_0),

and let `s_i = sqrt(n_i/lambda_i)`. Conditional on the draw's benchmark sufficient statistics, the
exact law is evaluated as the **one-dimensional outer integral**

    F(c; lambda) = E_U[ prod_{i=1}^{K} Phi(s_i (U + c d_i)) ],     U ~ N(0, 1/n_0),

equivalently `F(c; lambda) = int G(z; c, lambda) dPhi(z)` with `z = sqrt(n_0) U` and

    G(z; c, lambda) = prod_{i=1}^{K} Phi(s_i (z/sqrt(n_0) + c d_i)).

This is the form used by `02-bounded-rule/cert_sup.py` and by every certified object here. `Phi` is the standard
normal CDF. (The reduction is exact under the model; the draw's own sufficient statistics are held
fixed, which is legitimate because the region is a function of the draw.)

## 4. Ratio region and its coverage

    lambda_hat_i = S_i^2/S_0^2,   W = S_0^2/sigma_0^2,   Y_i = S_i^2/(lambda_i sigma_0^2)

with `W`, `Y_i` independent scaled chi-squares, and `lambda_hat_i = lambda_i Y_i / W`. The region is
the rectangle `Lambda(lambda_hat) = prod_i [lambda_lo_i, lambda_hi_i]` with

    lambda_lo_i = lambda_hat_i * w_lo / y_hi_i,     lambda_hi_i = lambda_hat_i * w_hi / y_lo_i,

where the four multipliers come from **separate equal-tail chi-square quantiles** (NOT from the
quantiles of the ratio `Y_i/W` - a literal reading of "quantiles of the ratio" under-covers by a
factor of about 17, as the P0 verifier round demonstrated). With `nu_i = n_i - 1`, benchmark share
`w_frac` (default `0.5`), level scale `s` (default: the tight scale, below) and `K` candidates:

    g_W = s * gamma * w_frac,     a_i = s * gamma * (1 - w_frac) / K,
    w_lo = chi2_{nu_0}.ppf(g_W/2) / nu_0,        w_hi = chi2_{nu_0}.ppf(1 - g_W/2) / nu_0,
    y_lo_i = chi2_{nu_i}.ppf(a_i/2) / nu_i,      y_hi_i = chi2_{nu_i}.ppf(1 - a_i/2) / nu_i,

exactly as implemented in `env_rule._levels` (two-sided `gamma/2` per factor, half of the budget to
the shared benchmark factor and half split over the candidates). The coverage identity is then

    P{lambda in Lambda} = int_0^inf prod_i [ F_{nu_i}(nu_i w a_hi_i) - F_{nu_i}(nu_i w a_lo_i) ]
                          * chi2_{nu_0}.pdf(nu_0 w) * nu_0 dw,     a_lo_i = y_lo_i/w_hi,
                          a_hi_i = y_hi_i/w_lo,

because `lambda_i in [lambda_hat_i w_lo/y_hi_i, lambda_hat_i w_hi/y_lo_i]` is equivalent to
`Y_i/W in [y_lo_i/w_hi, y_hi_i/w_lo]`, and given `W` the `Y_i` are independent. This is evaluated
deterministically by 1-D quadrature (`env_rule.region_coverage`, no Monte Carlo power problem).

**Level scale.** Coverage *decreases* in the scale `s`, so the tightest region uses the largest
feasible scale: `env_rule.tight_scale` bisects for the largest `s` whose exact quadrature coverage is
`>= 1 - gamma + margin` (`margin = 1e-6`), never using the union bound. For `n = [13, 13, 13]`,
`gamma = 5e-4` this gives `s = 35.877` and coverage `0.9995010000000001 = 1 - gamma + margin` exactly;
the union-bound allocation (`s = 1`) over-covers by roughly 17x in miss probability, i.e. it is far
more conservative than necessary. Independently measured on 240,000 draws
(`02-bounded-rule/deposits/cert/miss_probe.json`): pooled observed miss rate 4.96e-04 against the designed
5.00e-04, worst |z| 1.13.

**Decision-relevant identity**: `FWER <= (level used inside the region) + P{lambda not in Lambda}`.
The certified programme uses the internal level `alpha_0 = abar − gamma`, so

    FWER <= alpha_0 + gamma = abar.

## 5. The enclosure (Workstream B / A — the core of this programme)

For a box `Lambda = prod_i [lambda_lo_i, lambda_hi_i]` define the s-range
`s_lo_i = sqrt(n_i/lambda_hi_i)`, `s_hi_i = sqrt(n_i/lambda_lo_i)` and, writing
`x_i(z) = z/sqrt(n_0) + c d_i`,

    Psi_min,i(z) = Phi(s_lo,i x_i) if x_i >= 0 else Phi(s_hi,i x_i)
    Psi_max,i(z) = Phi(s_hi,i x_i) if x_i >= 0 else Phi(s_lo,i x_i)
    G_min(z) = prod_i Psi_min,i(z),      G_max(z) = prod_i Psi_max,i(z).

Both are continuous and **increasing in z**; `0 <= G_min <= G <= G_max <= 1` for every `lambda in Lambda`.
For any partition `z_0 < ... < z_M` with cell masses `dPhi_j = Phi(z_{j+1}) − Phi(z_j)`,

    F_lo(c; Lambda) = sum_j G_min(z_j) dPhi_j + G_min(z_M) (1 − Phi(z_M))
    F_hi(c; Lambda) = sum_j G_max(z_{j+1}) dPhi_j + G_max(z_0) Phi(z_0) + (1 − Phi(z_M))

and these are **rigorous**:

    F_lo(c; Lambda) <= min_{lambda in Lambda} F(c; lambda)
    max_{lambda in Lambda} F(c; lambda) <= F_hi(c; Lambda).

Proof sketch: an increasing integrand's left Riemann sum lower-bounds its Stieltjes integral and its
right sum upper-bounds it; on the lower tail `(−inf, z_0]` the integrand is at most `G_max(z_0)` for the
upper bound and at least `0` for the lower; on the upper tail `[z_M, inf)` it is at least `G_min(z_M)` and
at most `1`. A degenerate box (`s_lo = s_hi`) gives `G_min = G_max = G`, so the two bounds bracket
`F(c; lambda)`.

**No derivative bound, Lipschitz constant, mesh-covering radius, or vertex-attainment property is used.**

Arithmetic: every endpoint used in a *certified* claim is computed in Arb ball arithmetic
(`python-flint`, `arb`), with `Phi(a) = (1 + erf(a/sqrt(2)))/2`; the reported endpoint is the ball's
`.lower()`/`.upper()`. Float64→Arb conversion is exact, and `c`, `alpha`, `n`, `d`, `lambda` are floats
converted exactly, so a certified claim is about exactly the number that is reported. Cells whose mass
falls below `1e-18` are dropped and the omitted mass is added to the enclosure width (rigorous, since
`0 <= G <= 1`). **The mass-floor device belongs to this Riemann path only**: the GL certificate path
has no mass floor, and its `settings` block records `mass_floor: null` (verifier finding #4).

**Certificate quadrature (`enclose_gl.py`).** The monotone Riemann sum above is exact but converges
slowly (measured: `width ~ N^-2`, 160k nodes for 2.6e-06), so it serves as the fast search path and as
an independent cross-check. The certificate uses n-node Gauss-Legendre on a panel layout covering
`[-9, 9]`, with the remainder bounded by the Bernstein-ellipse theorem: for `f` analytic in `E_rho`
with `|f| <= M_rho` there, `|E_n| <= C * M_rho / (rho^{2n} (1 - rho^-2))` with
`C = (128/15) * GL_SAFETY` (`GL_SAFETY = 2`, i.e. twice Trefethen's constant). `M_rho` is obtained
rigorously **without sampling**: one `acb` ball evaluation on the rectangle inscribing the ellipse,
with each factor bounded by `max(|Phi(s_lo x)|, |Phi(s_hi x)|)`. Float64 node/weight rounding is
carried inside the interval arithmetic by the NODE_RAD / WEIGHT_RAD radii applied at the node
evaluations (the earlier `2 * GL_EPS * (b - a) * M_panel` claim was unfounded and removed, pass-3
finding #20c). The analyticity term is capped by the panel's own normal mass. Measured width at the
desk cell: `7.1e-14` in
probability, i.e. at the Arb/float64 resolution floor; `Delta p` is reported with that floor stated.

## 6. Certified critical values

| object | definition | property |
|---|---|---|
| `c_ub(Lambda)` | smallest float64 `c` for which `F_lo(c; Lambda) >= 1 − alpha_0` was **verified** | `c*(lambda) <= c_ub(Lambda)` for **every** `lambda in Lambda` |
| `c_lb(lambda_w)` | largest float64 `c` for which `F_hi(c; {lambda_w}) < 1 − alpha_0` was **verified** | `c*(lambda_w) > c_lb` at that **witness point** |
| `[c_L, c_U]` at a point | `[c_lb(lambda), c_ub({lambda})]` | `c*(lambda) in [c_L, c_U]`; width reported |

The **shipped** implementations are `enclose_gl.c_upper_box_gl` and `enclose_gl.c_lower_s_gl` (the GL
certificate path). `enclose.c_upper_box` and `enclose.c_lower_s` are the Riemann-path siblings used
for the float64 search and as an independent cross-check; they may be substituted only where their
wider enclosure width is acceptable and is reported (their width is ~1e5 times larger in probability
at the desk cell and floors at the partition budget rather than at the Arb resolution).

Note the directions: upper bounds on the supremum come from a **lower** enclosure of `min_lambda F`;
lower bounds on the supremum come from a certified evaluation at a **witness point**. (`F_hi` over a
non-degenerate box certifies bounds on `min_lambda c*`, not on the supremum.)

## 7. Certified continuous supremum (sub-box branch-and-bound)

Inputs: region `Lambda` (Section 4), drawn `d_row`, `n`, `alpha_0`, predeclared tolerance `tol = 1e-4`.

1. Queue = {`Lambda`}; witness set = vertices, face midpoints, centre, mixed corners suggested by the
   implicit-function derivative, plus quasi-random interior starts.
2. For each box: `c_ub(box)` (float64 search, then Arb verification for retained boxes) and
   `c_lb` at its best witness.
3. `L_lo = max` of certified witness lower bounds.
4. Prune a box when `c_ub(box) <= L_lo + tol` — the pruning proof is the certified comparison.
5. Otherwise split along the coordinate maximising `|dc*/d log lambda_i|` (tie-break: widest log-width).
6. Terminate when every box is pruned or the node budget is reached.
7. Report `c_L = L_lo` and **`c_U = max over ALL boxes of c_ub(box)`** (pruned boxes included; each
   pruned box has `c_ub <= L_lo + tol <= c_U`, so it can never dominate), with
   **gap = c_U - c_L** reported honestly. If the node budget is exhausted before every box is pruned
the gap may exceed `tol`; that is reported as a *failure to certify at the declared tolerance*
(region, achieved gap, cost) and never masked by substituting `L_lo + tol` for `c_U` (the P0
verifier round caught exactly that error in the first draft: substituting the prune threshold
returns a value *below* a box's true upper bound, i.e. a non-conservative threshold in the one case
where it matters). Together with the certificate tree (boxes, bounds, split order, pruning proofs,
termination reason, node and bisection counts, seconds).
8. **Shipped threshold**: `c_hat = c_U`. A region that cannot be certified within budget is *named*
   with its achieved gap, and falls back to `02-bounded-rule/cert_sup.c_rigorous`, labelled as such.

Because the statement is about the continuous supremum, it does **not** depend on vertex attainment.

## 8. Decision tiers (`tiers.py`, mechanically enforced)

| tier | meaning |
|---|---|
| `oracle_exact_law` | the fixed-ratio law / oracle critical value, exact to stated numerical tolerance |
| `plugin_measured` | plug-in rule measured on declared simulations; no finite-sample guarantee |
| `calibrated_measured` | calibrated/bootstrap rule measured under a declared budget; no theorem |
| `envelope_measured` | vertex-max envelope: empirically stress-tested, **not** certified over the continuum |
| `region_sup_certified` | certified continuous supremum over a coverage-certified region; carries the `alpha_0 + gamma` theorem |
| `global_sup_certified` | a formally valid conservative tier over the whole admissible domain, if implemented |

Rules: every deposit carries `tier`; every report states it; the gate rejects the words *guaranteed*,
*certified* and *controlled* in text attached to a measured-tier deposit. A calibrated rule must never
silently replace the certified threshold in a governance gate.

**Register policy (as implemented).** `tiers.gate_all` distinguishes *result reports* from
*documentation records*: a record is `SPEC.md`, `PROGRESS.md`, `ENVIRONMENT.md`, `README.md`,
`PRIOR-ART.md`, any `PHASE*`/`REPORT*` document, or any path with a `verify` component; detection is
path-form independent. For an UNDECLARED record the reserved words are reported as `record_hits` and
counted, not enforced - these documents must be able to name the certified tier and to quote the
claims the gate polices, and quoting an overclaim as a defect is not an overclaim. An explicit
`tier:` line **defeats** that demotion, so a document that declares a measured tier is gated exactly
like a result report. `run_modes.quick` prints the record-hit count and the manifest freshness, so
both are visible rather than silent.

Exemptions inside a clause are exactly two: (i) a reserved token lying **inside a reference phrase**
that names the certified tier (spans computed per clause), and (ii) a negation within 25 characters
before the token, where a comma, semicolon, full stop or exclamation mark inside that window breaks
the link (a new sub-clause). Two earlier versions are recorded here because re-verifiers broke them: a
clause-level "this clause is about the gate" exemption (six false negatives, e.g. "we certify the
envelope, a finding of this audit") and a 45-character negation lookback (which exempted "not that it
matters, we guarantee FWER control").

## 9. Error budget (per certificate)

`Delta p = F_hi - F_lo` at the certificate point is the reported probability-scale uncertainty, and
`c_equiv = Delta p / slope` (a certified local bound on `dF/dc`) is its critical-value equivalent.
The decomposition differs by path; `error_budget.py` (P3 step 16) assembles the per-certificate table:

* **GL certificate path** (the shipped path, `enclose_gl.py`): the analyticity remainder
  `C * M_rho / (rho^{2n}(1 - rho^-2))` per panel (`panel_error` minimises it over the t grid, with
  `M_rho` from the inscribing rectangle); float64 node/weight rounding
  carried by the `NODE_RAD`/`WEIGHT_RAD` radii (which set the floor - measured `7.1e-14` for the desk
  cell); the panel mass cap; the explicit tail terms; the Arb ball radii; and the bisection bracket width.
* **Riemann enclosure path** (`enclose.py`, the float64 search and the independent cross-check):
  partition width `sum_j dPhi_j (G_max(z_{j+1}) - G_min(z_j))` plus tail terms; the mass-floor
  omission (`dropped_up`, added to both bounds); Arb ball radii; bisection bracket width.

Both paths report a `settings` block (GL order, safety factor, `NODE_RAD`/`WEIGHT_RAD`, precision, panel edges, t
grid, `z_max`, and `mass_floor: null` for the GL path) so that a certificate can be reproduced from
its own deposit (verifier finding N2-b). **Predeclared**: `c_U - c_L <= 1e-4` absolute and
`Delta p <= 1e-6`, with the c-equivalent far below that. **Measured** at the desk cell:
`Delta p = 7.1e-14`, `c_equiv = 1.4e-12`, bracket width `6.9e-07`.

## 10. Review findings and their resolutions (honest record)

**P0 verifier round** (one read-only oracle on spec/code consistency, one code auditor; both ran
before any P1 result existed). Blocking findings, all now fixed:

| finding | resolution |
|---|---|
| SPEC section 4 described the region multipliers as "quantiles of `Y_i/W`"; the code uses ratios of separate quantiles. A literal reading under-covers (measured 0.99158, 17x the `gamma` budget) | section 4 rewritten above with the exact construction, the coverage identity and the level-scale rule |
| SPEC section 7 step 7 returned `min(c_ub, L_lo + tol)`, which is **not** an upper bound when the node budget is exhausted | corrected to `max` over all boxes with the gap reported honestly (above) |
| SPEC section 2 displayed a formula containing a literal ellipsis and an undefined `a_i` | replaced by the precise conditional/population identity (above) |
| `write_deposit` could leave a truncated deposit; the reader-side strict-JSON check accepted bare `NaN`/`Infinity` | atomic write plus a NaN sanitiser in `tiers.write_deposit`; `tiers.load_strict` rejects non-finite constants |
| the reserved-word gate had false negatives and false positives, and had **no call sites** | gate rewritten clause-wise with precise exemptions; `tiers.gate_all` runs over every deposit and report and is wired into `run_modes.quick` |
| the traceability gate accepted a missing field and never compared the claimed number with the deposited value | `manifest.check_map` now requires a concrete field and compares numeric claims against the deposit (non-comparable claims are listed) |
| `manifest.scan` silently dropped unreadable files; no staleness check existed; the stored manifest was already stale | failures recorded explicitly, cycle guard added, `manifest.verify()` reports changed/added/removed |
| `run_modes.quick` reported PASS when core checks were unavailable (fail-open) | fail-closed: an unavailable or failed check yields ATTENTION and a non-zero exit code |

**P1 falsifier F10** then caught a bug in the certificate error term itself: the float64 node/weight
rounding term used `M_rho` from the *largest* ellipse tried (up to `1e15` on far-tail panels), which
inflated a rigorous `1e-43` bound into `2.06` and made the enclosure useless at GL order 24. The term
now uses the integrand magnitude on the panel. Effect: enclosure width `5.3e-08 -> 7.1e-14` (and then, after the pass-3 radius fix, the floor is set by NODE_RAD/WEIGHT_RAD), point
bracket width `2.35e-06 -> 6.9e-07` at the desk cell, and the width is now monotone in the GL order
(ratio `4.5e6x` from order 6 to a mid order).

## 11. Required theorem statement (printed by every theorem-grade report)

> Under independent Gaussian samples across groups, with the stated sample sizes and the specified
> separate-variance Welch statistics, this decision uses internal level `alpha_0` and a ratio-region
> miss budget `gamma`. The certified region coverage is at least `1 − gamma`, and the reported threshold
> is a certified upper bound on the oracle critical value throughout that region. Therefore the
> familywise error rate is at most `alpha_0 + gamma = abar`.

Model exclusions (always printed): serial dependence; cross-group dependence beyond the model; factor
exposure; non-Gaussianity; rolling re-estimation; research-process multiplicity.

## 12. Deposits (field conventions)

Strict JSON only (no bare `NaN`/`Infinity`; non-applicable values are `null` with a `notes` entry).
Every deposit carries a `tier` (Section 8), enforced by `tiers.recheck_deposits`, and every deposit
writer goes through `tiers.write_deposit`, which sanitises non-finite values to `null`, records their
paths in `notes`, and writes atomically.

**Certificate deposits** (`region_sup_certified`) additionally carry: `alpha0`, `gamma`, `n`, the
region (bounds or `lambda_hat` plus the level scale), the certified objects with their enclosure
widths, the error budget of Section 9, the `settings` block that produced them, a certificate
identifier, `seconds`, and a `note`. **Measured-tier deposits** carry whatever their claim needs plus
the tier; the field list above is recommended rather than enforced for them. (Verifier finding #5:
the phase-1 self-test deposits predate this convention - they carry `alpha` rather than `alpha0` and
no region block, which is why the convention is stated per tier rather than as a blanket rule.)

## 13. Known failures, rejected approaches, and field-by-field deposits

### 13.1 Known failures and rejected approaches (recorded, not hidden)

| approach | verdict | reason |
|---|---|---|
| **Lipschitz-grid certificate** (plan section 2.1-2.3) | **REJECTED as a certificate** | the predeclared materiality test failed: the covering-radius correction bottoms out at `1.34e-01` against the predeclared `tol = 1e-04` on every feasible grid (`deposits/deriv_diag.json`). The derivative itself is correct (finite-difference agreement `5e-06`), so the route survives only as a split-direction diagnostic. |
| **Vertex-attainment theorem** | **NOT ATTEMPTED** | the plan allows either a vertex theorem or a certified continuous supremum; the certified supremum is sufficient and is what is shipped. Vertex attainment is measured (12/12 desk cells, 14/14 in the headline comparison) and reported as evidence, never as a theorem. |
| **Default-panel reference solver** | **REJECTED for the falsifiers** | `core.cond_noreject` with its default 40 panels is wrong by `-5.05e-07` at sharp configurations (`s ~ 114`); `npanels >= 2000` agrees with a 40-digit mpmath evaluation to `~6e-16`. The falsifiers use the refined reference. |
| **Clause-level "about the gate" exemption** in the tier gate | **REJECTED** | it produced six false negatives (e.g. "we certify the envelope, a finding of this audit"). Only two exemptions remain: a reserved token inside a reference phrase, and a negation within 25 characters with sub-clause breaks. |
| **The unfounded `eps * M` rounding term** | **REJECTED and removed** | it had no Lipschitz bound behind it and breached the panel bound by ~13x at a sharp-corner panel. The float64 node/weight perturbation is now carried by Arb radii inside the interval arithmetic. |
| **Comparing the conditional `c_U` against population quantities** | **REJECTED** | three separate instances were found and removed (the headline, the dominance test, the desk-map price). The rule is conditional; every comparison must be like-for-like. |
| **Two-sided levels** | **NOT COVERED** | the exact-law engine is the one-sided max statistic; the two-sided analogue needs the joint law of `|T_i|`, which is not implemented. It is not approximated. |
| **Homoscedastic comparison against published Dunnett tables** | **NOT INCLUDED** | `scipy` has no Dunnett routine; the population solver is instead validated against an independent transformation/panel architecture (02-bounded-rule kernel regression, max absolute difference `0.000e+00`). |

### 13.2 Deposit fields (field-by-field)

Every deposit is strict JSON (no bare `NaN`/`Infinity`; non-finite values are `null` with a `notes`
entry) and carries:

| field | meaning |
|---|---|
| `tier` | one of the six tier labels of section 8; written by `tiers.write_deposit`, which also records `tier_claim` |
| `notes` | honest qualifications, including the paths of any value replaced by `null` |
| `alpha0`, `gamma`, `abar` | the levels actually used (the certified tier always uses `alpha0 = abar - gamma`) |
| `tol`, `dp_target` | the predeclared materiality targets (`1e-4` on `c`, `1e-6` on `Delta p`) |
| `n`, `lam` / `lambda_hat`, `d_row` | the design and the draw |
| `region_lo`, `region_hi`, `coverage_exact`, `scale` | the ratio region and its exact quadrature coverage |
| `c_L`, `c_U`, `gap`, `certified_at_tol` | the certified bracket for the region supremum and the achieved gap |
| `settings` | the certificate-affecting constants: `gl_n`, `gl_safety`, `node_rad`, `weight_rad`, `prec`, `z_max`, `n_panels`, the panel edges, the `t_grid`, `mass_floor` (null for the GL path) |
| `witnesses`, `argmax_witness` | the certified witness lower bounds and which one attained `c_L` |
| `tree` | the certificate tree: pruned leaves with their `c_ub` and `pruned_at_L`, and the split records (axis, value, level) |
| `seconds`, `termination` | cost and why the run stopped ("all boxes pruned" / "split budget exhausted") |
| `verdict`, `verdict_scope` | the falsifier batteries' outcomes and what the composite verdict includes |

A certificate deposit is reproducible from its own `settings` block: the same `gl_n`,
`node_rad`/`weight_rad`, `prec`, panel layout and `t_grid` reproduce the same enclosure.

## 14. Reproducibility: the draw recipe and every underspecified constant

The blind replication (P5) reproduced the draw-INDEPENDENT half of this specification exactly - the
region multipliers, the exact coverage quadrature and the tight-scale root to 10+ digits on every
cell, and all twelve sealed deposit hashes - but could not reproduce any draw-DEPENDENT quantity,
because the following were underspecified. They are now pinned here; each is recorded as a defect the
replication exposed.

### 14.1 The draw recipe (defect D1/D8)

    import numpy as np
    rng = np.random.default_rng(seed)                      # numpy PCG64, not legacy RandomState
    s0_sq = rng.chisquare(n[0] - 1) / (n[0] - 1)           # sigma_0 = 1
    si_sq = [rng.chisquare(n[i] - 1) / (n[i] - 1) for i in 1..K]
    d_i   = sqrt(si_sq[i-1] / n[i] + s0_sq / n[0])         # the ESTIMATED SE of contrast i

Groups are drawn in increasing index order, one call per group, `sigma_0 = 1` (all variances are in
units of the benchmark's), and the seed is recorded in every deposit. `d_row` is always the
*estimated* standard error of the realised draw - never the population value
`sqrt(sigma_i^2/n_i + sigma_0^2/n_0)`. That distinction is the single most common source of
mismatches: the conditional law is evaluated at the former.

### 14.2 Certificate points (defect D2)

A certificate point in the kit is an `(n, lambda, c)` triple at which the *enclosure* is evaluated:
`F(c; lambda)` for the conditional law at the draw implied by the seed. The points are illustrative
and need not appear in any deposit; the deposit rows are the twelve desk cells plus the designs listed
in the branch-and-bound and root-certificate deposits. The conditional critical value at a point comes
from bisecting `F` at that point, not from reading a deposit.

### 14.3 Bisection, Delta p and the enclosure floor (defect D3/D4)

* Bisection stops when the bracket is at most the `tol` argument (default `1e-6`); a reported bracket
  width of `5e-07 ... 1.6e-06` is the stopping rule, not the `1e-15` float64 resolution that
  "smallest float64 c" might suggest.
* `Delta p` is `F_hi - F_lo` at the RETURNED `c`, on the GL certificate path. For a degenerate box the
  mathematical width is zero; the reported `7e-14`-scale value is the path's floor, set by
  `2 * GL_SUM_SLACK = 2e-14` plus the `NODE_RAD`/`WEIGHT_RAD` radii (1e-15 / 4e-15) that carry the
  float64 node/weight perturbation inside the interval arithmetic. A replica using plain float64 will
  report a smaller floor; that is a path difference, not a disagreement.
* The tight-scale search stops when the exact coverage is within `1e-9` of `1 - gamma + margin`
  (`margin = 1e-6`), or when the bracket collapses (defect D7).

### 14.4 Falsifier designs (defect D6)

`replicate/FALSIFIERS.md` carries the explicit `(n, lambda, c, seed)` for every falsifier, so a
replica can run them without inventing designs.

### 14.5 What is draw-independent (defect D5, benign)

The region, its coverage and the tight scale do NOT depend on the draw: `lambda_hat` cancels exactly
in the coverage integral, so they are functions of `(n, gamma, w_frac, scale)` alone. That is why the
replication matched them to ten digits while failing on everything downstream of the draw.
