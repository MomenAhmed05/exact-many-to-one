# REPORT, the efficiency frontier, and the desk procedure that comes out of it

Work area `04-efficient-rule/`.  Every number below is read from a strict-JSON deposit under
`deposits/`; the tier of each claim is the tier of the deposit it comes from.

---

## 1. The question, and the answer

The programme that preceded this one (`03-certified-supremum`) settled **validity**: a finite-sample FWER-valid rule
for unknown variance ratios exists, and its supremum over the nuisance region is certified rather
than measured.  What it did not settle is whether that rule is any use.  Measured against the
oracle it is not: at a 36-month benchmark its threshold is about three times the oracle constant,
and the only power measurement in that programme (`03-certified-supremum/deposits/paired_sim.json`) gives a maximum
power of **2.35%** across seventeen scenarios, zero in eleven of them.

This programme asks whether that is the price of validity or the price of one particular
construction.  It is the latter.

**The result.**  A rule with the same hard guarantee, costing **1.5 to 2.6 percentage points of
power** at two, three and five candidates once the smallest group has seven or more degrees of
freedom:

| | region rule (`03-certified-supremum`) | this rule |
|---|---|---|
| guarantee | `FWER <= abar`, certified | `FWER <= abar`; certified at K = 1, measured above |
| threshold at a 36-month benchmark | ~3x the oracle | the oracle constant at the estimated ratios, at a level ~2% below nominal |
| power | `<= 2.35%` across 17 scenarios, 0 in 11 | 1.5-2.6 pp below the oracle at K = 2, 3, 5 |
| decision-time cost | branch-and-bound, tens of seconds | one exact-law solve, ~0.2 s |

## 2. Why the region rule was paying too much

The certified region rule protects, on **every** draw, against **every** ratio in an admission
region.  The guarantee it has to meet is an average over draws:

    sup_lambda  E_S[ P_lambda( max_i T_i > c(S) | S ) ]  <=  abar .

Protecting conditionally on every draw is strictly more than that requires, and the surplus is not
small: the level a benign draw does not spend is simply thrown away instead of being made available
to a hostile one.

The rule here spends the budget where it is worth something.  It is the plug-in the paper already
studies - the exact fixed-ratio constant evaluated at the estimated ratios - run at an **internal
level** `alpha0` chosen so that the displayed constraint holds exactly.  One scalar per design.  No
region, no per-draw supremum, no simulation.

## 3. The ledger (Workstream A) - `deposits/wa_ledger.json`

Four objects on one design set, all conditional on the same data and all targeting the same
rejection event, which is the plan's first rule.  At K = 1 every one of them is closed form, so the
ledger carries no Monte-Carlo error at all.  The region rule is included exactly rather than
approximately: at K = 1 the conditional oracle constant is increasing in the ratio, so the supremum
over the admitted interval sits at its upper endpoint and needs no branch-and-bound.

At equal variances (`lambda = 1`), mean threshold, its ratio to the oracle, the realised
unconditional FWER, and the power lost at the effect where the oracle reaches 80%:

| design | object | mean c | x oracle | FWER | power lost |
|---|---|---|---|---|---|
| df3 | oracle_conditional | 2.253 | 1.00 | 0.02500 | 0.0 pp |
| | plugin_measured | 2.416 | 1.20 | **0.03108** | 15.3 |
| | region_sup_conditional | 6.951 | **2.72** | 0.00076 | **64.3** |
| | **uncond_certified** | 2.731 | 1.36 | 0.02165 | **23.9** |
| df11 | oracle_conditional | 2.027 | 1.00 | 0.02500 | 0.0 |
| | plugin_measured | 2.078 | 1.05 | 0.02460 | 3.8 |
| | region_sup_conditional | 4.843 | **2.38** | 0.00035 | **68.2** |
| | **uncond_certified** | 2.099 | **1.06** | 0.02354 | **4.4** |
| df35 | oracle_conditional | 2.030 | 1.00 | 0.02500 | 0.0 |
| | plugin_measured | 2.071 | 1.05 | **0.02551** | 3.7 |
| | region_sup_conditional | 4.903 | **2.48** | 0.00018 | **72.1** |
| | **uncond_certified** | 2.095 | **1.06** | 0.02428 | **4.5** |

Three readings, none of them available before this ledger existed:

1. **The region rule spends almost none of its budget.**  At a 36-month benchmark its realised FWER
   is 0.00018 against a ceiling of 0.025 - seven tenths of one % of what it is allowed - and
   it gives up 72 of the 80 percentage points of power the oracle has.  Its threshold is 2.4 to 2.7
   times the oracle's at every design here, which reproduces the `03-certified-supremum` programme's own finding
   from an independent K = 1 implementation.
2. **The plug-in at the full ceiling is not valid**, and the ledger shows where: 0.03108 at three
   degrees of freedom and 0.02551 even at thirty-five.  It is close to the oracle in power precisely
   because it is over-rejecting.
3. **The new rule buys validity for 4.4 to 4.5 pp** at eleven degrees of freedom and above, against
   the region rule's 68 to 72.  That is the whole result in one line.

## 4. The rule, and what it costs (Workstream B)

### 4.1 K = 1, certified over the continuum, `deposits/wb_certified_rule.json`

The deployed object is a step function `c(lambda_hat)` on 60 bins spanning three decades.  Its
unconditional FWER is bounded **above** by a rigorous enclosure over the continuum of ratios, not on
a grid: the state space is partitioned on the probability scale so every cell's mass is an exact
difference of chi-square distribution functions, the integrand is replaced on each cell by a
monotone corner substitution, and the ratio interval is decided by branch-and-bound.

| design | `alpha0` | certified | sup upper bound | power loss at the oracle's 80% point (median / max) |
|---|---|---|---|---|
| df3 | 0.01585 | yes | 0.02500 | 23.9 / 40.3 pp |
| df5 | 0.02160 | yes | 0.02500 | 11.2 / 20.7 |
| df7 | 0.02328 | yes | 0.02500 | 7.3 / 13.3 |
| df11 | 0.02393 | yes | 0.02500 | 5.9 / 7.8 |
| df19 | 0.02389 | yes | 0.02500 | 4.4 / 7.3 |
| df35 | 0.02377 | yes | 0.02500 | 4.5 / 7.6 |
| df59 | 0.02365 | yes | 0.02500 | 5.4 / 7.8 |
| df35, 4-df candidate | 0.01856 | yes | 0.02500 | 21.5 / 30.5 |
| balanced df11 | 0.02389 | yes | 0.02500 | 5.8 / 7.8 |

The reference here is the **known-lambda conditional oracle**, `c = z_{1-abar} sigma / d`, which is
the strongest reference available: it is the test a desk would run if it knew the variance ratios
exactly.  Nothing that estimates them can beat it.

### 4.2 K = 2, 3, 5, measured, `deposits/wb_general.json`

| design | plug-in sup FWER at the ceiling | `alpha0` | power loss (median / max) |
|---|---|---|---|
| K2 df3 | 0.03434 (+37%) | 0.01573 | 10.6 / 12.6 pp |
| K2 df5 | 0.02832 | 0.02130 | 3.8 / 4.6 |
| K2 df7 | 0.02594 | 0.02381 | 2.2 / 2.4 |
| K2 df11 | 0.02642 | 0.02352 | 2.2 / 2.5 |
| K2 df19 | 0.02567 | 0.02427 | 1.5 / 1.9 |
| K2 df35 | 0.02571 | 0.02445 | 1.5 / 1.8 |
| K3 df11 | 0.02658 | 0.02369 | 1.7 / 2.2 |
| K3 df35 | 0.02697 | 0.02329 | 1.8 / 2.6 |
| K5 df11 | 0.02628 | 0.02338 | 1.6 / 2.5 |
| K5 df35 | 0.02611 | 0.02403 | 1.7 / 2.2 |

**The two tables are not on the same scale and must not be read as one.**  Above K = 1 the reference
is the *population* oracle constant at the true ratios - the paper's own oracle - which is weaker
than the conditional one used at K = 1, so the losses are smaller by construction.  The column that
is comparable across both is `alpha0`, and it tells the same story: the level price is 2 to 6 per
cent once the smallest group has seven degrees of freedom, and 35 to 40% at three.

The first column is worth its own line.  The plug-in run at the full ceiling is **not valid**: its
worst-case FWER over the ratio domain is 0.0257 to 0.0343 against a nominal 0.025.  The paper
reports +13.8% at one particular design; the supremum over the domain is +37%.

## 5. Is the residual loss removable?  (Workstream C)

`deposits/wc_frontier.json` answers this by optimising over the whole class of threshold functions
of the estimated ratio - piecewise constant on 36 bins - subject to the same exact uniform
constraint.  Because any rule that sees only `lambda_hat` is of that form up to bin resolution, the
optimum bounds what **any** such rule can achieve.

| design | calibrated plug-in loss | best in class | gain from optimising |
|---|---|---|---|
| df3 | 25.2 pp | 20.9 pp | **+4.24** |
| df5 | 12.8 | 11.5 | +1.30 |
| df7 | 8.6 | 8.1 | +0.49 |
| df11 | 5.7 | 5.8 | -0.02 |
| df19 | 4.6 | 4.6 | +0.01 |
| df35 | 4.4 | 4.4 | +0.06 |

**From eleven degrees of freedom upward the calibrated plug-in is optimal in its class**: a full
numerical optimisation over every threshold function of the estimated ratio finds nothing.  The
4 to 6 pp that remains is therefore an information constraint, not a procedural failure - the
programme's terminal conclusion B, with numbers.

Below eleven there is removable slack, and the optimiser removes it: 4.2 pp at three degrees of
freedom.  The optimised step rule is deposited per design and is the better rule there.

### Why the loss cannot be driven to zero

Write `rho = lambda n_0/n_1`.  The conditional rejection probability of threshold `c` is
`Phi(-c d / sigma)`, and

    sigma^2 / d^2 = (rho + 1) / (rho u_1 + u_0) .

A rule attains the oracle exactly if and only if it can hold `c d / sigma` at the constant
`z_{1-abar}`, which requires knowing `(rho u_1 + u_0)/(rho + 1)`.  The data supply only
`rho_hat = rho u_1/u_0`: one equation, two unknowns.  No function of the estimated ratio can
equalise the conditional level across states, so some states are necessarily over-protected and the
uniform constraint binds at the least-favourable ratio.  That is the structural reason the frontier
sits where it does, and it is why more optimisation does not help above eleven degrees of freedom.

## 6. The planning surface (Workstream D), `deposits/wd_designmap.json`

A 7 x 5 grid over benchmark and candidate degrees of freedom, every cell calibrated and priced.

**The binding quantity is the smaller group, not the benchmark.**  Grouping the 35 cells by
`min(nu_0, nu_1)` leaves a spread in the worst-case loss of 13.1, 4.8, 1.1, 0.83, 0.12, 0.10 and
0.08 pp as that minimum runs 3, 5, 7, 11, 19, 23, 35.  Grouping the same cells by the benchmark's
degrees of freedom alone leaves 11.6, 17.0, 23.6, 32.3, 38.1, 42.0 and 43.4 - the spread *grows*.
The design map is symmetric to the last digit under swapping the two groups, which is the cleanest
possible statement of it: a 35-df benchmark with a 3-df candidate loses 31.1/44.1 pp, exactly what a
3-df benchmark with a 35-df candidate loses.

**This corrects a claim in the paper.**  Section 10.3 places the usable boundary in the benchmark's
degrees of freedom, because the sweep behind it varied only the benchmark's sample size.  It is the
smaller group that binds, and a long benchmark does not rescue a short candidate.

**The planning answer a desk can act on:**

| power budget | minimum degrees of freedom on the *smaller* group |
|---|---|
| lose at most 5 pp | 19 (about 20 monthly observations) |
| lose at most 10 pp | 11 (about 12) |
| lose at most 25 pp | 5 (about 6) |
| lose at most 50 pp | 3 |

## 7. A general-K class that can be certified, `deposits/wb_percand.json`

The joint plug-in couples all the coordinates, so its unconditional FWER is an expectation over
K+1 chi-squares with no tractable enclosure: it is measured above K = 1 and stays measured.

A **per-candidate** threshold, `reject i iff T_i > c(lambda_hat_i)`, restores the paper's own
factorisation for the FWER itself.  Conditioning on the benchmark pair, the candidates are
independent and each factor depends on one coordinate, so the unconditional FWER becomes a
two-dimensional outer rule over K one-dimensional inner rules - the same structure and the same cost
as the null law, at any K, with no draws.  That is what makes a rigorous enclosure affordable above
K = 1, and it is the route to certifying the general-K rule.

It is not free.  The restriction costs real power, and the measurement is reported rather than
buried:

| design | per-candidate `alpha0` | its worst loss | joint rule's `alpha0` | joint rule's worst loss |
|---|---|---|---|---|
| K2 df3 | 0.00999 | 24.0 pp | 0.01573 | 12.6 pp |
| K2 df7 | 0.01215 | 11.9 | 0.02381 | 2.4 |
| K2 df11 | 0.01222 | 9.9 | 0.02352 | 2.5 |
| K2 df35 | 0.01224 | 6.9 | 0.02445 | 1.8 |
| K3 df11 | 0.00807 | 16.2 | 0.02369 | 2.2 |
| K3 df35 | 0.00808 | 11.2 | 0.02329 | 2.6 |

The per-candidate rule needs roughly half the internal level and gives up three to four times the
power.  The reason is visible in the level: with per-candidate thresholds the rule cannot use one
candidate's variance to sharpen another's, so it must protect each comparison separately and the
multiplicity is paid in full.  The reading that follows
is that the joint rule is the one to run and the per-candidate rule is the one that can be proven;
closing that gap is the main open item.

## 7a. What an ex ante ratio restriction buys, `deposits/wd_restricted.json`

The obvious remedy for a uniform constraint is to make the domain smaller.  A desk that runs
volatility targeting, or has a contractual risk limit, can often justify a much narrower ratio range
before seeing any performance data.  The plan's section 6.2 asks what that is worth, and the answer
is: **much less than it sounds.**

| design | unrestricted (3 decades) | factor of 4 | factor of 2 | vol-targeted (1.5x) | best recovery |
|---|---|---|---|---|---|
| df3 | 23.7 pp | 23.7 | 23.0 | 22.1 | **1.6 pp** |
| df5 | 11.1 | 11.1 | 10.6 | 10.2 | **0.9** |
| df11 | 5.8 | 4.8 | 3.8 | 3.8 | **2.1** |
| df35 | 4.1 | 4.1 | 4.1 | 4.1 | **0.0** |

The reason is in the deposits already: the least-favourable ratio is **not** far out in the domain.
The uniform constraint binds near `lambda ~ 0.3` at three degrees of freedom, which is inside even a
factor-of-two restriction, so removing the distant ratios removes states the rule was not protecting
against anyway.

That is worth saying plainly because it closes a plausible-sounding route.  Narrowing the declared
domain is not the lever; observations on the smaller group are (Section 6).  It also means the rule
loses nothing by declaring a generous domain, which is why the shipped one is three decades either
side.

One condition travels with every row of that table: a restriction is worth its power only if it is
fixed before the performance data are seen.  A ratio bound inferred from the same returns being
tested is circular and voids the guarantee; the deposit carries that requirement in a `governance`
field rather than in a footnote.

## 7b. Two-sided (M6), `deposits/twosided.json`

The plan asks for the two-sided problem as a separate programme with its own factorisation, its own
enclosures and its own falsifiers.  At K = 1 it needs none of them, because of an identity that is
exact rather than asymptotic.  Conditioning on the variance estimates the numerator is symmetric
normal, so

    P(|T_1| > c | S) = 2 Phi(-c d / sigma) = 2 * P(T_1 > c | S)

at every draw, and the factor survives the expectation: `FWER_2s(lambda) = 2 * FWER_1s(lambda)` for
every lambda.  **The two-sided rule at ceiling abar is the one-sided rule at ceiling abar/2, and
every certificate already deposited transfers.**

The identity is falsified rather than assumed: over 300 random designs, ratios and thresholds the
conditional difference is **0.0** and the unconditional difference is **0.0**.  A third falsifier
exists to show the reduction is a K = 1 property and is not being smuggled upward: at K = 2 the ratio
is 1.9989, not 2, and that falsifier is required to fire.

| design | two-sided `alpha0` | sup FWER | power loss (median / max) |
|---|---|---|---|
| df3 | 0.00484 | 0.02451 | 35.6 / 55.8 pp |
| df5 | 0.00921 | 0.02452 | 17.6 / 31.9 |
| df7 | 0.01098 | 0.02500 | 11.1 / 20.1 |
| df11 | 0.01196 | 0.02500 | 10.1 / 11.1 |
| df19 | 0.01169 | 0.02446 | 6.5 / 10.8 |
| df35 | 0.01159 | 0.02450 | 6.5 / 11.3 |
| balanced df35 | 0.01244 | 0.02500 | **3.1 / 3.1** |

Above K = 1 the reduction fails and the module says so rather than extending it.  The two-sided law
is implemented exactly and the level is calibrated against it, but the **enclosure is not available
by the same route**: each factor `Phi(s(x+cd)) - Phi(s(x-cd))` is non-monotone in `s`, tending to one
inside the band and to zero outside it, so the corner substitution behind Section 4.1 does not apply.
That obstruction is named in the deposit rather than left for a reader to discover.

## 7c. Serial dependence (M6), `deposits/serial.json`

The preceding programmes measured this boundary and stopped: under AR(1) errors no rule in the tested
set survives `phi ~ 0.1` at n = 36.  That is true of rules built for independent data, and it is not
the end of the matter.

**With a declared coefficient it is not a boundary at all.**  The Prais-Winsten transform maps a
stationary Gaussian AR(1) to exactly independent Gaussians, so the whitened problem is the iid
problem this repository already solves and every rule and certificate transfers unchanged.

| `phi` | naive rule | whitened with the correct `phi` |
|---|---|---|
| 0.00 | 0.0254 | 0.0254 |
| 0.05 | 0.0308 | 0.0249 |
| 0.10 | 0.0394 | 0.0252 |
| 0.20 | 0.0546 | 0.0248 |
| 0.30 | 0.0748 | 0.0243 |
| 0.40 | 0.1000 | 0.0249 |

At 40,000 draws the whitened column runs from 0.0243 to 0.0254, every value within 0.9 standard errors of
the ceiling (z from -0.91 to +0.54).  The 4,000-draw run of the first pass put two of those cells 2.3
standard errors above the ceiling, which reads as a breach and was sampling noise; the precision has
to separate the two before the column means anything.

Two riders travel with it.  **Misspecification is not free**: declaring 0 when the truth is 0.3 gives
0.0697, and declaring 0.3 when the truth is 0 gives 0.0019 - valid, and throwing the budget away.  And
**a declared bound is the deployable version**: lowering the internal level until the worst case over
`phi <= phi_max` is under the ceiling gives a rule valid for every coefficient in the bound, at
`alpha0` 0.02494 for `phi <= 0.10` against 0.02494 for the iid rule, which is no reduction at all, because
whitening at the top of the bound is conservative everywhere inside it.  The worst realised rates
under the three bounds tested are 0.0247, 0.0244 and 0.0235.


### 7c (continued). The condition the whitening carries, and the transform that removes it

The claim above, declared coefficient, so the iid certificates transfer, was tested at one cell:
`n0 = n1 = 36` and a mean of zero.  Both of those matter, and neither was visible from inside that
cell.

Prais-Winsten rescales the first observation rather than dropping it, so the transformed mean is

    mu * k(phi),   k = (sqrt(1 - phi^2), 1 - phi, ..., 1 - phi),

which is **not constant in t**.  The deployed statistic takes the plain mean of the whitened series,
so it estimates `mu * mean(k)`, and `mean(k)` depends on the group's LENGTH.  With `n0 = n1` the two
groups share it and the difference of means is unbiased for zero under the null whatever `mu` is.
With `n0 != n1` the numerator carries a bias proportional to `mu`.

**Measured** (`deposits/serial_falsify.json`, 40,000 draws, each design calibrated at its own `(n0, n1)`
so a deviation cannot be an artefact of running one design's table on another's data):

| design | phi | mu | rescaled (Prais-Winsten) | dropped (Cochrane-Orcutt) |
|---|---|---|---|---|
| 36 vs 13 | 0.2 | 0.0 | 0.02480 | 0.02475 |
| 36 vs 13 | 0.2 | 0.5 | 0.02562 | 0.02475 |
| 36 vs 13 | 0.2 | 2.0 | 0.02772 | 0.02475 |
| 36 vs 13 | 0.4 | 0.0 | 0.02552 | 0.02562 |
| 36 vs 13 | 0.4 | 0.5 | 0.02678 | 0.02562 |
| 36 vs 13 | 0.4 | 2.0 | 0.03030 | 0.02562 |

Rescaling drifts with the mean and reaches **0.03030**, 6.2 standard errors above the ceiling.
Dropping the first observation is flat: identical to five decimals across every `mu`, because the
transformed mean vector is then constant at `1 - phi` and the statistic cannot see `mu` at all.

**The repair is a change of transform, not a change of estimator.**  Cochrane-Orcutt gives
`y_t = x_t - phi x_{t-1}` for `t >= 2` and discards `y_1`.  The mean vector is constant, so the
transformed sample is a location model on `n - 1` observations and lands in the Sections 3/4 family
with `n` replaced by `n - 1`.  No regression, no effective sample size, one observation per group.
Worst deviation at the unbalanced design **0.00062 (0.8 se)**; `repair_works` is deposited as a
boolean, not asserted here.

**What this changes upstream.**  The serial contribution was written as unconditional ("every
certificate transfers unchanged").  It is unconditional again, but only with the transform stated:
the paper now says which whitening to use and why, and the falsifier is what distinguishes them.
Table 25 remains correct and remains balanced and mean-zero, which is why it could not have found
this.

## 7d. Factor exposure (M6), `deposits/factor.json`

Exactly the same shape as serial dependence, which is the point.  With **declared** loadings the
factor-residual series is iid Gaussian, the estimand becomes factor alpha, and the certified rule
applies unchanged.  40,000 draws at n = 36, nominal 0.025:

| exposure | raw returns | declared loadings | estimated in sample, df corrected | estimated, df NOT corrected |
|---|---|---|---|---|
| none | 0.02388 | 0.02388 | 0.02513 | 0.02697 |
| same sign (1.0, 0.8) | **0.00477** (z = -59) | 0.02492 | 0.02645 | 0.02825 |
| opposite signs (1.0, -0.8) | **0.05065** (z = +23) | 0.02390 | 0.02595 | 0.02810 |
| benchmark unexposed | 0.02545 | 0.02455 | 0.02612 | 0.02788 |

The raw rule is catastrophically conservative when the loadings agree and twice nominal when they
oppose - which reproduces the preceding programme's finding and explains it: the factor inflates or
deflates the numerator's variance without touching the denominator's estimate.  Residualising with
declared loadings removes both, in every pattern.

Estimating the loadings in sample is a **plug-in, not a theorem**, and the deposit measures the size
of the one mistake a practitioner will actually make: the regression costs a degree of freedom, and
an arm that does not give it back to the rule runs at z = +2.4 to +3.9 rather than +0.2 to +1.8.
Loadings estimated on a pre-period that is not reused restore the exact statement: 0.02420, z = -1.04.

**The unifying result.**  Serial dependence and factor exposure are not two boundaries.  They are the
same one: a declared nuisance is a pre-processing step and the certificate survives it; an estimated
nuisance is a new parameter and the certificate does not.  Nothing in either case required a new
theorem, and neither is claimed to have one.

## 7e. Independent cross-check (M7), `deposits/crosscheck.json`

Every decisive quantity re-derived by a route sharing no code path with the one that produced it.

| check | result |
|---|---|
| X-A deposited thresholds vs `HeteroMaxT`, the paper's own implementation | worst 1.1e-03 over 42 bins in 6 designs |
| X-B FWER: Monte Carlo with the exact inner, vs the tensor quadrature | worst \|z\| = 1.36 over 18 comparisons |
| X-C power, the same way | worst \|z\| = 2.50 over 12 |
| X-D the certified bound must dominate a direct sample | 42 ratio points, **0 violations beyond 3 se** |

**This is a cross-check, not a blind replication**, and the deposit says so in a `scope` field: it was
written by the same author in the same session as the code it checks.  What it establishes is that
the quadrature, the tabulation and the enclosure are not all wrong in the same direction, which is
the failure mode a single implementation cannot detect.  The blind version is `replicate/`, which is
written and sealed - thirteen deposit hashes and ten canonical designs - and **has not been executed**.

The kit is built around the defect that broke the preceding programme's replication: that one failed
its draw-dependent half because the specification never pinned the draw recipe, the definition of
`d`, or the semantics of a certificate point.  **This construction has no draws in the decision path
at all**, so the entire class of ambiguity is absent by design; Monte Carlo appears only in falsifier
R3, where the seed does not need to match.

## 7f. Certified at ANY K, one- and two-sided, `deposits/certify_anyk.json`

Two things were reported as open in an earlier draft of this report: the two-sided enclosure above
K = 1, and the certification of any rule above K = 1.  Both are closed, and neither needed a new
idea - they needed one line of algebra rather than a bigger enclosure.

For a per-candidate threshold the rejection events are one per candidate, so pointwise in the
benchmark pair `1 - prod_i g_i <= sum_i (1 - g_i)`, and taking expectations

    FWER(lambda_1 .. lambda_K)  <=  sum_i FWER_1(lambda_i) ,

where `FWER_1` is the ONE-CANDIDATE error rate.  The right side separates, so the supremum over the
whole ratio cube is at most `K` times the one-candidate supremum - and that is exactly the quantity
Section 4.1's enclosure already certifies tightly over the continuum.  Certify it at `abar/K` and the
K-candidate rule is certified at any K; at `abar/(2K)` for two-sided, using the exact factor of two
of Section 7b.  **No K-dimensional enclosure is required.**

| design | sided | K | `alpha0` | one-candidate bound | union total | slack vs the exact law |
|---|---|---|---|---|---|---|
| K2_df35 | one | 2 | 0.01155 | 0.012500 | 0.025000 | 1.35e-03 |
| K2_df35 | two | 2 | 0.00556 | 0.006250 | 0.024999 | 1.39e-03 |
| K3_df35 | one | 3 | 0.00753 | 0.008333 | 0.025000 | 1.63e-03 |
| K3_df35 | two | 3 | 0.00359 | 0.004167 | 0.025000 | 1.73e-03 |
| K5_df35 | one | 5 | 0.00438 | 0.005000 | 0.025000 | 1.88e-03 |
| K5_df35 | two | 5 | 0.00204 | 0.002500 | 0.024999 | 2.25e-03 |
| K2_df11 | one | 2 | 0.01177 | 0.012500 | 0.025000 | 8.99e-04 |
| K2_df11 | two | 2 | 0.00571 | 0.006250 | 0.025000 | 9.15e-04 |
| K3_df11 | one | 3 | 0.00772 | 0.008333 | 0.024998 | 1.03e-03 |
| K3_df11 | two | 3 | 0.00371 | 0.004167 | 0.025000 | 1.10e-03 |

**Why this route and not the direct one.**  A three-dimensional enclosure over `(z, u_0, u_i)` was
written first, and it works: its substitutions are falsified at 0 violations in 1,200 points, one-
and two-sided, and it dominates the exact law.  It is also **loose by a factor of 2.7** - at K = 2 it
returns 0.0655 where the true supremum is 0.0246 - and tightening it costs cells in three axes at
once.  The union step costs a fiftieth of that: slack `5.4e-05` to `2.25e-03` at the ratios that
actually bind.  The losing construction is kept in `certify_percand.py` and reported here rather
than deleted, because which route was chosen and why is part of the result.

**What it costs.**  At K = 2 with a 36-month benchmark the certified level is 0.01155 against 0.01224
for the same rule merely measured, so the certificate itself costs about 6% of the level.
Against the *joint* plug-in's 0.02445 it is less than half - but that gap is the per-candidate
restriction of Section 7, not the enclosure: the joint rule cannot be certified above K = 1 by any
route in this repository, and the per-candidate rule can be certified at every K.

## 7g. The JOINT rule, certified, `deposits/certify_joint.json`

The joint plug-in was the last thing here reported as open.  It is the more
powerful of the two rules, so leaving it measured had a real cost, and the
obstruction was genuine: its threshold reads the whole estimated ratio vector,
so neither the union step of Section 7f nor the per-candidate factorisation
applies to it.

**What did not work.**  The conditional no-rejection probability is increasing in
the threshold, so replacing the threshold on a product cell by its minimum over
the ratios that cell can reach is a valid lower bound - and with the threshold
constant on the cell, the product factorises again.  That construction is
implemented and it dominates.  It is also loose by about a factor of two and
converges like `m^-1` in four axes at once: 0.0612 at 32 cells per axis, 0.0558
at 80, against a true value near 0.025.

**What did.**  Deploy the joint rule as a step table on a grid of estimated
ratios - which is what a desk runs anyway - and project it:
`psi_i(t) = min over the other coordinates of the deployed table`.  Because the
deployed rule IS the table, that minimum is exact rather than an approximation
of something else.  Rejection of candidate i then requires
`T_i > c(lambda_hat) >= psi_i(lambda_hat_i)`, so the union step of Section 7f
applies to `psi` and each term is certified by the tight two-dimensional
enclosure.  One projection handles the entire coupling.

| design | K | `alpha0` | union total | vs the certified per-candidate rule | projection gap |
|---|---|---|---|---|---|
| K2_df35 | 2 | 0.01332 | 0.024958 | 1.15x | 5.5% |
| K2_df11 | 2 | 0.01296 | 0.024990 | 1.10x | 3.7% |
| K3_df35 | 3 | 0.00928 | 0.024964 | 1.23x | 11.2% |
| K3_df11 | 3 | 0.00883 | 0.024983 | 1.14x | 6.1% |

So the joint rule certifies.  The `alpha0` column is NOT a comparison with the
per-candidate rule, and the earlier version of this section read it as one.  The
two levels are levels on different statistics - the joint level on the maximum
over candidates, the per-candidate level on one comparison - and at any common
level the constant on the maximum is the larger.  A higher joint `alpha0` can
therefore still mean a higher threshold, and here it does.

**The comparison in the quantity a desk applies** (`deposits/rule_compare.json`),
at equal estimated ratios:

| design | joint, certified | joint, measured | per-candidate, certified |
|---|---|---|---|
| K2_df35 | 2.6560 | 2.3828 | 2.4384 |
| K2_df11 | 2.6647 | 2.3856 | 2.4266 |
| K3_df35 | 2.9494 | 2.5558 | 2.6314 |
| K3_df11 | 2.9914 | 2.5308 | 2.6173 |

The certified per-candidate threshold is BELOW the certified joint threshold at
every design and every ratio pair tested.  The reason is the cost of the two
certificates: certifying the joint rule costs 44.9 to 62.7% of its own
measured level, against 3.7 to 6.8% for the per-candidate rule.  The
projection handles the coupling; the coupling is expensive to prove and cheap to
give up.  The joint rule remains the better rule when its level is merely
measured (2.3828 against 2.4384 at K2_df35), so this is a result about proof and not
about the rules.

The projection gap - how far the joint threshold falls when the other
candidates' ratios are pushed to their worst - grows with K as it must: 3.7 to
5.5% at two candidates, 6.1 to 11.2 at three.

**How conservative the certified joint rule actually is.**  The bound is what is
proven; the rule's realised error rate is what a desk pays for.  Measuring the
deployed step rule's own unconditional FWER, Rao-Blackwellised:

| design | union bound (proven) | realised FWER | budget spent |
|---|---|---|---|
| K2_df35 | 0.024958 | 0.013688 | 55% |
| K2_df11 | 0.024990 | 0.013409 | 54% |
| K3_df35 | 0.024964 | 0.009586 | 38% |
| K3_df11 | 0.024983 | 0.008948 | 36% |

So it spends 36 to 55% of its error budget, against the region rule's
0.7% (Section 3).  The slack is the union step plus the projection gap,
and it is what certifying a rule whose threshold couples the
candidates.

**A gate caught an error in this section before a reader could.**  The first
version of the deposit filled its comparison field by evaluating a PER-CANDIDATE
plug-in table at the joint rule's level - a different rule, which at the same
level has a lower threshold and therefore a higher error rate.  Its FWER sat
above the bound, which looks exactly like a bound failing to dominate, and gate
X11 failed on all four designs.  The certificates were never affected, because
they come from the projection and the enclosure and never used that field.  It
is recorded here rather than quietly corrected, because the interesting part is
that comparing a bound against the wrong object is a failure mode a gate can
catch and a reader cannot.

## 7h. Speed, `deposits/fast.json`

Two changes, both made after profiling rather than before it.  The profile put
**96%** of the time in the inner product-over-candidates kernel and only
4% in interpreter overhead, which ruled out the usual remedies.

1. **The distribution function does not depend on the level; only the root-find
   does.**  Every bisection over a level in this folder had been re-solving the
   same law from scratch at each step.  Tabulating the law once on a grid of `c`
   and inverting serves every level.  The inversion interpolates
   `log(rejection)` against `c`, which is near-linear because the tail is
   near-exponential.  Both scales are MEASURED on the same grid rather than
   asserted: **6.3e-05** in the constant on the log scale against **9.1e-04** on
   the probability scale, a factor of **15**, at identical cost.  Gate X10
   allows 1e-3, so the probability scale would have passed it with almost
   nothing to spare.
2. **Twelve cores were present and one was in use.**  The sweeps are
   embarrassingly parallel over table cells; the reduction is by index, so the
   deposit does not depend on scheduling.

Measured on a 144-cell, 8-level sweep: **255 s to 9.9 s, 26x**.  On the joint
certification the effect was 1,476 s to 261 s per design, **5.7x**, with
bit-identical output - `alpha0` 0.01332 and 0.01296 before and after.

Gate X10 exists because an optimisation that changes an answer is not an
optimisation: the parallel map must agree with the serial one **exactly**, and
the tabulated inversion with a direct root-find to the tolerance its node count
supports.

## 8. Verification, `deposits/gates.json`

Eight gates, each able to fail for the right reason:

| gate | what it falsifies |
|---|---|
| X1 | the K = 1 closed form against the general z-quadrature engine |
| X2 | the monotone corner substitutions behind the rigorous enclosure, on random cells and ratio boxes |
| X3 | the enclosure bracket must contain the quadrature value at every ratio tested |
| X4 | the per-candidate engine with a constant threshold must reproduce the paper's exact fixed-ratio law |
| X5 | the per-candidate engine at K = 1 must equal the K = 1 quadrature |
| X6 | every deposited rule re-checked on a grid three times finer than the one it was calibrated on |
| X7 | tier discipline: no deposit may claim a certified tier without a certified supremum |
| X8 | the certified rule may never be more liberal than its own certificate |

**One rigour defect was found and closed by this suite rather than by inspection.**  The enclosure
partitions the state space over `[1e-9, 1 - 1e-9]` on each axis, which silently dropped 4e-9 of
probability mass from the upper bound.  The integrand is a probability, so that mass can contribute
its own size to the truth; the tightest certificate margin in the deposit was 1.6e-8, so the omission
was within an order of magnitude of mattering.  `wb_certify.enclose_fwer_box` now adds it to every
upper bound, and `recertify.py` re-ran all ten certificates against the corrected version: every one
still certifies, the bounds moved by exactly 4e-9, and no design needed re-levelling.

## 8a. The audit programme, and what it cost

Three adversarial passes ran against the built document, plus one independent pass.
They are recorded because their yield is the argument for stopping, not a formality.

| pass | scope | confirmed | rejected on the evidence |
|---|---|---|---|
| 1 | numbers, figures, structure | 21 | 0 |
| external | numbers | 1 | 1 (a Rao-Blackwell standard error read as a binomial one) |
| 2 | fixes, plus what pass 1 missed | 12 | 3 |
| 3 | mathematics and prose logic only | 6 | 0 |

Pass 3 found the only defect in the programme that changed a result: the whitening condition in 7c
above.  Passes 1 and 2 found, between them, three fixes that pass 1 had reported as applied and
which a half-applied patch script had silently dropped, the reason patch scripts in this folder now
report per edit instead of asserting out.

**Two of the four passes were partly wrong**, and the rejections matter as much as the confirmations:

* the deposited `exact_sup_mcse` (0.000177 at a rate of 0.013688 over 4,000 draws) was called
  impossible because a binomial standard error would be 0.00184.  There are no Bernoulli draws:
  `fix_joint_compare.py` draws the states and evaluates the conditional rejection probability in
  closed form.  `Var(r) + E[r(1-r)] = 1.350064e-02` against `p(1-p) = 1.350064e-02`, to seven
  figures, the variance decomposition closes, and a 10x reduction is what conditioning buys;
* the `A_df3` cell was said to be reported as failing with a gap below tolerance.  It certifies
  (3.81e-05 against 1e-4); all five failures do exceed the tolerance;
* Table 16's region row was said to be the 9.8 certified tier.  `wa_ledger.py` builds an exact
  `1 - gamma` interval and takes the supremum of the conditional oracle constant over it, which is
  the 9.3 envelope; the 2.4785 ratio reproduces from first principles at `nu1 = 12`, where the
  interval's upper end is 7.9x the estimate.

## 8b. What the checks actually cover, `coverage_frontier.py`

"%d checks pass" does not say how many of the paper's claims a check ever looks at.  A claim no
check formats passes unconditionally, and counting passes cannot see it.  `coverage_frontier.py`
instruments the verifier to record every string it looks for, extracts every numeric token from the
sections this programme added, and reports the ones nothing covers.

First run: **86%%, 50 uncovered**, every one of them a quantity the prose DERIVES from a deposit
rather than quoting.  After adding checks for them: **100%% of %d tokens**, with six uncovered for a
stated reason (three of them figures the preceding programme owns, one a rounded characterisation,
one a definition, one Section 9's own price in the constant).

The checks were also shown to bind, which is a separate question:

* against `BACKUP-pre-frontier.docx` the verifier reports 347 of 396 failed, so the checks are not
  passing on presence alone;
* perturbing one deposit field by 1e-4 fails exactly the matching check;
* seven corruption probes, four falsified derived numbers, two intra-table transpositions, one
  draw count, each fail on the right check.  Before the row-scoped and derived-quantity checks were
  added, all seven passed.

That last point is the one worth keeping: substring containment cannot see a transposition, because
both values are still somewhere in the document.  Table rows are now matched whole.

## 9. Limitations

1. **The certificate above K = 1 is for the PER-CANDIDATE rule only.**  Section 7f certifies that
   rule over the continuum at every K and on either side.  The joint plug-in, which is the more
   powerful of the two, remains measured above K = 1: its threshold couples all the coordinates, so
   neither the union step nor the enclosure applies to it, and no route in this repository
   certifies it.  A desk that needs a certificate above one candidate runs the per-candidate rule
   and pays Section 7's power gap for it.
2. **The declared ratio domain is three decades either side.**  Validity outside `[1e-3, 1e3]` is
   not claimed.  Every desk and stress design in this repository is inside it, and the domain is
   stated rather than assumed away.
3. **Three degrees of freedom stays expensive.**  Even the in-class optimum loses 20.9 pp there.
   The design map says what to do about it: get more observations on the smaller group, or accept
   the loss knowingly.
4. **The K = 1 and K > 1 loss tables use different oracles** and are not comparable; the level price
   is the column that is.
5. **Two-sided is certified wherever one-sided is.**  The exact factor of two transfers the
   one-sided certificates at half the level at K = 1, and the union step of Section 7f carries it
   to every K at `abar/(2K)`.  The direct two-sided enclosure remains loose and is not the route
   used.
6. **Serial dependence and factor exposure are handled only for a DECLARED nuisance.**  Whitening
   with a declared coefficient, or residualising with declared loadings, makes the problem exactly
   iid and the certificate transfers.  An estimated autocorrelation or an estimated loading is a new
   nuisance parameter, is labelled a plug-in, and is not certified.  The `phi ~ 0.1` boundary of the
   preceding programmes stands for rules that ignore the dependence - which is what it always was.
7. **The blind replication has not been run.**  `replicate/` is written and sealed, thirteen
   deposit hashes and ten canonical designs, and no independent party has executed it.  Until one
   has, "independently reproduced" is not available; what is available is the cross-check of
   Section 7e.
8. **The optimiser in Workstream C is a local solver** on a smooth problem with analytic gradients.
   Its output is verified feasible on a grid three times finer than the one it optimised against,
   and shrunk until it is; it is not proven to be the global optimum of its own program.

## 10. Reproduce

    cd 04-efficient-rule
    python wb_k1.py             # the uniform calibration at K = 1
    python wc_frontier.py       # the in-class optimum, and the gain over the plug-in
    python wb_certified_rule.py # the certified step rules and their certificates
    python wb_general.py        # K = 2, 3, 5, measured
    python wb_percand.py        # the certifiable per-candidate class
    python wd_designmap.py      # the planning surface
    python wa_ledger.py         # the four-object ledger
    python twosided.py          # the exact factor of two, and its falsifiers
    python serial.py            # whitening with a declared coefficient
    python factor.py            # residualising with declared loadings
    python certify_anyk.py      # the per-candidate certificate at any K
    python certify_joint.py     # the joint certificate at K = 2, 3
    python fix_joint_compare.py # the deployed joint rule's own FWER
    python crosscheck.py        # the independent numerical cross-check
    python fast.py              # the speed layer and its falsifier
    python tests_frontier.py    # the gates
    python verify_frontier.py   # every number in this report against a deposit
    python deskrule.py          # three worked decisions, end to end
