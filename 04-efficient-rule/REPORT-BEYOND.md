# Beyond the exact Welch-Dunnett law: W1 to W4

What the roadmap asked for, what it got, and what it did not.  Every number here is formatted from a deposit under `deposits/` at build time and checked by `verify_beyond.py`.  Each carries a label from the taxonomy in `PROMOTION-PROTOCOL.md`: **exact in model**, **certified bound**, **measured**, **asymptotic**, **diagnostic**.

**One of those labels is missing from everything below, and that is the headline.**  No result in this programme is a *certified bound*.  The preceding paper's branch-and-bound construction remains the only certificate in the repository.  W1 established what that certificate's class actually is and bounded the discretisation inside it; W2 produced an exact law for a strictly larger model but only a grid maximum over its nuisance; W3 produced measurements, two of them negative; W4 produced a protocol.  Anything here described as holding over a domain holds over a GRID, and that difference is the thing this repository exists to take seriously.

---

## W1. What information the rule actually has

### The class Section 10 searched is the invariant class

The null has three free parameters and the sufficient statistic has 4 coordinates.  The only transformations that leave the problem alone are `x -> a + b x` applied to **both** groups, which is the 2-parameter affine group; scaling the groups separately changes the ratio and so changes the problem.  A maximal invariant is therefore 2-dimensional, and `(T, lambda_hat)` is one.

Both halves are checked rather than asserted.  Applying random affine maps to raw data moves the invariant by at most **4.25e-12** across three designs.  For maximality, a second sufficient statistic is constructed directly on the same `(T, lambda_hat)` level set, not by applying a group element, and the group element that maps one to the other is recovered in closed form as `b = sqrt(w0^2/s0^2), a = Ybar_0 - b Xbar_0`; the worst reconstruction error over 200,000 pairs per design is **1.46e-11**.

So Section 10 did not search a convenient parametric family.  It searched the non-randomised invariant tests, up to monotonicity and bin resolution.  *(exact in model)*

### The route to a lower bound

The group is solvable, because translations form a normal abelian subgroup with abelian quotient.  Solvable groups are amenable, and for an amenable group the Hunt-Stein theorem gives that a minimax test may be taken invariant, so a bound over (T, lambda_hat)-measurable tests bounds all tests for a worst-case criterion.

That is the sentence the preceding programme could not write.  It is a **derivation, not a computation**, and three conditions are still to discharge:

* the loss (regret against the known-ratio oracle) must be invariant under the group, which it is because both power and oracle power are invariant functions of lambda;
* the usual measurability and topological conditions of the Hunt-Stein theorem for the affine group;
* the extension from K = 1 to K > 1, where the group is the same but the invariant is (T_1..T_K, lambda_hat_1..K).

### What 36 bins can hide

The discretisation is itself an approximation, so the frontier was re-solved at 36, 72 and 144 bins with the same shrink step applied to each so that valid rules are compared with valid rules.

| design | 36 bins | 72 bins | 144 bins | hidden |
|---|---|---|---|---|
| df35 | 4.36 pp | 4.35 pp | 4.36 pp | 0.004 pp |
| df3 | 20.95 pp | 20.93 pp | 20.93 pp | 0.017 pp |

At most **0.017 points** anywhere.  The 36-bin result is not an artefact of the grid.  *(exact in model)*

### The optimal rule is monotone and deterministic

Section 10 restricted the rule to an upper set in `T` and to a deterministic decision.  Neither is implied by validity.  Dropping both turns the problem into a linear program over every region in the discretised plane, because the mass each cell receives is a closed-form function of the state quadrature and both the objective and the uniform constraint are linear in the rejection probabilities.

Monotonicity is imposed as a linear constraint rather than by rendering a threshold rule onto the grid, because rendering drops the cell straddling the threshold and hands the unrestricted rule about a point of power for free, the same order as the effect being measured.  The first version of this module did exactly that and reported a spurious 1.47-point gain.

| design | monotone optimum | unrestricted optimum | deterministic floor | gain |
|---|---|---|---|---|
| df3 | 20.69 pp | 20.69 pp | 21.15 pp | +0.00 pp |
| df5 | 11.23 pp | 11.23 pp | 11.92 pp | -0.00 pp |
| df7 | 7.89 pp | 7.89 pp | 8.42 pp | -0.00 pp |
| df11 | 5.52 pp | 5.52 pp | 6.14 pp | -0.00 pp |
| df19 | 4.40 pp | 4.40 pp | 5.06 pp | -0.00 pp |
| df35 | 4.17 pp | 4.17 pp | 5.03 pp | +0.00 pp |

The gain is **-0.00 to +0.00 points**: zero to the resolution of the computation, at every design.  Allowing a non-monotone region and allowing randomisation buy nothing.  That is numerical evidence for the proposition the roadmap listed as a target, that the optimal rule in this class is monotone in `T`, and it means Section 10's restriction to threshold rules was not a simplification.  Rounding the fractional cells down costs 0.46 to 0.85 points, which is the price of an exactly deterministic boundary rather than of determinism itself.

The constraint is imposed at 97 ratio values and verified at 193, so the verification grid finds a small overshoot: at most 3.69e-04, 1.48% of the ceiling.  It is reported rather than shrunk away because both arms carry it equally and the comparison is between them.  T-grid sensitivity at most 0.034 points.  *(measured)*

### The low-degrees-of-freedom regime

Both group sizes move, because the cost is set by the smaller of the two and a one-axis design map cannot see that.  A gain counts as material at 0.50 points, declared before the run.

| min(nu) | mean gain from a design-specific threshold |
|---|---|
| 3 | 3.34 pp |
| 4 | 1.81 pp |
| 5 | 1.09 pp |
| 7 | 0.38 pp |
| 9 | 0.09 pp |
| 11 | -0.05 pp |

Material in **26 of 36 cells**, best **4.09 points**.  The shipped table is deployable unaltered when the SMALLER group has more than 5 degrees of freedom; at or below that a design-specific threshold function is worth 0.50 points or more.  *(measured)*

---

## W2. Aligned returns

### The roadmap's route was wrong, and the replacement keeps the architecture

As written, Phase 2A formed paired differences.  That gives a maximum of K studentised statistics with different studentisers and a Wishart-coupled denominator set, which nothing factorises: the K-dimensional integral returns.

Keeping the returns unpaired and carrying the dependence in a latent factor instead, `r_it = mu_i + beta_i f_t + eps_it, f latent standard normal, eps independent across i, sigma_f = sigma_0 = 1`, conditioning on **`fbar`, `Vf`, `Xbar_0`, `s_0^2`** makes the K candidates exactly conditionally independent.  The reason is that `eps_bar_i` is orthogonal to the centred `eps_i`, and that `(n-1)s_i^2/sigma_i^2` is noncentral chi-square whose noncentrality sees the factor path only through its sum of squares, by rotational symmetry.

The null law is then a **4-dimensional outer rule over one one-dimensional inner rule per candidate**, cost linear in K.  Two of the four are the two the existing paper already integrates.

| cell | brute force (4,000,000 draws) | exact evaluator | z | what the independent law says |
|---|---|---|---|---|
| no loadings, independent | 0.02884 | 0.02872 | -1.42 | 0.02872 |
| equal loadings | 0.00282 | 0.00281 | -0.16 | 0.02872 |
| unequal loadings | 0.00895 | 0.00897 | +0.38 | 0.02872 |
| unequal loadings and ratios | 0.01118 | 0.01123 | +0.82 | 0.03055 |
| thirteen observations | 0.01025 | 0.01028 | +0.53 | 0.02602 |
| strong loadings | 0.00935 | 0.00933 | -0.38 | 0.02872 |
| one negative loading | 0.03153 | 0.03156 | +0.42 | 0.02928 |
| five candidates | 0.00707 | 0.00713 | +1.58 | 0.02819 |

Worst |z| **1.58** over 8 cells.  The last column is the non-vacuity check: the independent-groups law is wrong here by up to **90.2%** in relative terms, and the evaluator tracks the truth across that range.  A reduction that quietly ignored the factor would pass every other cell and fail that column.  Refinement gap against a coarser grid 3.5e-04.  *(exact in model)*

Critical constants, for the size of the correction:

| design | c at 0.025 |
|---|---|
| equal loadings, 36 observations | 1.6053 |
| unequal loadings, 36 observations | 1.9251 |
| equal loadings, 13 observations | 1.6882 |
| no loadings, 36 observations | 2.2596 |

### What it costs to ignore the factor, and what a certificate over it costs

The mechanism decides where the danger is, and the first version of this module could not see it.  In the aligned model the numerator carries (beta_1 - beta_0)^2 and the denominator carries beta_1^2 + beta_0^2; equal only when one loading is zero, conservative when the loadings share a sign, liberal when they do not.

That first version set the benchmark loading to zero, on the argument that only the contrast is identified.  It is not, and zero is exactly the boundary between the two regimes, so the scan found nothing.  A long-short candidate against a directional benchmark is the opposite-sign case, which is a configuration a desk is likely to be in.

| ratio | beta_0 | beta_1 | realised FWER | z |
|---|---|---|---|---|
| 0.25 | 0.0 | 1.0 | 0.02496 | -0.1 |
| 0.25 | 1.0 | 1.0 | 0.00117 | -382.4 |
| 0.25 | 1.0 | 0.2 | 0.01601 | -39.2 |
| 0.25 | 1.0 | -0.5 | 0.04841 | +59.7 |
| 0.25 | 1.0 | -1.0 | 0.06103 | +82.4 |
| 0.25 | 1.5 | -1.5 | 0.07155 | +98.9 |
| 1.00 | 0.0 | 1.0 | 0.02475 | -0.9 |
| 1.00 | 1.0 | 1.0 | 0.00312 | -215.0 |
| 1.00 | 1.0 | 0.2 | 0.01790 | -29.3 |
| 1.00 | 1.0 | -0.5 | 0.04281 | +48.2 |
| 1.00 | 1.0 | -1.0 | 0.05420 | +70.6 |
| 1.00 | 1.5 | -1.5 | 0.06543 | +89.6 |
| 4.00 | 0.0 | 1.0 | 0.02540 | +1.4 |
| 4.00 | 1.0 | 1.0 | 0.01066 | -76.5 |
| 4.00 | 1.0 | 0.2 | 0.02092 | -15.6 |
| 4.00 | 1.0 | -0.5 | 0.03441 | +28.3 |
| 4.00 | 1.0 | -1.0 | 0.04161 | +45.6 |
| 4.00 | 1.5 | -1.5 | 0.05288 | +68.2 |

**9 of 18 cells are liberal**, and the split is exactly along the sign: 9 of 9 opposite-sign cells, 0 of 9 same-sign cells.  Worst familywise error rate **0.07155**, which is **2.86 times** the ceiling; the worst same-sign cell only reaches 0.02540.  So the aligned factor is not a conservative nuisance, and which side of the ceiling a desk lands on is decided by a sign it never estimates.  *(measured)*

The nuisance is therefore **three-dimensional per candidate**, (ratio, beta_0, beta_1), and only the ratio is estimable from what the rule sees.  Over a declared box of 175 points the constant that holds everywhere is **2.7366**, against known-parameter oracle constants of 0.9003 to 2.7367, so the certificate costs **1.00 times** the worst oracle constant and 3.04 times the best.  The supremum is a **grid maximum**, not a continuum certificate: the branch-and-bound construction applies in principle but the enclosure would be three-dimensional per candidate rather than one and has not been built.  *(measured)*

### Unequal histories: what the incumbent's prehistory buys

A promotion gate has an asymmetry the clinical model does not: the incumbent has a long record and the candidate a short one.  Three predeclared designs on an overlap of 24 observations.

| prehistory | full-history power gain over strict pairing |
|---|---|
| 12 | +7.84 pp |
| 36 | +15.06 pp |
| 120 | +22.12 pp |
| 360 | +25.24 pp |

The two-window design, which uses the prehistory for the benchmark's VARIANCE only and so needs constant variance but not a constant mean, buys **0.07 points**: essentially nothing.  That is the expected answer and the reason is structural, since the prehistory enters only through `s_0^2`, whose sampling error is second order once the overlap is long enough to matter.

Full history, which uses the prehistory for the mean as well, buys **25.24 points**, and it does not saturate by 360.

That is a real gain and it is bought entirely with an assumption about the MEAN.  Breaking it takes the worst familywise error rate to **0.57120**, **22.85 times** the ceiling, on the `full_history` design.  The direction matters: a prehistory running ABOVE the later mean makes the test powerless, one running BELOW it makes the test liberal, and only the second is dangerous.  So there is no free lunch in the incumbent's track record.  Either take strict pairing and give up those 25.24 points, or buy them with a stationarity assumption on the mean whose failure costs a factor of 23.  *(measured)*

---

## W3. Dependence as three validity tiers, and what breaks the rule

### The tiers

| nuisance | A, declared | B, calibration period (360) | C, same period |
|---|---|---|---|
| serial (AR(1), Cochrane-Orcutt) | 0.02533 | 0.04143 | 0.28297 |
| factor (observed loadings) | 0.02527 | 0.02521 | 0.03426 |

Worst realised familywise error rate in each tier, against a ceiling of 0.025.  Tier B is reported at three calibration lengths in the deposit so the length-versus-error trade is visible rather than asserted.  *(measured)*

| calibration length | serial worst | factor worst |
|---|---|---|
| 36 | 0.15830 | 0.02518 |
| 120 | 0.07358 | 0.02537 |
| 360 | 0.04143 | 0.02521 |

### Why tiers B and C fail for serial dependence, and the repair

Tier A is exact because Cochrane-Orcutt maps a mean `mu` to `mu(1 - phi)` and, with `phi` declared, both groups are multiplied by the SAME constant, so the difference of transformed means is still zero under the null whatever `mu` is.

With the coefficient estimated SEPARATELY IN EACH GROUP the two multipliers differ and the numerator acquires `mu (phi_hat_0 - phi_hat_1)`: zero in expectation, not zero in any sample, and spreading linearly in the common mean.  It is the Prais-Winsten defect's shape arriving by another route.  There the multiplier varied with the group's **length**; here with the group's **estimate**.

**The repair is to pool**: one coefficient from both groups, applied to both, so the multiplier is identical by construction.

| | separate estimates | pooled estimate |
|---|---|---|
| worst realised rate | 0.48409 | 0.07918 |
| times the ceiling | 19.4 | 3.2 |
| worst at a zero mean | 0.08955 | 0.07818 |
| worst at a mean of 5 | 0.48409 | 0.07868 |
| cells breaching | 36 of 36 | 36 of 36 |

The clearest single view is one design held fixed while the common mean moves.  At the worst cell of the grid:

| common mean | separate estimates | pooled estimate |
|---|---|---|
| 0 | 0.08955 | 0.07818 |
| 2 | 0.28162 | 0.07918 |
| 5 | 0.48409 | 0.07868 |

The separate-estimate column grows with the common mean exactly as the algebra predicts; the pooled column does not move, because the repaired statistic cannot see the common mean at all.

The signature is that the pooled rate is flat in the common mean: spread **0.00102** at tier B and **0.00100** at tier C, against a separate-estimate arm that moves from 0.08955 to 0.48409 as the mean goes from 0 to 5.  What pooling does not remove is imperfect whitening at an estimated coefficient, which biases the variance estimate and which the zero-mean rows isolate.  That residual is on its own enough to breach: every one of the 36 cells in both arms is above the ceiling.  **Pooling is an improvement, not a repair**, and the operational conclusion is that Section 11's serial result is a tier-A result: with the coefficient estimated rather than declared, at these sample sizes, there is no guarantee and a calibration period of 360 does not deliver one.  *(measured)*

### What the asymptotic and resampling alternatives deliver

The first row is a REFERENCE, not a competitor.  It is given the coefficient; the others are not, and its job is to confirm that the instrument reads correctly before the rest of the column is believed.  The contest is between the rows that all have to estimate the coefficient from the same data.  Target size 0.025, one-sided, 12 cells per method, common mean 2.0.

| method | worst size | times target | cells breaching |
|---|---|---|---|
| the exact rule, coefficient declared | 0.02545 | 1.02 | 0 of 12 |
| the exact rule, coefficient per group | 0.28242 | 11.30 | 12 of 12 |
| the exact rule, coefficient pooled | 0.07830 | 3.13 | 12 of 12 |
| Newey-West, normal critical value | 0.12676 | 5.07 | 12 of 12 |
| Newey-West, t critical value | 0.10459 | 4.18 | 12 of 12 |
| moving-block bootstrap | 0.11075 | 4.43 | 12 of 12 |
| stationary bootstrap | 0.12458 | 4.98 | 12 of 12 |

The reference row holds at 1.02 times the target, so the instrument reads correctly.  Among the rows that must estimate the coefficient, nothing holds: the best of them is `the exact rule, coefficient pooled` at **3.13 times** the target and the worst is `the exact rule, coefficient per group` at **11.30 times**.  *(measured)*

Two things follow.  First, the exact law buys its guarantee from the coefficient being DECLARED, not from being exact: the same rule with the coefficient estimated separately in each group is the WORST row in the table, worse than any bootstrap.  Second, once the estimate is pooled the exact rule becomes the best of the estimated-coefficient methods by a clear margin, 3.13 times the target against 4.18 for the closest alternative.  Being least bad is not a guarantee and is not offered as one.  *(measured)*

### The stress battery

Every cell runs the same deployed rule; only the data-generating process changes.  A breach is 3 standard errors above the ceiling, declared before the run.  12 of 54 cells breach.

**Breaks it:** ar1 phi=+0.2 raw, ar1 phi=+0.4 raw, short variance window.

**Survives it:** ar1 phi=+0.4 repaired, ar1 phi=+0.6 repaired, ar1 phi=-0.4 raw, garch(.10,.85), normal (control), regime x2, residual correlation 0.3, second undeclared factor, skew chi2_3, t3, t5, undeclared factor and residual correlation.

Worst case **0.09991**, which is **4.00 times** the ceiling, at `ar1 phi=+0.4 raw`.  The direction matters and is recorded per cell: untreated positive serial correlation is violently liberal, negative serial correlation and undeclared factor structure are strongly conservative, and heavy tails and skew are conservative rather than liberal, which is the good direction to fail in.  *(diagnostic)*

The one genuinely new finding is `short variance window`: estimating the variance on a trailing subwindow while the mean uses the whole sample, which is what a rolling risk model does, breaches on its own.  It is not in the paper's list of boundaries and it is a thing desks do by default.

---

## W4. The promotion protocol

`PROMOTION-PROTOCOL.md`, `manifest.schema.json` and `w4_manifest.py`.  The guarantee is a statement about a family fixed before the outcome was seen; nothing in the mathematics can check that, but a file can.

Three verdicts: `GUARANTEED`, `CONDITIONAL`, `SCREEN`.  The example manifest returns **CONDITIONAL**, and its conditions are printed rather than suppressed.

**8 of 8 refusal probes are refused on their own ground:**

* locked after the window closed;
* ratio domain from the tested returns;
* research window overlaps evaluation;
* calibration overlaps evaluation;
* prais-winsten whitening;
* internal level above the target;
* benchmark listed as a candidate;
* missing nuisance domain.

A nuisance estimated on the evaluation period returns **CONDITIONAL**, not a refusal, which is the point: that run is worth doing and its number is worth having, as long as the report does not call it a certificate.  *(diagnostic)*

### Whether these checks bind

`verify_beyond.py` formats every quoted number from a deposit and looks for it in the report, so neither file can be edited into agreement on its own.  That is necessary and not sufficient: a check that never fails proves nothing.  `probe_beyond.py` moves one deposit field at a time, reruns the verifier in a fresh process against the unchanged report, and requires a new failure.

**9 of 9 probes caught**, including **4 of 4 transpositions**.  The transpositions are the ones that matter: they swap two values inside a single table row, so every value is still present in the document and substring containment cannot see the change.  The preceding programme found that seven such probes all passed before table rows were matched whole.

| probe | caught |
|---|---|
| scalar: invariant hidden pp | yes |
| scalar: region max gain | yes |
| scalar: stress breach count | yes |
| scalar: manifest refusals | yes |
| transposition: factor falsifier brute/exact | yes |
| transposition: region monotone/det-floor | yes |
| transposition: serial repair separate/pooled | yes |
| transposition: asymptotic two methods | yes |
| nested: certificate box point count | yes |

Binding is one question and REACH is another.  A number no check formats passes unconditionally, and counting passes cannot see it.  `coverage_beyond.py` instruments the verifier to record every string it looks for, extracts every numeric token from this report, and reports what nothing covers: **202 of 206 distinct tokens covered by a check**, 4 excused with a stated reason, **0 uncovered**.

The excused ones, each named rather than silently dropped:

* `1.47`, the spurious gain the first region module reported, quoted as prose about a superseded run rather than as a result;
* `144`, a bin count in the resolution table, covered by the row check;
* `193`, the verification grid size, a declared parameter;
* `200,000`, the draw count of the invariance check, a declared parameter;

One probe was replaced rather than recorded as a miss: swapping the monotone and unrestricted losses in the region table is a no-op, because the result of that section is precisely that the two are equal.  A vacuous probe is worse than no probe, since it reports a pass.  *(diagnostic)*

---

## What is still open

* **No new certificate.**  Every guarantee in this programme is either inherited from the preceding one or is a grid maximum.  Converting either of the two grid maxima here into a continuum enclosure is the single piece of work that would change the repository's claims rather than its understanding of them.
* **The Hunt-Stein lower bound is a derivation, not a proof carried through.**  The maximal invariant is established and checked; the amenability argument is standard; the conditions listed above are not discharged and the K > 1 extension is not written.
* **The aligned-factor certificate is a grid maximum.**  Making it a continuum certificate needs a two-dimensional enclosure per candidate and the branch-and-bound code is one-dimensional.
* **The loading contrast cannot be estimated from the data the rule sees**, so unlike the variance ratio it can only be declared.  A desk that cannot declare it has no guarantee under aligned returns.
* **Residual cross-correlation after one factor is removed** is measured as conservative in the stress battery, not certified.
* **Tier C has no guarantee and is the realistic case.**  Nothing here changes that; the protocol only makes it impossible to report a tier-C run as a certificate by accident.

