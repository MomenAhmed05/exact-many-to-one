# Milestones against the plan, with status and evidence

Plan: `welch_dunnett_efficiency_frontier_plan.md`.  Every status below is the status of a deposit
under `deposits/`, not of an intention.

---

## M1, efficiency ledger complete

**DONE.**  `deposits/wa_ledger.json`, `deposits/wb_k1.json`, `deposits/wb_general.json`.

The gate the plan sets is that no headline claim may be made until numerical tightness and economic
efficiency are reported separately.  They are: the branch-and-bound gap of the preceding programme
(`~1e-6`, a numerical quantity) and the power loss (72 pp at a 36-month benchmark, an economic one)
are different numbers about different things, and the ledger reports both with their objects named.

**The classification the plan asks for.**  The existing certified region rule is
**economically unusable under the intended desk designs**: realised FWER 0.00018 against a 0.025
ceiling, threshold 2.4 to 2.7 times the oracle, 64 to 72 percentage points of power given up.

## M2, measured unconditional candidates

**DONE.**  `deposits/wb_k1.json`, `deposits/wb_general.json`, `deposits/wb_percand.json`.

The gate is that at least one candidate must materially improve desk-map power before expensive
certification work expands.  The uniformly calibrated plug-in improves it from ~2% power to
within 1.5-2.6 pp of the oracle at K = 2, 3, 5 with seven or more degrees of freedom on the smallest
group.  The gate is passed by a wide margin, which is why M3 was attempted.

Three rule classes were built, not one:

| class | valid? | tier | worst loss |
|---|---|---|---|
| joint plug-in at a calibrated level | yes | measured above K = 1, certified at K = 1 | 1.8-12.6 pp |
| binned multiplier on it (the in-class optimum) | yes | measured | 20.9 pp at df3, ties above df11 |
| per-candidate threshold | yes | measured, but **certifiable at any K** | 6.9-24.0 pp |

## M3, certified unconditional-valid rule

**DONE AT EVERY K.**  `deposits/wb_certified_rule.json`, `deposits/certify_anyk.json`,
`deposits/gates.json`.

Ten designs, all certified: the uniform constraint is proven over the **continuum** of ratios in
`[1e-3, 1e3]` by rigorous enclosure plus branch-and-bound, not checked on a grid.  The enclosure's
monotone corner substitutions are falsified on 1,800 random cells (gate X2, zero violations), its
bracket is verified to contain the quadrature value at every ratio tested (X3), and every deposited
rule is re-checked on a grid three times finer than the one it was calibrated against (X6).

Above K = 1 the PER-CANDIDATE rule is certified over the continuum at every K and on either side
(`deposits/certify_anyk.json`, ten certificates, all pass).  The route is a union step rather than a
bigger enclosure: `sup` over the ratio cube is at most `K` times the one-candidate supremum, which
the K = 1 enclosure already certifies tightly, so certifying at `abar/K` suffices.  The measured
slack of that step at the binding ratios is `5.4e-05` to `2.25e-03`.

The JOINT plug-in is certified too, at K = 2 and K = 3 (`deposits/certify_joint.json`), by
deploying it as a step table and projecting that table onto each candidate: rejection of candidate i
requires the statistic to exceed the table's minimum over the other candidates' ratios, so the union
step applies to the projection.  It buys 10 to 23% more level than the per-candidate
certificate.  A direct four-dimensional enclosure was built first, dominates, and was abandoned for
being loose by a factor of two.

## M4, minimax lower-bound analysis

**DONE AT K = 1.**  `deposits/wc_frontier.json`.

Rather than a discretised primal-dual programme, the K = 1 structure admits the frontier directly:
maximise power over piecewise-constant threshold functions of the estimated ratio, subject to the
exact uniform constraint, with analytic gradients.  Any rule that sees only `lambda_hat` is of that
form up to bin resolution, so the optimum bounds what any such rule can do.

**The decision the gate asks for: the conservatism is NOT removable above eleven degrees of
freedom.**  The optimiser gains 4.24 pp at three, 1.30 at five, 0.49 at seven, and -0.02, +0.01,
+0.06 at eleven, nineteen and thirty-five.  The 4 to 6 pp that remains is an information constraint.

The plan wanted the dual to explain *why* the threshold must be large.  The explanation is
algebraic and does not need a dual: with `rho = lambda n_0/n_1`,

    sigma^2/d^2 = (rho + 1)/(rho u_1 + u_0) ,

so matching the oracle requires `(rho u_1 + u_0)/(rho + 1)` while the data supply only
`rho_hat = rho u_1/u_0`.  One equation, two unknowns, and therefore unavoidable over-protection in
some states.

## M5, information-set design map

**DONE.**  `deposits/wd_designmap.json`.

A 7 x 5 grid over benchmark and candidate degrees of freedom.  Minimum degrees of freedom on the
**smaller** group for a given power budget: 19 for 5 pp, 11 for 10 pp, 5 for 25 pp, 3 for 50 pp.

The map also corrects a claim in the paper: the binding quantity is `min(nu_0, nu_1)`, not the
benchmark's degrees of freedom.  Grouping by the minimum leaves a spread of 13.1 pp falling to
0.08 pp; grouping by the benchmark's leaves 11.6 pp rising to 43.4 pp.  The surface is symmetric
under swapping the two groups to the last digit.

**Section 6.2, ex ante bounded ratio classes: DONE and negative.**  `deposits/wd_restricted.json`
recalibrates on six declared domains from three decades down to a 1.5x volatility target.  The
recovery is 0.0 to 2.1 pp, because the least-favourable ratio sits near `lambda ~ 0.3` - inside even
a factor-of-two box - so narrowing the domain removes states the rule was not protecting against.
Worth reporting precisely because it closes a plausible-sounding route.

**NOT DONE** from this workstream: auxiliary volatility data (§6.3) and alternative estimands (§6.4)
beyond the factor-alpha case handled in M6.  Both change the information set and each needs its own
coverage argument.

## M6, two-sided and dependence extensions

> **Amended by M9.**  The whitening transfer below is stated unconditionally and is not:
> it holds at equal group lengths, or with the first observation dropped rather than
> rescaled.  See M9.

**DONE at K = 1; scoped above it.**  `deposits/twosided.json`, `deposits/serial.json`,
`deposits/factor.json`.

The gate is that no extension may inherit theorem-grade language without its own certificate and
falsifier suite.  Each of the three is reported against that gate separately:

| extension | status | what carries it |
|---|---|---|
| **two-sided** | CERTIFIED at K = 1 | an exact identity, not a new certificate: `FWER_2s = 2 FWER_1s` at every draw, falsified to 0.0 over 300 designs, so the deposited one-sided certificates transfer at half the level.  MEASURED above K = 1, with the enclosure obstruction named (the two-sided factor is non-monotone in `s`) |
| **serial dependence** | CERTIFIED given a declared coefficient | Prais-Winsten whitening makes the problem exactly iid, so every certificate transfers unchanged.  Verified at 40,000 draws: z from -1.84 to -0.42 across `phi` in {0, .05, .1, .2, .4}.  Misspecification and a declared-bound rule are measured |
| **factor exposure** | CERTIFIED given declared loadings | residualising makes the problem exactly iid for the factor-alpha estimand.  Estimated loadings are a labelled plug-in, and the deposit measures the cost of the one mistake a practitioner will make - not returning the regression's degree of freedom to the rule |

**The unifying result**, which neither the plan nor I anticipated: serial dependence and factor
exposure are not two boundaries but one.  A **declared** nuisance is a pre-processing step and the
certificate survives it; an **estimated** nuisance is a new parameter and the certificate does not.
No new theorem was needed for either, and none is claimed.

The two-sided enclosure above K = 1, listed as open in an earlier draft, is closed by the same
union step at `abar/(2K)`.  What remains genuinely open: a certificate for the JOINT plug-in above
K = 1, and a certified treatment of estimated loadings or an estimated autocorrelation.

## M7, external reproducibility

**KIT WRITTEN AND SEALED; NOT EXECUTED.  Cross-check done.**  `replicate/`,
`deposits/crosscheck.json`.

The gate is that only after an external replication may the result be called independently
reproducible end to end.  **It may not be so called**, and the report says so in the same paragraph
that reports the cross-check.

What was done instead, and is weaker: every decisive quantity re-derived by a route sharing no code
path with the one that produced it - the thresholds against `HeteroMaxT` (worst 1.1e-03), the FWER
and power against Monte Carlo with the exact inner (worst |z| 1.36 and 2.50), and the certified bound
against a direct sample at 42 ratio points (**0 breaches**).  That cannot substitute for a blind run,
but it does exclude the failure mode a single implementation cannot see: all three routes being wrong
in the same direction.

The kit is built around the defect that broke the preceding programme's replication.  That one failed
its draw-dependent half because the specification never pinned the draw recipe, the definition of
`d`, or the semantics of a certificate point.  **This construction has no draws in the decision path
at all** - the rule, the constraint and the certificate are all deterministic - so that entire class
of ambiguity cannot arise.  Monte Carlo appears only in falsifier R3, where the seed need not match.

---

## M8, the paper carries the programme

**DONE.**  `paper_content_frontier.py`, `paper_build_frontier.py`,
`paper_figs_frontier.py`, `paper_verify_frontier.py`, `coverage_frontier.py`.

Section 9.11 (the certified supremum over the continuum), Section 10 (the efficient rule, 10.1 to
10.10), Section 11 (dependence as a declared nuisance, 11.1 to 11.3), equations 37 to 46, Tables 15
to 27 and A3, Figures 14 to 20, with the Abstract, section 1.4, Limitations and Conclusion extended
and old Sections 10 to 12 renumbered to 12 to 14.

The gate is that no number reaches the document except by being formatted from a deposit at build
time.  **540 checks pass**, and `coverage_frontier.py` reports **100%% of the new material's numeric
tokens covered by at least one check**, with six uncovered for a reason it names.  The checks bind:
347 of 396 fail against the pre-update document, one perturbed deposit field fails exactly its own
check, and seven corruption probes each fail on the right check.

## M9, adversarial audit, and the one result it changed

**DONE.**  Four passes; `REPORT-FRONTIER.md` section 8a carries the yield table and the three
findings that were rejected on the evidence.

The pass scoped to mathematics and prose logic found the only defect that changed a result.  M6
recorded the whitening transfer as unconditional.  Prais-Winsten rescales the first observation, so
the transformed mean is `mu * k(phi)` with `k` not constant in `t`, and the plain sample mean of a
whitened series estimates `mu * mean(k)`, which depends on the group's LENGTH.  At `n0 = n1` the two
groups share it and the transfer holds for any `mu`; at `n0 != n1` the numerator carries a bias in
`mu` and the rule drifts off its ceiling, **%.5f at 36 against 13, phi 0.4, mu 2, which is %.1f
standard errors above** (`deposits/serial_falsify.json`).  M6's own Table 25 was run at `n0 = n1`
and `mu = 0`, the one cell where both failure modes vanish.

**The repair is a change of transform.**  Cochrane-Orcutt drops the first observation instead of
rescaling it, the mean vector is then constant at `1 - phi`, and the transformed sample is a
location model on `n - 1` observations.  Worst deviation at the unbalanced design **%.5f (%.1f
se)**, and the realised rate is identical to five decimals across every `mu` tested, which is the
invariance a constant mean vector predicts.  M6's conclusion stands with the transform named.

**What this says about the programme.**  Every certificate in M3 and M7 was unaffected: the defect
was in a transfer argument, not in an enclosure.  But it was invisible to 540 numeric checks,
because the numbers were right about the cell that was run.  A deposit can only be as general as its
design grid, and the design grid is the thing no verifier reads.

## What the plan's final success criteria look like now

| # | criterion | status |
|---|---|---|
| 1 | a valid finite-sample one-sided iid Gaussian rule exists | **done** (preceding programme) |
| 2 | its numerical enclosure and continuous nuisance supremum are certified | **done** (preceding programme) |
| 3 | its practical power cost is measured against like-for-like references | **done** - and it is 64-72 pp, which is why this programme exists |
| 4 | a better valid rule is found, or a lower bound shows why one cannot exist | **both**: a rule losing 1.5-4.5 pp instead of 72, and a lower bound showing the residual is an information constraint above eleven degrees of freedom.  Both the per-candidate and the joint rule are certified over the continuum |
| 5 | a quantitative design map | **done** |
| 6 | two-sided, factor and serial-dependence claims separately solved or excluded | **solved**: two-sided certified at every K via the exact factor of two plus the union step; serial dependence and factor exposure by declared-nuisance pre-processing.  Estimated loadings and estimated autocorrelation remain labelled plug-ins |
| 7 | the full draw-dependent decision path independently reproduced | **kit sealed, not executed.**  There is no draw-dependent path in this construction, which removes the class of defect that broke the last replication; an independent numerical cross-check passes on every decisive quantity |

## Terminal conclusion

The plan offers two acceptable endings.  The evidence supports a compound of both, and the split is
at eleven degrees of freedom on the smaller group:

> Under independent Gaussian samples across groups, for one-sided separate-variance Welch
> comparisons of K candidates against a shared benchmark, a data-computable threshold with
> finite-sample familywise error at most `abar` is obtained by running the exact fixed-ratio
> constant at the estimated ratios at an internal level calibrated against the exact unconditional
> constraint.  At K = 1 that constraint is certified over the continuum of ratios.  The rule costs
> 4.4 to 4.5 percentage points of power at eleven or more degrees of freedom on the smaller group,
> against 64 to 72 for the previously certified region rule.  Above eleven degrees of freedom it is
> optimal within the class of all threshold functions of the estimated ratios, so the residual cost
> is a property of the information available and not of the procedure.  Below five degrees of
> freedom the cost is 21 to 40 percentage points, it is only partly removable, and the design map
> states how many observations buy it back.

## W1-W4: beyond the independent-groups model

Executed against `ROADMAP.md` with its Workstream 0 removed and three amendments recorded there.  Report: `REPORT-BEYOND.md`, checked by `verify_beyond.py`.

### W1: the remaining gap, characterised

**DONE.**  `w1_invariant.py`.  The maximal invariant of the one-candidate problem under the affine group is two-dimensional and `(T, lambda_hat)` is it: invariance drift 4.25e-12, maximality recovery error 1.46e-11 over 200,000 constructed pairs per design.  So Section 10 searched the invariant class, not a convenient family.  The affine group is solvable hence amenable, so Hunt-Stein gives that a minimax test may be taken invariant, which is the route to a lower bound over ALL tests.  That step is a **derivation, not a proof carried through**, and three conditions are listed unclosed.

The 36-bin discretisation hides at most **0.017 points**, measured by re-solving at 72 and 144 bins with the shrink step applied to each.

**DONE.**  `w1_region.py`.  Monotonicity and determinism priced by imposing monotonicity as a linear constraint, so both arms are optima of the same program on the same grid.  Cost of the two restrictions together: **-0.00 to +0.00 points**.  The first version of this module rendered the threshold rule onto the T grid instead and dropped the straddling cell, which handed the unrestricted rule about a point for free; that is why the comparison is now a constraint.

**DONE.**  `w1_lowdf.py`.  Both group sizes moved over (3, 4, 5, 7, 9, 11).  Material in 26 of 36 cells at a 0.50-point bar declared in advance, best **4.09 points**.  The shipped table is deployable unaltered when the SMALLER group has more than 5 degrees of freedom; at or below that a design-specific threshold function is worth 0.50 points or more.

### W2: aligned returns

**DONE, at the oracle tier.**  `w2_factor.py`.  The roadmap's paired-difference route was abandoned: it produces a maximum of K studentised statistics with different studentisers and a Wishart-coupled denominator set, which nothing factorises.  Keeping the statistic unpaired and carrying the dependence in a LATENT factor instead, conditioning on (fbar, Vf, Xbar_0, s_0^2) makes the K candidates exactly conditionally independent: a **4-dimensional outer rule over K one-dimensional inner rules**, cost linear in K, two of the four dimensions already present in the existing paper.

Worst |z| **1.58** over 8 cells against 4,000,000 brute-force draws.  Non-vacuity is built into the same table: the independent-groups law is wrong by up to 90.2% here and the evaluator tracks the truth across that range.

An earlier Monte Carlo probe of the same claim reported it FALSIFIED at 8.75 sigma.  That probe sampled one variable by picking a Gauss-Legendre NODE INDEX uniformly, which is not a draw from the law; corrected, it agrees.  Recorded because a falsifier that is itself wrong is the most expensive kind.

**PARTIAL.**  `w2_cert.py`.  The numerator of the aligned statistic carries the SQUARED CONTRAST of the loadings and the denominator carries their SUM OF SQUARES, so the iid rule is conservative when the loadings share a sign and liberal when they do not.  Measured: **9 of 18 cells liberal, and the split is exactly along the sign** (9 of 9 opposite-sign, 0 of 9 same-sign), worst **0.07155**, **2.86 times** the ceiling.  A long-short candidate against a directional benchmark is the opposite-sign case.

The first version of this module set the benchmark loading to zero, which is exactly the boundary between the two regimes, and found nothing.  The nuisance is **three-dimensional per candidate**, and only the variance ratio is estimable from what the rule sees.  A constant holding over the 175-point declared box costs 1.00 times the worst known-parameter oracle constant and 3.11 times the best, and that supremum is a **grid maximum, not a continuum certificate**.

**DONE.**  `w2_unequal.py`.  The benchmark's prehistory buys **0.07 points** through the two-window design, which needs constant variance only, and that is essentially nothing; and **25.24 points** through full history, which needs a constant mean too.  Breaking that assumption in the dangerous direction, a prehistory running BELOW the later mean, takes the worst rate to **0.57120**, 22.85 times the ceiling, on `full_history`.  There is no free lunch in the incumbent's track record.

### W3: validity tiers, a new defect, and its repair

**DONE, and it found the programme's second result-changing defect.**  `w3_tiers.py`, `w3_serial_repair.py`.

Tier A holds exactly, as Section 11 says.  Tier B and tier C do not, and not marginally: with the coefficient estimated on the evaluation period the worst realised rate is **0.28297** against a ceiling of 0.025.

The mechanism is that Cochrane-Orcutt maps a mean `mu` to `mu(1 - phi)`, so with the coefficient estimated SEPARATELY IN EACH GROUP the numerator acquires `mu (phi_hat_0 - phi_hat_1)`, zero in expectation but not in any sample, and spreading linearly in the common mean.  It is the Prais-Winsten defect's shape arriving by another route: there the multiplier varied with the group's LENGTH, here with the group's ESTIMATE.

**The repair is to pool**: one coefficient from both groups, applied to both, so the multiplier is identical by construction and the common mean cancels whatever the estimate is.  Separate worst **0.48409** (19.4 times the ceiling); pooled worst **0.07918** (3.2 times).  The signature is that the pooled rate is flat in the common mean: spread 0.00102 at tier B and 0.00100 at tier C.  What pooling does NOT remove is imperfect whitening at an estimated coefficient, which the mu = 0 rows isolate, and that residual breaches on its own: every cell in both arms is above the ceiling.  **Pooling is an improvement, not a repair.**  Section 11's serial result is a tier-A result.

**DONE.**  `w3_stress.py`.  12 of 54 cells breach a bar of 3 standard errors declared in advance.  Survives heavy tails, skew, GARCH, a variance regime shift and negative serial correlation, all conservatively.  Breaks on untreated positive serial correlation (worst 0.09991, 4.00 times the ceiling) and on a **short variance window** -- estimating the variance on a trailing subwindow while the mean uses the whole sample, which is what a rolling risk model does.  That last one is not in the paper's list of boundaries and is a thing desks do by default.

### W4: the promotion protocol

**DONE.**  `PROMOTION-PROTOCOL.md`, `manifest.schema.json`, `w4_manifest.py`.  Three verdicts (GUARANTEED, CONDITIONAL, SCREEN); 8 of 8 refusal probes refused on their own ground; a nuisance estimated on the evaluation period returns CONDITIONAL rather than a refusal, because that run is worth doing as long as the report does not call it a certificate.

**DONE.**  `w3_asymptotic.py`.  HAC and block bootstraps measured against the exact rule on the same draws.  With the coefficient DECLARED the exact rule holds at 1.02 times the target and nothing else is given that information.  Among the rows that must estimate it, none holds: the exact rule with a pooled estimate is best at 3.13 times, the closest alternative is 4.18 times, and the exact rule with separate estimates is the WORST row in the table.  The exact law buys its guarantee from the coefficient being declared, not from being exact.

### The checks were shown to bind

**DONE.**  `probe_beyond.py`.  9 of 9 probes caught, including 4 of 4 transpositions, each moving one deposit field and requiring a new verifier failure against the unchanged report.  One probe was replaced rather than recorded as a miss: swapping the monotone and unrestricted losses is a no-op, because the result of that section is that they are equal, and a vacuous probe reports a pass.

### M10 - the paper carries W1, W3 and the factor case

**DONE.**  `paper_content_beyond.py`, `paper_figs_beyond.py`, `paper_verify_beyond.py`, `coverage_paper_beyond.py`, with the insertion driver extended.

New subsections 10.8, 11.4, 11.5 and 11.6; equations 47 to 50; Tables 28 to 34; Figures 18, 22 and 23, with the three that 10.8 displaced renumbered to 19, 20 and 21 so the figures run in document order; the old 10.8 to 10.10 renumbered to 10.9 to 10.11; contributions (14) to (16); and additions to the abstract, the limitations and the conclusion.

**166 new checks pass** beside the existing 558, and `coverage_paper_beyond.py` reports **100% of the new material's numeric tokens covered**, with none uncovered.  Getting there moved two declared parameters out of prose and into their deposits, because a number typed into a sentence is a number no check can bind.

**What the build found that the modules had not.**  Four defects, none of them numerical:

* the first draft used nine capitalised lead-ins and forty capitalised words for emphasis, in a document whose only all-capital words are FWER, ANOVA and JSON, plus one markdown bold that would have printed as literal asterisks;
* inserting a figure into 10.8 put Figure 21 between 17 and 18, so the figures no longer ran in document order;
* Section 9.10 called the general factor case open three sentences before the amended text said 11.6 implements it;
* Section 11.3 said both nuisance repairs are available to a desk that plans ahead, which 11.4 then measures to be false for the serial coefficient.

The first was found by reading the rendered pages, which is the only way it could have been; the other three by cross-section consistency checks that now live in the verifier.

**What the paper does NOT now claim.**  Section 11.6 presents the reduction as the one Section 9.10 already names, carried out at more than one candidate, not as a new construction; the bound of 10.8 is over invariant tests with Hunt-Stein labelled a route; and the supremum over the loading box is labelled a grid maximum wherever it appears.

### What this programme did not close

* the Hunt-Stein lower bound, as a proof rather than a route;
* the aligned-factor certificate, as a continuum enclosure rather than a grid maximum;
* tier C, which remains without a guarantee and remains the realistic case;
* residual cross-correlation after one factor is removed, measured as conservative but not certified;
* **no new certificate at all.**  Every guarantee here is inherited from the preceding programme or is a grid maximum.

