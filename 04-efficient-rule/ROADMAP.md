# Roadmap Beyond Exact Welch–Dunnett

## Purpose

The current work has reached a coherent endpoint for its stated model:

- independent Gaussian groups;
- unknown and unequal variances;
- one shared benchmark/control;
- a fixed, predeclared set of candidates;
- known variance ratios: deterministic exact null evaluation;
- estimated variance ratios: a calibrated or certified upper-bound procedure rather than an exactly similar test.

The next programme should **not** attempt to erase the Behrens–Fisher boundary by wording or numerical force. The boundary is real: for useful non-randomised threshold rules, unknown unequal variances do not generally permit finite-sample exact similarity. The objective is instead to solve the next statistically distinct problems while retaining the present paper’s virtues:

1. separate exact distributional calculations from nuisance estimation;
2. attach guarantees only where their assumptions genuinely support them;
3. quantify the price of each guarantee in power;
4. publish a procedure a desk can pre-specify and run.

The programme has two main research tracks:

1. **Finish the independent-groups theory at the actual boundary.**
2. **Build a new dependent-returns procedure for a realistic desk panel.**

The second is the economically important sequel. The first should be completed first because it produces clean lower/upper bounds and prevents the new paper from overstating optimality.

---

## Workstream 1 — Characterise the remaining Behrens–Fisher gap

### Objective

Turn “the remaining loss appears informational” into the sharpest defensible statement possible.

The current numerical optimisation is already strong: it searches arbitrary piecewise-constant threshold functions of the estimated variance ratio, which approximates the full class

\[
\phi(T,\widehat\lambda)=\mathbf{1}\{T>h(\widehat\lambda)\}.
\]

At moderate degrees of freedom it finds almost no improvement over the calibrated rule. The next goal is to delimit exactly what this does—and does not—establish.

### Question 1: Is \((T,\widehat\lambda)\) the relevant invariant information set?

For the one-candidate normal model, derive the maximal invariant after removing:

- common location under the null;
- common scale;
- any other nuisance transformation admissible under the target testing problem.

Then establish whether a non-randomised invariant test can always be represented as a rejection region in \((T,\widehat\lambda)\), or whether useful information remains in the full sufficient statistic beyond the threshold representation.

**Outputs**

- Formal proposition identifying the statistic available to invariant procedures.
- Explicit distinction between:
  - all tests measurable in the full sufficient statistic;
  - all invariant tests;
  - all \((T,\widehat\lambda)\)-measurable tests;
  - monotone threshold rules \(T>h(\widehat\lambda)\).

### Question 2: What can a non-threshold region recover?

Optimise a general rejection region in the discretised \((T,\widehat\lambda)\) plane, not only an upper-threshold boundary.

For each ratio bin \(b\), let the decision be a binary or randomised acceptance function \(q_b(T)\in[0,1]\), subject to a uniform size constraint:

\[
\sup_{\lambda\in\Lambda}
\Pr_{\lambda,H_0}\{\text{reject}\}
\le \alpha.
\]

Optimise average or least-favourable power against a declared alternative design. Compare:

1. calibrated monotone threshold;
2. best deterministic non-monotone region;
3. best randomised region;
4. oracle known-ratio test.

This does not claim a universal optimum, but it measures precisely what monotonicity and deterministic deployment cost.

### Question 3: Low-degrees-of-freedom regime

The existing results show meaningful recoverable power at 3, 5, and, weakly, 7 df. Treat this as a first-class result rather than an exception.

- Map the best rule over \(\nu_0,\nu_1\in\{3,4,5,7,9,11\}\).
- Identify whether the gain comes from asymmetric treatment of extreme ratio estimates, non-monotonicity, or a different allocation of type-I error across \(\widehat\lambda\).
- Report the gain against both the calibrated rule and the known-ratio oracle.
- State a practical rule: below which minimum df should the current desk table not be used without a low-df enhancement?

### Desired theorem, if obtainable

A full minimax theorem is ambitious. A useful intermediate result would be any of:

- a lower bound on unavoidable oracle-power loss for every uniformly valid \((T,\widehat\lambda)\)-measurable test;
- an admissibility result for the calibrated rule within monotone threshold rules;
- a proof that the optimal deterministic rule is monotone in \(T\), reducing the search to \(h(\widehat\lambda)\);
- a quantified upper bound showing the 36-bin approximation cannot hide more than a stated number of power points.

### Gate

Do not call the current rule “optimal” unless the test class, objective, nuisance domain, alternative weighting, and randomisation policy are all explicit. “Numerically near-optimal among variance-ratio-measurable monotone threshold rules above 11 df” is strong and defensible.

---

## Workstream 2 — The real desk problem: aligned dependent returns

### Objective

Replace independent groups with the economically natural object: candidate and benchmark returns observed on overlapping dates with contemporaneous dependence.

This is not a cosmetic extension. With aligned returns,

\[
\operatorname{Var}(\bar r_i-\bar r_0)
=
\operatorname{Var}(\bar r_i)+
\operatorname{Var}(\bar r_0)-
2\operatorname{Cov}(\bar r_i,\bar r_0),
\]

and the unknown covariance terms alter the testing problem itself.

### Phase 2A — Fully aligned iid Gaussian return panel

Start with the cleanest genuinely new model:

\[
\mathbf r_t=(r_{0,t},r_{1,t},\ldots,r_{K,t})^\top
\overset{\mathrm{iid}}{\sim}
N_{K+1}(\boldsymbol\mu,\Sigma),
\qquad t=1,\ldots,n,
\]

with unrestricted unknown positive-definite \(\Sigma\), and hypotheses

\[
H_{0i}:\mu_i-\mu_0\le 0.
\]

Define paired differences:

\[
d_{i,t}=r_{i,t}-r_{0,t},
\qquad
\mathbf d_t=(d_{1,t},\ldots,d_{K,t})^\top.
\]

Then every marginal candidate comparison is a one-sample problem, but the vector \(\mathbf d_t\) has an unknown covariance matrix \(\Omega\). The problem becomes simultaneous inference for several correlated one-sample means with unknown correlation.

#### Subproblems

1. **Known covariance/correlation oracle**
   - Derive an exact deterministic null evaluator for \(\max_i \bar d_i/(s_i/\sqrt n)\) conditional on a known correlation matrix, or identify the best exact finite-sample representation.
   - Establish special cases: independence, equicorrelation, one-factor covariance, and perfect correlation.

2. **Unknown covariance but confidence-set guarantee**
   - Construct a confidence region \(\mathcal C(\Omega)\) for correlation/covariance from the Wishart sampling law.
   - Use a least-favourable critical value:

     \[
     c_{\mathrm{env}}
     =
     \sup_{\Omega\in\mathcal C(\Omega)} c_{\alpha'}(\Omega).
     \]

   - Allocate error explicitly between covariance-region miss probability and conditional test level, as in the current variance-ratio envelope.
   - Quantify the power price as a function of \(K\), \(n\), and true correlation geometry.

3. **Structured covariance models**
   - One-factor residual covariance;
   - equicorrelation band;
   - sparse precision matrix;
   - sector/block factor models.

   Each structure must be declared before testing. A narrower credible covariance class should buy power in the same way a narrower variance-ratio domain buys power.

4. **Certificate feasibility**
   - Avoid claiming a full \(K(K-1)/2\)-dimensional covariance enclosure unless a real numerical certificate exists.
   - Seek reductions via factor loadings, eigenvalue bounds, Slepian-type monotonicity where applicable, or a conservative projection/certificate analogous to the current per-candidate construction.

#### Primary deliverable

A paper with a clear title such as:

> Exact and certified many-to-one testing for correlated paired returns under unknown covariance.

Its honest endpoint may again be a bound, not equality. That is acceptable and expected.

### Phase 2B — Unequal histories and missing overlap

A realistic promotion gate has an incumbent with a long history and a new candidate with a shorter history. Pairing only overlap dates discards benchmark information.

Study three predeclared designs:

| Design | Data retained | Statistical cost |
|---|---|---|
| Strict common-window pairing | Only overlap dates | Exact paired structure, but wastes incumbent history |
| Two-window model | Common overlap plus benchmark-only prehistory | Requires a stationarity / change-point assumption |
| State-space or hierarchical model | All returns | Model dependence and estimation risk become central |

The priority is not to optimise a rich finance model immediately. It is to specify a minimal normal model under which the longer benchmark history genuinely provides test-relevant information and derive the correct null law.

Possible first model:

\[
r_{i,t}=\mu_i+\beta_i f_t+\varepsilon_{i,t},
\]

with a declared observed factor \(f_t\), Gaussian residuals, and candidate availability indicators. Work out exactly which factor design quantities and residual covariance parameters need to be known or estimated.

### Gate

Do not describe paired differences as “solving the original problem.” They remove Behrens–Fisher variance-ratio uncertainty but replace it with unknown cross-candidate correlation. The sequel must state that trade honestly.

---

## Workstream 3 — Serial dependence and factor estimation

### Objective

Turn the current declared-nuisance transformations into procedures with explicit validity tiers.

### Tier A — Parameters declared known

For an AR(1), GLS, or factor model with parameters fixed externally, retain exact transformations:

- whitening for known serial covariance;
- residualisation for known factor loadings/design;
- loss of observations, degrees of freedom, and effective sample size stated exactly.

This tier can inherit the existing exact/certified result after transformation if the transformed residual groups meet the model assumptions.

### Tier B — Parameters estimated on an independent calibration period

Split time into:

1. a calibration period estimating factor/serial parameters;
2. an evaluation period applying the fixed transform and the selection rule.

This preserves a cleaner conditional-inference interpretation. Quantify the calibration-length versus evaluation-length tradeoff.

### Tier C — Parameters estimated from the same evaluation period

This is the difficult realistic case.

- Build a confidence region for serial/factor parameters.
- Propagate it through a worst-case threshold or a sample-splitting/cross-fitting design.
- Compare to moving-block bootstrap, stationary bootstrap, HAC, and self-normalised alternatives, but label them as asymptotic/empirical unless finite-sample coverage is genuinely proved.

### Stress tests

Every claimed certificate should be attacked under:

- AR(1) values spanning negative to positive persistence;
- GARCH-like conditional heteroskedasticity;
- factor misspecification;
- residual correlation left after factor removal;
- skewed and heavy-tailed innovations;
- rolling versus expanding windows;
- regime shifts.

### Gate

No exactness label may survive an estimated transformation without specifying what is conditioned on and what error is newly introduced.

---

## Workstream 4 — Research protocol and multiplicity

### Objective

Make the testing procedure a credible final-stage promotion gate, not a claim to retrospectively repair an adaptive research process.

### Pre-registration / audit object

Before outcome evaluation, lock a machine-readable manifest containing:

- candidate identifiers and code hashes;
- benchmark identity;
- return construction, including costs/slippage and rebalance timing;
- sample start/end dates and missing-data policy;
- hypothesis direction and FWER target;
- dependence model and transformation tier;
- variance-ratio / covariance parameter domain;
- primary and sensitivity specifications;
- threshold-table version and seed-independent certificate hash.

### Research pipeline separation

Use at least three periods or data roles:

1. **Research/training:** generate hypotheses and tune models;
2. **Calibration:** estimate declared dependence nuisance parameters, if needed;
3. **Promotion evaluation:** apply exactly one locked testing family and threshold table.

If strategy discovery occurred using the evaluation data, report the test as a conditional screen rather than a discovery-valid familywise guarantee.

### Relation to other multiplicity methods

The procedure should be described as complementary to, not a replacement for:

- White’s Reality Check;
- Hansen’s SPA;
- Romano–Wolf stepdown;
- false-discovery-rate procedures;
- deflated Sharpe ratio / probability of backtest overfitting methods;
- walk-forward and embargoed cross-validation.

Each controls a different declared target under different resampling/model assumptions.

### Deliverables

- `PROMOTION-PROTOCOL.md`.
- A JSON schema for the manifest.
- A verifier that refuses to run if required fields or pre-evaluation hashes are absent.
- A report template that labels claims by validity tier: exact, bounded, estimated, asymptotic, or diagnostic.

---

## Validation standard

Every future workstream should preserve the current paper’s unusually good discipline: falsifiers before persuasive examples.

### Required validation categories

| Category | Requirement |
|---|---|
| Algebraic limits | Recover independent, perfect-correlation, one-candidate, and known-parameter special cases |
| Numerical refinement | Show convergence under independently implemented quadrature/integration rules |
| Brute-force simulation | Compare against large Monte Carlo simulations only where a simulation is an independent check, not the definition of truth |
| Coverage | Measure confidence-region miss probability against its claimed coverage |
| Certificate challenge | Search interior nuisance points more finely than the construction grid |
| Code independence | Maintain at least one independent checker that reads deployed tables rather than regenerating them through the same logic |
| Robustness | Deliberately violate normality, independence, factor structure, and stationarity; report where guarantees fail |
| Reproducibility | Pin environment, hash inputs, regenerate every table/figure, and preserve failed checks |

### Claim taxonomy

Every headline number must carry one of these labels:

- **Exact in model:** deterministic evaluation of the correct finite-sample quantity conditional on stated nuisance values;
- **Certified bound:** a proved upper/lower bound over a stated nuisance domain;
- **Measured:** simulation estimate with Monte Carlo uncertainty;
- **Asymptotic:** valid only in a stated limiting regime;
- **Diagnostic:** stress-test result, not a guarantee.

---

## Suggested sequence

### Phase I — 4 to 8 weeks

- Complete Workstream 1's invariant/statistic classification, and the Hunt-Stein route to a lower bound.
- Bound what the 36-bin discretisation can hide.
- Expand the low-df optimisation and write the resulting operational warning.
- Falsify the latent-factor reduction before building on it.

### Phase II — 2 to 4 months

- Solve the aligned iid Gaussian paired-return problem with known covariance.
- Identify whether a low-dimensional conditional representation exists under one-factor/equicorrelation structures.
- Build a covariance-confidence-set envelope prototype for small \(K\).

### Phase III — 4 to 8 months

- Handle unknown covariance under declared structure.
- Compare strict pairing with use of longer benchmark history under a minimal two-window model.
- Build the first certificate and power map.

### Phase IV — 3 to 6 months

- Add serial dependence / factor estimation as validity tiers.
- Integrate the research-protocol manifest and promotion-gate tooling.
- Run a realistic desk-style case study only after the statistical scope is fixed.

---

## Success criteria

The programme succeeds if it produces one of the following—not necessarily all of them:

1. A theorem or sharp numerical bound locating the remaining independent-groups information gap beyond \((T,\widehat\lambda)\)-threshold rules.
2. A certified many-to-one procedure for aligned correlated returns under a declared covariance class.
3. A quantified impossibility/price result showing how much power is lost when covariance is unknown.
4. A transparent protocol that separates a valid final promotion test from adaptive strategy discovery.

The programme fails if it relabels a simulation as exact, a measured grid maximum as a continuum certificate, an estimated nuisance transformation as known, or a retrospective candidate set as a predeclared family.

---

## Final positioning

The existing paper should be treated as the completed first module:

\[
\text{independent groups}
+
\text{unknown unequal variances}
+
\text{shared-control sampling dependence}.
\]

The next paper should not be sold as “adding realism” in a vague sense. It should name the changed object:

\[
\text{aligned dependent return vectors}
+
\text{unknown residual covariance}
+
\text{possibly unequal histories}.
\]

That is the real unresolved desk problem. It is harder for a reason, but it is also more valuable: it replaces an elegant clinical-design model with the dependence structure strategy selection actually has.

---

## Amendments adopted before execution

Three changes were made to this plan before any of it was run.  They are
recorded here rather than folded in silently.

**A1.  Workstream 2 keeps the returns unpaired.**  As originally written, Phase
2A formed paired differences `d_{i,t} = r_{i,t} - r_{0,t}` and faced a maximum
of K studentised statistics with different studentisers and a Wishart-coupled
denominator set.  No conditioning variable factorises that, so the K-dimensional
integral returns and the existing architecture is lost.

Instead the dependence is carried by a LATENT factor and the statistic is left
alone:

    r_{i,t} = mu_i + beta_i f_t + eps_{i,t},   i = 0..K,   t = 1..n,

with `f` latent Gaussian and `eps` independent across `i`.  Conditioning on FOUR
quantities -- `fbar`, `Vf = sum_t (f_t - fbar)^2`, `Xbar_0`, `s_0^2` -- makes the
K candidates exactly conditionally independent, because `eps_bar_i` is
orthogonal to the centred `eps_i`, and because `(n-1) s_i^2 / sigma_i^2` is
noncentral chi-square on `n-1` degrees of freedom with noncentrality
`beta_i^2 Vf / sigma_i^2`, which sees the factor path only through its sum of
squares.  Two of the four conditioning quantities are the two the existing paper
already integrates.  The desk's statistic does not change.

**A2.  The minimax lower bound goes first, not last.**  The maximal invariant of
the one-candidate problem under common location and common scale is
two-dimensional, and `(T, lambda_hat)` is it.  The location-scale group is
solvable, hence amenable, so Hunt-Stein applies and a minimax test may be taken
invariant.  A numerical bound over `(T, lambda_hat)`-measurable tests then
becomes a bound over ALL tests, which is success criterion 1.  This was listed
as "ambitious"; it is days of work and it strengthens the existing paper, so it
runs first.

**A3.  Inherit the programme's dead ends.**  A jointly designed rule and
certificate was already built and loses (45% of level against 6% for the
per-candidate construction); the confidence-region rule is superseded; one
factor falsifier was mis-specified and discarded rather than deposited.  Any
workstream that re-opens those must say why.

And one addition: every phase carries a falsifier before it carries a
deliverable.  The only defect in the preceding programme that changed a result
was found by a falsifier, not by a check.
