# The invariance reduction, proved

Section 10.8 of the paper bounds every rule measurable in `(T, lambda_hat)` by a
linear program and calls the extension to all tests a route with three
conditions open.  This closes it, for a criterion that has to be stated
precisely because the one Section 10.8 currently reports is not the one the
argument covers.

Everything numerical here is deposited and checked: `w1_maximal_k.json`,
`w1_minimax.json`, `w1_noninvariant.json`.

---

## 1. Setting

`K + 1` independent normal samples, the benchmark indexed 0:

    X_{i,1} .. X_{i,n_i}  iid N(mu_i, sigma_i^2),   i = 0 .. K,

all mutually independent.  The parameter is

    theta = (mu_0 .. mu_K, sigma_0^2 .. sigma_K^2)
          in  Theta = R^{K+1} x (0, inf)^{K+1} .

The sufficient statistic is

    S = (Xbar_0 .. Xbar_K, s_0^2 .. s_K^2) ,

of dimension `2K + 2`.  Write `lambda_i = sigma_i^2 / sigma_0^2` and let
`Lambda` be the declared ratio domain.  Write

    Theta_0 = { theta : mu_0 = mu_1 = .. = mu_K, lambda in Lambda }

for the boundary of the global null inside the declared domain, of dimension
`K + 2`, and let `A` be the declared alternative region, a set specified in
terms of the standardised shifts and the ratios.

The class under consideration is

    Phi_abar = { phi : sup over Theta_0 of E_theta[phi] <= abar } .

---

## 2. The group

    G = { g_{a,b} : a in R, b > 0 },     g_{a,b} x = a + b x ,

applied to **every sample at once**, with
`g_{a,b} g_{a',b'} = g_{a + b a', b b'}`.  This is the affine group of the line.

It acts on the sufficient statistic and on the parameter by

    (Xbar_i, s_i^2)  ->  (a + b Xbar_i,  b^2 s_i^2) ,
    (mu_i, sigma_i^2) ->  (a + b mu_i,   b^2 sigma_i^2) .

Nothing smaller is available and nothing larger is admissible.  Scaling the
samples by different factors would change the ratios and so change the problem;
`b > 0` is needed because `b < 0` would reverse the one-sided alternative.

**The problem is invariant.**  `Theta_0` maps into itself because equal means
stay equal and the ratios are untouched, so `Lambda` is preserved.  `A` maps
into itself because the standardised shifts and the ratios are both preserved.

---

## 3. Lemma 1: the maximal invariant, at any K

> **Lemma 1.**  `V(S) = (T_1 .. T_K, lambda_hat_1 .. lambda_hat_K)`, with
> `T_i = (Xbar_i - Xbar_0) / sqrt(s_i^2/n_i + s_0^2/n_0)` and
> `lambda_hat_i = s_i^2 / s_0^2`, is a maximal invariant for `G` acting on the
> sufficient statistic.

**Invariance.**  Under `g_{a,b}` the numerator of `T_i` becomes
`b (Xbar_i - Xbar_0)` and the denominator becomes
`sqrt(b^2 s_i^2/n_i + b^2 s_0^2/n_0) = b sqrt(...)`, so `T_i` is unchanged;
`lambda_hat_i` is a ratio of two quantities that both scale by `b^2`.

**Maximality.**  Suppose `V(S) = V(S')`.  Put

    b = sqrt(s_0'^2 / s_0^2) > 0 ,      a = Xbar_0' - b Xbar_0 .

Then, in order and with nothing else fitted:

* `s_0'^2 = b^2 s_0^2` and `Xbar_0' = a + b Xbar_0`, by construction;
* `lambda_hat_i' = lambda_hat_i` gives `s_i'^2 / s_0'^2 = s_i^2 / s_0^2`, hence
  `s_i'^2 = b^2 s_i^2`;
* therefore `d_i' = b d_i`, where `d_i` is the denominator of `T_i`;
* `T_i' = T_i` then gives `Xbar_i' - Xbar_0' = b (Xbar_i - Xbar_0)`, and
  adding `Xbar_0' = a + b Xbar_0` gives `Xbar_i' = a + b Xbar_i`.

So `S' = g_{a,b} S`.  ∎

**The dimension count agrees.**  `2K + 2` coordinates, a two-dimensional group,
so a maximal invariant has `2K` coordinates, which is what `V` has.

**The check is not a tautology.**  The group element is fitted on the benchmark
pair alone and then verified against all `2K + 2` coordinates, leaving `2K`
free.  `w1_maximal_k.py` runs it at `K = 1, 2, 3, 5` on 200,000 independently
constructed pairs per case, with the second point built from the invariant
rather than by applying a group element to the first.  Worst invariance drift
7.22e-12, worst reconstruction error 2.91e-11.

The corresponding maximal invariant on the parameter is
`(eta_1 .. eta_K, lambda_1 .. lambda_K)` with
`eta_i = (mu_i - mu_0)/sqrt(sigma_i^2/n_i + sigma_0^2/n_0)`, by the same
argument applied to `Theta`.  Orbits in `Theta` are exactly the level sets of
that map.

---

## 4. The criterion, and why it had to change

This is where the work actually was.

Section 10.8 reports the value of a program that maximises **weighted average
power** over five declared alternatives subject to the uniform size constraint.
The invariance argument does not cover that criterion, and the reason is worth
stating because it is the trap.

For an amenable group the device is to average a test over the group.  The
resulting invariant test has, at a parameter point, the **mean of the original
test's power over that point's orbit**.  For a supremum criterion a mean over
the orbit is at least the infimum over the orbit and the argument closes.  For a
weighted sum of powers evaluated at chosen representatives it does not: the
orbit mean can sit below the value at the representative, and nothing follows.

So the criterion is restated as **worst-case regret against the known-ratio
oracle**.  Let `beta*(theta)` be the power of the test that knows
`lambda(theta)` and uses the exact fixed-ratio constant at level `abar`.  Define

    R(theta, phi) = beta*(theta) - E_theta[phi] ,   theta in A ,
    Rsup(phi)     = sup over A of R(theta, phi) .

`beta*` depends on `theta` only through `(eta, lambda)`, so `beta*(g theta) =
beta*(theta)`: **the loss is invariant**.  That is the property the whole
argument turns on, and it is the reason regret works where weighted average
power does not.

For an invariant test the power is constant on orbits and `Rsup` reduces to a
supremum over the invariant.  For a general test the supremum runs over every
point of every orbit, so it picks up the worst representative.  That is the
quantity the reduction equates.

---

## 5. Lemma 2: the group is amenable

> **Lemma 2.**  `G` is amenable.

`N = { g_{a,1} }` is a closed normal subgroup isomorphic to `(R, +)`, which is
abelian; `G / N` is isomorphic to `(R_{>0}, x)`, which is abelian.  So `G` is
solvable.  Abelian locally compact groups are amenable, and an extension of an
amenable group by an amenable group is amenable, so `G` is amenable.  ∎

This is the substantive condition and it is the one that fails in the classical
counterexamples: the full linear group `GL(p)` for `p >= 2` is not amenable, and
the Hunt-Stein conclusion is false for it.  The affine group of the line is on
the right side of that line.

---

## 6. Proposition

> **Proposition.**  Under conditions (C1) to (C4) below,
>
>     inf over Phi_abar of Rsup  =  inf over invariant members of Phi_abar of Rsup .

**Proof.**  By Lemma 2 there is a left-invariant mean `m` on the bounded
measurable functions on `G`.  For `phi` in `Phi_abar` define

    phi*(x) = m_g [ phi(g x) ] .

`phi*` is almost invariant, and by (C4) may be taken invariant.

*It stays in the class.*  For `theta` in `Theta_0`,

    E_theta[phi*] = m_g [ E_{g^{-1} theta} phi ]
                 <= sup over Theta_0 of E[phi]  <=  abar ,

because `Theta_0` is `G`-invariant, so `g^{-1} theta` is in `Theta_0` for every
`g` and the mean of a family of numbers is at most their supremum.

*Its worst regret is no larger.*  For `theta` in `A`, using the invariance of
`beta*` to take it inside the mean,

    R(theta, phi*) = beta*(theta) - m_g [ E_{g^{-1} theta} phi ]
                   = m_g [ beta*(g^{-1} theta) - E_{g^{-1} theta} phi ]
                   = m_g [ R(g^{-1} theta, phi) ]
                  <= sup over A of R(., phi)  =  Rsup(phi) ,

because `A` is `G`-invariant.  Taking the supremum over `theta` in `A` gives
`Rsup(phi*) <= Rsup(phi)`.

So every member of `Phi_abar` is matched by an invariant member with worst
regret no larger, and the two infima coincide.  ∎

> **Corollary.**  By Lemma 1 an invariant test is a function of
> `(T_1..T_K, lambda_hat_1..lambda_hat_K)`.  So the linear program of
> `w1_minimax.py`, which minimises worst-case regret over every rule measurable
> in that pair, computes the best worst-case regret achievable by **any** test
> in `Phi_abar`, up to its discretisation.

---

## 7. The conditions, discharged

**(C1) The problem is invariant.**  `Theta_0` and `A` are `G`-invariant, shown
in Section 2, and `beta*` is a `G`-invariant function on `Theta` because it
depends on `theta` only through the maximal invariant of Lemma 1's parameter
version.  *Discharged here.*

**(C2) `G` is amenable.**  Lemma 2.  *Discharged here.*

**(C3) Domination and measurability.**  The family is dominated by Lebesgue
measure on `R^{sum n_i}`, the action `(g, x) -> g x` is jointly continuous hence
jointly measurable, and `G` is a locally compact second countable group.
*Discharged here, routine.*

**(C4) Almost invariance implies invariance.**  The averaged test `phi*` is
almost invariant, and the step to an invariant test is the standard lemma for
groups that are locally compact, σ-compact and act measurably with a σ-finite
invariant measure on the quotient.  The affine group satisfies these.

*This is the one condition taken from the literature rather than verified from
first principles here*, and it is the place a reader should look first.  The
references are Bondar and Milnes (1981), which surveys exactly which groups
satisfy the Hunt-Stein conditions and places the amenable ones on the right
side, and Lehmann and Romano's treatment of invariance and minimaxity.  Neither
is currently in the paper's bibliography.

---

## 8. What this does not claim

* **Not weighted average power.**  Section 4 above says why.  The number
  Section 10.8 currently reports is a different and also useful quantity, and it
  is bounded only over invariant tests.  `w1_minimax.py` computes the quantity
  the proposition covers, and it is a different number.
* **Not a claim about randomised versus deterministic.**  The class
  `Phi_abar` already contains randomised tests and the linear program optimises
  over them; that the optimum turns out to be deterministic is a separate
  measured fact, not part of this argument.
* **Not a statement outside the declared domain.**  `Lambda` enters through
  `Theta_0` and `A`; outside it nothing here applies, as elsewhere in the
  programme.
* **Not an existence claim for a minimax test.**  The proposition equates two
  infima.  Whether either is attained is not needed and is not shown.

---

## 8a. What the covered criterion gives

`w1_minimax.py` computes the same program in epigraph form, minimising the worst
regret over the declared alternatives instead of maximising their weighted
average.  By the Proposition its value is the best worst-case regret achievable
by any test in `Phi_abar`, not only by an invariant one.

| Benchmark | minimax-regret floor | the deployed rule | gap |
|---|---|---|---|
| 3 df | 30.72 pp | 40.02 pp | +9.30 pp |
| 5 df | 17.75 pp | 20.51 pp | +2.76 pp |
| 7 df | 12.14 pp | 13.18 pp | +1.04 pp |
| 11 df | 7.22 pp | 7.43 pp | +0.21 pp |
| 19 df | 6.68 pp | 6.96 pp | +0.28 pp |
| 35 df | 6.80 pp | 7.23 pp | +0.43 pp |

Above eleven degrees of freedom the deployed rule sits **0.21 to 0.43
percentage points** above the floor that no test can beat.

Three things are worth reading off this table rather than the one in Section
10.8.

* **The floor itself is higher.**  At 35 degrees of freedom the worst-case
  regret floor is 6.80 points against a weighted-average loss of about four,
  because the worst declared alternative is at an extreme ratio and the average
  hides it.  The two numbers price the same rule on different bases and are not
  to be read across.
* **The gap is wider and it is the one the theorem covers.**  0.21 to 0.43
  points against 0.21 to 0.24 on the average basis.  Reporting the narrower
  number under a minimax theorem would be exactly the conflation this document
  exists to prevent.
* **Monotonicity is free here too.**  The unrestricted and monotone programs
  agree to 0.000 points at every design, independently of the result on the
  average basis, so the finding is not an artefact of one objective.


---

## 9. The falsifier

A theorem cannot be falsified by simulation; its conclusion can.

The linear-program design described next was ABANDONED after three
failures.  What was built instead is in "Where it stands" below,
and it checks the hypotheses directly and tests the conclusion with
rules written down rather than optimised.  The design is kept here
because the reasoning about what makes a probe non-vacuous carried
over to the instrument that replaced it.

`w1_noninvariant.py` hands a linear program a coordinate that an invariant rule
throws away, one for each generator of the group.  The sufficient statistic is
in bijection with `(T, lambda_hat, Xbar_0, s_0^2)`, so `s_0^2` is the scale
coordinate and `Xbar_0` the location one.  Each probe lets the rule read its
coordinate and then requires it to hold its size across the nuisance that
coordinate is sensitive to, which under the null is a free parameter.

If the proposition is being misapplied, the program should find a rule that uses
the coordinate, meets the size constraint and beats the bound.  The comparison
is against the **same program** with equality across that axis imposed as a
linear constraint, so no rendering difference sits between the two arms.

**The first version of this probe was wrong, and how it was wrong is the
point.**  It maximised weighted average power at ONE value of the nuisance and
reported gains of up to 13.9 percentage points, with the sanity arm returning
exactly the same numbers as the real probe.  That identity was the tell: a
constraint over a range of the nuisance was doing nothing, because the
objective sat at a single point of the orbit.  Those gains are not evidence
against the Proposition.  They are a demonstration of Section 4: at one point
of an orbit the discarded coordinate really does carry information, and average
power at a chosen representative is not a criterion the invariance argument
covers.  The probe now takes its objective over the orbit as well as its
constraint, which is the statement the Proposition makes.

Three things make the probe non-vacuous, and the first draft had none of them:

* `Xbar_0` is **not** independent of `T`, since it sits inside `T`'s numerator.
  Treating them as independent would have guaranteed the probe found nothing.
  The joint law is built from `A = Xbar_0 - mu` and `B = Xbar_1 - mu - delta`,
  independent given the state, so that conditionally on the state and on `A` the
  event `T <= t` is a single normal tail.
* the **objective is a supremum over the orbit**, not a value at one point of
  it, because that is what the Proposition is about;
* a **sanity arm** keeps a single nuisance value, where the orbit supremum
  degenerates to a point and the coordinate genuinely is informative, and there
  the free arm must beat the flat one.  Without it, finding no gain in the real
  probes would be indistinguishable from a mesh that cannot express a
  non-invariant rule at all.

### Where it stands

**The instrument is built and the conclusion survives it.**  Not by the route
described above: the linear program was abandoned after three failures, and what
replaced it is cheaper and sharper.  The history is kept at the end of this
section because the failures are informative.

**Why the hypotheses are the right target.**  A theorem's conclusion does not
need numerical testing; what needs testing is whether it has been pointed at the
right problem.  Three of its hypotheses are specific to this setting rather than
textbook, and all three are checkable.

| what | how it is checked | result |
| --- | --- | --- |
| H1, the maximal invariant is `(T, lambda_hat)` | closed in exact algebra by `sympy`, at `K = 1` and `K = 2` | 4 and 6 identities, all identically zero |
| H2, the null and the alternative are invariant | random group elements applied to random parameters | invariant drifts by `1.3e-13`, both sets preserved |
| H3, the loss is invariant | `beta*` computed from the RAW parameter, then the group applied | moved by `5.4e-15` over 48 draws |

H1 is now algebra rather than a tolerance, which is as far as a proof checker
would usefully go here: it is the only part of the argument specific to this
problem, and Lemma 2 and the averaging step are textbook.

H3 is the load-bearing one, since it is what lets the invariant mean pass
through the regret and it is exactly the step that fails for weighted average
power.  A first version of that check was vacuous: drawing the mean freely put
`beta*` at 0 or 1 on almost every draw, and a quantity pinned at an endpoint is
unmoved by anything, so it would have passed whether or not the hypothesis held.
The shift is now drawn to place `beta*` across `0.058` to `0.927`, with all 48
draws interior.

**The conclusion, tested without a mesh.**  Rather than optimise, the probes
write rules down: each reads `s_0^2` alongside `(T, lambda_hat)`, is calibrated
against the uniform constraint over `sigma_0^2` by the same exact quadrature the
paper uses, and is scored on the same instrument as the deployed rule.  Since
nothing is discretised, none of the three failure modes below can apply.

| family | rules | shape | best result |
| --- | --- | --- | --- |
| bump, step, linear, band | 24, of which 20 can be liberal somewhere | generic tilts in `log s_0^2` | `-1.33` points |
| deproj | 12, all liberal, all calibrated | uses `s_0^2` to de-bias `lambda_hat` | `-2.05` points |

Zero of 36 beat the floor.  The `deproj` family is the directed probe and the
rule a sceptic would actually propose: since `s_0^2 = sigma_0^2 u_0` and
`lambda_hat = lambda u_1 / u_0`, a rule that guesses `sigma_0^2` recovers an
estimate of `u_0` and multiplies the ratio back by it.  At its blend parameter
`gamma = 0` it IS the deployed rule, so that arm is a sanity check that has to
reproduce the baseline, and it does.  Losses grow monotonically in `gamma`, and
every calibrating shift is positive: buying accuracy where the guess is right is
paid for by a uniformly more conservative threshold, and the trade is bad at
every setting.

**The one real threat, measured.**  The size constraint runs over a continuum
and the probes sweep it finitely.  Two rules differing only in their guess
`v0` are the same rule at relocated `sigma_0^2`, so their scores must agree once
the sweep contains both, which makes the gap between them a direct read on
convergence.  On the original sweep they disagreed by `0.391` points, the same
order as the gap the Proposition concerns.  Widening the sweep from `e^6` to
`e^18`:

| `gamma` | gap at `e^6` | gap at `e^18` | monotone | converged |
| --- | --- | --- | --- | --- |
| 0.25 | `0.391` | `0.098` | yes | no |
| 1.00 | `0.174` | `0.004` | yes | yes |

At `gamma = 0.25` the gap is still above the `0.02` threshold set in advance, so
the deposit records `sup_has_converged: false` and no figure from that arm is
quoted.  **This does not threaten the finding, and the reason is directional.**
Truncating a supremum understates a rule's worst size, so it undercalibrates and
leaves the rule too liberal, which FLATTERS it.  Widening therefore makes every
non-invariant rule worse, and it does: at `gamma = 0.25` the measured regret
rises from `9.67` to `12.01` points across the four sweeps, so the loss against
the floor grows from `2.45` to `4.79`.  An unconverged supremum means the
competitor is still being flattered, not that the verdict could turn.

**And the baseline cannot move.**  The deployed rule reads `lambda_hat` alone,
and the law of `lambda_hat` does not involve `sigma_0^2`, so widening the sweep
cannot touch it.  Measured across all four sweeps, spanning 16 orders of
magnitude in `sigma_0^2`, it returns worst FWER `0.0249960512` and worst regret
`7.227473` points every time, with spread exactly zero.  That is what makes the
rising losses a comparison rather than a common drift, and it is the mechanism
of the Proposition visible as arithmetic: widening the sweep punishes every rule
that reads the scale and leaves the invariant rule untouched.

**What is not settled.**  Three things, so they are not mistaken for strength:

1. of 36 rows there are about **three informative shapes**.  The linear family
   duplicated the step family exactly, all three band windows duplicated each
   other, and the four bump rules were incapable of being liberal anywhere.
   Both duplications have the same cause, that a supremum over `sigma_0^2`
   absorbs a translation in `log s_0^2`, which is the invariance itself showing
   up in the harness;
2. the guess `v0` was predicted to be a non-probe for that same reason and was
   **not**.  It is a truncation artefact, and that failed prediction is how the
   sweep defect surfaced at all;
3. `gamma = 0.25` has not converged, so its regrets are directional evidence
   and not quotable values.

**The three failed attempts**, kept because each one failed differently:

1. the first design scored average power at one point of the orbit and found
   gains of 13.9 points.  Not a refutation, a demonstration of Section 4.  The
   tell was that its sanity arm returned the same numbers as its real probe, so
   the nuisance range was doing nothing;
2. the redesign onto the covered criterion completed one probe in 602 seconds
   and then failed numerically on the flat arm, where about 17,000 equality rows
   on a 19,520-variable mesh left the solver with no solution;
3. the refinement study solved the flat arm on the two-dimensional mesh instead.
   That is exact in exact arithmetic and wrong on a grid: the two arms then
   carry different discretisation errors, and their difference came out near
   `-62` points and did not move as the nuisance grid went from 3 values to 15.
   Flatness in that grid is the signature, because the quantity being measured
   is supposed to depend on exactly it.

Both of those deposits carry `VALID: false` and say why.  The lesson they paid
for is in the design that replaced them: construct rules and evaluate them
exactly, rather than optimise two arms on meshes whose errors do not cancel.

---

## 10. What changes in the paper, if this is adopted

Nothing is being edited now.  The shape of the change would be:

* Section 10.8's three open conditions are replaced by Lemma 1, Lemma 2 and the
  Proposition, with (C4) cited rather than proved;
* the sentence "over rules measurable in the estimated ratio" becomes "over
  every test", **for the minimax-regret criterion only**;
* a second table reports the minimax-regret floor beside the weighted-average
  one, because the two criteria give different numbers and reporting one under
  the other's name would be the same error this document exists to avoid;
* Bondar and Milnes (1981) and Lehmann and Romano join the bibliography.
