# The procedure

What a desk runs, what it is guaranteed, and where the guarantee stops.  Every number here is in a
deposit under `deposits/` and is checked by `paper_verify_frontier.py`.

---

## The setting

One benchmark (group 0) and `K` candidate strategies.  For each candidate the question is whether
its mean return beats the benchmark's, relative to the noise in that difference.  The statistic is
the separate-variance Welch one,

    T_i = (Xbar_i - Xbar_0) / sqrt(s_i^2/n_i + s_0^2/n_0) ,   i = 1..K ,

and the discipline is familywise: at most `abar` chance that ANY candidate looks like a winner when
none is.  `abar = 0.025` one-sided throughout below.

---

## Before you see the performance data

**1. Fix the two sample sizes.**  The cost is set by the SMALLER of the two groups, not by the
benchmark.  From `deposits/wd_designmap.json`:

| power you are willing to give up | degrees of freedom needed on the smaller group |
|---|---|
| under 5 points | 19 |
| under 10 points | 11 |
| under 25 points | 5 |
| under 50 points | 3 |

On monthly data, 19 degrees of freedom is 20 observations.

**2. Declare the ratio domain.**  The guarantee holds for every variance ratio
`lambda_i = sigma_i^2 / sigma_0^2` inside a domain you declare in advance.  The shipped one is
`[0.001, 1000]` in the variance ratio, three decades in volatility, which admits a candidate about
thirty times as volatile as the benchmark or the reverse.

This declaration is the one place a desk can void the guarantee.  Admissible sources: a published
volatility target, a contractual risk limit, a construction constraint, or a volatility history from
a period not reused for selection.  **A bound inferred from the same returns you are about to test
is circular and voids the guarantee.**  Narrowing it buys very little, at most 2.1 points at eleven
degrees of freedom, nothing at 35, so the default is usually right.

**3. Declare any dependence you intend to remove.**  Serial correlation and factor exposure are
handled by transforming the returns before the statistic is formed, not by changing the rule.  Both
require the nuisance to be declared in advance; estimating it from the same returns makes the result
a plug-in and the guarantee no longer applies.  See "Dependence" below.

---

## Computing the threshold

**4. Look up the internal level `alpha0` for your design.**  This is a single number per
`(n_0, n_1, ..., K)`, computed offline by quadrature.  It does not depend on your data.

One candidate (`deposits/wb_certified_rule.json`, candidate at 12 degrees of freedom):

| benchmark | `alpha0` | certified bound on the error rate | power given up |
|---|---|---|---|
| 3 df | 0.01585 | 0.024999 | 23.9 points |
| 5 df | 0.02160 | 0.025000 | 11.2 points |
| 7 df | 0.02328 | 0.025000 | 7.3 points |
| 11 df | 0.02393 | 0.024998 | 5.9 points |
| 19 df | 0.02389 | 0.024998 | 4.4 points |
| 35 df | 0.02377 | 0.025000 | 4.5 points |
| 59 df | 0.02365 | 0.024996 | 5.4 points |

More than one candidate, per-candidate rule (`deposits/certify_anyk.json`, one-sided; halve for
two-sided):

| benchmark | K | `alpha0` | the rule's own error rate |
|---|---|---|---|
| 35 df | 2 | 0.01155 | 0.02365 |
| 35 df | 3 | 0.00753 | 0.02337 |
| 35 df | 5 | 0.00438 | 0.02312 |
| 11 df | 2 | 0.01177 | 0.02410 |
| 11 df | 3 | 0.00772 | 0.02397 |

**5. On the data, compute the estimated ratios and read off the threshold.**

    lambda_hat_i = s_i^2 / s_0^2
    c_i          = c*( lambda_hat_i ; alpha0 )        # the exact fixed-ratio constant
    reject i     iff  T_i > c_i

`c*` is the exact constant of Sections 3 and 4 of the paper: a two-dimensional outer quadrature over
the benchmark's two sufficient statistics and one one-dimensional inner rule per candidate.  It is
deterministic, carries no seed, and a desk precomputes it once per design as a table in
`lambda_hat`.

That is the whole procedure.  Three numbers in, one threshold out, no simulation anywhere in the
decision path.

---

## What is guaranteed

**At one candidate:** the familywise error rate is at most 0.025 at every ratio in the declared
domain.  Certified over the continuum by branch and bound with interval arithmetic, not by a grid
search, 12 of 17 cells at a tolerance of 1e-4 declared in advance, the other five reported with the
gap they reached (worst 3.17e-03) and still valid as upper bounds.

**At any number of candidates, one- or two-sided:** the same, by a union step that costs no
high-dimensional enclosure.  Ten certificates, all passing, with the deployed rule spending 91% to
96% of the ceiling it is allowed.

**What it costs:** 4.5 points of power at a 36-month benchmark against an oracle that knows the true
ratios.  The previously certified region rule costs 72.1 points on the same design.  Above eleven
degrees of freedom the residual is not removable by any threshold function of the estimated ratios
that a numerical search over the class could find.

---

## Dependence

**Serial correlation.**  Declare the coefficient `phi`, then transform each group before forming the
statistic:

    y_t = x_t - phi x_{t-1}   for t >= 2 ,   and DISCARD y_1 .

Dropping the first observation rather than rescaling it is what makes the transformed mean constant,
so the problem lands back in the independent family with `n` replaced by `n - 1` and every
certificate transfers.  The cost is one observation per group.  Rescaling it instead (Prais-Winsten)
leaves a mean that depends on the group's length, and with unequal lengths and a non-zero mean the
rule drifts off its ceiling, 0.03030 against 0.025 in the deposited falsifier.

If you can only declare a BOUND `phi <= phi_max`, whiten at the top of the bound: that is
conservative everywhere inside it and costs no reduction in the internal level.

**A common factor.**  Declare the loadings and residualise: `r_tilde = r - beta f`.  The hypothesis
tested is then about the factor alpha rather than the raw mean, which is a change of question and
has to be declared as one.  Declared loadings return every exposure pattern to the ceiling.
Loadings fitted on the same returns are a plug-in; fitted on a disjoint pre-period they are not.

---

## Where the guarantee stops

* **Outside the declared ratio domain there is no guarantee at all.**  This is the trade the rule
  makes: a ceiling of exactly 0.025 rather than 0.025 plus a miss probability, bought with a domain
  of validity declared in advance instead of one built from the data.
* **Below about five degrees of freedom** the cost is 21 to 40 points and the 2023 state of the art
  should be preferred.
* **Above one candidate the certified level's POWER is not reported.**  The 1.5 to 2.2 point figures
  are the rule at its measured level; the certified per-candidate level is lower and its power cost
  is not measured anywhere.
* **Normality and independence across candidates** are assumed throughout, with variance estimates
  of the usual chi-square form.
* **An estimated autocorrelation or estimated factor loadings** make the result a plug-in, not a
  certificate.
* **The blind replication kit is sealed and has not been executed by a second party.**

---

## Does this solve Dunnett's problem?

The distribution theory, yes.  The joint null law of the maximum Welch statistic against a shared
benchmark under unequal variances is computable exactly, deterministically, and at a cost linear in
the number of candidates.  That law had been treated as intractable since 1955 and the current state
of the art approximates it in four places and then simulates the result.

The testing problem, no, and not for want of computation.  With unknown unequal variances this is a
Behrens-Fisher situation.  Exactly similar tests of it are not literally non-existent, but every
known construction buys similarity with something a desk cannot use: randomisation, dependence on
the parity of the sample sizes, or a discontinuous acceptance region.  None is a threshold you can
compute and act on, which is the object this procedure is.  What exists here is the next
best thing and it is what a desk can act on: a threshold computed from data, with a familywise error
rate **bounded rather than estimated**, certified over a continuum of nuisance values that the desk
declares in advance, at a power cost of about four and a half points rather than the seventy-two the
previously certified rule was paying.
