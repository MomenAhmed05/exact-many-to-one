# METHOD — everything a replica needs, and nothing it must not see

This kit exists because the preceding programme's blind replication failed on its draw-dependent
half, and the cause was specification under-determination rather than a mathematical disagreement.
Eight ambiguities were found there.  This specification is written to be immune to all eight: there
are **no draws anywhere in the decision path**, so there is no draw recipe to get wrong.

Read this file and `inputs/designs.json`.  Do not open `targets/expected.json` until your numbers
exist.  Do not read the programme source.

---

## 1. The model

`K + 1` independent Gaussian samples.  Group 0 is the benchmark, group `i` a candidate.  Sample
sizes `n_i`, unknown variances `sigma_i^2`, and

    lambda_i = sigma_i^2 / sigma_0^2 ,     T_i = (Xbar_i - Xbar_0) / sqrt(S_i^2/n_i + S_0^2/n_0) .

Everything is scale invariant, so set `sigma_0 = 1` throughout.  This kit covers `K = 1`.

## 2. The conditional law, in closed form

Condition on the realised sample variances `S^2`.  Under the global null the numerator is normal
with variance `sigma^2_diff = lambda/n_1 + 1/n_0`, and with `d = sqrt(S_1^2/n_1 + S_0^2/n_0)`,

    P( T_1 > c | S )      = Phi( -c d / sigma_diff ) ,
    P( T_1 > c | S, shift delta ) = Phi( (delta - c d) / sigma_diff ) .

There is no quadrature in this step and no approximation.  Derive it rather than trusting it: it is
the identity `E_z[Phi(a + b z)] = Phi(a / sqrt(1 + b^2))` applied once.

## 3. The state law

With `u_i = S_i^2 / sigma_i^2`, each `u_i ~ chi2_{n_i - 1} / (n_i - 1)`, independent.  The estimate a
desk sees is

    lambda_hat = S_1^2 / S_0^2 = lambda * u_1 / u_0 .

## 4. The rule

A step function of `lambda_hat` on the bins in `inputs/designs.json`:

    reject iff T_1 > c( lambda_hat ) .

The bin values are the exact fixed-ratio critical value at the bin's geometric midpoint, evaluated
at an internal level `alpha0`.  The exact fixed-ratio critical value at ratio `L` is the root in `c`
of

    E_{u_0, u_1} [ Phi( -c d / sigma_diff ) ] = alpha0      (with lambda = L throughout) .

`alpha0` is the largest level for which the rule satisfies section 5.  Everything is deterministic:
bisect.

## 5. What must be true

    sup over lambda in [1e-3, 1e3] of    E_{u_0, u_1}[ Phi( -c(lambda u_1/u_0) d / sigma_diff ) ]   <=   abar .

The expectation is over the two scaled chi-squares; evaluate it by any means you can defend.  The
reference implementation uses a tensor Gauss-Legendre rule on the probability scale with the panel
and order counts in `inputs/designs.json`, but the quantity is defined independently of how it is
computed, and a replica that reaches it another way is the stronger check.

## 6. The certificate

The supremum in section 5 is over a continuum, so a grid does not establish it.  The reference route:

1. Partition each state axis on the PROBABILITY scale, so each cell's mass is an exact difference of
   chi-square distribution functions and the masses are known without quadrature.
2. On a cell, and over a ratio box `[L_lo, L_hi]`, bound the integrand above by substituting
   **`c` at its minimum** over the cell's reachable `lambda_hat` range, **`d` at its minimum**
   (`L_lo, u_0 lo, u_1 lo`), and **`sigma_diff` at its maximum** (`L_hi`).  Each substitution is
   monotone in the right direction; prove that before using it.
3. The mass left outside the partition must be ADDED to the upper bound, because the integrand is a
   probability and can be 1 there.  With the span in `inputs/designs.json` it is 4e-9.
4. Branch and bound over the ratio interval: accept a box when its upper bound is at or below
   `abar`, split it otherwise, and report failure if the budget runs out.

## 7. Deliverables

Produce, for each design in `inputs/designs.json`:

* `alpha0`;
* the certified upper bound on the supremum in section 5;
* whether every box was decided;
* the power at the effect where the **known-lambda** oracle reaches 0.8, at each alternative ratio.
  The known-lambda oracle rejects when `Xbar_1 - Xbar_0 > z_{1-abar} sigma_diff`, so its power is
  `Phi(delta/sigma_diff - z_{1-abar})` exactly and the effect solves that for 0.8.

## 8. Falsifiers you must run before comparing

* **R1** the closed form of section 2 against a direct quadrature of
  `E_z[Phi(s(z/sqrt(n_0) + c d))]`, `s = sqrt(n_1/lambda)`.  They must agree to float64.  Refine the
  reference and confirm the residual falls; if it does not, the error is yours and not the identity's.
* **R2** the monotone substitutions of section 6.2, on random cells and random ratio boxes, against
  direct evaluation at interior points.  Zero violations, or the bound is not a bound.
* **R3** your certified bound must dominate a direct Monte Carlo of the same quantity at every ratio
  you test.  A breach beyond three standard errors refutes the construction.
* **R4** two-sided: `P(|T_1| > c | S)` must equal exactly twice `P(T_1 > c | S)`.  If it does not,
  one of your two laws is wrong.

## 9. Then, and only then

Open `targets/expected.json`.  It carries the canonical values and a SHA-256 for every deposit as of
the freeze.  Report agreements and disagreements both; a disagreement is information, not a failure.
