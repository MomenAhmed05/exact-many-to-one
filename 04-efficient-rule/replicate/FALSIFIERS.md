# The four falsifiers, and what each one can refute

A check that cannot fail is not a check.  Each of these can.

## R1 — the closed form

`P(T_1 > c | S) = Phi(-c d / sigma_diff)` against a direct quadrature of
`E_z[Phi(s(z/sqrt(n_0) + c d))]` with `s = sqrt(n_1/lambda)`.

Refutes: an algebra error in the conditional law, which would make every number downstream wrong.

Run it at extreme ratios (`1e-3`, `1e3`) as well as at 1, and **refine your reference**: a residual
that does not fall when the reference is refined belongs to the identity, not to the quadrature.  In
the reference implementation the residual is 1.7e-12 at a coarse reference and 1.3e-15 at a fine one.

## R2 — the monotone substitutions

On random state cells and random ratio boxes, evaluate the integrand at interior points and check it
lies inside the corner bounds used by the enclosure.

Refutes: a wrong direction in the enclosure, which would make the "upper bound" not an upper bound.
This is the single most dangerous failure mode in the whole construction, because it produces a
certificate that certifies nothing.  Zero violations is the only passing result.

## R3 — the bound must dominate a sample

Monte Carlo the unconditional FWER directly and compare against your certified upper bound at a
spread of ratios.

Refutes: the enclosure as a whole.  A breach beyond three standard errors is a refutation, not a
discrepancy.  Note the asymmetry: the bound sitting far ABOVE the sample is expected and is not a
defect; sitting below it even once is fatal.

## R4 — the two-sided factor of two

`P(|T_1| > c | S)` must equal exactly `2 * P(T_1 > c | S)`, conditionally at every draw and hence
unconditionally.

Refutes: either of the two laws.  It is also the reason the one-sided certificate transfers to the
two-sided problem at half the level, so a replica that cannot reproduce it cannot reproduce that
claim either.
