# The four falsifiers a replica must run

1. **Enclosure containment.** For random boxes, check `F_lo <= min over corners F` and
   `max over corners F <= F_hi` against your own reference evaluation of
   `F(c, lambda) = E_U[prod_i Phi(s_i (U + c d_i))]`.  Expect zero violations and a width of order
   `1e-13` or better.
2. **Bracket against an independent route.** At a certificate point, compute `c*` by your own solver
   and check `c_L <= c* <= c_U`.  Then check the defining property directly: the conditional
   no-rejection probability at `c_U` must be at or above `1 - alpha0`, and at `c_L` below it.
3. **Draw-average identity.** Average your conditional no-rejection probability at the production
   population oracle over several thousand draws; expect `1 - alpha` within Monte-Carlo error (the
   03-certified-supremum measurement is `0.975123` against `0.9755`, `z = -1.19`).
4. **Sharp-corner stress.** At `lambda ~ 1e-3` (so `s ~ 114`, where the integrand is a step of width
   `~1/s`), check that your enclosure still contains a high-precision evaluation
   (03-certified-supremum uses `mpmath` at 30-50 digits and finds margins of `~1.1e-14`).

Report the outcome of each, with the numbers, before opening `targets/expected.json`.


## Explicit designs (added after the blind replication flagged their absence)

| falsifier | n | lambda | c | seed | notes |
|---|---|---|---|---|---|
| 1 containment (point) | [13, 13, 13] | [1, 25] | 2.0, 2.5, 2.75, 3.0 | 20260918 | point box |
| 1 containment (box) | random K in 1..5, df0 in 2..39 | random | 2.0 + 1.5 U | 20260919 | width factor 1.0 + 0.9 U |
| 2 bracket | [13, 13, 13] | [1, 25] | returned by bisection | 20260918 | K = 2 |
| 2 bracket | [7, 7, 9, 11] | [0.25, 1, 4] | returned by bisection | 20260918 | K = 3 |
| 3 draw-average | [13, 13, 13] | [1, 25] | the population oracle at alpha = 0.0245 | 20260919 | 3000-4000 draws |
| 4 sharp corner | [4, 13, 13] | the region's lower corner (~1e-3) | 5.0 | 20260918 | mpmath at 30-50 digits |

The draw recipe for every falsifier is SPEC section 14.1: `numpy.random.default_rng(seed)`,
chi-square draws in group order, `sigma_0 = 1`, and `d_i` the ESTIMATED standard error of the
realised draw.
