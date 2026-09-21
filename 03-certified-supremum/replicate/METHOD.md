# Method summary (the specification is `../SPEC.md`; this is the shortest correct route)

1. `alpha0 = abar - gamma` with `abar = 0.025`, `gamma = 0.0005` (so `alpha0 = 0.0245`).
2. From the draw's variances: `lambda_hat_i = S_i^2 / S_0^2`.
3. Region: `lambda_lo_i = lambda_hat_i * w_lo / y_hi_i`, `lambda_hi_i = lambda_hat_i * w_hi / y_lo_i`
   with the four multipliers from separate equal-tail chi-square quantiles (`w_frac = 0.5`), at the
   largest feasible scale whose EXACT 1-D quadrature coverage is at least `1 - gamma + 1e-6`.
4. Enclosure at a box: `s_lo_i = sqrt(n_i/lambda_hi_i)`, `s_hi_i = sqrt(n_i/lambda_lo_i)`;
   `Psi_min,i = Phi(s_lo x)` for `x >= 0` else `Phi(s_hi x)`; `G_min = prod Psi_min` (and `G_max`
   with the branches swapped).  `G_min` and `G_max` are increasing in `z`, so a left Riemann sum with
   `G_min` lower-bounds `min_lambda F` and a right sum with `G_max` upper-bounds `max_lambda F`, with
   explicit tail terms.
5. Certificate quadrature: Gauss-Legendre per panel with the Bernstein-ellipse remainder
   `C * M_rho / (rho^{2n}(1 - rho^-2))`, `C = (128/15) * 2`; `M_rho` from ONE complex-ball evaluation
   on the rectangle inscribing the ellipse; each node/weight carries a small Arb radius so the float64
   perturbation is inside the interval arithmetic; panel edges include an edge at every kink
   `z*_i = -c d_i sqrt(n0)` and subdivide panels wider than `8/s_max`.
6. Certified values: `c_ub(box)` = smallest `c` with `F_lo(c; box) >= 1 - alpha0` (an UPPER bound on
   the box supremum); `c_lb(lambda_w)` = largest `c` with `F_hi(c; {lambda_w}) < 1 - alpha0` at a
   WITNESS (a LOWER bound).
7. Branch-and-bound: prune when `c_ub <= L_lo + tol`; split along the largest
   `|d c*/d log lambda_i|`; report `c_L = L_lo`, `c_U = max over ALL leaf boxes of c_ub`, and the gap.
8. Decide: promote candidate `i` iff `T_i > c_U`.
