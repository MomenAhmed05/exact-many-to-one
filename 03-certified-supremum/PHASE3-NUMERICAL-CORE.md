# PHASE 3 — the interval-certified numerical core (plan M4: section 3, B)

tier: region_sup_certified

**Status: complete; P3 verifier round fired after the P2 fix-verify pass.** Every number is read from a
strict-JSON deposit under `deposits/`.

---

## 1. What P3 delivered

| file | role |
|---|---|
| `error_budget.py` | the per-certificate error budget with **measured** components (SPEC §9) |
| `enclose_iv.py` | the independent high-precision cross-check (the third implementation class) |
| `impl_compare.py` | the analytic limits and cross-implementation comparisons, with widths |
| deposits | `error_budget.json`, `enclose_iv.json`, `impl_compare.json` |

## 2. The error budget — `deposits/error_budget.json`

For four certificate points (two desk cells, a sharp corner, a K = 5 design) the enclosure width
`Delta p = F_hi - F_lo` is decomposed into components that are **measured**, not asserted — each by
differencing runs with the relevant setting overridden:

| component | how it is measured | value |
|---|---|---|
| analyticity (Bernstein-ellipse) remainder | `diag["err_total"]` | 3.9e-20 … 8.3e-18 |
| node/weight radius contribution | re-run with `NODE_RAD = WEIGHT_RAD = 0` and difference | 1.4e-15 … 1.8e-15 |
| Arb ball radii | re-run at `prec = 32` and difference (an upper bound: the width grows as precision falls) | ≤ 2.3e-08 |
| summation slack | `2 * GL_SUM_SLACK` | 2.0e-14 |
| tail masses | the explicit tail terms | 1.1e-19 each |
| bisection width | the search's own bracket | `tol`-limited (1e-6) |
| **`Delta p`** | `F_hi - F_lo` at the certificate point | **max 8.38e-14** (target 1e-6) |
| **c-equivalent** | `Delta p / slope`, local `dF/dc` bound | 1.4e-12 … 3.1e-12 |

So the numerical uncertainty is ~7 orders of magnitude below the predeclared `Delta p` target, and the
c-equivalent is ~8 orders below the `1e-4` critical-value target. The node/weight radius term is now
carried *inside* the interval arithmetic — it replaced the earlier, unfounded `2*2*eps*M` claim
(pass-3 finding #20c).

## 3. The independent cross-check — `deposits/enclose_iv.json`

Five certificate points (including the sharp `s ≈ 114` corner and a K = 5 design) are re-evaluated with
`mpmath` at 30/40/50 digits, splitting the range at every kink so each piece is smooth, and the result
is checked for containment inside the Arb enclosure:

| point | enclosure width | mpmath value | contained | margins |
|---|---|---|---|---|
| desk_B_df35 | 2.20e-14 | 0.989537451559 | **yes** | +1.10e-14 / +1.10e-14 |
| desk_A_df3 | 2.16e-14 | 0.796819153200 | **yes** | +1.09e-14 / +1.08e-14 |
| sharp_A_df3_lo | 2.20e-14 | 0.963715279918 | **yes** | +1.10e-14 / +1.10e-14 |
| K3_mixed | 2.20e-14 | 0.960373360229 | **yes** | +1.11e-14 / +1.09e-14 |
| K5_balanced | 2.20e-14 | 0.969806888799 | **yes** | +1.10e-14 / +1.10e-14 |

Convergence between 30 and 50 digits: **0.0** (the high-precision values agree exactly at float64
resolution). This is the falsifier class the plan describes: it certifies nothing, but a shared
conceptual bug would show up here immediately.

## 4. The three implementation classes and the analytic limits — `deposits/impl_compare.json`

| comparison | result |
|---|---|
| **K = 1 closed form.** For one candidate the conditional statistic is a scaled normal, so `c* = z_{1-alpha} * sqrt(lambda_1/n_1 + 1/n_0) / sqrt(S_1^2/n_1 + S_0^2/n_0)` exactly. Four configurations (including a thin benchmark and `lambda_1 = 0.04`) | the enclosure's bracket contains the exact value in **4/4**; max bracket width `7.95e-07` |
| **Population law, independent route.** The production solver `core.oracle_c` vs the draw-average identity (solve for the `c` whose average conditional no-rejection probability equals `1 - alpha`) | agrees to **6.6e-04** (bracket width 9.1e-14 over 3000 draws) |
| **Vanishing-ratio limit.** As `lambda_1 -> 0` the candidate factor becomes a step, so `c* = z_{1-alpha} / (d_1 sqrt(n_0))` exactly | contained in **2/2** configurations (this exercises the kink handling directly) |

**Honest gap:** a homoscedastic comparison against published Dunnett tables is **not** included —
`scipy` has no Dunnett routine, and the 02-bounded-rule kernel regression already validates the population solver
against an independent transformation/panel architecture (max absolute difference 0.000e+00 over 12
cells). The limit is recorded rather than approximated.

## 5. Root certification (SPEC §3.4)

Every certified value is reported as a **bracket** with its width: `c_L` from the witness lower bounds
and `c_U` from the box upper bounds, both rigorous, with the bisection width recorded. Across the 91
point certificates of P2 the widths are `6.9e-07 … 1.6e-06`, and the branch-and-bound gap — the
supremum enclosure itself — is `9.5e-07` on most cells and at most `3.2e-03` on the widest family-A
cells, where it is reported as a failure to certify at the declared tolerance rather than hidden.

## 6. Verifier round P3

Fired at the close of this phase (after the P2 fix-verify pass completes, one verifier at a time):
an independent re-derivation of the error-budget components on fresh points, a number audit of this
report against the deposits, and a code audit of `enclose_iv.py` and `impl_compare.py` (including the
question of whether the K = 1 closed form is derived correctly and whether the draw-average route is
genuinely independent of `core.oracle_c`).

## 7. Open items carried forward

1. The Arb-radii component is only bounded by precision differencing (≤ 2.3e-08); a direct
   interval-arithmetic measurement of the ball radii would tighten the report but cannot change the
   total, which is dominated by `2 * GL_SUM_SLACK = 2e-14`.
2. The kink-edge sliver (float64 `sqrt(n0)` vs Arb `sqrt(n0)`, ≤ ~9e-16) is documented in the code; it
   is far below `GL_SUM_SLACK`.
3. `mpmath` is used as a falsifier, never as a certificate: the certified path is Arb.

## 8. Reproduction

    python error_budget.py    # the per-certificate budget -> error_budget.json
    python enclose_iv.py      # the mpmath cross-check    -> enclose_iv.json
    python impl_compare.py    # limits and cross-checks   -> impl_compare.json
