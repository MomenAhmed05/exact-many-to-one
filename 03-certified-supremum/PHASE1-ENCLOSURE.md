# PHASE 1 — the enclosure and the root certificate (plan M2: sections 2.1–2.3, A.1–A.3)

tier: region_sup_certified

**Status: complete; P1 verifier round run, its findings fixed, fix-verify round fired.**
The enclosure and the bracket are `region_sup_certified`; the derivative diagnostics are
`envelope_measured`. Every number below is read from a strict-JSON deposit under `deposits/`.

Work area `03-certified-supremum/`, plan `plans/welch-dunnett-completion.md`.

---

## 1. What P1 delivered

| file | role |
|---|---|
| `enclose.py` | two-sided monotone Riemann enclosure (float64 search + Arb certificate path), adaptive partition, explicit tail terms, mass floor, **directed** endpoint rounding |
| `enclose_gl.py` | **the certificate quadrature**: Gauss-Legendre per panel with a rigorous Bernstein-ellipse remainder, kink-aware and s-aware panel layout, `c_upper_box_gl` / `c_lower_s_gl` / `bracket_point_gl`, the `settings` block |
| `test_enclose.py` | falsifier battery F1–F10 |
| `derivs.py`, `deriv_diag.py` | corrected implicit-function derivative (+ 02-bounded-rule-defect evidence), finite-difference validation, the predeclared Lipschitz test, argmax geometry |
| `rootcert.py` | the point root-certificate battery over the design suite (P2 step 10) |
| deposits | `enclose_falsifiers.json`, `enclose_selftest.json`, `enclose_gl_selftest.json`, `deriv_diag.json`, `rootcert.json` |

## 2. The enclosure, in one page

For a box `Lambda = prod_i [lam_lo_i, lam_hi_i]`, `s_lo_i = sqrt(n_i/lam_hi_i)`,
`s_hi_i = sqrt(n_i/lam_lo_i)`, `x_i(z) = z/sqrt(n_0) + c d_i`:

    Psi_min,i = Phi(s_lo,i x_i) if x_i >= 0 else Phi(s_hi,i x_i)
    Psi_max,i = Phi(s_hi,i x_i) if x_i >= 0 else Phi(s_lo,i x_i)
    G_min = prod_i Psi_min,i,   G_max = prod_i Psi_max,i

Both are continuous and increasing in `z`, and `0 <= G_min <= G <= G_max <= 1` on the box.

* **Riemann path** (`enclose.py`, float64 search + independent cross-check): for any partition the
  left sum with `G_min` lower-bounds `min_lambda F` and the right sum with `G_max` upper-bounds
  `max_lambda F`, with explicit tail terms. No derivative bound, Lipschitz constant or vertex
  property is used.
* **GL certificate path** (`enclose_gl.py`): Gauss-Legendre per panel with the remainder
  `|E_n| <= C M_rho / (rho^{2n}(1 - rho^{-2}))`, `C = (128/15) * GL_SAFETY` — i.e. **4x** Trefethen's
  `64/15`, since `GL_SAFETY = 2` (verifier finding: the earlier text said "2x"; the code was always
  the conservative one). `M_rho` is obtained **rigorously without sampling**: one `acb` ball
  evaluation on the rectangle inscribing the ellipse, each factor bounded by
  `max(|Phi(s_lo x)|, |Phi(s_hi x)|)`.

**Panel layout (`layer_edges`) is part of the certificate.** Two properties matter:

1. **Kink edges.** For a box the integrand is only *piecewise* analytic: `Psi_min`/`Psi_max` switch
   branch at `x_i = 0`, i.e. at `z*_i = -c d_i sqrt(n0)` (verifier finding: the "entire integrand"
   premise is false for boxes, so the ellipse theorem did not apply verbatim). An **even**
   `per_layer` places an edge exactly at each `z*_i`, so every panel lies wholly on one side of the
   kink and *is* analytic on it. `plant_layer=True` is the default, so no caller can omit it.
2. **S-aware subdivision.** The ellipse bound is tight only while `s * h * sinh(t)` is moderate; a
   panel wider than about `8/s` cannot be bounded tightly for *any* `t`. Panels wider than that are
   subdivided where the mass lives.

**Directions.** Upper bounds on the supremum come from a **lower** enclosure of `min_lambda F`
(`c_upper_box_gl` requires `F_lo(c) >= 1 - alpha`); lower bounds come from a **witness point**
(`c_lower_s_gl` requires `F_hi(c; {lambda}) < 1 - alpha`). An upper bound on `max_lambda F`
certifies bounds on `min c*`, not on the supremum, and is not used for the sup.

## 3. Falsifiers (`deposits/enclose_falsifiers.json`), all PASS

| falsifier | result |
|---|---|
| F1 point containment: Riemann-float, Riemann-Arb, GL-Arb vs the **independent** solver | PASS; GL width `2.5e-14` |
| F2 random-box containment, K = 1..5, df0 = 2..39, unbalanced n, **two-sided** spreads `1e-3..1e3` (small lambda included) | **24/24 contained** |
| F3 box collapse | monotone narrowing |
| F4 refinement: Riemann width vs nodes; GL width | Riemann `2.7e-03 → 5.4e-06` (200 → 100k nodes); GL `7.1e-14`; ratio `2.5e8` |
| F5 mass floor widens, never narrows | PASS |
| F6 cross-path agreement (GL-Arb vs Riemann-Arb) | all overlap |
| F7 **decisive check** `cond_noreject(c_U) >= 1-alpha` and `cond_noreject(c_L) < 1-alpha`, K = 1,2,3 | PASS; widths `6.9e-07 … 1.35e-06` |
| F8 population/conditional identity (4000 draws) | `0.975123` vs `0.975500`, z = −1.19 |
| F9 cost by K (deposited run) | K=1 `0.91 s`, K=2 `1.93`, K=3 `3.22`, K=5 `6.73`, K=8 `13.52` — all under the 60 s budget |
| F10 error-estimate consistency | `1.4e-08 -> 2.5e-14 -> 2.2e-14`, plateau from index 3 (the node-radius floor), ratio `6.3e5`; this falsifier **found defect #9** |

The composite verdict now also requires F4 (GL must beat the finest Riemann width) and F9 (every
bracket inside the runtime budget), and records `verdict_scope`; the self-test deposits carry a
`verdict` field (verifier finding: they previously had none and were written even on failure).

## 4. Certificate performance

Desk cell `n = [13,13,13]`, `lambda = (1, 25)` (point bracket):

| quantity | value | predeclared target |
|---|---|---|
| enclosure width (probability) | `7.1e-14` (the NODE_RAD/WEIGHT_RAD floor, provable) | — |
| point bracket `[c_L, c_U]` | `[2.44757627533, 2.44757696800]`, width `6.93e-07` | `<= 1e-4` ✓ (144x inside) |
| `Delta p` | `7.1e-14` | `<= 1e-6` ✓ |
| c-equivalent of `Delta p` | `4.36e-13` | `<< 1e-4` ✓ |
| time per bracket | `1.93 s` (K=2) ... `13.52 s` (K=8) | `<= 60 s` ✓ |

**Root-certificate battery (`deposits/rootcert.json`, P2 step 10):** 17 cells (12 desk + K = 1, 3, 5,
extreme spread, quiet benchmark), **91 brackets**, each verified against the refined independent
solver:

| quantity | value | target |
|---|---|---|
| independent check (`c_L` below, `c_U` at/above the crossing) | **91/91 pass** | all |
| bracket width | max `1.60e-06`, median `8.79e-07` | `<= 1e-4` ✓ |
| `Delta p` | max `8.03e-11`, median `2.47e-14` | `<= 1e-6` ✓ |

**Structural box slack (verifier finding N-A; now deposited as `deposits/box_slack.json`).** For a *box*, `int min_lambda G dPhi <= min_lambda int G dPhi` is strict because the pointwise minimiser moves with `z`. In **c-space** (the quantity the branch-and-bound must drive below `tol`), measured on fixed draws at fixed `c`:

| box width factor | 3.0 | 2.0 | 1.5 | 1.1 | 1.01 | 1.001 |
|---|---|---|---|---|---|---|
| slack (c-units) | 2.88e-02 | 1.91e-02 | 1.56e-02 | 4.81e-03 | 5.33e-04 | 5.34e-05 |

(The probability-space analogue is necessarily negative - `F_lo(box) <= min over corners F` by construction - so c-space is the meaningful unit; an earlier version of `box_slack.py` measured it in probability space and reported negative values.) **Consequence for P2**: reaching `tol = 1e-4` needs boxes about 0.2% wide, which the branch-and-bound cannot afford on wide regions - see `PHASE2-CERTIFIED-SUP.md` for the measured gaps and the honest boundary.

## 5. Derivative diagnostics (`deposits/deriv_diag.json`)

1. **Implicit-function derivative vs finite differences — PASS.** Max relative deviation **5e-06**
   over the 12 desk cells (`A df3` analytic `[0.019488, 0.019581]` vs FD identical; `B df3`
   `[0.077852, 0.920035]`). *Correction (verifier finding): the earlier version of this report quoted
   `[0.587625, 0.121009]` for `B df3`; those are the **02-bounded-rule-defect synthetic draw** values from
   `deriv_diag.w26_derivative_defect.fixed_values`, not the desk cell.*
2. **Predeclared Lipschitz materiality test — FAIL, recorded.** Minimum covering-radius correction
   `1.340e-01` against the predeclared `TOL = 1e-04`; `feasible_and_within_tol = 0`. The Lipschitz
   route of plan sections 2.1–2.3 is therefore not a certificate (nor a useful cross-check at this
   tolerance) and is retained only as a split-direction diagnostic. The certificate path is the
   branch-and-bound of plan section 2.4, which does not need a Lipschitz constant.
3. **Argmax geometry — 12/12 vertices.** The population critical value attains its grid maximum at a
   vertex in all 12 desk cells (`n_vertex = 12`, `n_mixed = 0`), seeding P2's witness set.

## 6. Defects found and fixed (P1)

| # | defect | found by | fix and effect |
|---|---|---|---|
| 9 | the float64 rounding term used `M_rho` from the **largest** ellipse (up to `1e15`), inflating a rigorous `1e-43` bound to `2.06` and breaking monotonicity in GL order | falsifier F10 | the term uses the panel's own magnitude; width `5.3e-08 → 2.5e-14`, monotone in order |
| 10 | the first bracket falsifier compared a *conditional* bracket with the *population* oracle | F7 as first written | F7 now checks `cond_noreject` at `c_L`/`c_U`; F8 checks the draw-average identity |
| 11 | `02-bounded-rule/cert_sup.dcdloglam` paired a 201-node grid with 200 weights, so `np.dot` always raised and a bare `except` returned `[]`: **every row of `02-bounded-rule/deposits/cert/rigorous.json` has an empty `dcdloglam_at_oracle` column (0 of 12)**, and the Stage-1 derivative diagnostic was silently empty. No published number is affected | `derivs.w26_version_defect()` | `derivs.py` implements it correctly on the probability-scale nodes (`w26_empty: true`, `fixed_nonempty: true`); `02-bounded-rule/` untouched |
| 12 | `deriv_diag.c_cond` compared the *rejection* tail against `1 - alpha` (never true) → NaN → silently empty derivative columns while reporting "FAIL with zero deviation" | its own output | bisection condition corrected; a failed row records an explicit `error` field instead of an empty list |
| 13 | **`core.cond_noreject` at its default 40 panels is inaccurate at sharp configurations** (`s ~ 109`, lambda ~ 1e-3): `0.9754993407` vs the kink-split truth `0.9755000050`, an error of `−6.6e-07` | P1 verifier (coverage gap) + the rootcert failures | the falsifiers and `rootcert.py` use `npanels = REF_PANELS = 8000`, which agrees with the independent kink-split truth to `4e-14` |
| 14 | F2 only spread lambda **upward** (`lambda >= 0.65`), so the large-`s` regime was untested | P1 verifier | F2 now uses two-sided spreads `1e-3..1e3`; 24/24 still contained |
| 15 | at tiny lambda the integrand is a step of width `1/s ~ 0.009` while the finest base panel is `0.14`, so the ellipse bound exploded (nine brackets with width up to 25 and `Delta p` up to 4) | rootcert battery | `layer_edges` now subdivides panels wider than `8/s` where the mass lives; the sharp corner's width fell `4.97e-03 → 2.4e-14` (216 edges, 109 ms/call) |
| 16 | **undirected float endpoints** (MAJOR): `_f(x.lower())` is round-to-nearest, so a mass used as a lower bound exceeded the true mass in 24/400 cells and the final endpoint conversion leaked ~`5.5e-17`; `up_tail_up = 0` while the true upper tail is `1.13e-19` | P1 verifier | `_lo`/`_hi` do directed rounding (`np.nextafter`) for masses, tails and final endpoints; the tail term is computed in Arb. Effect: no observed containment failure before or after, but the rigour leak is closed |
| 17 | the rounding term was called "capped" when only `min(best, h)` applied, and the floor rationale was wrong | P1 verifier | the rounding term now carries its **proof** in the docstring; the cap is the panel's *normal mass* `Phi(b) - Phi(a)` (tighter than `h`, and the reason the far-tail panels stopped contributing); the width floor is `2 * err_total`, the rounding term, **not** a float64 floor |
| 18 | dead `FALLBACK[0]` counter; the mass-cap `except` swallowed silently | P1 verifier | counter removed; the cap failure is recorded in `ERRLOG` |

## 7. Verifier rounds touching this phase

* **P0 rounds** (spec + code, four passes; reports in `verify/`): 7 blocking defects found and fixed,
  then 2 more (six gate false negatives from a clause-level exemption; the analytic suite still
  fail-open) plus a stale manifest — all fixed and re-verified in pass 5, which closed N1/N2/N3/N4
  with 0 false negatives and 0 false positives over a 60-sentence battery plus 19 new multi-clause
  probes.
* **P1 round** (`verify/P1-VERIFY-CODE.md`; oracle report summarised here): independent containment
  **30 configurations / 342 corners, 0 violations**, enclosure tracking the corner bracket to 8–10
  decimals; the GL error bound verified point by point; 5 MAJOR and several MINOR defects, all fixed
  above (13–18); the number audit found 1 mis-attribution and 1 timing drift, both corrected here.
* **P1 fix-verify round**: fired after these fixes.

## 8. Open items carried into P2

1. `cert_bnb.py`: the branch-and-bound with `c_U = max` over ALL boxes (SPEC §7 step 7), certificate
   trees, witness seeding from `deriv_diag`, derivative-guided splitting, and the structural box
   slack of section 4 as the driving quantity.
2. `error_budget.py` (P3 step 16) assembles the per-certificate decomposition SPEC §9 specifies.
3. The `plant_layer` layout makes sharp configurations ~4x more expensive per evaluation; the B&B
   will need the search/final split of layouts noted in section 2.

## 9. Reproduction

    python enclose.py         # Riemann self-test        -> enclose_selftest.json
    python enclose_gl.py      # GL certificate self-test -> enclose_gl_selftest.json
    python test_enclose.py    # falsifiers F1-F10        -> enclose_falsifiers.json
    python deriv_diag.py      # FD check, Lipschitz, argmax -> deriv_diag.json
    python rootcert.py        # point root certificates  -> rootcert.json
    python run_modes.py quick # gates: strict JSON, tiers, traceability, reserved words

---

## 10. P1 fix-verify round (pass 2) — findings and resolutions

A second, targeted verifier round re-checked the pass-1 fixes. It **confirmed the load-bearing
claim with the strongest evidence in the programme so far** and found two MAJOR defects, both fixed.

**Containment — CLOSED.** 20 fresh configurations (λ down to 1e-4, `s` up to 475), point boxes and
boxes of width factor 1.01–2: 19 apparent violations, all ≤ 6e-14, while the verifier's own float64
reference self-converges only to 2.4e-13 — i.e. reference noise, not enclosure failures. The decisive
check: at the hardest configuration (λ = 1e-3, `s ≈ 114`, `c` = 2.0 and 3.0), an independent
**40-digit mpmath** evaluation gives `0.9884211170787645`, and the enclosure is
`[0.9884211170787522, 0.9884211170787766]` — the truth lies **inside with ~1.2e-14 margins on both
sides**.

**Defect #13 confirmed, with the right magnitude.** Against that 40-digit truth: the default-40-panel
`cond_noreject` errs by **−5.05e-07** (c = 2.0) and **−5.28e-08** (c = 3.0); `npanels = 2000` and
`8000` agree with the truth to ~6e-16, and the GL enclosure contains it. The refined reference is the
correct one, and the falsifiers use it.

| # | defect | severity | resolution |
|---|---|---|---|
| 19 | **kink edges were never planted**: `per_layer + (1 - per_layer % 2)` forces *odd* (10 → 11), but the planted offsets `(2k/per_layer - 1)` hit zero only for *even* — measured `min |edge - z*|` = half a layer spacing, so the box integrand was **not** analytic per panel and the ellipse remainder was formally unjustified for box evaluations | MAJOR | `+ (per_layer % 2)`; verified: in-range kinks now have an edge exactly at `z*` (exact hits 1/1 for in-range, the other lying outside `[-9, 9]` by construction) |
| 20 | **the float64 node/weight rounding term was unfounded**: the `2*2*eps*M` "proof" cannot yield an eps-scaled term (no Lipschitz bound exists), and the verifier measured a **12.9× breach** of the panel bound on a kink-free sharp-corner panel (returned 6.78e-21 vs actual 8.74e-20) | MAJOR | the term is **eliminated**, not patched: every node and weight now carries a tiny Arb radius (`NODE_RAD`/`WEIGHT_RAD`), so the interval arithmetic propagates the perturbation rigorously. The width now has a provable floor (`7.1e-14`), and F10's monotonicity check was made plateau-aware |
| 21 | F4/F9 were claimed in `verdict_scope` but missing from the composite `ok` (my earlier patch silently failed to apply) | MEDIUM | both conditions are now in the composite verdict (verified by reading the live function) |
| 22 | `test_enclose.py`'s main guard was `contextlib.suppress(Exception)` — any failure exited 0 with no RESULT | MEDIUM | explicit failure path: prints `BATTERY ERROR` and exits 2, or exits 1 on a FAIL verdict |
| 23 | `bounds_arb` applied the `dropped_up` adjustment *after* the directed conversion with plain `+=`/`-=` (leak ~6.6e-18) | MINOR | the adjustment is folded into the Arb sums *before* the conversion; verified finite and mass-floor-monotone |
| 24 | the mass cap converted with round-to-nearest `_f` (cap < true mass in 1505/3000 random panels) | MINOR | converted with `_hi` (directed) |
| 25 | the on-disk self-test deposits still carried `verdict: None` | MINOR | both regenerated with the fixes (verdicts present) |
| 26 | F10 demanded strict monotone decay, which the new provable floor breaks | LOW | plateau-aware: large shrinkage factor required, increases beyond a 1e-14 plateau tolerance rejected, plateau index reported |
| 27 | F1/F2 use float64 references (~1e-13) that are coarser than the enclosure width (~2e-14), so they cannot detect a 1e-14-level failure | LOW | recorded; the decisive checks use the refined reference, and a high-precision (mpmath) reference belongs to P3 (`enclose_iv.py`, step 17) |

**Also corrected here** (verifier number audit): the section-4 figures are now read from the current
deposits (width `7.1e-14`, `Delta p` `7.1e-14`, c-equivalent `4.36e-13`, bracket
`[2.44757627533, 2.44757696800]`, timings `0.91/1.93/3.22/6.73/13.52 s`), and the structural box slack
— previously report-only — is now deposited as `deposits/box_slack.json`.

**Status after pass 2:** the falsifier battery is PASS (F1–F10) and both self-tests PASS with the
corrected code; the pass-2 fixes are subject to the pass-3 re-verification round.
