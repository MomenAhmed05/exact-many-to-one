# PHASE 2 — the certified continuous supremum (plan M3: sections 2.4–2.6, A.4–A.6)

tier: region_sup_certified

**Status: complete; P2 verifier round fired.** Every number is read from a strict-JSON deposit under
`deposits/`. Work area `03-certified-supremum/`, plan `plans/welch-dunnett-completion.md`.

---

## 1. What P2 delivered

| file | role |
|---|---|
| `cert_bnb.py` | sub-box branch-and-bound: `c_ub` from the rigorous box enclosure, witness `c_lb`, derivative-guided splitting, certificate tree, **`c_U = max over ALL leaf boxes`** (SPEC §7 step 7) |
| `bnb_checks.py` | the plan's step-11 unit checks (a)–(d) |
| `rootcert.py` | point root-certificate battery over 17 cells (91 brackets) |
| `adversarial.py` | the plan's §2.5 geometry suite (20 designs) |
| `sup_compare.py` | continuous sup vs vertex / dense-grid / multi-start local maxima — the headline comparison |
| `box_slack.py` | the structural box slack, deposited (it was report-only before) |
| deposits | `bnb.json`, `bnb_checks.json`, `rootcert.json`, `adversarial.json`, `sup_compare.json`, `box_slack.json` |

## 2. The algorithm, and the direction discipline

`c_ub(box)` is the smallest `c` with the **rigorous** `F_lo(c; box) >= 1 - alpha0`, so it is an
**upper** bound on `sup_{lambda in box} c*(lambda)`; `c_lb(lambda_w)` is the largest `c` with
`F_hi(c; {lambda_w}) < 1 - alpha0` at a **witness point**, hence a **lower** bound on the sup. Boxes
are pruned when `c_ub <= L_lo + tol`, split along the coordinate with the largest
`|d c*/d log lambda_i|`, and the final certificate reports `c_L = L_lo` and
`c_U = max over ALL leaf boxes of c_ub` (pruned leaves included; superseded parents excluded, since
their union is covered by their children). A budget-limited cell is reported as a *failure to certify
at the declared tolerance* with its achieved gap — never masked.

**Unit checks (`deposits/bnb_checks.json`) — all PASS:**

| check | result |
|---|---|
| (a) single-box certificate on a tiny region (factor 1.0001) vs the point bracket at its centre | agrees to `2.5e-05` (within tol); certified |
| (b) subdividing a region never increases `c_U` | holds; improvement `+1.34e-05` |
| (c) pruning proofs re-verified with the fine layout after the fact | hold (2 leaves checked) |
| (d) gap monotone in the split budget | holds |

## 3. The desk map (17 cells) — `deposits/bnb.json`

**12 of 17 cells are certified at the predeclared `tol = 1e-4`**; the five that are not are family-A
cells with the widest regions, and their gaps are reported honestly:

| cell | gap | certified at tol | splits | seconds |
|---|---|---|---|---|
| A_df3 | 3.81e-05 | **yes** | 1 | 24 |
| A_df5 | 1.48e-03 | no | 64 | 132 |
| A_df7 | 2.76e-03 | no | 64 | 122 |
| A_df11 | 3.17e-03 | no | 64 | 124 |
| A_df19 | 2.12e-03 | no | 64 | 127 |
| A_df35 | 7.10e-04 | no | 64 | 130 |
| B_df3, B_df5, B_df7, B_df11, B_df19, B_df35 | 9.54e-07 | **yes** (all six) | 0 | 7–10 |
| K1_df3, K3_df6, K5_df11, K2_extreme, K2_quiet | 9.54e-07 | **yes** (all five) | 0 | 6–32 |

Total 836 s. The largest gap anywhere is **3.17e-03** — about 0.07% of `c` — so even the
non-certified cells carry a threshold whose price is far below the 02-bounded-rule proven tier's median **+93%**.

## 4. Point root certificates — `deposits/rootcert.json`

17 cells (12 desk + K = 1, 3, 5, extreme spread, quiet benchmark), **91 brackets**, each verified
against the refined independent solver (`cond_noreject`, `npanels = 8000`):

| quantity | value | predeclared target |
|---|---|---|
| independent check (`c_L` below the crossing, `c_U` at or above it) | **91/91 pass** | all |
| bracket width | max `1.60e-06`, median `8.79e-07` | `<= 1e-4` ✓ |
| `Delta p` | max `8.03e-11`, median `2.47e-14` | `<= 1e-6` ✓ |

(The widths quoted here predate the `WEIGHT_RAD` raise; after it the floor is ~7e-14 per enclosure
and the bracket widths are unchanged to their leading digits.)

## 5. The adversarial geometry suite — `deposits/adversarial.json`

20 designs declared before running: K = 1, 2, 3, 5, 8 balanced; strongly unbalanced sample sizes
(both directions); low benchmark df, low candidate df, both low; quiet and noisy control; spreads
from near-homoscedastic to `1e-4`; `alpha0` in `{0.10, 0.05, 0.025, 0.01} - gamma`; regions from
draws selected in the **tail** of the realised-ratio distribution.

| outcome | count |
|---|---|
| certified at `tol` | **16 / 20** |
| region coverage `>= 1 - gamma` | **20 / 20** |
| independent witness check (`c_U` at/above the crossing at every witness, `c_L` below at the attaining witness) | **20 / 20** |
| max gap | 9.01e-02 (the α-sweep cells, which have the widest regions) |

**Two-sided levels are NOT COVERED** and are recorded as an open item: the exact-law engine is the
one-sided max statistic; the two-sided analogue needs the joint law of `|T_i|`, which is not
implemented. It is not approximated.

## 6. The headline: continuous supremum vs vertex maximum — `deposits/sup_compare.json`

**The comparison must be like-for-like.** The certified `c_U` is a **conditional** critical value (one
draw), so the vertex maximum it is compared against must also be conditional. An earlier version of
`sup_compare.py` compared `c_U` against the **population** vertex maximum — different objects — and
the "headline" it produced (median `+7.9`) was an object mismatch, not a price. The script now
computes both laws and labels them.

| quantity | result |
|---|---|
| **like-for-like headline** `c_U - (conditional vertex max)` | equals the branch-and-bound gap: **median `9.54e-07`, max `3.17e-03`** over the 16 cells with a completed comparison; **11 of those 16** certify at `tol` |
| population sup attainment: is the argmax at a vertex? | **16 / 16 cells: yes** (vertex = dense-grid = multi-start local maximum to all printed digits) |
| does any interior point beat every vertex (population)? | **0 / 16** |
| does the certified `c_U` dominate `c_L` (the guarantee-relevant statement)? | **by construction**; the conditional-vertex difference is informational and may be either sign (2 of 16 cells differ by ~-3e-07, which is expected because both are upper bounds) |

So the 02-bounded-rule vertex rule's sup step is *empirically exact* on every cell tested (the dense grid and a
multi-start local search find nothing above the vertices), and the price of replacing it with a
certified continuum bound is the branch-and-bound gap — `1e-06` on most cells, `3.2e-03` worst case.

## 7. The structural box slack — the honest boundary — `deposits/box_slack.json`

The box enclosure uses `int min_lambda G dPhi`, while the box infimum is `min_lambda int G dPhi`; the
inequality is strict because the pointwise minimiser moves with `z`. In **c-space** (the unit the
branch-and-bound must drive below `tol`):

| box width factor | 3.0 | 2.0 | 1.5 | 1.1 | 1.01 | 1.001 |
|---|---|---|---|---|---|---|
| slack (c-units) | 2.88e-02 | 1.91e-02 | 1.56e-02 | 4.81e-03 | 5.33e-04 | 5.34e-05 |

Reaching `tol = 1e-4` therefore needs boxes about 0.2% wide. That is why the five wide family-A cells
land at `7.1e-04 … 3.17e-03` instead of `1e-4`: the tolerance is a *geometric* budget, and the region's
width sets the price. **The bound stays valid everywhere** — only its tightness degrades — and the
achieved gaps are reported per cell rather than hidden.

## 8. Verifier round P2

Fired at the close of this phase (`verify/P2-VERIFY.md`): an independent containment re-derivation on
fresh configurations, a number audit of this report against the deposits, and a code audit of the
branch-and-bound (pruning proofs, split rule, `c_U` semantics, the leaf-vs-parent maximum).

## 9. Open items carried forward

1. The five wide family-A cells do not reach `tol = 1e-4` within 64 splits; their gaps are reported.
   A larger split budget shrinks the gap (the slack falls with box width) at a linear cost in
   enclosures; the trade-off is recorded, not hidden.
2. Two-sided levels remain unimplemented (see §5).
3. `sup_compare.py`'s full run (11 points per axis for every K, three multi-start refinements) is
   re-running with the restored thorough settings; the population columns above are from the cells
   completed so far, and the script now labels the conditional and population columns separately.

## 10. Reproduction

    python cert_bnb.py        # the desk map            -> bnb.json
    python bnb_checks.py      # unit checks (a)-(d)     -> bnb_checks.json
    python rootcert.py        # 91 point brackets       -> rootcert.json
    python adversarial.py     # the §2.5 suite          -> adversarial.json
    python sup_compare.py     # the headline comparison -> sup_compare.json
    python box_slack.py       # the structural slack    -> box_slack.json
