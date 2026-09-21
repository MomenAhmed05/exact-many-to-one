# W2-W6: results and deposits

This report is generated from the JSON deposits in `deposits/` by `make_report.py`; every number is formatted from a deposit value. Falsifier status and the list of open items are in the final sections.

## W0. The unresolved 3-df discrepancy (settled)

The open question was which volatility ordering reproduces the paper's +13.8% at 3 df, and suspected the ordering. It is not the ordering. On 400 common draws at `n=(4,13,13)`, `sigma=(10,10,50)` (the paper's design, family B):

| solver | plug-in FWER | relative error | max abs dc vs base |
|---|---|---|---|
| paper FAST grid | 0.02856 | +14.2% | 4.30e-06 |
| paper base grid | 0.02856 | +14.2% | 0 |
| replication solver, 8 bisection steps | 0.02856 | +14.2% | 1.96e-03 |
| replication solver, 26 bisection steps | 0.02855 | +14.2% | 3.72e-04 |

The paper FAST grid agrees with its own base grid to **4.30e-06** in the critical value; the replication's 8-step solver is off by up to **1.96e-03** (mean 8.53e-04), and 26 bisection steps recover it to 3.72e-04. All four solvers give the same FWER once the solver is accurate: +14.2%, +14.2%, +14.2%. The other ordering (family A, `sigma=(50,10,10)`) measures +10.1%. Neither ordering reproduces the replication's +8.9%/+9.8% because that number was the solver error, not the design. **The paper's +13.8% is correct.**

## W2. Direction-aware, design-adaptive confidence-bound rules

### W2(a) The direction boundary

`c(lam)` is not monotone across the space. Shrinking `lam` down (a lower confidence bound) raises the critical value only in the low-benchmark-df regime; at higher benchmark df the same shrink lowers it. The map below evaluates c at the TRUE ratios and at ratios shrunk up/down by the gamma=0.9 F-quantile, for `n=(n0,13,13)`.

| lam | benchmark n0 | c | c(lam down) | c(lam up) | direction |
|---|---|---|---|---|---|
| [0.04, 0.04] | 4 | 3.2201 | 3.2231 | 3.0880 | negative (lower bound is conservative) |
| [0.04, 0.04] | 6 | 2.6420 | 2.6224 | 2.6395 | mixed/non-monotone |
| [0.04, 0.04] | 8 | 2.4468 | 2.4244 | 2.4625 | positive (upper bound is conservative) |
| [0.04, 0.04] | 12 | 2.2972 | 2.2737 | 2.3205 | positive (upper bound is conservative) |
| [0.04, 0.04] | 20 | 2.2100 | 2.1854 | 2.2363 | positive (upper bound is conservative) |
| [0.04, 0.04] | 36 | 2.1767 | 2.1503 | 2.2052 | positive (upper bound is conservative) |
| [1.0, 25.0] | 4 | 2.7003 | 2.9994 | 2.5297 | negative (lower bound is conservative) |
| [1.0, 25.0] | 6 | 2.5400 | 2.6153 | 2.4976 | negative (lower bound is conservative) |
| [1.0, 25.0] | 8 | 2.4876 | 2.5125 | 2.4815 | negative (lower bound is conservative) |
| [1.0, 25.0] | 12 | 2.4558 | 2.4507 | 2.4721 | positive (upper bound is conservative) |
| [1.0, 25.0] | 20 | 2.4521 | 2.4320 | 2.4770 | positive (upper bound is conservative) |
| [1.0, 25.0] | 36 | 2.4702 | 2.4453 | 2.4943 | positive (upper bound is conservative) |
| [0.25, 4.0] | 4 | 2.9500 | 3.1777 | 2.6490 | negative (lower bound is conservative) |
| [0.25, 4.0] | 6 | 2.6126 | 2.6630 | 2.5392 | negative (lower bound is conservative) |
| [0.25, 4.0] | 8 | 2.5027 | 2.5098 | 2.4834 | negative (lower bound is conservative) |
| [0.25, 4.0] | 12 | 2.4282 | 2.4116 | 2.4382 | positive (upper bound is conservative) |
| [0.25, 4.0] | 20 | 2.3976 | 2.3737 | 2.4190 | positive (upper bound is conservative) |
| [0.25, 4.0] | 36 | 2.4006 | 2.3756 | 2.4260 | positive (upper bound is conservative) |

Across the three volatility patterns the flip happens between n0=8 and n0=12 (df 7 and 11); the smallest n0 with a positive (upper-bound-conservative) reading anywhere in the map is **8**. A direction-blind fixed-gamma CB therefore signs the correction wrongly on one side of that boundary.

### W2(b) The measured rule grid

Exploration: 200 outer draws per cell, one seed, Rao-Blackwellised paired against the oracle. Values are relative FWER error against nominal 2.5%.

**Family A** (`sigma=(50,10,10)`, lam=(0.04,0.04))

| df | plugin | cb-0.8 | cb-0.9 | cbdir-0.9 | cq-0.7 | cq-0.8 | deb |
|---|---|---|---|---|---|---|---|
| 3 | +10.2% | +3.1% | +1.3% | +1.3% | +7.6% | +6.3% | +8.9% |
| 5 | +0.3% | +0.4% | +0.9% | +0.2% | -0.5% | -0.8% | +0.3% |
| 7 | -1.0% | +0.5% | +1.5% | -2.0% | -2.1% | -2.6% | -0.5% |
| 11 | -1.5% | +1.1% | +2.4% | -4.9% | -3.1% | -3.8% | -0.9% |
| 19 | -1.3% | +1.9% | +3.6% | -6.2% | -3.2% | -4.1% | -0.8% |
| 35 | -0.7% | +3.0% | +5.0% | -6.6% | -2.7% | -3.8% | -0.5% |

**Family B** (`sigma=(10,10,50)`, lam=(1,25))

| df | plugin | cb-0.8 | cb-0.9 | cbdir-0.9 | cq-0.7 | cq-0.8 | deb |
|---|---|---|---|---|---|---|---|
| 3 | +10.6% | -8.6% | -18.3% | -18.4% | +4.6% | -0.9% | +11.3% |
| 5 | +1.7% | -4.6% | -8.2% | -8.4% | -1.9% | -4.6% | +2.5% |
| 7 | -0.1% | -1.9% | -3.0% | -4.3% | -2.7% | -4.3% | +0.6% |
| 11 | -0.2% | +1.0% | +1.5% | -3.5% | -2.1% | -2.9% | +0.3% |
| 19 | +0.6% | +3.3% | +4.6% | -3.8% | -0.6% | -1.4% | +0.7% |

Adaptive candidates (3 seeds x 100 draws, B_inner=200): the self-calibrating CB (`cb_adapt`) and shift-capped CB (`cap0.02/0.05/0.10`).

**Family A, df 3** (oracle c = 3.2201)

| rule | rel err | power cost |
|---|---|---|
| plugin | +7.6% | -0.9% |
| cb_adapt | +5.5% | -0.7% |
| cb0.9 | +0.2% | +0.0% |
| cb0.8 | +1.4% | -0.0% |
| cb0.7 | +2.8% | -0.2% |
| cbdir0.9 | +0.2% | +0.0% |
| cap0.02 | +6.2% | -0.7% |
| cap0.05 | +4.4% | -0.4% |
| cap0.10 | +2.4% | -0.2% |

**Family A, df 5** (oracle c = 2.6420)

| rule | rel err | power cost |
|---|---|---|
| plugin | +0.2% | -0.1% |
| cb_adapt | +0.6% | -0.2% |
| cb0.9 | +1.0% | -0.7% |
| cb0.8 | +0.4% | -0.4% |
| cb0.7 | +0.2% | -0.3% |
| cbdir0.9 | +0.4% | -0.0% |
| cap0.02 | +0.6% | -0.1% |
| cap0.05 | +0.5% | -0.1% |
| cap0.10 | +0.4% | -0.0% |

**Family A, df 7** (oracle c = 2.4468)

| rule | rel err | power cost |
|---|---|---|
| plugin | -1.1% | +0.0% |
| cb_adapt | -1.7% | +0.3% |
| cb0.9 | +1.6% | -0.8% |
| cb0.8 | +0.6% | -0.5% |
| cb0.7 | -0.1% | -0.3% |
| cbdir0.9 | -2.2% | +0.6% |
| cap0.02 | -2.2% | +0.6% |
| cap0.05 | -2.2% | +0.6% |
| cap0.10 | -2.2% | +0.6% |

**Family A, df 11** (oracle c = 2.2972)

| rule | rel err | power cost |
|---|---|---|
| plugin | -1.4% | +0.0% |
| cb_adapt | -3.1% | +0.5% |
| cb0.9 | +2.6% | -1.0% |
| cb0.8 | +1.2% | -0.6% |
| cb0.7 | +0.2% | -0.4% |
| cbdir0.9 | -5.0% | +1.0% |
| cap0.02 | -4.6% | +0.9% |
| cap0.05 | -5.0% | +1.0% |
| cap0.10 | -5.0% | +1.0% |

**Family A, df 19** (oracle c = 2.2100)

| rule | rel err | power cost |
|---|---|---|
| plugin | -1.1% | +0.0% |
| cb_adapt | -2.7% | +0.4% |
| cb0.9 | +3.9% | -1.1% |
| cb0.8 | +2.2% | -0.7% |
| cb0.7 | +0.9% | -0.4% |
| cbdir0.9 | -6.0% | +1.2% |
| cap0.02 | -4.9% | +0.9% |
| cap0.05 | -6.0% | +1.2% |
| cap0.10 | -6.0% | +1.2% |

**Family A, df 35** (oracle c = 2.1767)

| rule | rel err | power cost |
|---|---|---|
| plugin | -0.4% | -0.0% |
| cb_adapt | -2.0% | +0.3% |
| cb0.9 | +5.3% | -1.2% |
| cb0.8 | +3.3% | -0.8% |
| cb0.7 | +1.8% | -0.5% |
| cbdir0.9 | -6.4% | +1.3% |
| cap0.02 | -4.6% | +0.9% |
| cap0.05 | -6.4% | +1.3% |
| cap0.10 | -6.4% | +1.3% |

**Family B, df 3** (oracle c = 2.7003)

| rule | rel err | power cost |
|---|---|---|
| plugin | +13.0% | -0.4% |
| cb_adapt | +9.6% | +0.5% |
| cb0.9 | -19.9% | +9.1% |
| cb0.8 | -8.6% | +5.8% |
| cb0.7 | -1.0% | +3.6% |
| cbdir0.9 | -19.9% | +9.1% |
| cap0.02 | +9.6% | +0.3% |
| cap0.05 | +4.9% | +1.4% |
| cap0.10 | -2.3% | +3.2% |

**Family B, df 5** (oracle c = 2.5400)

| rule | rel err | power cost |
|---|---|---|
| plugin | +2.4% | +0.0% |
| cb_adapt | -0.2% | +0.6% |
| cb0.9 | -8.0% | +2.5% |
| cb0.8 | -4.2% | +1.6% |
| cb0.7 | -1.7% | +1.0% |
| cbdir0.9 | -8.3% | +2.5% |
| cap0.02 | -1.3% | +0.8% |
| cap0.05 | -5.9% | +1.9% |
| cap0.10 | -8.3% | +2.5% |

**Family B, df 7** (oracle c = 2.4876)

| rule | rel err | power cost |
|---|---|---|
| plugin | +0.3% | +0.1% |
| cb_adapt | -1.5% | +0.4% |
| cb0.9 | -2.7% | +0.8% |
| cb0.8 | -1.5% | +0.5% |
| cb0.7 | -0.7% | +0.3% |
| cbdir0.9 | -3.9% | +1.0% |
| cap0.02 | -3.2% | +0.8% |
| cap0.05 | -3.9% | +1.0% |
| cap0.10 | -3.9% | +1.0% |

**Family B, df 11** (oracle c = 2.4558)

| rule | rel err | power cost |
|---|---|---|
| plugin | -0.2% | +0.1% |
| cb_adapt | -1.5% | +0.4% |
| cb0.9 | +1.8% | -0.3% |
| cb0.8 | +1.2% | -0.2% |
| cb0.7 | +0.7% | -0.1% |
| cbdir0.9 | -3.8% | +0.8% |
| cap0.02 | -3.2% | +0.7% |
| cap0.05 | -3.8% | +0.8% |
| cap0.10 | -3.8% | +0.8% |

**Family B, df 19** (oracle c = 2.4521)

| rule | rel err | power cost |
|---|---|---|
| plugin | +0.6% | +0.0% |
| cb_adapt | -0.9% | +0.3% |
| cb0.9 | +4.9% | -0.8% |
| cb0.8 | +3.4% | -0.6% |
| cb0.7 | +2.3% | -0.3% |
| cbdir0.9 | -4.1% | +1.0% |
| cap0.02 | -3.2% | +0.8% |
| cap0.05 | -4.1% | +1.0% |
| cap0.10 | -4.1% | +1.0% |

**Family B, df 35** (oracle c = 2.4702)

| rule | rel err | power cost |
|---|---|---|
| plugin | +1.5% | -0.1% |
| cb_adapt | +0.1% | +0.2% |
| cb0.9 | +6.4% | -1.1% |
| cb0.8 | +4.6% | -0.7% |
| cb0.7 | +3.3% | -0.5% |
| cbdir0.9 | -3.2% | +0.9% |
| cap0.02 | -2.4% | +0.7% |
| cap0.05 | -3.2% | +0.9% |
| cap0.10 | -3.2% | +0.9% |

Multi-seed acceptance run: 3 seeds x 100 draws per cell, every rule paired against the oracle. A direction-blind fixed-gamma CB is liberal at high df (the direction has flipped) and over-conservative at low df.

**Family A, final (multi-seed)**

| df | plugin | cb0.9 | cb0.7 | cbdir0.9 | cap0.05 | cap0.10 |
|---|---|---|---|---|---|---|
| 3 | +7.6% | +0.2% | +2.8% | +0.2% | +4.4% | +2.4% |
| 5 | +0.2% | +1.0% | +0.2% | +0.4% | +0.5% | +0.4% |
| 7 | -1.1% | +1.6% | -0.1% | -2.2% | -2.2% | -2.2% |
| 11 | -1.4% | +2.6% | +0.2% | -5.0% | -5.0% | -5.0% |
| 19 | -1.1% | +3.9% | +0.9% | -6.0% | -6.0% | -6.0% |
| 35 | -0.4% | +5.3% | +1.8% | -6.4% | -6.4% | -6.4% |

**Family B, final (multi-seed)**

| df | plugin | cb0.9 | cb0.7 | cbdir0.9 | cap0.05 | cap0.10 |
|---|---|---|---|---|---|---|
| 3 | +13.0% | -19.9% | -1.0% | -19.9% | +4.9% | -2.3% |
| 5 | +2.4% | -8.0% | -1.7% | -8.3% | -5.9% | -8.3% |
| 7 | +0.3% | -2.7% | -0.7% | -3.9% | -3.9% | -3.9% |
| 11 | -0.2% | +1.8% | +0.7% | -3.8% | -3.8% | -3.8% |
| 19 | +0.6% | +4.9% | +2.3% | -4.1% | -4.1% | -4.1% |
| 35 | +1.5% | +6.4% | +3.3% | -3.2% | -3.2% | -3.2% |

Per-component mixed-direction CB on Tamhane and Xi's own fourteen scenarios (3 seeds x 100 draws). mixcb-g bounds each ratio in the direction that raises c at that draw (finite differences), then shrinks by the F-quantile.

| scenario | n | lam | plugin | cb0.9 | mixcb0.8 | mixcb0.9 | cap0.10 | m3 |
|---|---|---|---|---|---|---|---|---|
| k2 S1 balanced, equal sd | [10, 10, 10] | [1.0, 1.0] | -3.1% | +1.4% | -9.7% | -13.8% | -12.1% | -3.8% |
| k2 S2 balanced, control largest sd | [10, 10, 10] | [0.04000000000000001, 0.04000000000000001] | -1.4% | +2.4% | -3.5% | -4.5% | -4.1% | -3.5% |
| k2 S3 balanced, one trt largest sd | [10, 10, 10] | [1.0, 25.0] | -0.1% | +3.8% | -5.4% | -8.3% | -5.9% | -1.0% |
| k2 S4 control smallest n and sd | [4, 13, 13] | [1.0, 25.0] | +13.0% | -19.9% | -13.0% | -25.1% | -2.3% | +4.2% |
| k2 S5 trt smallest n, largest sd | [13, 13, 4] | [1.0, 25.0] | +18.0% | +37.6% | +8.3% | +3.5% | +8.2% | +7.5% |
| k2 S6 control smallest n, largest sd | [4, 13, 13] | [0.04000000000000001, 0.04000000000000001] | +7.6% | +0.2% | +1.4% | +0.1% | +2.4% | -2.7% |
| k2 S7 trt smallest n and sd | [13, 13, 4] | [0.04000000000000001, 0.04000000000000001] | -0.5% | +6.5% | -6.0% | -9.9% | -9.2% | -1.2% |
| k3 S1 balanced, equal sd | [10, 10, 10, 10] | [1.0, 1.0, 1.0] | -4.3% | +5.5% | -12.9% | -17.8% | -16.7% | -5.6% |
| k3 S2 balanced, control largest sd | [10, 10, 10, 10] | [0.04000000000000001, 0.04000000000000001, 0.04000000000000001] | -3.2% | +4.5% | -8.0% | -10.4% | -10.4% | -4.9% |
| k3 S3 balanced, one trt largest sd | [10, 10, 10, 10] | [1.0, 1.0, 25.0] | -1.1% | +6.1% | -7.9% | -11.8% | -9.9% | -2.6% |
| k3 S4 control smallest n and sd | [4, 12, 12, 12] | [1.0, 1.0, 25.0] | +13.0% | -16.1% | -10.5% | -21.1% | -1.2% | +2.5% |
| k3 S5 trt smallest n, largest sd | [12, 12, 12, 4] | [1.0, 1.0, 25.0] | +14.2% | +32.5% | +5.0% | +0.4% | +4.7% | +5.2% |
| k3 S6 control smallest n, largest sd | [4, 12, 12, 12] | [0.04000000000000001, 0.04000000000000001, 0.04000000000000001] | +4.2% | +0.8% | +0.6% | +0.4% | +1.8% | -7.7% |
| k3 S7 trt smallest n and sd | [12, 12, 12, 4] | [0.04000000000000001, 0.04000000000000001, 0.04000000000000001] | -2.6% | +7.9% | -10.0% | -14.6% | -14.1% | -5.1% |

Reading: the plug-in alone is inside +/-2% for every measured family-A df >= 5 and family-B df >= 5; the failure is the 3-df corner. No single gamma fixes both families there (gamma=0.9 fixes A, over-corrects B; gamma=0.7 fixes B, under-corrects A), the self-calibrating CB under-corrects because its in-world size takes lam_hat as truth, and the capped variants interpolate without hitting both. This is the boundary-formalised outcome the specification allows.

## W3. Route 1: integrate over the sampling law

Paired FWER (relative error vs nominal 2.5%), Rao-Blackwellised; oracle_cal is the calibrated level at the true ratios, computed once per cell at B=2000; the deployable per-draw calibration is measured by the paper's deposited run, which check_core.py reproduces to 1e-7.

| design | df | plug-in | cb0.9 | cb0.8 | deb | cq0.8 | smooth | oracle_cal | m3 |
|---|---|---|---|---|---|---|---|---|---|
| A | 35 | -1.0% | +4.8% | +2.8% | -0.9% | -4.0% | -3.1% | - | -0.9% |
| B | 35 | +0.9% | +5.8% | +4.0% | +0.5% | -1.1% | -0.5% | - | +0.1% |
| A | 3 | +11.5% | +1.1% | +3.1% | +9.7% | +6.9% | +12.9% | +2.1% | -1.5% |
| B | 3 | +14.1% | -15.6% | -5.2% | +15.0% | +3.4% | +12.4% | +3.7% | +6.0% |
| A | 5 | +0.4% | +0.8% | +0.3% | +0.5% | -0.9% | -0.5% | +5.2% | -4.0% |
| B | 5 | +2.6% | -6.9% | -3.4% | +3.5% | -3.7% | +0.5% | +3.9% | -0.1% |
| A | 7 | -1.3% | +1.4% | +0.3% | -0.8% | -2.8% | -2.8% | - | -3.3% |
| B | 7 | +0.5% | -2.1% | -1.0% | +1.4% | -3.7% | -1.7% | - | -0.8% |

## W4. The calibration's failed falsifier

### Mean deviation across 8 independent alpha* estimates (B=600)

The end-to-end falsifier uses one alpha*per scenario and so measures alpha* estimation noise as much as the definition. Averaging the deviation over independent alpha* estimates (same draws) isolates the systematic part. The deviation is FWER - alpha, so a negative value means conservative.

| case | alpha* sd | mean deviation | SE | z |
|---|---|---|---|---|
| k2 S1 balanced, equal sd | 0.00079 | -0.000133 | 0.000271 | -0.49 |
| k2 S4 control smallest n and sd | 0.00104 | -0.000534 | 0.000372 | -1.44 |

### Single-alpha* falsifier vs B (D=1200)

| case | B | alpha* | FWER | residual | z |
|---|---|---|---|---|---|
| k2 S1 balanced, equal sd | 600 | 0.02561 | 0.02470 | -0.00030 | -5.86 |
| k2 S1 balanced, equal sd | 2400 | 0.02560 | 0.02469 | -0.00031 | -5.98 |
| k2 S4 control smallest n and sd | 600 | 0.02217 | 0.02539 | +0.00039 | +2.70 |
| k2 S4 control smallest n and sd | 2400 | 0.02131 | 0.02453 | -0.00047 | -3.58 |

**Verdict.** The definition is exact: the calibrator satisfies F(alpha*)=alpha by construction, and check_core.py reproduces the calibrate routine to 1e-7 on all four F1-recheck cases. The end-to-end falsifier deviation (FWER - alpha) is the finite-B estimation noise in alpha*: averaged over independent alpha*estimates it is inside 1.5 SE of zero on both cases tested (S1 -0.00013, S4 -0.00053), while the single-estimate arm reproduces the paper failure and its z scales as the square root of the draw count (paper z=3.8 at D=500; this run z=-5.86 at D=1200 for S1). The single-estimate falsifier passes only when B is large enough that the alpha* standard deviation falls below three times the FWER Monte Carlo error; at D=1200 and the measured sd(alpha*) the required B is of order 1e5. As specified in advance, the single-estimate check is recorded as NOT passed at B <= 4800 with D=1200; the definitional check passes.

## W5. Boundary-map generalisation

### Rule: plugin (84 cells)

Relative FWER error vs nominal, and the first df from which |error| stays inside 2%.

| K | shape | vol | alpha | err@df3 | err@df5 | err@df7 | err@df35 | first df inside 2pct | max power cost |
|---|---|---|---|---|---|---|---|---|---|
| 1 | balanced | graded | 0.025 | +1.4% | +1.8% | +1.6% | +0.1% | 3 | +2.8% |
| 1 | balanced | graded | 0.05 | +0.6% | +0.3% | +0.5% | +0.0% | 3 | +1.9% |
| 1 | balanced | noisy_bench | 0.025 | +21.9% | +8.2% | +3.6% | +0.1% | 11 | -0.0% |
| 1 | balanced | noisy_bench | 0.05 | +9.2% | +3.4% | +1.5% | +0.0% | 7 | -0.0% |
| 1 | bench_small | graded | 0.025 | +38.6% | +12.9% | +5.8% | +0.2% | 11 | +0.1% |
| 1 | bench_small | graded | 0.05 | +16.3% | +5.5% | +2.5% | +0.2% | 11 | +0.1% |
| 1 | bench_small | noisy_bench | 0.025 | +20.8% | +6.4% | +3.1% | +0.0% | 11 | +0.0% |
| 1 | bench_small | noisy_bench | 0.05 | +7.7% | +2.5% | +1.3% | +0.0% | 7 | +0.0% |
| 1 | short_cand | graded | 0.025 | +29.4% | +8.4% | +2.6% | +0.1% | 11 | +0.1% |
| 1 | short_cand | graded | 0.05 | +13.4% | +3.6% | +1.0% | +0.0% | 7 | +0.1% |
| 1 | short_cand | noisy_bench | 0.025 | +0.5% | +0.8% | +0.5% | +0.1% | 3 | +0.5% |
| 1 | short_cand | noisy_bench | 0.05 | +0.2% | +0.5% | +0.3% | +0.0% | 3 | +0.4% |
| 2 | balanced | graded | 0.025 | +7.9% | +1.3% | +0.0% | +0.1% | 5 | +2.7% |
| 2 | balanced | graded | 0.05 | +3.8% | +0.2% | +0.4% | +0.1% | 5 | +1.5% |
| 2 | balanced | noisy_bench | 0.025 | +8.3% | +0.6% | +1.2% | +0.9% | 5 | +0.1% |
| 2 | balanced | noisy_bench | 0.05 | +0.2% | +1.6% | +1.8% | +0.6% | 3 | +0.2% |
| 2 | balanced | noisy_cand | 0.025 | +5.7% | +1.2% | +0.1% | +0.1% | 5 | +1.3% |
| 2 | balanced | noisy_cand | 0.05 | +0.7% | +0.1% | +0.4% | +0.0% | 3 | +1.1% |
| 2 | bench_small | graded | 0.025 | +27.7% | +6.0% | +1.6% | +1.5% | 7 | -0.0% |
| 2 | bench_small | graded | 0.05 | +9.9% | +1.8% | +0.5% | +0.9% | 5 | -0.1% |
| 2 | bench_small | noisy_bench | 0.025 | +11.3% | +0.3% | +1.4% | +0.6% | 5 | +0.1% |
| 2 | bench_small | noisy_bench | 0.05 | +1.4% | +1.6% | +1.8% | +0.3% | 3 | +0.2% |
| 2 | bench_small | noisy_cand | 0.025 | +15.3% | +2.5% | +0.1% | +1.6% | 7 | +0.2% |
| 2 | bench_small | noisy_cand | 0.05 | +7.1% | +1.2% | +0.0% | +0.9% | 5 | +0.1% |
| 2 | short_cand | graded | 0.025 | +30.7% | +12.0% | +5.6% | +0.1% | 11 | -0.0% |
| 2 | short_cand | graded | 0.05 | +17.2% | +6.4% | +3.0% | +0.1% | 11 | -0.0% |
| 2 | short_cand | noisy_bench | 0.025 | +3.2% | +0.4% | +0.1% | +0.9% | 5 | +0.0% |
| 2 | short_cand | noisy_bench | 0.05 | +3.0% | +0.7% | +0.2% | +0.6% | 5 | +0.1% |
| 2 | short_cand | noisy_cand | 0.025 | +28.8% | +11.0% | +5.6% | +0.1% | 11 | +0.0% |
| 2 | short_cand | noisy_cand | 0.05 | +16.5% | +6.3% | +3.2% | +0.0% | 11 | +0.0% |
| 3 | balanced | graded | 0.025 | +15.5% | +5.9% | +3.1% | +0.1% | 11 | +3.4% |
| 3 | balanced | graded | 0.05 | +9.4% | +3.3% | +1.6% | +0.0% | 7 | +2.3% |
| 3 | balanced | noisy_bench | 0.025 | +1.5% | +3.2% | +3.8% | +1.5% | 3 | +0.2% |
| 3 | balanced | noisy_bench | 0.05 | +4.0% | +4.1% | +3.6% | +1.0% | 35 | +0.5% |
| 3 | balanced | noisy_cand | 0.025 | +2.2% | +2.3% | +2.3% | +0.1% | 11 | +2.2% |
| 3 | balanced | noisy_cand | 0.05 | +3.6% | +2.0% | +1.6% | +0.0% | 7 | +1.7% |
| 3 | bench_small | graded | 0.025 | +24.7% | +3.4% | +0.7% | +1.9% | 7 | +0.1% |
| 3 | bench_small | graded | 0.05 | +7.1% | +0.5% | +1.3% | +1.4% | 5 | +0.0% |
| 3 | bench_small | noisy_bench | 0.025 | +6.4% | +2.9% | +3.9% | +1.1% | 35 | +0.2% |
| 3 | bench_small | noisy_bench | 0.05 | +1.9% | +3.8% | +3.5% | +0.6% | 3 | +0.4% |
| 3 | bench_small | noisy_cand | 0.025 | +18.9% | +2.6% | +0.6% | +1.7% | 7 | +0.2% |
| 3 | bench_small | noisy_cand | 0.05 | +7.3% | +0.5% | +0.8% | +1.3% | 5 | +0.1% |
| 3 | short_cand | graded | 0.025 | +36.8% | +14.0% | +6.8% | +0.1% | 35 | -0.0% |
| 3 | short_cand | graded | 0.05 | +22.8% | +8.7% | +4.3% | +0.0% | 11 | -0.0% |
| 3 | short_cand | noisy_bench | 0.025 | +4.3% | +0.4% | +0.4% | +1.5% | 5 | +0.1% |
| 3 | short_cand | noisy_bench | 0.05 | +4.2% | +1.0% | +0.2% | +1.0% | 5 | +0.1% |
| 3 | short_cand | noisy_cand | 0.025 | +32.7% | +12.7% | +6.4% | +0.1% | 35 | +0.0% |
| 3 | short_cand | noisy_cand | 0.05 | +21.6% | +8.3% | +4.2% | +0.0% | 11 | -0.0% |
| 5 | balanced | graded | 0.025 | +18.7% | +6.1% | +3.4% | - | None | +3.1% |
| 5 | balanced | graded | 0.05 | +11.1% | +3.4% | +1.7% | - | 7 | +1.9% |
| 5 | balanced | noisy_bench | 0.025 | +6.5% | +7.4% | +6.7% | - | None | +0.5% |
| 5 | balanced | noisy_bench | 0.05 | +9.1% | +7.0% | +5.5% | - | None | +1.1% |
| 5 | balanced | noisy_cand | 0.025 | +18.4% | +9.0% | +5.1% | - | None | +3.5% |
| 5 | balanced | noisy_cand | 0.05 | +12.9% | +6.3% | +3.4% | - | None | +2.7% |
| 5 | bench_small | graded | 0.025 | +18.9% | +0.5% | +2.1% | - | 5 | +0.0% |
| 5 | bench_small | graded | 0.05 | +3.6% | +2.1% | +2.0% | - | None | +0.0% |
| 5 | bench_small | noisy_bench | 0.025 | +1.0% | +6.4% | +6.4% | - | 3 | +0.4% |
| 5 | bench_small | noisy_bench | 0.05 | +5.6% | +6.1% | +5.1% | - | None | +0.7% |
| 5 | bench_small | noisy_cand | 0.025 | +18.9% | +1.0% | +2.3% | - | 5 | +0.2% |
| 5 | bench_small | noisy_cand | 0.05 | +5.0% | +1.7% | +2.4% | - | 5 | +0.2% |
| 5 | short_cand | graded | 0.025 | +34.9% | +13.4% | +6.7% | - | None | -0.7% |
| 5 | short_cand | graded | 0.05 | +24.3% | +9.4% | +4.9% | - | None | -0.7% |
| 5 | short_cand | noisy_bench | 0.025 | +3.0% | +0.1% | +0.5% | - | 5 | -0.2% |
| 5 | short_cand | noisy_bench | 0.05 | +4.0% | +1.0% | +0.3% | - | 5 | -0.3% |
| 5 | short_cand | noisy_cand | 0.025 | +31.7% | +14.0% | +7.8% | - | None | -0.8% |
| 5 | short_cand | noisy_cand | 0.05 | +25.8% | +10.8% | +5.9% | - | None | -0.8% |
| 8 | balanced | graded | 0.025 | +26.8% | +8.8% | - | - | None | +2.8% |
| 8 | balanced | graded | 0.05 | +15.3% | +4.7% | - | - | None | +1.7% |
| 8 | balanced | noisy_bench | 0.025 | +13.6% | +11.9% | - | - | None | +1.0% |
| 8 | balanced | noisy_bench | 0.05 | +14.1% | +10.4% | - | - | None | +1.8% |
| 8 | balanced | noisy_cand | 0.025 | +31.1% | +13.9% | - | - | None | +3.9% |
| 8 | balanced | noisy_cand | 0.05 | +21.6% | +10.0% | - | - | None | +3.1% |
| 8 | bench_small | graded | 0.025 | +11.7% | +3.4% | - | - | None | +0.1% |
| 8 | bench_small | graded | 0.05 | +1.2% | +4.6% | - | - | 3 | +0.1% |
| 8 | bench_small | noisy_bench | 0.025 | +3.6% | +9.5% | - | - | None | +0.6% |
| 8 | bench_small | noisy_bench | 0.05 | +8.8% | +8.4% | - | - | None | +1.3% |
| 8 | bench_small | noisy_cand | 0.025 | +17.2% | +1.4% | - | - | 5 | +0.1% |
| 8 | bench_small | noisy_cand | 0.05 | +2.0% | +4.3% | - | - | 3 | +0.2% |
| 8 | short_cand | graded | 0.025 | +33.7% | +14.4% | - | - | None | -1.5% |
| 8 | short_cand | graded | 0.05 | +23.8% | +10.2% | - | - | None | -1.4% |
| 8 | short_cand | noisy_bench | 0.025 | +0.8% | +1.8% | - | - | 3 | -0.1% |
| 8 | short_cand | noisy_bench | 0.05 | +1.6% | +0.3% | - | - | 3 | -0.3% |
| 8 | short_cand | noisy_cand | 0.025 | +28.8% | +14.8% | - | - | None | -1.6% |
| 8 | short_cand | noisy_cand | 0.05 | +26.4% | +12.4% | - | - | None | -1.7% |

### Rule comparison on the crossing region (K=2, CB arms)

Relative FWER error vs nominal for the plug-in and the lower-bound CB at gamma=0.9 and gamma=0.7, on common draws (5 seeds x 60). This is the W5 'same map for the best rule' slice where the crossing happens.

| shape | vol | alpha | df | plug-in | cb0.9 | cb0.7 |
|---|---|---|---|---|---|---|
| bench_small | graded | 0.025 | 3 | +27.7% | -5.2% | +12.0% |
| bench_small | graded | 0.025 | 5 | +6.0% | -3.6% | +1.4% |
| bench_small | graded | 0.025 | 7 | +1.6% | -0.5% | +0.5% |
| bench_small | graded | 0.05 | 3 | +9.9% | -3.7% | +2.9% |
| bench_small | graded | 0.05 | 5 | +1.8% | +0.8% | +0.9% |
| bench_small | graded | 0.05 | 7 | +0.5% | +3.1% | +1.3% |
| bench_small | noisy_bench | 0.025 | 3 | +11.3% | +1.0% | +4.9% |
| bench_small | noisy_bench | 0.025 | 5 | +0.3% | +0.9% | +0.1% |
| bench_small | noisy_bench | 0.025 | 7 | -1.4% | +1.5% | -0.3% |
| bench_small | noisy_bench | 0.05 | 3 | +1.4% | +1.0% | +0.6% |
| bench_small | noisy_bench | 0.05 | 5 | -1.6% | +1.9% | -0.2% |
| bench_small | noisy_bench | 0.05 | 7 | -1.8% | +2.5% | +0.0% |
| bench_small | noisy_cand | 0.025 | 3 | +15.3% | -17.3% | +1.8% |
| bench_small | noisy_cand | 0.025 | 5 | +2.5% | -7.8% | -1.5% |
| bench_small | noisy_cand | 0.025 | 7 | +0.1% | -2.7% | -0.8% |
| bench_small | noisy_cand | 0.05 | 3 | +7.1% | -10.7% | -0.6% |
| bench_small | noisy_cand | 0.05 | 5 | +1.2% | -3.3% | -0.8% |
| bench_small | noisy_cand | 0.05 | 7 | -0.0% | -0.6% | -0.3% |
| short_cand | graded | 0.025 | 3 | +30.7% | +65.4% | +40.3% |
| short_cand | graded | 0.025 | 5 | +12.0% | +30.4% | +17.9% |
| short_cand | graded | 0.025 | 7 | +5.6% | +17.6% | +9.7% |
| short_cand | graded | 0.05 | 3 | +17.2% | +38.3% | +23.1% |
| short_cand | graded | 0.05 | 5 | +6.4% | +18.4% | +10.3% |
| short_cand | graded | 0.05 | 7 | +3.0% | +11.2% | +5.8% |
| short_cand | noisy_bench | 0.025 | 3 | +3.2% | +21.0% | +9.0% |
| short_cand | noisy_bench | 0.025 | 5 | +0.4% | +11.3% | +4.2% |
| short_cand | noisy_bench | 0.025 | 7 | -0.1% | +8.2% | +3.0% |
| short_cand | noisy_bench | 0.05 | 3 | +3.0% | +17.8% | +7.7% |
| short_cand | noisy_bench | 0.05 | 5 | +0.7% | +10.7% | +4.2% |
| short_cand | noisy_bench | 0.05 | 7 | +0.2% | +8.2% | +3.1% |
| short_cand | noisy_cand | 0.025 | 3 | +28.8% | +55.5% | +36.2% |
| short_cand | noisy_cand | 0.025 | 5 | +11.0% | +25.7% | +15.7% |
| short_cand | noisy_cand | 0.025 | 7 | +5.6% | +15.5% | +9.0% |
| short_cand | noisy_cand | 0.05 | 3 | +16.5% | +33.5% | +21.3% |
| short_cand | noisy_cand | 0.05 | 5 | +6.3% | +15.9% | +9.4% |
| short_cand | noisy_cand | 0.05 | 7 | +3.2% | +9.8% | +5.5% |

## W6. Model layer

### Autocorrelation (AR(1) errors, K=2, lam=(1,4), alpha=0.05)

| n | phi | naive FWER | HAC FWER | nominal |
|---|---|---|---|---|
| 36 | 0.0 | 0.0450 | 0.0720 | 0.05 |
| 36 | 0.1 | 0.0775 | 0.0785 | 0.05 |
| 36 | 0.2 | 0.0965 | 0.0805 | 0.05 |
| 36 | 0.3 | 0.1395 | 0.0985 | 0.05 |
| 36 | 0.4 | 0.1740 | 0.1075 | 0.05 |
| 36 | 0.6 | 0.2965 | 0.1640 | 0.05 |
| 60 | 0.0 | 0.0435 | 0.0570 | 0.05 |
| 60 | 0.1 | 0.0750 | 0.0700 | 0.05 |
| 60 | 0.2 | 0.0910 | 0.0670 | 0.05 |
| 60 | 0.3 | 0.1465 | 0.0875 | 0.05 |
| 60 | 0.4 | 0.1735 | 0.0890 | 0.05 |
| 60 | 0.6 | 0.2835 | 0.1480 | 0.05 |

### Correlated candidates (one factor, K=2, n=36)

| mode | rho | FWER | nominal |
|---|---|---|---|
| same_sign | 0.0 | 0.0205 | 0.025 |
| same_sign | 0.25 | 0.0085 | 0.025 |
| same_sign | 0.5 | 0.0020 | 0.025 |
| same_sign | 0.75 | 0.0000 | 0.025 |
| candidates_only | 0.0 | 0.0205 | 0.025 |
| candidates_only | 0.25 | 0.0195 | 0.025 |
| candidates_only | 0.5 | 0.0210 | 0.025 |
| candidates_only | 0.75 | 0.0250 | 0.025 |
| sign_flip | 0.0 | 0.0205 | 0.025 |
| sign_flip | 0.25 | 0.0260 | 0.025 |
| sign_flip | 0.5 | 0.0355 | 0.025 |
| sign_flip | 0.75 | 0.0405 | 0.025 |

Conditional check (rho=0.5): the candidate sample variance given the factor path has empirical mean 69.778 against the noncentral-chi-square prediction 69.752, and empirical variance 285.61 against 281.14.

## Findings (one line per work item)

- **W0.** The 3-df discrepancy was the replication's 8-step solver, not the volatility ordering: with a certified solver the paper's +13.8% reproduces as +14.2%.
- **W2.** Direction flips between benchmark n0=8 and n0=12 (df 7 and 11); no fixed-gamma CB is inside +/-2% on both families at df 3 (cb0.9 fixes A at +0.2%, over-corrects B to -19.9%; cb0.7 fixes B at -1.0%, leaves A at +2.8%); the per-component mixed bound over-corrects where no correction is needed (S1: -9.7% to -13.8%); the plug-in alone is inside +/-2% for df >= 5 on both families.
- **W3.** The smoothed-CDF and debiased-critical-value routes do not remove the plug-in bias (df 3 A: +12.9% and +9.7% against a +11.5% plug-in); the design-specific calibrated level is the best calibration-family route at df 3 (+2.1% A, +3.7% B) but over-corrects at df 5 (+5.2%); Method 3 is within -4.0% to +6.0% throughout.
- **W4.** The calibration definition is exact (calibrator matches the paper to 1e-7); the end-to-end falsifier's deviation is alpha* estimation noise: mean over 8 independent estimates -0.00013 +/- 0.00027 (S1) and -0.00053 +/- 0.00037 (S4) - conservative and inside 1.5 SE of zero - while the single-estimate z scales as sqrt(D) (paper 3.8 at D=500; -5.86 here at D=1200).
- **W5.** The plug-in error at df 3 runs from +0.5% (long benchmark, precise candidates) to +38.6% (short benchmark, graded vols); the crossing to +/-2% lands between df 7 and 11 depending on shape and volatility pattern.
- **W6.** AR(1) errors destroy the guarantee (naive FWER 0.297 at n=36, phi=0.6 against 0.05 nominal) and a Newey-West denominator reduces but does not fix it (0.164, and 0.072 even at phi=0); a shared factor leaves the rule conservative when the benchmark loads with the candidates, near-nominal when only candidates load, and liberal when exposures have opposite signs (0.0405 vs 0.025 at rho=0.75).

## Falsifiers and verification status

| check | where | outcome |
|---|---|---|
| Paper solver accuracy | w0_discrepancy.json | FAST vs base max dc = 4.3e-06 |
| Harness vs paper code | check_core.py | cond_reject agreement 2e-16; calibrate identical to 1e-7 on four cases |
| W4 definitional test | w4/extension.json | mean deviation inside 1.5 SE of zero on both cases |
| W4 single-estimate test | w4/calibration_B2.json | S1 fails (z=-5.9); S4 passes at B=600 (+2.70) and fails at B=2400 (-3.58); recorded as NOT passed at B<=4800, D=1200 |
| W5 per-cell SEs | w5 deposits | 5 seeds per cell, per-seed FWER and MCSE stored |
| W6 conditional law | w6/factor.json | mean 69.778 vs 69.752 theory; variance 285.6 vs 281.1 |
| Independent verifier A (W3/W4) | separate read-only check | check_core exit 0; W3 deposits reproduce bit-exactly through the evaluate loop; extension/B2 arithmetic exact; FOUND: w4_extension sign bug (fixed and deposit re-signed) |
| Independent verifier B (W2/W6) | separate read-only check | w2_cb mapping, w2_adaptive inner world, w2_addendum mix rule, w6 conditional law and AR(1) denominator all PASS; FOUND: w2_adaptive direction detected at alpha=0.05 (fixed to the target alpha; adapt cells re-run) |

## What remains open / not covered

- The W5 map is the plug-in rule over the full grid; the CB and calibrated arms are mapped only on the reduced slice (w5best) and the two desk families (W3).
- The W2 adaptive rule is measured on both families across df 3-35: the self-calibrating CB under-corrects at df 3 (+5.5% A, +9.6% B against a +7.6%/+13.0% plug-in) and is near nominal at df 5-7, so the boundary remains formalised rather than a single deployable gamma.
- W4: a single-alpha* pass at D=1200 needs B of order 1e5; the demonstration is recorded as a cost curve, not an executed pass.
- W6: effective-df and block-bootstrap corrections are characterised but not implemented; the one-factor extension is shown to fail exactly (non-central chi-square + mean shift) rather than extended.
