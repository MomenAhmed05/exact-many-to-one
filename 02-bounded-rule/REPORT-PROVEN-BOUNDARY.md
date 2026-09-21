# REPORT-PROVEN-BOUNDARY

**The desk objective, closed with a stated boundary.** A reproducible, data-computable critical value
that controls FWER across broad many-to-one designs, delivered with a *proven* tier and a *measured*
tier, each with its price, its scope, and its failure modes named. This document consolidates Stages
1–5; each stage's detail lives in its own report (`STAGE1-CERTIFICATE.md`, `STAGE2-MAP.md`,
`STAGE3-INSTRUMENTS.md`, `STAGE4-FRONTIER.md`, `STAGE5-PUBLICATION.md`) and every number below is read
from a deposit under `deposits/`.

---

## 1. The boundary statement (the headline deliverable)

| tier | rule | what is proven | what it costs |
|---|---|---|---|
| **theorem-grade** | `c_rigorous` — smallest `c` whose monotone-step lower bound on `min_{lambda in Lambda} F(c,lambda)` reaches `1-alpha` | conditional rejection `<= alpha` for **every** `lambda` in the region and **every** draw; the **sup step carries no `gamma`** (no candidate set, no Lipschitz constant, no quadrature-error
assumption); the region is the same `gamma`-budgeted one, so the unconditional bound stays `alpha + gamma` | `c` inflated **+14.1…+33.5%** (family A) and **+133.0…+324.4%** (family B); realised FWER 0.0018–0.0048 (A) |
| **measured-grade** | `c_env` — population exact critical value maximised over the tight region's candidate vertices (the rule actually shipped) | sup attained at a vertex in **12/12** cells on an 11×11 per-axis log-grid (worst deficit **0.0000**); guarantee `FWER <= alpha + P(lambda not in Lambda)`, with the region's coverage certified by exact 1-D quadrature | `c` inflated **+1.58…+22.99%** (mean +6.01%) on the desk cells, versus +1.78…+23.64% (mean +7.24%) for the union-bound allocation |

**The guarantee as shipped, in full:**

    FWER  <=  alpha  +  P(lambda not in Lambda),      coverage(Lambda) >= 1 - gamma  certified by quadrature,
    alpha = 0.025, gamma = 0.0005, Lambda = the tight chi-square region at the largest feasible level scale.

Measured on the desk sweep: **12/12 dominance** (`rej(c_env) <= rej(c_oracle)` per draw, excess
-1.9e-13…-8.3e-4), realised FWER **0.0091–0.0239 <= 0.0255**, **0** region escapes in the sweep.

## 2. The sup step, and the one violation, explained

The measured tier needs exactly one fact — that the region's supremum of `c*_pop` is attained on the
candidate set. That fact is measured, not assumed: corners attain the grid sup in **12/12** cells
(worst corner deficit 0.0000), and the certified rule's margin above the grid sup is +0.41…+0.78 (A)
and +3.48…+10.33 (B).

The desk sweep's single substantive dominance violation (`+0.01180`, family B at df 3) was dissected
draw by draw and is a **region miss — a `gamma`-budget event, not a missed supremum**:

- on that draw the corner set attains the population sup exactly (corner max = fine-grid sup = 2.576260,
  deficit 0.0000, argmax at a corner);
- `lambda_true = (1.0, 25.0)` lay **outside** that draw's region `[3.1874, 23291.66] x [216.21, 1579907.90]`;
- the observed excess decomposes as +0.00503 implied at the argmax plus the draw's remainder.

**Rate, stated without smoothing:** 1 event in the sweep's ~1200 draws = 8.3e-4 against a designed
`gamma = 5e-4` (the denominator belongs to the sweep; the probe slice itself is 200 draws)
(~1.6x). At the 200-draw slice level a single event has ~9% probability even if the true rate were
`gamma`, so the observation is *consistent* with the design — but that slice alone cannot prove the
designed rate, and this report does not claim it does. **Resolved by direct probe** (20,000 draws per cell, 240,000 total, membership tests only): the pooled observed miss rate is **4.96e-04** against the designed `gamma = 5.00e-04`, worst |observed - designed| = 2.01e-04, max |z| = 1.13 (`deposits/cert/miss_probe.json`) - the region machinery reproduces its own quadrature certificate in every cell, so the designed rate is certified as well as empirically confirmed and the lone sweep event is no longer the only evidence.

**K=5/K=8 rigour (45 rows, full plan scope, accurate grid):** **0 bound violations** in both
(`fwer_bound` -0.0012…0.0226 at K=5; -0.0038…0.0228 at K=8), coverage escapes **0.000 in all 45 rows**,
and one known quadrature-level pointwise artifact (+1.595e-06) carried forward rather than suppressed.
This evidence is **screening-grade**: at these draw counts the MCSE rules out violations larger than
~1e-2, not the 5e-4 slack — a limitation fixed in the plan *before* the run.

## 3. Price in power, and the scope of the guarantee

**Power price at the desk** (effect size `nu80` where the *oracle* reaches 80% power):
loss `= 0.80 - power_env(nu80)` runs **0.017 to 0.185**. The cost is concentrated exactly where the
region is wide and the benchmark is thin: family B at df 3 loses **18.5 percentage points**, df 5 loses
8.8 pp, and every other cell loses at most 5.5 pp.

**Correlation scope (S2.1/S2.2, R = 2000).** The guarantee is an *independence* guarantee, and the
boundary is much tighter than the pilot suggested:

- AR(1) at n = 36, `phi = 0.1`: plug-in FWER **0.0500** — exactly twice nominal, ~5 sigma;
- `phi = 0` is clean (0.0250 / 0.0245), so the inflation is a correlation effect, not an artefact;
- the envelope does **not** help (0.0480), because its region absorbs *variance-ratio* uncertainty, not
  dependence;
- the moving-block bootstrap is **never** within 5% at any `phi` (and is worse than the plug-in at
  small `phi`);
- the **effective-df plug-in does not repair it**: 2/10 cells within 5%, `max phi = 0.00`, and both
  variants move the excess by only 1–3% of its size.

**Practical reading for the desk:** use the envelope with its priced guarantee when observations are
independent; **do not** rely on any of the three rules once serial correlation is plausible
(`phi >= ~0.1` at these sample sizes).

## 4. Instruments: what is trustworthy, and what is not

| instrument | verdict | evidence |
|---|---|---|
| **Calibration gate** (pre-registered) | **PASS at `B >= 1500`** — 9/10 aggregate cells satisfy the criterion frozen before the run; the single failure is at the smallest `B` on the hardest case and is *over*-conservative | `deposits/calib_vr/gate.json`; VR residual shrinks with `B` (S5: -0.0017 -> -0.0001); plain arm fails 5 cells |
| **Selector v2** | **BAR MISSED (0.788 vs 0.85) — not used as a rule** | `deposits/selector/selector2.json`; six features incl. K, bench/candidate df, spread, power cost; routing stays a default-with-fallback heuristic |
| **Region shape (S4.1)** | **At the family floor** — every split/shape tried changes `c` by at most **+0.39%** at equal certified coverage; the HPD shape can *raise* `c` while narrowing width | `deposits/frontier/shapes.json`; all shapes at coverage 0.999500 |
| **D2 factor (S4.2)** | **Boundary, quantified** — with known loadings the numerator-scale correction is *exact* for `cand_only` (`w = 1.000`), fixes drastic over-conservatism for `same` (0.0000 -> 0.025), and **overshoots** `sign_flip`; the exact K=1 law requires a **2.644x** inflation at `rho = 0.75` where the closed form supplies 1.323x, both from
`deposits/cert/k1_reference.json` | `deposits/w6/factor_d2.json`; exact K=1 `c` = 2.5746 / 3.4537 / 5.2733 at `rho` = 0.25 / 0.5 / 0.75 |

## 5. What is new, stated against the prior art

Prior art divides on whether the variance ratios are known. With **known** ratios an exact
non-simulation critical value exists — Li (2009 dissertation Ch. 2, published as Li 2012, JSTA 11(2)) —
but it is not applicable here, and its pooling step is legitimate *only* because the ratios are known.
With **estimated** ratios Li's own treatment (dissertation Ch. 6.3) falls back on Bonferroni/Slepian
**conservative inequalities** (the `z_i` correlation depends on the unknown variances), and the recent
literature resorts to **simulation** (Tamhane & Xi 2023, Biom. J.; JAS 2023 parametric bootstrap).

**Novel here:** an exact non-simulation critical value for the *estimated*-ratio problem, obtained by
admitting a coverage-certified region and paying an explicit additive `gamma` instead of an inequality
penalty; a two-tier boundary with both prices measured; the sup step measured on a dense grid with the
one violation attributed to a `gamma` event; the desk power-price table; the AR(1) boundary with its
negative effective-df result; and a pre-registered calibration gate. **Not claimed:** the device itself
(maximise over a confidence region for the nuisance) is Berger & Boos (1994), and the conditioning
algebra is Li's — both cited as such.

## 6. Verification chain (what makes any of this load-bearing)

- **Falsifier suite** (`tests_falsifiers.py`): 9/9 gates, including the widened T2 gate that enforces
  strict JSON across the whole deposit tree (262 deposits, 0 non-strict).
- **Core checks** (`check_core.py`): ALL CORE CHECKS PASSED, including the envelope identity checks.
- **Kernel regression**: the lint-only sweep of the shared verified kernel was proven not to move a
  number — A/B against a pre-sweep backup gives max |diff| **0.000e+00** for `cdf`/`critical_value` at
  K = 1,2,5,8 on both grids, and `oracle_c` reproduces all 12 deposited desk cells exactly.
- **Certificate falsifier**: the certified rule's worst conditional rejection anywhere on the dense
  grid is **0.024964 <= alpha**, in 12/12 cells.
- **Independent replication kit** (`replicate/`): method spec + inputs with the *actual recorded seeds*,
  outputs withheld, four falsifiers the replicator must run, and the paired-dominance statistic named as
  the decisive quantity.
- **Independent replication (run, not just shipped)**: a separate implementation written from the kit alone reproduced every design-level quantity (family B to four decimals in all six cells; family A to the third-fourth decimal, worst case 9e-4 at the hardest corner A df3), passed all four falsifiers, kept paired dominance `<= 0` in 12/12, and reproduced the B-df3 region-miss rate at 1/200 - while exposing two kit specification defects (conditional-vs-unconditional law; unpinned draw order), both now fixed. Detail in `STAGE5-PUBLICATION.md` §2.
- **Stage verifier rounds**: one read-only audit per stage, fired after the stage closed (see the
  stage reports for findings and their resolution).

## 7. Honest limitations

1. **Sup attainment is measured, not proven.** The certificate is the proven object; the corner
   envelope's sup step is measured on a 121-points-per-cell grid and stays measured-grade by
   construction.
2. **The `gamma`-miss rate rests on one event in ~1200 draws** — consistent with the design, not
   verified at the designed rate.
3. **K >= 5 evidence is screening-grade** (MCSE ~1e-2); it does not resolve the 5e-4 slack.
4. **AR(1) scope**: `phi <= 0.4`, n in {36, 60}; the boundary is named, not fully mapped, and no rule
   in the tested set survives.
5. **Effective-df results carry a documented 0.25-wide binning** of `nu_eff` (`max_nu_bin_err = 0.125`),
   which moves `c` by O(1e-3) — an order of magnitude below the reported effects.
6. **D2's composite law for K >= 2 is open**: only the K = 1 exact quadrature and the closed-form
   boundary are delivered; the full factor law needs its own nuisance and therefore its own envelope.
7. **Selector v2 missed its pre-registered bar** and is therefore not used as a rule.
8. **The two tiers cannot currently be collapsed**: the certificate proves the continuum sup is at most `c_rigorous`, but not that the vertex set attains it; coordinate monotonicity of `lambda -> c*_pop` is the obvious route and is false here (mixed-corner argmaxes, e.g. A df3 selects `lambda* = [4.43, 9.5e-05]`), so attainment needs an argument specific to the equal-tail construction.

## 8. Reproduction

    # Stage 1 (certificate, sup falsifier, leak probe, certified rule, K=5/K=8 rigour)
    python cert_sup.py; python cert_falsify.py; python cert_leak.py; python cert_rule.py
    W26_ENV_KS=5 ... python w5_env.py ; W26_ENV_KS=8 ... python w5_env.py
    # Stage 2 (AR(1), effective df, power)
    W26_AR_R=2000 python w6_ar2.py ; W26_EFF_R=2000 python ar1_effdf.py ; W26_POW_D=40 python power.py
    # Stage 3 (instruments)
    W26_POOL=4 python calib_gate.py ; python selector2.py
    # Stage 4 (frontier, factor)
    python frontier.py ; W26_D2_B=2000 python factor_d2.py
    # Stage 5 (kit)
    cd replicate && cat README.md METHOD.md inputs/designs.json
