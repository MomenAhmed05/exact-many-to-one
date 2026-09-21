# Stage 2 — Map: correlation, effective df, and the power price

**Question this stage answers.** The guarantee `FWER <= alpha + P(lambda not in Lambda)` assumes
independent normal observations. Two things follow: (i) what is the *price in power* the envelope
charges on the desk designs, and (ii) where does the guarantee stop holding once observations are
serially correlated? Stage 2 answers both, and settles whether an effective-degrees-of-freedom
plug-in can rescue the correlated case.

Status: **all four steps complete** (S2.1 AR(1) at R>=2000, S2.2 effective-df plug-in, S2.3 shifted
quadrature, S2.4 power table). Every number below is read from `deposits/ar2/ar1.json`,
`deposits/ar2/effdf.json` and `deposits/power/power_table.json`.

---

## 1. S2.1 — the AR(1) boundary is far tighter than the pilot suggested (R = 2000)

Design: 3 groups, n = 36 or 60, `lambda = (1, 4)`, AR(1) innovations with lag-1 correlation `phi`;
R = 2000 draws per cell; moving-block bootstrap with block length `round(n^(1/3))` (3 at n=36, 4 at
n=60). MCSE at `alpha = 0.025` is 0.0035-0.0077; at the cell carrying the boundary claim (n=36, `phi`=0.1) the doubling is 5.1 sigma, while at the largest MCSE in the table it is 3.2 sigma — the claim rests on the former, not on every cell.

| n | phi | plug-in | envelope | bootstrap |
|---|---|---|---|---|
| 36 | 0.0 | 0.0250 (+0%) | 0.0245 (-2%) | 0.0380 (+52%) |
| 36 | 0.1 | **0.0500 (+100%)** | 0.0480 (+92%) | 0.0520 (+108%) |
| 36 | 0.2 | 0.0570 (+128%) | 0.0545 (+118%) | 0.0510 (+104%) |
| 36 | 0.3 | 0.1010 (+304%) | 0.0980 (+292%) | 0.0650 (+160%) |
| 36 | 0.4 | 0.1245 (+398%) | 0.1200 (+380%) | 0.0765 (+206%) |
| 60 | 0.0 | 0.0245 (-2%) | 0.0225 (-10%) | 0.0335 (+34%) |
| 60 | 0.1 | 0.0330 (+32%) | 0.0320 (+28%) | 0.0360 (+44%) |
| 60 | 0.2 | 0.0690 (+176%) | 0.0670 (+168%) | 0.0495 (+98%) |
| 60 | 0.3 | 0.0910 (+264%) | 0.0885 (+254%) | 0.0545 (+118%) |
| 60 | 0.4 | 0.1360 (+444%) | 0.1320 (+428%) | 0.0665 (+166%) |

**Findings, stated plainly.**

1. **The failure boundary is `phi ~ 0.1`, not 0.2-0.4.** At n = 36, `phi = 0.1` already gives FWER
   0.0500 — exactly twice nominal, 5 sigma from 0.025. The earlier R = 150 pilot reported
   "phi <= 0.2 unresolved"; that was a power problem, and 13x the replicates resolve it against us.
2. **`phi = 0` is clean** (0.0250 and 0.0245), so the harness and the exact law agree where they
   must; the inflation is a correlation effect, not a coding artefact.
3. **The envelope does not help under serial correlation** (0.0480 vs 0.0500 at n=36, phi=0.1): its
   region absorbs *variance-ratio* uncertainty, not *dependence*. This is a scope statement for the
   desk guarantee, not a defect in the envelope.
4. **The bootstrap is worse at small phi and better at large phi** (+52% at phi=0 to +206% at 0.4,
   versus the plug-in's +0% to +398%). Within the tested grid (n in {36,60}, `phi` in {0,0.1,...,0.4}) no cell puts it within 5%
   (`bootstrap_boundary_phi_max: null`; that field means *no qualifying cell in this grid* — and it is also the value the writer emits on exception — so the statement is about the grid, not "any phi"; the bootstrap also carries an undisclosed extra Monte-Carlo layer, `R_boot = 200` per draw). It does not rescue the guarantee:
   it simply fails differently. 27 of the 30 rule-cells (10 cells x 3 rules) fall outside +/-5%.

## 2. S2.2 — the effective-df plug-in does not repair the inflation

Two variants were implemented, both plugging an estimate of the effective degrees of freedom into
the *critical value only* (the observed statistics keep their real n), with every rule seeing the
same replicates so the comparison is paired:

- (a) sample-size route: `nu_eff = (n-1) * mean_i (1 - phi_i)/(1 + phi_i)` from lag-1 autocorrelations;
- (b) Satterthwaite route: `nu_eff = (tr AV)^2 / tr((AV)^2)`, exact for `Var(s^2)` under AR(1).

Result: **`effdf_satt` is within 5% in only 2 of 10 cells, and `max phi = 0.00`** — i.e. it holds only
where the naive rule already held (phi = 0). The corrections move the excess by 0-3% of its size (Satterthwaite route; the sample-size route reaches 7% at n=60, phi=0.4):

| n | phi | naive | effdf (Satterthwaite) |
|---|---|---|---|
| 36 | 0.1 | 0.0500 (+100%) | 0.0495 (+98%) |
| 36 | 0.2 | 0.0570 (+128%) | 0.0570 (+128%) |
| 36 | 0.3 | 0.1010 (+304%) | 0.0995 (+298%) |
| 36 | 0.4 | 0.1245 (+398%) | 0.1215 (+386%) |
| 60 | 0.2 | 0.0690 (+176%) | 0.0680 (+172%) |
| 60 | 0.4 | 0.1360 (+444%) | 0.1340 (+436%) |

**Conclusion: a negative result, and a clean one.** Shrinking the degrees of freedom changes the
denominator's *spread*, but the inflation under AR(1) is driven by the *dependence between numerator
and denominator* across comparisons — a mechanism no df substitution touches. Practical reading:
for serially correlated data the desk should not rely on any of these three rules; the guarantee's
scope statement must name `phi` explicitly.

Implementation note (auditability): per-row `nu_eff` cannot be fed to the law one row at a time, so
rows are grouped into bins of width `bin_w = 0.25` in `nu`; the deposit records the maximum rounding
error per cell (`max_nu_bin_err_satt = 0.125`), which moves `c` by O(1e-3) — an order of magnitude
below the effects reported above.

## 3. S2.3 — shifted quadrature, and how it was validated

The power question needs the *alternative* law, which the harness did not expose. The exact
alternative rejection probability at a mean shift is

    P(reject | s^2, delta) = E_U [ prod_i Phi( s_i (U + c d_i) - s_i nu_i SE_i^pop ) ],
    SE_i^pop = sqrt( lambda_i / n_i + 1 / n_0 ),    nu_i = (mu_i - mu_0) / SE_i^pop,

evaluated with the same 40-panel x order-10 rule as the null law. This was validated against a
brute-force simulation of the whole procedure (draw data, compute the rule, compare rejection
rates), on the **same** `s^2` draws so the comparison is paired:

- unconditional differences (fresh draws on both sides): max 0.071 against a Monte-Carlo SE of ~0.019
  (600 reps), i.e. a few SE — consistent with noise rather than a quadrature error;
- **paired** differences (identical `s^2`): max **0.0028** — this is the number that validates the
  quadrature, and it is the one reported as evidence.

Bugs found and fixed during this validation, recorded because they were real: a scalar-`nu`
broadcast, a dot-shape mismatch, a units error (the shift is `s_i * Delta_i`, not `s_i * d_i * nu_i`),
and an endpoint-rule bias (switched to the harness's own node rule).

## 4. S2.4 — the power price at the desk

For each desk cell: the effect `nu80` at which the **oracle** rule reaches 80% power, and the
envelope's power there. `loss = 0.80 - power_env(nu80)`.

| cell | c_oracle | c_env | nu80 | power_env@nu80 | **loss** |
|---|---|---|---|---|---|
| A df3 | 3.2201 | 3.2702 | 4.000 | 0.745 | 0.055 |
| B df3 | 2.7003 | 3.3224 | 2.917 | 0.615 | **0.185** |
| A df5 | 2.6420 | 2.7114 | 3.587 | 0.783 | 0.017 |
| B df5 | 2.5400 | 2.8253 | 2.777 | 0.712 | 0.088 |
| A df7 | 2.4468 | 2.5264 | 3.257 | 0.780 | 0.020 |
| B df7 | 2.4876 | 2.6577 | 2.679 | 0.746 | 0.054 |
| A df11 | 2.2972 | 2.3811 | 3.253 | 0.777 | 0.023 |
| B df11 | 2.4558 | 2.5417 | 2.720 | 0.774 | 0.026 |
| A df19 | 2.2100 | 2.3087 | 2.962 | 0.772 | 0.028 |
| B df19 | 2.4521 | 2.5293 | 2.658 | 0.775 | 0.025 |
| A df35 | 2.1767 | 2.2892 | 2.929 | 0.767 | 0.033 |
| B df35 | 2.4702 | 2.5341 | 2.525 | 0.780 | 0.020 |

**Reading.** The envelope costs **1.7 to 18.5 percentage points of power** at the effect size where
the oracle would have had 80%. The cost is concentrated exactly where the admission region is wide:
family B at low df (18.5 pp at df 3, 8.8 pp at df 5), while every other cell loses at most 5.5 pp.
A note on metric discipline: an earlier draft of this table labelled `1 - power` as "loss", which
inflated the headline to 0.385; the corrected definition above is the one in the deposit.

**One cell needs a caveat.** A df3's `nu80` sits at the grid endpoint (4.000), where the oracle's
power is 0.7561 — not 0.80 — so 0.0439 of that cell's 0.0546 "loss" is oracle shortfall rather
than envelope cost (its interior gap is 0.0107). The worst (0.185, B df3) and the minimum (0.017,
A df5) are interior and sound, so the A-family range is quoted from interior cells only.

**One cell needs a caveat.** A df3's `nu80` sits at the grid endpoint (4.000), where the oracle's
power is 0.7561 — not 0.80 — so 0.0439 of that cell's 0.0546 "loss" is oracle shortfall rather
than envelope cost (its interior gap is 0.0107). The worst (0.185, B df3) and the minimum (0.017,
A df5) are interior and sound, so the A-family range is quoted from interior cells only.

**One cell needs a caveat.** A df3's `nu80` sits at the grid endpoint (4.000), where the oracle's
power is 0.7561 — not 0.80 — so 0.0439 of that cell's 0.0546 "loss" is oracle shortfall rather
than envelope cost (its interior gap is 0.0107). The worst (0.185, B df3) and the minimum (0.017,
A df5) are interior and sound, so the A-family range is quoted from interior cells only.

## 5. What Stage 2 adds to the boundary statement

- **Power price, measured**: the envelope is cheap for family A (+1.7 to 5.5 pp) and expensive only
  at wide regions with few benchmark degrees of freedom (B df3: 18.5 pp).
- **Correlation scope**: the guarantee is an *independence* guarantee. Under AR(1) it fails from
  `phi ~ 0.1` at n=36, and no member of {plug-in, envelope, moving-block bootstrap} survives, nor
  does the effective-df plug-in (S2.2).
- **Honest limits**: the AR(1) study covers `phi <= 0.4` and n in {36, 60}; the power table is
  conditional on a single alternative direction per cell; the effdf result carries a documented
  0.25-wide binning of `nu_eff`.

## 6. Reproduction

    W26_AR_R=2000 W26_AR_RB=200 python w6_ar2.py     # S2.1  (~18.5k s)
    W26_EFF_R=2000 python ar1_effdf.py               # S2.2  (~8.7k s)
    W26_POW_D=40 python power.py                     # S2.3/S2.4 (~811 s)
