# S5.1 — Prior art, and the exact boundary of the novelty claim

**Purpose.** Before the terminal report claims novelty, the field has to be read, not guessed.
This file records what each relevant work computes, under which assumption about the variance
ratios, with what guarantee — and what it therefore cannot deliver. Every claim below is either
(a) read in full text, or (b) marked as abstract-level with the retrieval failure
recorded. Machine-readable version: `deposits/prior_art/prior_art.json`.

---

## 1. The dividing line

The desk problem is *many-to-one comparisons with a control, unequal variances, ratios unknown*.
Prior art splits on one axis — **are the ratios known?** — and the answer to that question decides
whether an exact non-simulation critical value exists at all.

| work | ratios | method | guarantee | what it cannot deliver |
|---|---|---|---|---|
| Dunnett (1955) | equal (`λᵢ=1` known) | exact multivariate-t | exact under equal variances | anything when variances differ |
| **Li (2009 dissertation Ch. 2; published as Li 2012 JSTA 11(2))** | **KNOWN `λ`** | exact: condition on control mean, and pool `σ̂₀² = (1/v)Σᵢ Sᵢ²(nᵢ−1)/λᵢ`, `v = Σ(nᵢ−1)`; critical value from a double integral over `z` and `r = σ̂₀/σ₀ ~ √(χ²ᵥ/v)` | exact **given** the ratios | nothing for estimated ratios — pooling is only legitimate because `λ` is known |
| **Li (2009 dissertation Ch. 6.3)** — the closest prior art | **UNKNOWN `λ`** | plug in the UMVU ratio estimator `λ̂ᵢ = c·Sᵢ²/S₀²`, `c = (n₀−3)/(n₀−1)`, with `σ̂₀ = S₀` (no pooling); then, because the correlation of the `zᵢ` depends on the unknown variances and the k-fold integral is infeasible, fall back on two **conservative inequalities**: Bonferroni, and Slepian (`P(Zᵢ ≤ cᵢ ∀i) ≥ ∏P(Zᵢ ≤ cᵢ)`, valid as `cov(zᵢ,zⱼ) ≥ 0`) | **conservative only** — and *measured* as such: §6.4 reports 0.0370 / 0.0408 / 0.0396 (n = (20,10,10)) and 0.0452 / 0.0416 / 0.0446 (n = (30,20,20)) against `alpha = 0.05` at `lambda` = 2/3/4 — below nominal, "slightly conservative" in the author's own words — while Tamhane's Welch-based method *exceeds* nominal (0.0442 / **0.0560** / 0.0502 and 0.0482 / 0.0510 / **0.0532**) | an exact, non-conservative critical value for estimated ratios; any coverage accounting for the `λ̂ − λ` gap; any statement of the power price paid |
| Tamhane & Xi (2023/24, Biom. J., doi:10.1002/bimj.202200300) | unknown | approximates the denominator law by correlated χ² and **generates it by simulation** | simulation-based | determinism, exactness, a proof |
| JAS 2023 (doi:10.1080/02664763.2023.2245179) | unknown, unbalanced | **parametric bootstrap** | simulation-based | determinism, exactness, a proof |
| classical conservative/approximate family (Tamhane 1977/1979; Games–Howell 1976; Hochberg 1976) | unknown | Welch df + Sidak/Bonferroni step-down, studentized maximum modulus | conservative or approximate (some documented liberal in corners) | exactness; a certificate |
| **this work** | **estimated `λ̂`** | tight χ² admission region + population corner envelope + `γ` budget; theorem-grade monotone-step certificate (encloses the continuum — no candidate set, no Lipschitz constant, no `γ`) | certified rule: conditional rejection `≤ α` for **every** λ in the region and **every** draw; envelope rule: `FWER ≤ α + P(λ ∉ Λ)` with coverage certified by exact quadrature and the sup step measured (corners attain the grid sup 12/12, worst deficit 0.0000) | the price of the proof is real (+14…+34% family A, +133…+324% family B in `c`), and the envelope's `γ` term is a genuine weak point, evidenced by one event in ~1200 draws |

## 2. Why the boundary sits where it does

Li's Ch. 6.3 shows precisely *why* the estimated-ratio problem is hard, and the reason is
structural, not technical: once `λ` is estimated, the standardised statistics `tᵢ` are correlated
through a correlation matrix **that is itself a function of the unknown variances**, so the joint
law cannot be evaluated exactly by the same conditioning that makes the known-ratio case exact.
The published responses to that obstruction are (i) an inequality (Bonferroni/Slepian), which is
valid but pays in power and restricts the applicable correlation sign, or (ii) simulation, which is
valid but pays in determinism.

This work's move is different in kind: **do not try to evaluate the joint law at `λ̂`; admit a
region of `λ` values with certified coverage `1−γ` and maximise the exact critical value over that
region.** The gap between `λ̂` and `λ` then appears as an explicit additive term
`P(λ ∉ Λ) ≤ γ` in the guarantee rather than as an inequality penalty — and the region is built
tight (largest level scale whose exact quadrature coverage is still `≥ 1−γ`), which *spends* the `γ` budget instead of
over-covering it. The direction matters and an earlier draft had it backwards: at the A df3 design the
union-bound allocation's miss probability is ~17× *smaller* than the tight region's (2.9e-5 versus
5.0e-4) — it over-covers — and that surplus coverage is what makes its region wider and its critical
value larger (union +7.24% in `c` versus tight +6.01% on the desk sweep). The tight scale takes the
coverage at the target and spends the difference on width and power.

## 3. What is claimed as novel (precisely)

1. **An exact, non-simulation critical value for the estimated-ratio problem**, with the
   `λ̂ − λ` gap given a coverage-certified region instead of an inequality.
2. **A two-tier boundary with prices**: a theorem-grade rule (monotone-step certificate; no
   candidate set, no Lipschitz constant, no `γ`) whose cost is measured (+14…+34% / +133…+324% in
   `c`), alongside the measured-grade envelope actually shipped.
3. **The sup step measured against a dense grid**: corners attain the population sup in 12/12 cells
   (worst deficit 0.0000), and the single observed dominance violation is attributed to a region
   miss (a `γ` event), not to a missed supremum.
4. **Measured price table on desk designs** (12-cell map, dual-arm), plus the **AR(1) boundary**
   (φ ≥ 0.4 all rules liberal) and the **effective-df plug-in**.
5. **A pre-registered calibration gate** for the plug-in's residual noise (criteria frozen before
   the run), so the "no residual bias" claim is falsifiable by construction.

**Not claimed:** the device itself (condition on the nuisance, maximise over a confidence region
for it) is Berger & Boos (1994) and is cited as such. The conditioning algebra for this statistic
family is Li's. What is new is the exact non-simulation treatment of the estimated-ratio case, the
certificate, and the priced boundary.

## 4. Retrieval log (so the record shows what was and was not read)

| source | status |
|---|---|
| Li (2009), BGSU dissertation, OhioLINK `bgsu1254160726` | **full text retrieved** — Ch. 2 §2.3 (exact method), §2.6 (5000-run error-rate study), §2.7 (logistic/lognormal robustness), Ch. 6 §6.1–6.3 (unknown ratios, Bonferroni/Slepian fallback) |
| Li (2012), JSTA 11(2):165–181 | **Citation independently confirmed** (author, year, title, volume/issue, pages) by an unrelated 2015 bibliography — IARJSET 2(12), entry 13-14: "Hong,Li.(2012): A Multiple Comparison Procedure for Populations with Unequal Variances,11(2),165-181". **Full text still not retrievable**: the `mscs.mu.edu` JSTA host times out on every route (direct, plain HTTP, reader proxy 422, archive snapshot). The **content identity** with dissertation Ch. 2 (same author, topic and title family) remains a *labeled inference*, not page-verified; Ch. 2 serves as the full-text proxy |
| IARJSET 2(12), 2015 (Sarmah & Gogoi), "Multiple Comparison Procedures under Equal and unequal population Variances" | **Background only, not prior art**: an empirical (Monte-Carlo) comparison of the classical approximation family (Dunnett, Tamhane T2/T3, Dunn, Spjøtvoll-Stoline, …); no exact critical value, no known-ratio construction. Caveat recorded: its printed abstract is about testing exponentiality, mismatching its own title — a venue-quality flag |
| Tamhane & Xi (2023/24, Biom. J.) | **full text read** (`prior_art_src/tamhane_xi_2023.txt`, 728 lines). Verified against the paper's own claims: the abstract describes the denominator law as "approximated by correlated chi-squares"; their line 278 says of their nested construction "This formula does not have a product structure like (8), and so, it is not possible to use it to derive expressions for [the degrees of freedom]. Instead, we propose to use..."; their Methods 0-3 are numbered as the head-to-head assumes; and **their reference list contains no entry for the known-ratio exact method** (zero occurrences of "Li", "Hong", 2009 or 2012 anywhere in the paper) |
| JAS (2023) parametric bootstrap | abstract verified (T&F full text HTTP 403) |
| classical conservative family | contextual, via the dissertation's review chapter and search results |

## 5. Consequence for the terminal report

The novelty sentence in `REPORT-PROVEN-BOUNDARY.md` can be stated constructively rather than
defensively: *with known ratios the exact answer is Li's; with estimated ratios Li's own treatment
is inequality-conservative and the recent literature simulates; this work gives the estimated-ratio
case an exact non-simulation critical value plus a theorem-grade certificate, and prices both.*
