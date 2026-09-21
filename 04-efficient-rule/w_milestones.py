"""Append the W1-W4 milestones to MILESTONES.md, formatted from deposits.

Idempotent: it removes any previous "## W1-W4" block before writing, so it can
be rerun after any module is rerun.
"""

from __future__ import annotations

import io
import os
import re

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, "MILESTONES.md")
MARK = "## W1-W4: beyond the independent-groups model"


def dep(name):
    p = os.path.join(F.DEP, name)
    if not os.path.exists(p):
        return None
    import json
    with io.open(p, encoding="utf-8") as fh:
        return json.load(fh)


def main():
    inv, reg, low = dep("w1_invariant.json"), dep("w1_region.json"), \
        dep("w1_lowdf.json")
    fac, cer, une = dep("w2_factor.json"), dep("w2_cert.json"), \
        dep("w2_unequal.json")
    tie, rep, strs = dep("w3_tiers.json"), dep("w3_serial_repair.json"), \
        dep("w3_stress.json")
    man = dep("w4_manifest.json")
    prb = dep("probe_beyond.json")
    asy = dep("w3_asymptotic.json")

    L = [MARK, ""]
    A = L.append
    A("Executed against `ROADMAP.md` with its Workstream 0 removed and three "
      "amendments recorded there.  Report: `REPORT-BEYOND.md`, checked by "
      "`verify_beyond.py`.")
    A("")

    A("### W1: the remaining gap, characterised")
    A("")
    if inv:
        A("**DONE.**  `w1_invariant.py`.  The maximal invariant of the "
          "one-candidate problem under the affine group is two-dimensional "
          "and `(T, lambda_hat)` is it: invariance drift %.2e, maximality "
          "recovery error %.2e over 200,000 constructed pairs per design.  So "
          "Section 10 searched the invariant class, not a convenient family.  "
          "The affine group is solvable hence amenable, so Hunt-Stein gives "
          "that a minimax test may be taken invariant, which is the route to "
          "a lower bound over ALL tests.  That step is a **derivation, not a "
          "proof carried through**, and three conditions are listed unclosed."
          % (inv["worst_invariance_drift"], inv["worst_recovery_error"]))
        A("")
        A("The 36-bin discretisation hides at most **%.3f points**, measured "
          "by re-solving at 72 and 144 bins with the shrink step applied to "
          "each." % inv["max_hidden_pp"])
        A("")
    if reg:
        A("**DONE.**  `w1_region.py`.  Monotonicity and determinism priced by "
          "imposing monotonicity as a linear constraint, so both arms are "
          "optima of the same program on the same grid.  Cost of the two "
          "restrictions together: **%+.2f to %+.2f points**.  The first "
          "version of this module rendered the threshold rule onto the T grid "
          "instead and dropped the straddling cell, which handed the "
          "unrestricted rule about a point for free; that is why the "
          "comparison is now a constraint."
          % (reg["min_gain_pp_from_dropping_monotonicity"],
             reg["max_gain_pp_from_dropping_monotonicity"]))
        A("")
    if low:
        A("**DONE.**  `w1_lowdf.py`.  Both group sizes moved over %s.  "
          "Material in %d of %d cells at a %.2f-point bar declared in "
          "advance, best **%.2f points**.  %s."
          % (str(tuple(low["df_grid"])), low["n_material_cells"],
             low["n_cells"], low["material_threshold_pp"], low["max_gain_pp"],
             low["operational_rule"][0].upper()
             + low["operational_rule"][1:]))
        A("")

    A("### W2: aligned returns")
    A("")
    if fac:
        A("**DONE, at the oracle tier.**  `w2_factor.py`.  The roadmap's "
          "paired-difference route was abandoned: it produces a maximum of K "
          "studentised statistics with different studentisers and a "
          "Wishart-coupled denominator set, which nothing factorises.  "
          "Keeping the statistic unpaired and carrying the dependence in a "
          "LATENT factor instead, conditioning on (fbar, Vf, Xbar_0, s_0^2) "
          "makes the K candidates exactly conditionally independent: a "
          "**%d-dimensional outer rule over K one-dimensional inner rules**, "
          "cost linear in K, two of the four dimensions already present in "
          "the existing paper."
          % fac["outer_dimension"])
        A("")
        A("Worst |z| **%.2f** over %d cells against %s brute-force draws.  "
          "Non-vacuity is built into the same table: the independent-groups "
          "law is wrong by up to %.1f%% here and the evaluator tracks the "
          "truth across that range."
          % (fac["worst_abs_z"], len(fac["falsifiers"]),
             "{:,}".format(fac["brute_draws"]),
             100 * fac["max_relative_departure_from_independent_law"]))
        A("")
        A("An earlier Monte Carlo probe of the same claim reported it "
          "FALSIFIED at 8.75 sigma.  That probe sampled one variable by "
          "picking a Gauss-Legendre NODE INDEX uniformly, which is not a draw "
          "from the law; corrected, it agrees.  Recorded because a falsifier "
          "that is itself wrong is the most expensive kind.")
        A("")
    if cer:
        A("**PARTIAL.**  `w2_cert.py`.  The numerator of the aligned "
          "statistic carries the SQUARED CONTRAST of the loadings and the "
          "denominator carries their SUM OF SQUARES, so the iid rule is "
          "conservative when the loadings share a sign and liberal when they "
          "do not.  Measured: **%d of %d cells liberal, and the split is "
          "exactly along the sign** (%d of %d opposite-sign, %d of %d "
          "same-sign), worst **%.5f**, **%.2f times** the ceiling.  A "
          "long-short candidate against a directional benchmark is the "
          "opposite-sign case."
          % (cer["n_liberal_cells"], cer["n_cells"],
             cer["n_opposite_sign_liberal"], cer["n_opposite_sign_cells"],
             cer["n_same_sign_liberal"],
             cer["n_cells"] - cer["n_opposite_sign_cells"],
             cer["worst_fwer_iid_rule_on_aligned_data"],
             cer["worst_ratio_to_abar"]))
        A("")
        A("The first version of this module set the benchmark loading to "
          "zero, which is exactly the boundary between the two regimes, and "
          "found nothing.  The nuisance is **three-dimensional per "
          "candidate**, and only the variance ratio is estimable from what "
          "the rule sees.  A constant holding over the %d-point declared box "
          "costs %.2f times the worst known-parameter oracle constant and "
          "%.2f times the best, and that supremum is a **grid maximum, not a "
          "continuum certificate**."
          % (cer["box"]["points"],
             cer["certificate_price_vs_worst_oracle_x"],
             cer["certificate_price_vs_best_oracle_x"]))
        A("")
    if une:
        A("**DONE.**  `w2_unequal.py`.  The benchmark's prehistory buys "
          "**%.2f points** through the two-window design, which needs "
          "constant variance only, and that is essentially nothing; and "
          "**%.2f points** through full history, which needs a constant mean "
          "too.  Breaking that assumption in the dangerous direction, a "
          "prehistory running BELOW the later mean, takes the worst rate to "
          "**%.5f**, %.2f times the ceiling, on `%s`.  There is no free lunch "
          "in the incumbent's track record."
          % (une["two_window_max_power_gain_pp"],
             une["full_history_max_power_gain_pp"],
             une["broken_stationarity_worst_fwer"],
             une["broken_ratio_to_abar"],
             une["broken_stationarity_worst_design"]))
        A("")

    A("### W3: validity tiers, a new defect, and its repair")
    A("")
    if tie and rep:
        A("**DONE, and it found the programme's second result-changing "
          "defect.**  `w3_tiers.py`, `w3_serial_repair.py`.")
        A("")
        A("Tier A holds exactly, as Section 11 says.  Tier B and tier C do "
          "not, and not marginally: with the coefficient estimated on the "
          "evaluation period the worst realised rate is **%.5f** against a "
          "ceiling of %.3f."
          % (tie["serial_tierC_worst_fwer"], tie["abar"]))
        A("")
        A("The mechanism is that Cochrane-Orcutt maps a mean `mu` to "
          "`mu(1 - phi)`, so with the coefficient estimated SEPARATELY IN "
          "EACH GROUP the numerator acquires `mu (phi_hat_0 - phi_hat_1)`, "
          "zero in expectation but not in any sample, and spreading linearly "
          "in the common mean.  It is the Prais-Winsten defect's shape "
          "arriving by another route: there the multiplier varied with the "
          "group's LENGTH, here with the group's ESTIMATE.")
        A("")
        A("**The repair is to pool**: one coefficient from both groups, "
          "applied to both, so the multiplier is identical by construction "
          "and the common mean cancels whatever the estimate is.  Separate "
          "worst **%.5f** (%.1f times the ceiling); pooled worst **%.5f** "
          "(%.1f times).  The signature is that the pooled rate is flat in "
          "the common mean: spread %.5f at tier B and %.5f at tier C.  What "
          "pooling does NOT remove is imperfect whitening at an estimated "
          "coefficient, which the mu = 0 rows isolate, and that residual "
          "breaches on its own: every cell in both arms is above the ceiling. "
          " **Pooling is an improvement, not a repair.**  Section 11's serial "
          "result is a tier-A result."
          % (rep["separate_worst_fwer"], rep["separate_worst_ratio"],
             rep["pooled_worst_fwer"], rep["pooled_worst_ratio"],
             rep["pooled_spread_over_mu_tierB"],
             rep["pooled_spread_over_mu_tierC"]))
        A("")
    if strs:
        A("**DONE.**  `w3_stress.py`.  %d of %d cells breach a bar of %.0f "
          "standard errors declared in advance.  Survives heavy tails, skew, "
          "GARCH, a variance regime shift and negative serial correlation, "
          "all conservatively.  Breaks on untreated positive serial "
          "correlation (worst %.5f, %.2f times the ceiling) and on a **short "
          "variance window** -- estimating the variance on a trailing "
          "subwindow while the mean uses the whole sample, which is what a "
          "rolling risk model does.  That last one is not in the paper's list "
          "of boundaries and is a thing desks do by default."
          % (strs["n_breaches"], strs["n_cells"], strs["breach_se"],
             strs["worst_fwer"], strs["worst_ratio_to_abar"]))
        A("")

    A("### W4: the promotion protocol")
    A("")
    if man:
        A("**DONE.**  `PROMOTION-PROTOCOL.md`, `manifest.schema.json`, "
          "`w4_manifest.py`.  Three verdicts (%s); %d of %d refusal probes "
          "refused on their own ground; a nuisance estimated on the "
          "evaluation period returns %s rather than a refusal, because that "
          "run is worth doing as long as the report does not call it a "
          "certificate."
          % (", ".join(man["verdict_labels"]), man["n_refused"],
             man["n_probes"], man["tier_c_verdict"]))
        A("")

    if asy:
        A("**DONE.**  `w3_asymptotic.py`.  HAC and block bootstraps measured "
          "against the exact rule on the same draws.  With the coefficient "
          "DECLARED the exact rule holds at %.2f times the target and nothing "
          "else is given that information.  Among the rows that must estimate "
          "it, none holds: the exact rule with a pooled estimate is best at "
          "%.2f times, the closest alternative is %.2f times, and the exact "
          "rule with separate estimates is the WORST row in the table.  The "
          "exact law buys its guarantee from the coefficient being declared, "
          "not from being exact."
          % (asy["worst_ratio_to_abar"]["the exact rule, coefficient declared"],
             asy["worst_ratio_to_abar"]["the exact rule, coefficient pooled"],
             min(asy["worst_ratio_to_abar"][m] for m in asy["methods"]
                 if not m.startswith("the exact rule"))))
        A("")

    if prb:
        A("### The checks were shown to bind")
        A("")
        A("**DONE.**  `probe_beyond.py`.  %d of %d probes caught, including "
          "%d of %d transpositions, each moving one deposit field and "
          "requiring a new verifier failure against the unchanged report.  "
          "One probe was replaced rather than recorded as a miss: swapping "
          "the monotone and unrestricted losses is a no-op, because the "
          "result of that section is that they are equal, and a vacuous probe "
          "reports a pass."
          % (prb["n_caught"], prb["n_live"], prb["n_transposition_caught"],
             prb["n_transposition_probes"]))
        A("")

    A("### M10 - the paper carries W1, W3 and the factor case")
    A("")
    A("**DONE.**  `paper_content_beyond.py`, `paper_figs_beyond.py`, "
      "`paper_verify_beyond.py`, `coverage_paper_beyond.py`, with the "
      "insertion driver extended.")
    A("")
    A("New subsections 10.8, 11.4, 11.5 and 11.6; equations 47 to 50; Tables "
      "28 to 34; Figures 18, 22 and 23, with the three that 10.8 displaced "
      "renumbered to 19, 20 and 21 so the figures run in document order; the "
      "old 10.8 to 10.10 renumbered to 10.9 to 10.11; contributions (14) to "
      "(16); and additions to the abstract, the limitations and the "
      "conclusion.")
    A("")
    A("**166 new checks pass** beside the existing 558, and "
      "`coverage_paper_beyond.py` reports **100% of the new material's "
      "numeric tokens covered**, with none uncovered.  Getting there moved "
      "two declared parameters out of prose and into their deposits, because "
      "a number typed into a sentence is a number no check can bind.")
    A("")
    A("**What the build found that the modules had not.**  Four defects, none "
      "of them numerical:")
    A("")
    A("* the first draft used nine capitalised lead-ins and forty capitalised "
      "words for emphasis, in a document whose only all-capital words are "
      "FWER, ANOVA and JSON, plus one markdown bold that would have printed "
      "as literal asterisks;")
    A("* inserting a figure into 10.8 put Figure 21 between 17 and 18, so the "
      "figures no longer ran in document order;")
    A("* Section 9.10 called the general factor case open three sentences "
      "before the amended text said 11.6 implements it;")
    A("* Section 11.3 said both nuisance repairs are available to a desk that "
      "plans ahead, which 11.4 then measures to be false for the serial "
      "coefficient.")
    A("")
    A("The first was found by reading the rendered pages, which is the only "
      "way it could have been; the other three by cross-section consistency "
      "checks that now live in the verifier.")
    A("")
    A("**What the paper does NOT now claim.**  Section 11.6 presents the "
      "reduction as the one Section 9.10 already names, carried out at more "
      "than one candidate, not as a new construction; the bound of 10.8 is "
      "over invariant tests with Hunt-Stein labelled a route; and the "
      "supremum over the loading box is labelled a grid maximum wherever it "
      "appears.")
    A("")

    A("### What this programme did not close")
    A("")
    A("* the Hunt-Stein lower bound, as a proof rather than a route;")
    A("* the aligned-factor certificate, as a continuum enclosure rather than "
      "a grid maximum;")
    A("* tier C, which remains without a guarantee and remains the realistic "
      "case;")
    A("* residual cross-correlation after one factor is removed, measured as "
      "conservative but not certified;")
    A("* **no new certificate at all.**  Every guarantee here is inherited "
      "from the preceding programme or is a grid maximum.")
    A("")

    s = io.open(P, encoding="utf-8").read()
    s = re.sub(re.escape(MARK) + r".*?(?=\n## |\Z)", "", s, flags=re.S)
    s = s.rstrip() + "\n\n" + "\n".join(L) + "\n"
    io.open(P, "w", encoding="utf-8").write(s)
    print("MILESTONES.md updated with the W1-W4 block (%d lines)" % len(L))


if __name__ == "__main__":
    main()
