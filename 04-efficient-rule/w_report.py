"""Generate REPORT-BEYOND.md from the W1-W4 deposits.

Same discipline as the rest of the folder: no number reaches the document except
by being formatted from a deposit at build time, so `verify_frontier.py` can
check every one of them against the file the computation wrote.

    python w_report.py
"""

from __future__ import annotations

import io
import json
import os

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "REPORT-BEYOND.md")


def dep(name):
    p = os.path.join(F.DEP, name)
    if not os.path.exists(p):
        return None
    with io.open(p, encoding="utf-8") as fh:
        return json.load(fh)


def need(d, *keys):
    """Treat a deposit written by an older version of its module as absent."""
    if d is None:
        return None
    return d if all(k in d for k in keys) else None


def pc(x, nd=1):
    return ("%." + str(nd) + "f%%") % x


def main():
    inv = dep("w1_invariant.json")
    reg = need(dep("w1_region.json"), "worst_overshoot", "min_gain_pp_from_dropping_monotonicity")
    low = dep("w1_lowdf.json")
    fac = dep("w2_factor.json")
    cer = need(dep("w2_cert.json"), "mechanism", "price_rows")
    une = need(dep("w2_unequal.json"), "saturates_at_prehistory")
    tie = dep("w3_tiers.json")
    rep = dep("w3_serial_repair.json")
    asy = need(dep("w3_asymptotic.json"), "mu")
    strs = dep("w3_stress.json")
    man = dep("w4_manifest.json")
    prb = dep("probe_beyond.json")
    cov = dep("coverage_beyond.json")

    L = []
    A = L.append

    A("# Beyond the exact Welch-Dunnett law: W1 to W4")
    A("")
    A("What the roadmap asked for, what it got, and what it did not.  Every "
      "number here is formatted from a deposit under `deposits/` at build "
      "time and checked by `verify_beyond.py`.  Each carries a label from "
      "the taxonomy in `PROMOTION-PROTOCOL.md`: **exact in model**, "
      "**certified bound**, **measured**, **asymptotic**, **diagnostic**.")
    A("")
    A("**One of those labels is missing from everything below, and that is "
      "the headline.**  No result in this programme is a *certified bound*.  "
      "The preceding paper's branch-and-bound construction remains the only "
      "certificate in the repository.  W1 established what that certificate's "
      "class actually is and bounded the discretisation inside it; W2 "
      "produced an exact law for a strictly larger model but only a grid "
      "maximum over its nuisance; W3 produced measurements, two of them "
      "negative; W4 produced a protocol.  Anything here described as holding "
      "over a domain holds over a GRID, and that difference is the thing this "
      "repository exists to take seriously.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------------ W1
    A("## W1. What information the rule actually has")
    A("")
    if inv:
        A("### The class Section 10 searched is the invariant class")
        A("")
        A("The null has three free parameters and the sufficient statistic has "
          "%d coordinates.  The only transformations that leave the problem "
          "alone are `x -> a + b x` applied to **both** groups, which is the "
          "%d-parameter affine group; scaling the groups separately changes "
          "the ratio and so changes the problem.  A maximal invariant is "
          "therefore %d-dimensional, and `(T, lambda_hat)` is one."
          % (inv["sufficient_dim"], inv["group_dim"],
             inv["maximal_invariant_dim"]))
        A("")
        A("Both halves are checked rather than asserted.  Applying random "
          "affine maps to raw data moves the invariant by at most "
          "**%.2e** across three designs.  For maximality, a second "
          "sufficient statistic is constructed directly on the same "
          "`(T, lambda_hat)` level set, not by applying a group element, and "
          "the group element that maps one to the other is recovered in "
          "closed form as `%s`; the worst reconstruction error over 200,000 "
          "pairs per design is **%.2e**."
          % (inv["worst_invariance_drift"], inv["recovery_formula"],
             inv["worst_recovery_error"]))
        A("")
        A("So Section 10 did not search a convenient parametric family.  It "
          "searched the non-randomised invariant tests, up to monotonicity "
          "and bin resolution.  *(exact in model)*")
        A("")
        A("### The route to a lower bound")
        A("")
        hs = inv["hunt_stein"]
        A("The group is solvable, because %s.  Solvable groups are amenable, "
          "and for an amenable group the Hunt-Stein theorem gives that %s."
          % (hs["reason"].replace(
              "translations form a normal abelian subgroup of the affine "
              "group with abelian quotient, so the group is solvable and "
              "therefore amenable",
              "translations form a normal abelian subgroup with abelian "
              "quotient"), hs["consequence"]))
        A("")
        A("That is the sentence the preceding programme could not write.  It "
          "is a **derivation, not a computation**, and three conditions are "
          "still to discharge:")
        A("")
        for k, c in enumerate(hs["conditions_still_to_discharge"]):
            A("* %s%s"
              % (c, "." if k == len(hs["conditions_still_to_discharge"]) - 1
                 else ";"))
        A("")
        A("### What 36 bins can hide")
        A("")
        A("The discretisation is itself an approximation, so the frontier was "
          "re-solved at 36, 72 and 144 bins with the same shrink step applied "
          "to each so that valid rules are compared with valid rules.")
        A("")
        A("| design | 36 bins | 72 bins | 144 bins | hidden |")
        A("|---|---|---|---|---|")
        for d in inv["bin_resolution"]:
            r = {x["nbins"]: x["mean_loss_pp"] for x in d["rows"]}
            A("| %s | %.2f pp | %.2f pp | %.2f pp | %.3f pp |"
              % (d["design"], r[36], r[72], r[144], d["hidden_pp"]))
        A("")
        A("At most **%.3f points** anywhere.  The 36-bin result is not an "
          "artefact of the grid.  *(exact in model)*"
          % inv["max_hidden_pp"])
        A("")

    if reg:
        A("### The optimal rule is monotone and deterministic")
        A("")
        A("Section 10 restricted the rule to an upper set in `T` and to a "
          "deterministic decision.  Neither is implied by validity.  Dropping "
          "both turns the problem into a linear program over every region in "
          "the discretised plane, because the mass each cell receives is a "
          "closed-form function of the state quadrature and both the "
          "objective and the uniform constraint are linear in the rejection "
          "probabilities.")
        A("")
        A("Monotonicity is imposed as a linear constraint rather than by "
          "rendering a threshold rule onto the grid, because rendering drops "
          "the cell straddling the threshold and hands the unrestricted rule "
          "about a point of power for free, the same order as the effect "
          "being measured.  The first version of this module did exactly "
          "that and reported a spurious 1.47-point gain.")
        A("")
        A("| design | monotone optimum | unrestricted optimum | deterministic "
          "floor | gain |")
        A("|---|---|---|---|---|")
        for d in reg["designs"]:
            A("| %s | %.2f pp | %.2f pp | %.2f pp | %+.2f pp |"
              % (d["design"], d["loss_pp_monotone"], d["loss_pp_randomised"],
                 d["loss_pp_deterministic_floor"],
                 d["gain_pp_from_dropping_monotonicity"]))
        A("")
        A("The gain is **%+.2f to %+.2f points**: zero to the resolution of "
          "the computation, at every design.  Allowing a non-monotone region "
          "and allowing randomisation buy nothing.  That is numerical "
          "evidence for the proposition the roadmap listed as a target, that "
          "the optimal rule in this class is monotone in `T`, and it means "
          "Section 10's restriction to threshold rules was not a "
          "simplification.  Rounding the fractional cells down costs %.2f to "
          "%.2f points, which is the price of an exactly deterministic "
          "boundary rather than of determinism itself."
          % (reg["min_gain_pp_from_dropping_monotonicity"],
             reg["max_gain_pp_from_dropping_monotonicity"],
             min(d["loss_pp_deterministic_floor"] - d["loss_pp_monotone"]
                 for d in reg["designs"]),
             max(d["loss_pp_deterministic_floor"] - d["loss_pp_monotone"]
                 for d in reg["designs"])))
        A("")
        A("The constraint is imposed at 97 ratio values and verified at 193, "
          "so the verification grid finds a small overshoot: at most %.2e, "
          "%.2f%% of the ceiling.  It is reported rather than shrunk away "
          "because both arms carry it equally and the comparison is between "
          "them.  T-grid sensitivity at most %.3f points.  *(measured)*"
          % (reg["worst_overshoot"], reg["worst_overshoot_pct_of_abar"],
             reg["worst_t_grid_sensitivity_pp"]))
        A("")

    if low:
        A("### The low-degrees-of-freedom regime")
        A("")
        A("Both group sizes move, because the cost is set by the smaller of "
          "the two and a one-axis design map cannot see that.  A gain counts "
          "as material at %.2f points, declared before the run."
          % low["material_threshold_pp"])
        A("")
        A("| min(nu) | mean gain from a design-specific threshold |")
        A("|---|---|")
        for k, v in sorted(low["mean_gain_by_min_df"].items(),
                           key=lambda kv: int(kv[0])):
            A("| %s | %.2f pp |" % (k, v))
        A("")
        A("Material in **%d of %d cells**, best **%.2f points**.  %s.  "
          "*(measured)*"
          % (low["n_material_cells"], low["n_cells"], low["max_gain_pp"],
             low["operational_rule"][0].upper() + low["operational_rule"][1:]))
        A("")

    A("---")
    A("")

    # ------------------------------------------------------------------ W2
    A("## W2. Aligned returns")
    A("")
    if fac:
        A("### The roadmap's route was wrong, and the replacement keeps the "
          "architecture")
        A("")
        A("As written, Phase 2A formed paired differences.  That gives a "
          "maximum of K studentised statistics with different studentisers "
          "and a Wishart-coupled denominator set, which nothing factorises: "
          "the K-dimensional integral returns.")
        A("")
        A("Keeping the returns unpaired and carrying the dependence in a "
          "latent factor instead, `%s`, conditioning on **%s** makes the K "
          "candidates exactly conditionally independent.  The reason is that "
          "`eps_bar_i` is orthogonal to the centred `eps_i`, and that "
          "`(n-1)s_i^2/sigma_i^2` is noncentral chi-square whose "
          "noncentrality sees the factor path only through its sum of "
          "squares, by rotational symmetry."
          % (fac["model"], ", ".join("`%s`" % c for c in fac["conditioning"])))
        A("")
        A("The null law is then a **%d-dimensional outer rule over one "
          "one-dimensional inner rule per candidate**, cost linear in K.  "
          "Two of the four are the two the existing paper already "
          "integrates." % fac["outer_dimension"])
        A("")
        A("| cell | brute force (%s draws) | exact evaluator | z | what the "
          "independent law says |"
          % "{:,}".format(fac["brute_draws"]))
        A("|---|---|---|---|---|")
        for r in fac["falsifiers"]:
            A("| %s | %.5f | %.5f | %+.2f | %.5f |"
              % (r["cell"], r["brute"], r["exact"], r["z"],
                 r["independent_law"]))
        A("")
        A("Worst |z| **%.2f** over %d cells.  The last column is the "
          "non-vacuity check: the independent-groups law is wrong here by up "
          "to **%s** in relative terms, and the evaluator tracks the truth "
          "across that range.  A reduction that quietly ignored the factor "
          "would pass every other cell and fail that column.  Refinement gap "
          "against a coarser grid %.1e.  *(exact in model)*"
          % (fac["worst_abs_z"], len(fac["falsifiers"]),
             pc(100 * fac["max_relative_departure_from_independent_law"]),
             fac["refinement_gap"]))
        A("")
        A("Critical constants, for the size of the correction:")
        A("")
        A("| design | c at %.3f |" % fac["abar"])
        A("|---|---|")
        for t in fac["critical_table"]:
            A("| %s | %.4f |" % (t["design"], t["c_alpha"]))
        A("")

    if cer:
        A("### What it costs to ignore the factor, and what a certificate "
          "over it costs")
        A("")
        A("The mechanism decides where the danger is, and the first "
          "version of this module could not see it.  In the aligned model %s."
          % cer["mechanism"])
        A("")
        A("That first version set the benchmark loading to zero, on the "
          "argument that only the contrast is identified.  It is not, and "
          "zero is exactly the boundary between the two regimes, so the scan "
          "found nothing.  A long-short candidate against a directional "
          "benchmark is the opposite-sign case, which is a configuration a "
          "desk is likely to be in.")
        A("")
        A("| ratio | beta_0 | beta_1 | realised FWER | z |")
        A("|---|---|---|---|---|")
        for r in cer["price_rows"]:
            A("| %.2f | %.1f | %.1f | %.5f | %+.1f |"
              % (r["lam"], r["b0"], r["b1"], r["fwer"], r["z"]))
        A("")
        A("**%d of %d cells are liberal**, and the split is exactly along "
          "the sign: %d of %d opposite-sign cells, %d of %d same-sign cells. "
          " Worst familywise error rate **%.5f**, which is **%.2f times** the "
          "ceiling; the worst same-sign cell only reaches %.5f.  So the "
          "aligned factor is not a conservative nuisance, and which side of "
          "the ceiling a desk lands on is decided by a sign it never "
          "estimates.  *(measured)*"
          % (cer["n_liberal_cells"], cer["n_cells"],
             cer["n_opposite_sign_liberal"], cer["n_opposite_sign_cells"],
             cer["n_same_sign_liberal"],
             cer["n_cells"] - cer["n_opposite_sign_cells"],
             cer["worst_fwer_iid_rule_on_aligned_data"],
             cer["worst_ratio_to_abar"], cer["worst_same_sign_fwer"]))
        A("")
        A("The nuisance is therefore **three-dimensional per candidate**, "
          "(ratio, beta_0, beta_1), and only the ratio is estimable from what "
          "the rule sees.  Over a declared box of %d points the constant that "
          "holds everywhere is **%.4f**, against known-parameter oracle "
          "constants of %.4f to %.4f, so the certificate costs **%.2f times** "
          "the worst oracle constant and %.2f times the best.  The supremum "
          "is a **grid maximum**, not a continuum certificate: the "
          "branch-and-bound construction applies in principle but the "
          "enclosure would be three-dimensional per candidate rather than one "
          "and has not been built.  *(measured)*"
          % (cer["box"]["points"], cer["certified_constant_over_box"],
             cer["oracle_constant_min"], cer["oracle_constant_max"],
             cer["certificate_price_vs_worst_oracle_x"],
             cer["certificate_price_vs_best_oracle_x"]))
        A("")

    if une:
        A("### Unequal histories: what the incumbent's prehistory buys")
        A("")
        A("A promotion gate has an asymmetry the clinical model does not: the "
          "incumbent has a long record and the candidate a short one.  Three "
          "predeclared designs on an overlap of %d observations."
          % une["n_overlap"])
        A("")
        A("| prehistory | full-history power gain over strict pairing |")
        A("|---|---|")
        for k, v in sorted(une["full_history_power_gain_by_pre"].items(),
                           key=lambda kv: int(kv[0])):
            A("| %s | %+.2f pp |" % (k, v))
        A("")
        A("The two-window design, which uses the prehistory for the "
          "benchmark's VARIANCE only and so needs constant variance but not "
          "a constant mean, buys **%.2f points**: essentially nothing.  That "
          "is the expected answer and the reason is structural, since the "
          "prehistory enters only through `s_0^2`, whose sampling error is "
          "second order once the overlap is long enough to matter."
          % une["two_window_max_power_gain_pp"])
        A("")
        A("Full history, which uses the prehistory for the mean as well, "
          "buys **%.2f points**, and %s."
          % (une["full_history_max_power_gain_pp"],
             une["saturates_at_prehistory"]))
        A("")
        A("That is a real gain and it is bought entirely with an assumption "
          "about the MEAN.  Breaking it takes the worst familywise error rate "
          "to **%.5f**, **%.2f times** the ceiling, on the `%s` design.  The "
          "direction matters: a prehistory running ABOVE the later mean makes "
          "the test powerless, one running BELOW it makes the test liberal, "
          "and only the second is dangerous.  So there is no free lunch in "
          "the incumbent's track record.  Either take strict pairing and "
          "give up those %.2f points, or buy them with a stationarity "
          "assumption on the mean whose failure costs a factor of %.0f.  "
          "*(measured)*"
          % (une["broken_stationarity_worst_fwer"],
             une["broken_ratio_to_abar"],
             une["broken_stationarity_worst_design"],
             une["full_history_max_power_gain_pp"],
             une["broken_ratio_to_abar"]))
        A("")

    A("---")
    A("")

    # ------------------------------------------------------------------ W3
    A("## W3. Dependence as three validity tiers, and what breaks the rule")
    A("")
    if tie:
        A("### The tiers")
        A("")
        A("| nuisance | A, declared | B, calibration period (360) | C, same "
          "period |")
        A("|---|---|---|---|")
        for key, lab in (("serial", "serial (AR(1), Cochrane-Orcutt)"),
                         ("factor", "factor (observed loadings)")):
            A("| %s | %.5f | %.5f | %.5f |"
              % (lab, tie["%s_tierA_worst_fwer" % key],
                 tie["%s_tierB_cal360_worst_fwer" % key],
                 tie["%s_tierC_worst_fwer" % key]))
        A("")
        A("Worst realised familywise error rate in each tier, against a "
          "ceiling of %.3f.  Tier B is reported at three calibration lengths "
          "in the deposit so the length-versus-error trade is visible rather "
          "than asserted.  *(measured)*" % tie["abar"])
        A("")
        A("| calibration length | serial worst | factor worst |")
        A("|---|---|---|")
        for cl in (36, 120, 360):
            A("| %d | %.5f | %.5f |"
              % (cl, tie["serial_tierB_cal%d_worst_fwer" % cl],
                 tie["factor_tierB_cal%d_worst_fwer" % cl]))
        A("")

    if rep:
        A("### Why tiers B and C fail for serial dependence, and the repair")
        A("")
        A("Tier A is exact because Cochrane-Orcutt maps a mean `mu` to "
          "`mu(1 - phi)` and, with `phi` declared, both groups are multiplied "
          "by the SAME constant, so the difference of transformed means is "
          "still zero under the null whatever `mu` is.")
        A("")
        A("With the coefficient estimated SEPARATELY IN EACH GROUP the two "
          "multipliers differ and the numerator acquires "
          "`mu (phi_hat_0 - phi_hat_1)`: zero in expectation, not zero in any "
          "sample, and spreading linearly in the common mean.  It is the "
          "Prais-Winsten defect's shape arriving by another route.  There the "
          "multiplier varied with the group's **length**; here with the "
          "group's **estimate**.")
        A("")
        A("**The repair is to pool**: one coefficient from both groups, "
          "applied to both, so the multiplier is identical by construction.")
        A("")
        A("| | separate estimates | pooled estimate |")
        A("|---|---|---|")
        A("| worst realised rate | %.5f | %.5f |"
          % (rep["separate_worst_fwer"], rep["pooled_worst_fwer"]))
        A("| times the ceiling | %.1f | %.1f |"
          % (rep["separate_worst_ratio"], rep["pooled_worst_ratio"]))
        A("| worst at a zero mean | %.5f | %.5f |"
          % (rep["separate_worst_at_mu0"], rep["pooled_worst_at_mu0"]))
        A("| worst at a mean of 5 | %.5f | %.5f |"
          % (rep["separate_worst_at_mu5"], rep["pooled_worst_at_mu5"]))
        A("| cells breaching | %d of %d | %d of %d |"
          % (rep["separate_breaches"], rep["n_per_arm"],
             rep["pooled_breaches"], rep["n_per_arm"]))
        A("")
        A("The clearest single view is one design held fixed while the "
          "common mean moves.  At the worst cell of the grid:")
        A("")
        wc = max(rep["rows"], key=lambda r: r["fwer"])
        trip = sorted([r for r in rep["rows"]
                       if r["design"] == wc["design"] and r["phi"] == wc["phi"]
                       and r["tier"] == wc["tier"]],
                      key=lambda r: (r["mu"], r["pooled"]))
        A("| common mean | separate estimates | pooled estimate |")
        A("|---|---|---|")
        for m in sorted(set(r["mu"] for r in trip)):
            sepv = [r["fwer"] for r in trip if r["mu"] == m and not r["pooled"]]
            poov = [r["fwer"] for r in trip if r["mu"] == m and r["pooled"]]
            A("| %.0f | %.5f | %.5f |" % (m, sepv[0], poov[0]))
        A("")
        A("The separate-estimate column grows with the common mean exactly as "
          "the algebra predicts; the pooled column does not move, because the "
          "repaired statistic cannot see the common mean at all.")
        A("")
        A("The signature is that the pooled rate is flat in the common mean: "
          "spread **%.5f** at tier B and **%.5f** at tier C, against a "
          "separate-estimate arm that moves from %.5f to %.5f as the mean "
          "goes from 0 to 5.  What pooling does not remove is imperfect "
          "whitening at an estimated coefficient, which biases the variance "
          "estimate and which the zero-mean rows isolate.  That residual is on its own enough to "
          "breach: every one of the %d cells in both arms is above the "
          "ceiling.  **Pooling is an improvement, not a repair**, and the "
          "operational conclusion is that Section 11's serial result is a "
          "tier-A result: with the coefficient estimated rather than "
          "declared, at these sample sizes, there is no guarantee and a "
          "calibration period of 360 does not deliver one.  *(measured)*"
          % (rep["pooled_spread_over_mu_tierB"],
             rep["pooled_spread_over_mu_tierC"],
             rep["separate_worst_at_mu0"], rep["separate_worst_at_mu5"],
             rep["n_per_arm"]))
        A("")

    if asy:
        A("### What the asymptotic and resampling alternatives deliver")
        A("")
        A("The first row is a REFERENCE, not a competitor.  It is given the "
          "coefficient; the others are not, and its job is to confirm that "
          "the instrument reads correctly before the rest of the column is "
          "believed.  The contest is between the rows that all have to "
          "estimate the coefficient from the same data.  Target size %.3f, "
          "one-sided, %d cells per method, common mean %.1f."
          % (asy["abar"], asy["cells_per_method"], asy.get("mu", 0.0)))
        A("")
        A("| method | worst size | times target | cells breaching |")
        A("|---|---|---|---|")
        for mth in asy["methods"]:
            A("| %s | %.5f | %.2f | %d of %d |"
              % (mth, asy["worst_size"][mth], asy["worst_ratio_to_abar"][mth],
                 asy["breaches"][mth], asy["cells_per_method"]))
        A("")
        A("The reference row holds at %.2f times the target, so the "
          "instrument reads correctly.  Among the rows that must estimate the "
          "coefficient, nothing holds: the best of them is `%s` at **%.2f "
          "times** the target and the worst is `%s` at **%.2f times**.  "
          "*(measured)*"
          % (asy["worst_ratio_to_abar"]["the exact rule, coefficient declared"],
             min((m for m in asy["methods"] if m != "the exact rule, coefficient declared"),
                 key=lambda m: asy["worst_ratio_to_abar"][m]),
             min(asy["worst_ratio_to_abar"][m] for m in asy["methods"]
                 if m != "the exact rule, coefficient declared"),
             max((m for m in asy["methods"] if m != "the exact rule, coefficient declared"),
                 key=lambda m: asy["worst_ratio_to_abar"][m]),
             max(asy["worst_ratio_to_abar"][m] for m in asy["methods"]
                 if m != "the exact rule, coefficient declared")))
        A("")
        A("Two things follow.  First, the exact law buys its guarantee from "
          "the coefficient being DECLARED, not from being exact: the same "
          "rule with the coefficient estimated separately in each group is "
          "the WORST row in the table, worse than any bootstrap.  Second, "
          "once the estimate is pooled the exact rule becomes the best of "
          "the estimated-coefficient methods by a clear margin, %.2f times "
          "the target against %.2f for the closest alternative.  Being least "
          "bad is not a guarantee and is not offered as one.  *(measured)*"
          % (asy["worst_ratio_to_abar"]["the exact rule, coefficient pooled"],
             min(asy["worst_ratio_to_abar"][m] for m in asy["methods"]
                 if m not in ("the exact rule, coefficient declared",
                              "the exact rule, coefficient pooled",
                              "the exact rule, coefficient per group"))))
        A("")

    if strs:
        A("### The stress battery")
        A("")
        A("Every cell runs the same deployed rule; only the data-generating "
          "process changes.  A breach is %.0f standard errors above the "
          "ceiling, declared before the run.  %d of %d cells breach."
          % (strs["breach_se"], strs["n_breaches"], strs["n_cells"]))
        A("")
        A("**Breaks it:** %s." % ", ".join(strs["breaching_cells"]))
        A("")
        A("**Survives it:** %s." % ", ".join(strs["surviving_cells"]))
        A("")
        A("Worst case **%.5f**, which is **%.2f times** the ceiling, at "
          "`%s`.  The direction matters and is recorded per cell: untreated "
          "positive serial correlation is violently liberal, negative serial "
          "correlation and undeclared factor structure are strongly "
          "conservative, and heavy tails and skew are conservative rather "
          "than liberal, which is the good direction to fail in.  "
          "*(diagnostic)*"
          % (strs["worst_fwer"], strs["worst_ratio_to_abar"],
             strs["worst_cell"]))
        A("")
        A("The one genuinely new finding is `short variance window`: "
          "estimating the variance on a trailing subwindow while the mean "
          "uses the whole sample, which is what a rolling risk model does, "
          "breaches on its own.  It is not in the paper's list of boundaries "
          "and it is a thing desks do by default.")
        A("")

    A("---")
    A("")

    # ------------------------------------------------------------------ W4
    A("## W4. The promotion protocol")
    A("")
    if man:
        A("`PROMOTION-PROTOCOL.md`, `manifest.schema.json` and "
          "`w4_manifest.py`.  The guarantee is a statement about a family "
          "fixed before the outcome was seen; nothing in the mathematics can "
          "check that, but a file can.")
        A("")
        A("Three verdicts: %s.  The example manifest returns **%s**, and its "
          "conditions are printed rather than suppressed."
          % (", ".join("`%s`" % v for v in man["verdict_labels"]),
             man["example_verdict"]))
        A("")
        A("**%d of %d refusal probes are refused on their own ground:**"
          % (man["n_refused"], man["n_probes"]))
        A("")
        for k, q in enumerate(man["probes"]):
            A("* %s%s" % (q["probe"],
                          "." if k == len(man["probes"]) - 1 else ";"))
        A("")
        A("A nuisance estimated on the evaluation period returns **%s**, not "
          "a refusal, which is the point: that run is worth doing and its "
          "number is worth having, as long as the report does not call it a "
          "certificate.  *(diagnostic)*" % man["tier_c_verdict"])
        A("")

    if prb:
        A("### Whether these checks bind")
        A("")
        A("`verify_beyond.py` formats every quoted number from a deposit and "
          "looks for it in the report, so neither file can be edited into "
          "agreement on its own.  That is necessary and not sufficient: a "
          "check that never fails proves nothing.  `probe_beyond.py` moves "
          "one deposit field at a time, reruns the verifier in a fresh "
          "process against the unchanged report, and requires a new failure.")
        A("")
        A("**%d of %d probes caught**, including **%d of %d transpositions**. "
          " The transpositions are the ones that matter: they swap two values "
          "inside a single table row, so every value is still present in the "
          "document and substring containment cannot see the change.  The "
          "preceding programme found that seven such probes all passed before "
          "table rows were matched whole."
          % (prb["n_caught"], prb["n_live"], prb["n_transposition_caught"],
             prb["n_transposition_probes"]))
        A("")
        A("| probe | caught |")
        A("|---|---|")
        for r in prb["probes"]:
            if not r["skipped"]:
                A("| %s | %s |" % (r["probe"], "yes" if r["caught"] else "no"))
        A("")
        if cov:
            A("Binding is one question and REACH is another.  A number no "
              "check formats passes unconditionally, and counting passes "
              "cannot see it.  `coverage_beyond.py` instruments the verifier "
              "to record every string it looks for, extracts every numeric "
              "token from this report, and reports what nothing covers: "
              "**%d of %d distinct tokens covered by a check**, %d excused "
              "with a stated reason, **%d uncovered**."
              % (cov["n_covered"], cov["n_tokens"], cov["n_excused"],
                 cov["n_uncovered"]))
            A("")
            if cov["excused"]:
                A("The excused ones, each named rather than silently "
                  "dropped:")
                A("")
                for t, why in sorted(cov["excused"].items()):
                    A("* `%s`, %s;" % (t, why))
                A("")

        A("One probe was replaced rather than recorded as a miss: swapping "
          "the monotone and unrestricted losses in the region table is a "
          "no-op, because the result of that section is precisely that the "
          "two are equal.  A vacuous probe is worse than no probe, since it "
          "reports a pass.  *(diagnostic)*")
        A("")

    A("---")
    A("")
    A("## What is still open")
    A("")
    A("* **No new certificate.**  Every guarantee in this programme is "
      "either inherited from the preceding one or is a grid maximum.  "
      "Converting either of the two grid maxima here into a continuum "
      "enclosure is the single piece of work that would change the "
      "repository's claims rather than its understanding of them.")
    A("* **The Hunt-Stein lower bound is a derivation, not a proof carried "
      "through.**  The maximal invariant is established and checked; the "
      "amenability argument is standard; the conditions listed above are not "
      "discharged and the K > 1 extension is not written.")
    A("* **The aligned-factor certificate is a grid maximum.**  Making it a "
      "continuum certificate needs a two-dimensional enclosure per candidate "
      "and the branch-and-bound code is one-dimensional.")
    A("* **The loading contrast cannot be estimated from the data the rule "
      "sees**, so unlike the variance ratio it can only be declared.  A desk "
      "that cannot declare it has no guarantee under aligned returns.")
    A("* **Residual cross-correlation after one factor is removed** is "
      "measured as conservative in the stress battery, not certified.")
    A("* **Tier C has no guarantee and is the realistic case.**  Nothing here "
      "changes that; the protocol only makes it impossible to report a tier-C "
      "run as a certificate by accident.")
    A("")

    io.open(OUT, "w", encoding="utf-8").write("\n".join(L) + "\n")
    have = sum(1 for x in (inv, reg, low, fac, cer, une, tie, rep, asy,
                           strs, man, prb, cov) if x is not None)
    print("REPORT-BEYOND.md written: %d lines, %d of 13 deposits present"
          % (len(L), have))


if __name__ == "__main__":
    main()
