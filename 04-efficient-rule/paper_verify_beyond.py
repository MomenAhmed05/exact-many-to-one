"""Check every number the W1/W3 sections quote against its deposit.

Same discipline as `paper_verify_frontier.py`, and the same helpers: each
expected string is FORMATTED FROM A DEPOSIT and then looked for in the built
document, so editing the paper cannot make a check pass and editing a deposit
cannot make one pass either.

Three kinds of check run:

  numeric        a quantity formatted from a deposit and found in the text;
  whole row      a table row matched end to end, because substring containment
                 cannot see a transposition inside a row;
  structural     the headings, the renumbering, the equation and table numbers,
                 and the typographic conventions the document follows, which is
                 where the first draft of this material failed: it used nine
                 capitalised lead-ins and a markdown bold in a document that
                 uses neither.

    python paper_verify_beyond.py [path-to-docx]
Env: DUMP_EXPECTED=file writes every expected string, for coverage_beyond.
"""

from __future__ import annotations

import io
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WF = os.path.join(HERE, "deposits")
DEFAULT = os.path.join(HERE, "paper",
                       "Exact-Many-to-One-Under-Estimated-Volatilities.docx")

sys.path.insert(0, HERE)
from paper_verify_frontier import doc_text, table_rows  # noqa: E402

ok = True
npass = nfail = 0
FAILED = []
EXPECTED = []


def load(name):
    with open(os.path.join(WF, name), encoding="utf-8") as fh:
        return json.load(fh)


def chk(label, cond, detail=""):
    global ok, npass, nfail
    cond = bool(cond)
    ok &= cond
    npass += cond
    nfail += not cond
    if not cond:
        FAILED.append(label)
    print(("  OK   " if cond else "  FAIL ") + label
          + ("   " + detail if detail else ""))


def main(path=None):
    path = path or DEFAULT
    TXT = doc_text(path)
    ROWS = table_rows(path)

    def has(label, value, fmt="%s"):
        v = fmt % value
        EXPECTED.append(v)
        chk("%s -> %s" % (label, v), v in TXT)

    def row(label, cells):
        cells = [str(c) for c in cells]
        EXPECTED.append(" | ".join(cells))
        hit = any(list(cells) == list(r[:len(cells)]) for r in ROWS)
        chk("row: %s" % label, hit, "" if hit else " | ".join(cells)[:95])

    inv = load("w1_invariant.json")
    reg = load("w1_region.json")
    low = load("w1_lowdf.json")
    fro = load("wc_frontier.json")
    tie = load("w3_tiers.json")
    rep = load("w3_serial_repair.json")
    asy = load("w3_asymptotic.json")
    strs = load("w3_stress.json")
    fac = load("w2_factor.json")
    cer = load("w2_cert.json")
    une = load("w2_unequal.json")
    idn = load("w2_identify.json")

    fd = {d["design"]: d for d in fro["designs"]}
    rd = {d["design"]: d for d in reg["designs"]}
    DF = ["df3", "df5", "df7", "df11", "df19", "df35"]
    LAB = {"df3": "3 df", "df5": "5 df", "df7": "7 df", "df11": "11 df",
           "df19": "19 df", "df35": "35 df"}

    def mean_pp(key, d):
        v = fd[d][key]
        return sum(v) / len(v)

    def gap(d):
        return mean_pp("plugin_loss_pp", d) - rd[d]["loss_pp_monotone"]

    # ------------------------------------------------------------ structure
    print("STRUCTURE")
    for h in ("10.8 The floor the class reaches",
              "10.9 Which group binds, and a map for planning",
              "10.10 What a declared ratio domain is worth",
              "10.11 Cost, and what checks the numbers",
              "11.4 Three validity tiers, and the effect of estimating the "
              "coefficient",
              "11.5 What the alternatives deliver at these sample sizes",
              "11.6 The factor case without residualising"):
        chk("heading: %s" % h[:44], h in TXT)
    for n in (28, 29, 30, 31, 32, 33, 34):
        chk("Table %d present" % n, "Table %d." % n in TXT)
    for n in (18, 22, 23):
        chk("Figure %d present" % n, "Figure %d." % n in TXT)
    # the three that Section 10.8 displaced must have moved, not vanished
    for n in (19, 20, 21):
        chk("displaced Figure %d present" % n, "Figure %d." % n in TXT)
    for n in (47, 48, 49, 50):
        chk("equation (%d) present" % n, "(%d)" % n in TXT)
    chk("old 10.8 heading is gone",
        "10.8 Which group binds" not in TXT)
    chk("no stale reference to the old 10.8 map",
        "Section 10.8 shows that the binding quantity" not in TXT)

    print("\nREGISTER")
    caps = re.findall(r"(?<![A-Za-z])[A-Z]{4,}(?![A-Za-z])", TXT)
    stray = sorted(set(c for c in caps
                       if c not in ("FWER", "ANOVA", "JSON", "METHOD",
                                    "GARCH", "SLSQP", "HTML", "AR")))
    chk("no capitalised emphasis", not stray, ", ".join(stray[:8]))
    chk("no markdown bold", "**" not in TXT)
    chk("no em dash", "—" not in TXT)
    chk("no spelled-out per cent", "per cent" not in TXT.lower())
    for w in ("caveat", "the twist", "uncomfortable truth"):
        chk("no filler: %s" % w, w not in TXT.lower())

    # ------------------------------------------------------------- 10.8
    print("\n10.8 THE FLOOR")
    has("invariance drift", inv["worst_invariance_drift"], "%.2e")
    has("recovery error", inv["worst_recovery_error"], "%.2e")
    # the paper spells small numbers, so the check has to as well, and it has
    # to fail if the deposit's dimension count ever stops being four
    WORD = {2: "two", 3: "three", 4: "four"}
    chk("sufficient statistic dimension agrees with the deposit",
        ("%s coordinates" % WORD[inv["sufficient_dim"]]) in TXT)
    chk("group dimension agrees with the deposit",
        ("%s parameters" % WORD[inv["group_dim"]]) in TXT)
    chk("invariant dimension agrees with the deposit",
        ("%s-dimensional" % WORD[inv["maximal_invariant_dim"]]) in TXT)
    has("max hidden", inv["max_hidden_pp"], "%.3f percentage points")
    has("gap low", min(gap(d) for d in ("df11", "df19", "df35")),
        "%.2f to")
    has("gap high", max(gap(d) for d in ("df11", "df19", "df35")),
        "%.2f percentage points")
    has("overshoot", reg["worst_overshoot"], "%.2e")
    has("overshoot share", reg["worst_overshoot_pct_of_abar"], "%.2f%%")
    has("constraint grid", reg["n_constraint_grid"], "at %d ratio values")
    has("verify grid", reg["n_verify_grid"], "verified at %d")
    has("pairs per design", "{:,}".format(inv["pairs_per_design"]),
        "over %s pairs per design")
    dfl = [x["loss_pp_deterministic_floor"] - x["loss_pp_monotone"]
           for x in reg["designs"]]
    has("det floor low", min(dfl), "gives up %.2f to")
    has("det floor high", max(dfl), "%.2f percentage points")
    for d in ("df35", "df3"):
        r = next(x for x in inv["bin_resolution"] if x["design"] == d)
        for q in r["rows"]:
            has("bins %s %d" % (d, q["nbins"]), q["mean_loss_pp"], "%.2f")
        has("bin count %s" % d, r["rows"][-1]["nbins"], "and %d with")
    has("t-grid", reg["worst_t_grid_sensitivity_pp"],
        "%.3f percentage points")
    for d in DF:
        row("Table 28 %s" % d,
            [LAB[d], "%.2f pp" % mean_pp("plugin_loss_pp", d),
             "%.2f pp" % mean_pp("optimal_loss_pp", d),
             "%.2f pp" % rd[d]["loss_pp_monotone"],
             "%.2f pp" % rd[d]["loss_pp_randomised"],
             "%.2f pp" % gap(d)])
    chk("monotonicity is not binding, to six decimals",
        all(abs(x["gain_pp_from_dropping_monotonicity"]) < 1e-6
            for x in reg["designs"]))
    chk("the paper says so rather than printing a signed zero",
        "agree to six decimal places" in TXT
        and "+0.00 percentage points" not in TXT)
    # The extension is now taken, so these guard the claim actually made.
    # The trap the section has to avoid is letting the theorem appear to
    # cover the mean criterion Table 28 reports, which it does not.
    mmx = load("w1_minimax.json")
    hyp = load("w1_hypotheses.json")
    bfx = load("w1_baseline_fixed.json")
    chk("the step is taken rather than deferred",
        "That step is taken here" in TXT)
    chk("the criterion is named as worst-case regret",
        "Worst-case regret against the known-ratio oracle" in TXT)
    chk("the average criterion is excluded by name",
        "not the quantity the theorem covers" in TXT)
    chk("the reason the average criterion is excluded is given",
        "bounds a supremum but not a weighted sum" in TXT)
    chk("C4 is flagged as taken from the literature",
        "taken from the literature rather than proved here" in TXT)
    chk("its references are named",
        "Bondar and Milnes 1981" in TXT and "Lehmann and Romano 2005" in TXT)
    chk("the open extension is still declared",
        "extension above one candidate is still unwritten" in TXT)
    gp = {r["design"]: r["gap_deployed_above_floor_pp"]
          for r in mmx["designs"]}
    for d, w in (("df11", "above the floor at eleven"),
                 ("df19", "at nineteen"), ("df35", "at thirty-five")):
        has("minimax gap %s" % d, gp[d], "%.2f percentage points " + w
            if d == "df11" else "%.2f " + w)
    has("probe count", hyp["n_rules"] + 12, "%d rules that read")
    has("best probe shortfall", abs(hyp["worst_gain_pp"]),
        "falls %.2f percentage points short")
    has("fixed baseline under widening", bfx["worst_regret_pp"],
        "worst regret at %.6f")
    chk("the baseline is stated as unchanged in every digit",
        "unchanged in every digit" in TXT)
    for k, v in low["mean_gain_by_min_df"].items():
        row("Table 34 min nu %s" % k, [k + " df", "%.2f pp" % v])
    has("material cells", low["n_material_cells"], "material in %d of the")
    has("max lowdf gain", low["max_gain_pp"], "%.2f percentage points")
    has("material bar", low["material_threshold_pp"],
        "%.2f percentage points, declared before the run")
    chk("Tables 28 and 34 are declared non-comparable",
        "not to be read across" in TXT)
    has("df7 gap on the two bases", gap("df7"), "the two give %.2f and")
    has("df7 lowdf mean", low["mean_gain_by_min_df"]["7"],
        "%.2f percentage points. The difference")

    # ------------------------------------------------------------- 11.4
    print("\n11.4 VALIDITY TIERS")
    for key, lab in (("serial", "AR(1) coefficient"),
                     ("factor", "factor loadings")):
        row("Table 29 %s" % key,
            [lab, "%.4f" % tie["%s_tierA_worst_fwer" % key]]
            + ["%.4f" % tie["%s_tierB_cal%d_worst_fwer" % (key, c)]
               for c in (36, 120, 360)]
            + ["%.4f" % tie["%s_tierC_worst_fwer" % key]])
    has("serial tier C ratio", tie["serial_tierC_worst_fwer"] / tie["abar"],
        "%.1f times the ceiling")
    has("tier draws", "{:,}".format(tie["draws"]), "%s draws per cell")
    has("repair draws", "{:,}".format(rep["draws"]), "%s draws per cell")
    has("calibration length", 360, "at %d observations")
    has("pooled ratio", rep["pooled_worst_ratio"], "to %.1f, and every")
    wc = max(rep["rows"], key=lambda r: r["fwer"])
    trip = [r for r in rep["rows"] if r["design"] == wc["design"]
            and r["phi"] == wc["phi"] and r["tier"] == wc["tier"]]
    for m in sorted(set(r["mu"] for r in trip)):
        row("Table 30 mean %.0f" % m,
            ["%.0f" % m,
             "%.4f" % next(r["fwer"] for r in trip
                           if r["mu"] == m and not r["pooled"]),
             "%.4f" % next(r["fwer"] for r in trip
                           if r["mu"] == m and r["pooled"])])
    has("separate worst", rep["separate_worst_fwer"], "%.4f")
    has("pooled worst", rep["pooled_worst_fwer"], "%.4f")
    has("separate ratio", rep["separate_worst_ratio"], "%.1f times the")
    has("pooled spread", max(rep["pooled_spread_over_mu_tierB"],
                             rep["pooled_spread_over_mu_tierC"]), "%.4f")
    chk("pooling not called a repair",
        "an improvement and not a repair" in TXT)
    chk("the narrowing of 11.1 is stated", "is a tier-A result" in TXT)

    # ------------------------------------------------------------- 11.5
    print("\n11.5 THE ALTERNATIVES")
    for m in asy["methods"]:
        row("Table 31 %s" % m[:28],
            [m, "%.4f" % asy["worst_size"][m],
             "%.2f" % asy["worst_ratio_to_abar"][m],
             "%d of %d" % (asy["breaches"][m], asy["cells_per_method"])])
    has("reference row", asy["worst_ratio_to_abar"][
        "the exact rule, coefficient declared"], "%.2f times the target")
    chk("the reference is labelled as one",
        "a reference and not a competitor" in TXT)
    chk("no ASCII phi reaches the page", "phi" not in TXT)
    has("closed-form draws", "{:,}".format(asy["draws_fast"]),
        "%s draws per cell")
    has("bootstrap draws", "{:,}".format(asy["draws_boot"]), "and %s with")
    has("bootstrap reps", asy["boot_reps"], "with %d resamples")
    has("common mean", asy["mu"], "common mean is %.1f")

    # ------------------------------------------------------------- 11.6
    print("\n11.6 THE FACTOR CASE")
    has("outer dimension", fac["outer_dimension"], "%d-dimensional")
    has("worst z", fac["worst_abs_z"], "%.2f standard errors")
    has("departure", 100 * fac["max_relative_departure_from_independent_law"],
        "%.1f%%")
    has("refinement", fac["refinement_gap"], "%.1e")
    for r in fac["falsifiers"]:
        row("Table 32 %s" % r["cell"][:26],
            [r["cell"], "%.5f" % r["brute"], "%.5f" % r["exact"],
             "%+.2f" % r["z"], "%.5f" % r["independent_law"]])
    for r in cer["price_rows"]:
        row("Table 33 %.2f %.1f %.1f" % (r["lam"], r["b0"], r["b1"]),
            ["%.2f" % r["lam"], "%.1f" % r["b0"], "%.1f" % r["b1"],
             "%.5f" % r["fwer"], "%+.1f" % r["z"]])
    has("liberal cells", cer["n_liberal_cells"], "%d of the")
    has("worst aligned", cer["worst_fwer_iid_rule_on_aligned_data"], "%.5f")
    has("worst ratio", cer["worst_ratio_to_abar"], "%.2f times the ceiling")
    has("same-sign worst", cer["worst_same_sign_fwer"], "%.5f")
    has("box points", cer["box"]["points"], "box of %d points")
    has("brute draws", "{:,}".format(fac["brute_draws"]), "%s draws per cell")
    has("price vs best", cer["certificate_price_vs_best_oracle_x"],
        "%.2f times the")
    has("certified constant", cer["certified_constant_over_box"], "%.4f")
    has("oracle min", cer["oracle_constant_min"], "%.4f")
    has("oracle max", cer["oracle_constant_max"], "%.4f")
    chk("grid maximum labelled", "a grid maximum and not a continuum" in TXT)
    # the loadings being outside the RULE must not be written as the loadings
    # being unknowable: they are identified above one candidate
    has("sampling error", idn["sampling"][0]["median_abs_error"],
        "median absolute error of %.3f")
    has("at n", idn["sampling"][0]["n"], "%d observations")
    has("parameters", idn["counts"][0]["parameters"], "%d parameters against")
    has("entries", idn["counts"][0]["covariance_entries"],
        "%d distinct covariance entries")
    chk("the claim is about the rule, not about the data",
        "a function of what the rule sees" in TXT
        and "the loadings are recoverable from the covariance" in TXT)
    chk("the K=1 boundary is identification, not method",
        "the model is under-identified" in TXT)
    chk("no claim that the loadings can only be declared",
        "can only be declared" not in TXT
        and "declared and bounded and nothing more" not in TXT)
    chk("9.10's open item is closed, not re-claimed",
        "Section 9.10 names the reduction" in TXT)
    chk("9.10 no longer calls the general case open",
        "the general case remains open" not in TXT)
    chk("9.10 and 11.6 agree on what 9.10 reported",
        "reports it as implemented at one candidate only" in TXT
        and "carries the construction out for any number of candidates" in TXT)
    chk("11.6 does not re-claim the reduction as new",
        "names the reduction that supplies it" in TXT)
    has("two-window gain", une["two_window_max_power_gain_pp"],
        "%.2f percentage points")
    has("full-history gain", une["full_history_max_power_gain_pp"],
        "%.2f percentage points")
    has("broken worst", une["broken_stationarity_worst_fwer"], "%.5f")
    has("broken ratio", une["broken_ratio_to_abar"], "%.0f times the ceiling")

    # ------------------------------------------------------------ additions
    print("\nABSTRACT, CONTRIBUTIONS, LIMITATIONS, CONCLUSION")
    chk("contribution (14)", "(14) The efficiency result placed in its class"
        in TXT)
    chk("contribution (15)", "(15) The declared nuisance shown to be "
        "load-bearing" in TXT)
    chk("contribution (16)", "(16) The general factor reduction implemented"
        in TXT)
    has("abstract tier B", tie["serial_tierB_cal360_worst_fwer"], "%.4f")
    sv = max(r["fwer"] for d in strs["designs"] for r in d["rows"]
             if r["cell"] == "short variance window")
    has("rolling window breach", sv, "%.4f")
    has("stress cells", strs["n_cells"], "%d cells in the stress battery")
    has("stress breaches", strs["n_breaches"], "%d breach")
    chk("the limitation of no certificate is stated",
        "supplies no certificate" in TXT)

    dump = os.environ.get("DUMP_EXPECTED")
    if dump:
        with io.open(os.path.join(HERE, dump), "w", encoding="utf-8") as fh:
            for e in EXPECTED:
                fh.write(e + chr(10))
        print("\nwrote %d expected strings to %s" % (len(EXPECTED), dump))

    print("\n%d checks passed, %d failed" % (npass, nfail))
    if FAILED:
        print("failures:")
        for f in FAILED[:30]:
            print("   " + f)
    print("PAPER UPDATE VERIFIED" if ok else "PAPER UPDATE HAS FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
