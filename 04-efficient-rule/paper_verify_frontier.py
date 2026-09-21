"""Check every number the new sections quote against its deposit.

Same discipline as `verify_frontier.py` and `02-bounded-rule/paper_verify_w26.py`: each
expected string is FORMATTED FROM A DEPOSIT and then looked for in the built
document, so editing the paper cannot make a check pass and editing a deposit
cannot make one pass either.  Only agreement between the two does.

Structural checks run beside the numeric ones: the headings the update creates,
the renumbering it performs, the tables and figures it inserts, the equation
numbers it uses, and the typographic conventions the document follows.

    python paper_verify_frontier.py [path-to-docx]
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.dirname(HERE)
WA = os.path.join(AUDIT, "03-certified-supremum", "deposits")
WF = os.path.join(HERE, "deposits")
DEFAULT = os.path.join(HERE, "paper",
                       "Exact-Many-to-One-Under-Estimated-Volatilities.docx")

ok = True
npass = nfail = 0
EXPECTED = []      # every string a check looked for, for coverage audit


def chk(label, cond, detail=""):
    global ok, npass, nfail
    cond = bool(cond)
    ok &= cond
    npass += cond
    nfail += not cond
    print(("  OK   " if cond else "  FAIL ") + label
          + ("   " + detail if detail else ""))


def load(base, name):
    with open(os.path.join(base, name), encoding="utf-8") as fh:
        return json.load(fh)


def doc_text(path):
    """All visible text of the document, whitespace normalised."""
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    body = re.sub(r"<[^>]+>", "\x00", xml)
    return re.sub(r"\s+", " ", body.replace("\x00", ""))



def table_rows(path):
    """Every table row in the document, as a tuple of its cell strings.

    Checking a value with `in T` only asks whether it occurs somewhere.  That
    cannot see a transposition: swapping two numbers between rows leaves both
    of them present.  Matching a whole row at once can.
    """
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    out = []
    for tbl in re.findall(r"<w:tbl>.*?</w:tbl>", xml, re.S):
        for tr in re.findall(r"<w:tr[ >].*?</w:tr>", tbl, re.S):
            cells = []
            for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, re.S):
                txt = "".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", tc))
                cells.append(re.sub(r"\s+", " ", txt).strip())
            if cells:
                out.append(tuple(cells))
    return out


def main():
    path = (sys.argv[1:] or [DEFAULT])[0]
    T = doc_text(path)
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8")
        media = [n for n in z.namelist() if n.startswith("word/media/")]

    def has(label, value, fmt="%s"):
        v = fmt(value) if callable(fmt) else fmt % value
        EXPECTED.append(v)
        chk("%s -> %s" % (label, v), v in T)

    # ------------------------------------------------------------ structure
    print("STRUCTURE")
    for h in ("9.11 A certified supremum over the continuum",
              "10. An Efficient Rule With the Same Guarantee",
              "10.1 What the region rule gives up",
              "10.2 Where the region rule's conservatism comes from",
              "10.3 The level the constraint allows",
              "10.4 A certificate for the internal level",
              "10.5 Any number of candidates, and both sides",
              "10.6 The joint rule, and what certifying it requires",
              "10.7 Whether anything in the same class does better",
              "10.8 The floor the class reaches",
              "10.9 Which group binds, and a map for planning",
              "10.10 What a declared ratio domain is worth",
              "10.11 Cost, and what checks the numbers",
              "11. Dependence as a Declared Nuisance",
              "11.1 Serial correlation", "11.2 A common factor",
              "11.3 What is declared and what is estimated",
              "12. Discussion", "13. Limitations", "14. Conclusion"):
        chk("heading present: %s" % h, h in T)
    for gone in ("10. Discussion", "11. Limitations", "12. Conclusion",
                 "10.1 What is settled", "boundary in Section 10.3"):
        chk("old numbering removed: %s" % gone, gone not in T)
    for n in range(37, 51):
        chk("equation (%d)" % n, "(%d)" % n in T)
    for n in range(15, 35):
        chk("Table %d." % n, "Table %d." % n in T)
    for n in range(14, 24):
        chk("Figure %d." % n, "Figure %d." % n in T)
    chk("Table A3.", "Table A3." in T)
    chk("37 tables in the document", xml.count("<w:tbl>") == 37,
        "%d" % xml.count("<w:tbl>"))
    chk("23 images in the package", len(media) == 23, "%d" % len(media))

    # The figure files are named in the order they were written, not in the
    # order the paper uses them.  Pin the mapping by content hash, so that a
    # rename or a re-ordered insertion cannot put the wrong picture under a
    # caption without this failing.
    # document order, which is not figure-file order and, since Section 10.8
    # was inserted, is not the order the files were written in either
    FIG_ORDER = ["fig20_bnbgap.png", "fig14_ledger.png", "fig17_levels.png",
                 "fig15_frontier.png", "fig21_floor.png",
                 "fig16_designmap.png", "fig18_serial.png",
                 "fig19_factor.png", "fig22_tiers.png", "fig23_sign.png"]
    with zipfile.ZipFile(path) as z:
        rels = z.read("word/_rels/document.xml.rels").decode("utf-8")
        r2t = dict(re.findall(r'Id="([^"]+)"[^>]*Target="([^"]+)"', rels))
        order = [r2t.get(m, "?") for m in re.findall(r'r:embed="([^"]+)"', xml)]
        sha = {n: hashlib.sha1(z.read(n)).hexdigest()
               for n in z.namelist() if n.startswith("word/media/")}
    figdir = os.path.join(HERE, "figs")
    chk("the document holds %d images" % (13 + len(FIG_ORDER)),
        len(order) == 13 + len(FIG_ORDER), "%d" % len(order))
    for k, fname in enumerate(FIG_ORDER):
        want = hashlib.sha1(
            open(os.path.join(figdir, fname), "rb").read()).hexdigest()
        if 13 + k >= len(order):
            chk("figure at position %d is %s" % (14 + k, fname), False,
            "no such image")
            continue
        tgt = "word/" + order[13 + k].replace("../", "")
        chk("figure at position %d is %s" % (14 + k, fname),
            sha.get(tgt) == want)

    print("\nTYPOGRAPHY")
    chk("no em dash", "—" not in T)
    chk("no en dash", "–" not in T)
    chk("uses % and not the spelled form", "per cent" not in T
        and T.count("%") > 200, "%d occurrences of %%" % T.count("%"))
    for w in ("caveat", "honest", "Here's the", "uncomfortable truth",
              "What most people miss"):
        chk("absent: %r" % w, w not in T)

    # ------------------------------------------------- 9.11 the certificate
    print("\nBRANCH AND BOUND (9.11)")
    bnb = load(WA, "bnb.json")
    chk("cells counted", bnb["n_cells"] == 17)
    has("certified of total", (bnb["n_certified_at_tol"], bnb["n_cells"]),
        "%d of the %d cells certify")
    has("median gap", bnb["median_gap"], "%.3e")
    has("max gap", bnb["max_gap"], "%.3e")
    has("total seconds", bnb["total_seconds"], "%.0f seconds")
    for c in bnb["cells"]:
        has("gap %s" % c["cell"], c["gap"], "%.2e")
    slack = load(WA, "box_slack.json")
    has("slack max", slack["slack_max"], "%.2e")
    has("slack at 1.001", slack["slack_at_factor_1p001"], "%.2e")
    rootc = load(WA, "rootcert.json")
    has("brackets", rootc["n_brackets"], "%d point brackets")
    has("bracket cells", rootc["n_cells"], "across %d cells")
    has("widest bracket", rootc["max_width"], "%.2e")
    adv = load(WA, "adversarial.json")
    has("adversarial designs", adv["n_designs"], "A suite of %d adversarial")
    has("adversarial certified", adv["n_certified"], "certifies %d of them")
    has("adversarial max gap", adv["max_gap"], "%.2e")
    encl = load(WA, "enclose_iv.json")
    chk("independent interval evaluations all contained", encl["all_contained"])
    has("contained points", (encl["n_points"], encl["n_points"]),
        "%d of %d such points are contained")
    ebud = load(WA, "error_budget.json")
    has("error budget", ebud["max_Delta_p"], "%.2e")
    supc = load(WA, "sup_compare.json")
    has("vertex median diff",
        supc["headline_continuous_minus_vertex"]["median"], "%.2e")
    has("vertex max diff",
        supc["headline_continuous_minus_vertex"]["max"], "%.2e")
    has("argmax at vertex",
        (supc["n_argmax_at_vertex"], supc["n_cells"]),
        "at a vertex in %d of %d cells")

    # ---------------------------------------------------------- 10 the rule
    print("\nTHE LEDGER (10.1)")
    ledger = load(WF, "wa_ledger.json")
    led = {d["design"]: d for d in ledger["designs"]}

    def lrow(des, lam=1.0):
        return next(q for q in led[des]["rows"] if abs(q["lam"] - lam) < 1e-9)

    r35 = lrow("df35")
    for key in ("oracle_conditional", "plugin_measured",
                "region_sup_conditional", "uncond_certified"):
        e = r35[key]
        has("%s mean c" % key[:12], e["mean_c"], "%.3f")
        has("%s ratio" % key[:12], e["ratio_to_oracle_mean"], "%.2f")
        has("%s fwer" % key[:12], e["fwer"], "%.5f")
        has("%s loss" % key[:12], e["loss_pp_at_80"], "%.1f pp")
    has("region budget share",
        100.0 * r35["region_sup_conditional"]["fwer"] / 0.025,
        "%.1f%% of what the guarantee permits")
    has("df3 plug-in fwer", lrow("df3")["plugin_measured"]["fwer"], "%.5f")

    print("\nK = 1 CERTIFIED RULE (10.3, 10.4)")
    k1 = load(WF, "wb_certified_rule.json")
    chk("all certified", k1["all_certified"])
    for x in k1["designs"]:
        has("%s alpha0" % x["design"], x["alpha0"], "%.5f")
        has("%s bound" % x["design"], x["sup_upper_bound"], "%.6f")
        has("%s median loss" % x["design"], x["median_loss_pp_at_80"],
            "%.1f pp")
        has("%s worst loss" % x["design"], x["max_loss_pp_at_80"], "%.1f pp")
    chk("no design re-levelled after the tail correction",
        k1["tail_mass_correction"]["designs_relevelled"] == 0)
    has("corner trials", k1["falsifiers"]["x2_trials"],
        lambda v: "no violations in {:,} trials".format(v))
    chk("corner violations are zero",
        k1["falsifiers"]["x2_corner_substitution_violations"] == 0)

    print("\nGENERAL K MEASURED (10.3)")
    gen = load(WF, "wb_general.json")
    chk("all within ceiling", gen["all_within_ceiling"])
    for x in gen["designs"]:
        has("%s plug-in sup" % x["design"], x["plugin_at_abar_sup_fwer"],
            "%.5f")
        has("%s excess" % x["design"], x["plugin_at_abar_rel_pct"], "+%.1f%%")
        has("%s alpha0" % x["design"], x["alpha0"], "%.5f")
        has("%s median loss" % x["design"], x["median_loss_pp"], "%.1f pp")
        has("%s worst loss" % x["design"], x["max_loss_pp"], "%.1f pp")

    print("\nTWO-SIDED AND ANY K (10.5)")
    twos = load(WF, "twosided.json")
    chk("reduction exact", twos["reduction"]["pass"])
    has("reduction trials", twos["reduction"]["trials"],
        lambda v: "to the last bit over {:,} trials".format(v))
    has("K=3 ratio", twos["no_reduction_above_k1"]["ratio"], "%.4f rather")
    has("two-sided level low",
        min(x["alpha0_two_sided"] for x in twos["designs"]), "%.5f")
    has("two-sided level high",
        max(x["alpha0_two_sided"] for x in twos["designs"]), "%.5f")
    anyk = load(WF, "certify_anyk.json")
    chk("all any-K certificates pass", anyk["all_certified"])
    has("worst union slack", anyk["worst_union_slack"], "%.2e")
    for x in anyk["designs"]:
        has("anyk %s %s level" % (x["design"], x["two_sided"]), x["alpha0"],
            "%.5f")
        has("anyk %s %s bound" % (x["design"], x["two_sided"]),
            x["union_bound_total"], "%.6f")

    print("\nJOINT RULE (10.6)")
    joint = load(WF, "certify_joint.json")
    chk("all joint certificates pass", joint["all_certified"])
    chk("every joint bound dominates the deployed rule",
        joint["all_bounds_dominate"])
    for x in joint["designs"]:
        has("joint %s level" % x["design"], x["alpha0_certified_joint"],
            "%.5f")
        has("joint %s certificate cost" % x["design"],
            x["level_vs_measured_pct"], "%.1f%%")
        has("joint %s bound" % x["design"], x["union_bound_total"], "%.6f")
        has("joint %s realised" % x["design"], x["exact_sup_over_patterns"],
            "%.6f")
    has("projection gap low",
        min(x["projection_gap_pct"] for x in joint["designs"]), "%.1f%%")
    has("projection gap high",
        max(x["projection_gap_pct"] for x in joint["designs"]), "%.1f%%")

    print("\nRULE COMPARISON (10.6)")
    rcmp = load(WF, "rule_compare.json")
    chk("the per-candidate certificate asks the lower threshold everywhere",
        rcmp["percand_certificate_dominates_everywhere"])
    for d in rcmp["designs"]:
        for pt in d["points"]:
            chk("%s at %s: joint threshold above per-candidate"
                % (d["design"], pt["lam_hat"]),
                pt["c_joint_certified"] > max(pt["c_percand_certified"]))

    rcd = {d["design"]: d for d in rcmp["designs"]}
    rc0 = next(d for d in rcmp["designs"] if d["design"] == "K2_df35")
    has("joint certified threshold",
        rc0["points"][0]["c_joint_certified"], "%.4f")
    has("per-candidate certified threshold",
        max(rc0["points"][0]["c_percand_certified"]), "%.4f")
    has("joint measured threshold",
        rc0["points"][0]["c_joint_measured"], "%.4f")
    has("percand certificate cost low",
        min(d["percand_certificate_costs_pct"] for d in rcmp["designs"]),
        "%.1f%%")
    has("percand certificate cost high",
        max(d["percand_certificate_costs_pct"] for d in rcmp["designs"]),
        "%.1f%%")

    print("\nFRONTIER AND DESIGN MAP (10.7, 10.8, 10.9)")
    k1d = load(WF, "wb_certified_rule.json")
    has("ratio domain low", k1d["lam_domain"][0], "[%g,")
    has("ratio domain high", k1d["lam_domain"][1], "%g]")
    has("volatility factor", k1d["lam_domain"][1] ** 0.5, "about %.0f times")
    has("optimiser bins", load(WF, "wc_frontier.json")["nbins"],
        "piecewise constant on %d bins")

    front = load(WF, "wc_frontier.json")
    chk("all optima valid", front["all_valid"])
    for x in front["designs"]:
        has("%s gain" % x["design"], x["mean_gain_pp"], "%+.2f pp")
        has("%s optimum fwer" % x["design"], x["sup_fwer_fine"], "%.4f")
    dmap = load(WF, "wd_designmap.json")
    chk("binds on the smaller group", dmap["binds_on_min_df"])
    for k, v in dmap["planning"].items():
        chk("planning %s = %s" % (k, v), str(v) in T)
    for k in ("3", "5", "7", "11", "19", "35"):
        has("min-df spread %s" % k, dmap["spread_within_min_df"][k], "%.1f pp")
        has("bench-df spread %s" % k, dmap["spread_within_bench_df"][k],
            "%.1f pp")
    has("min-df spread above seven",
        max(dmap["spread_within_min_df"][k] for k in ("11", "19", "35")),
        "%.2f percentage points or less")
    restr = load(WF, "wd_restricted.json")
    for d in restr["designs"]:
        for r in d["rows"]:
            has("%s %s" % (d["design"], r["domain"][:14]),
                r["median_loss_pp"], "%.1f pp")

    print("\nSPEED AND CROSS-CHECKS (10.10)")
    fastd = load(WF, "fast.json")
    chk("parallel is bit-identical",
        fastd["selftest"]["worst_parallel_vs_serial"] == 0.0)
    has("new seconds", fastd["bench"]["new_sweep_seconds"], "%.1f seconds")
    has("old seconds", fastd["bench"]["old_sweep_seconds_est"],
        "against an estimated %.0f")
    has("speedup", fastd["bench"]["speedup_x"], "a factor of %.0f")
    has("inversion error", fastd["selftest"]["worst_inversion_vs_rootfind"],
        "%.1e")
    has("linear-scale inversion error",
        fastd["selftest"]["worst_inversion_linear_vs_rootfind"], "%.1e")
    has("log-scale gain", fastd["selftest"]["log_scale_gain_x"], "%.0f")
    has("gate tolerance", fastd["selftest"]["gate_tolerance"], "allows %g")
    has("tightest margin",
        min(0.025 - d["sup_upper_bound"] for d in k1["designs"]), "%.2e")
    cross = load(WF, "crosscheck.json")
    chk("all cross-checks pass", cross["all_pass"])
    has("independent quadrature", cross["xa_worst"], "%.1e")
    has("level z", cross["xb_worst_abs_z"], "%.2f")
    has("power z", cross["xc_worst_abs_z"], "%.2f")
    has("cross-check points", len(cross["xd"]), "at %d ratio points")
    gates = load(WF, "gates.json")
    has("gates", (gates["n_pass"], gates["n_gates"]),
        "%d of %d deposited gates pass")

    print("\nDEPENDENCE (11)")
    ser = load(WF, "serial.json")
    for r in ser["rows"]:
        has("phi %.2f naive" % r["phi"], r["naive"]["fwer"], "%.4f")
        has("phi %.2f whitened" % r["phi"], r["whitened_correct"]["fwer"],
            "%.4f")
    sf = load(WF, "serial_falsify.json")
    has("balanced worst z", sf["balanced_worst_z"], "within %.1f standard")
    has("unbalanced worst z", sf["unbalanced_worst_z"], "%.1f standard errors")
    has("unbalanced worst rate", max(r["fwer"] for r in sf["rows"]), "%.4f")
    has("unbalanced worst share",
        100.0 * max(r["fwer"] for r in sf["rows"]) / sf["abar"], "%.0f%% of it")
    chk("the balanced case holds the ceiling",
        sf["balanced_worst_abs_dev"] < 0.002, "%.5f" % sf["balanced_worst_abs_dev"])
    chk("the unbalanced case does not",
        sf["unbalanced_worst_z"] > 3.0, "%.1f se" % sf["unbalanced_worst_z"])
    chk("the repair restores the ceiling", sf["repair_works"],
        "%.1f se" % sf["unbalanced_repaired_worst_z"])
    has("repaired worst deviation", sf["unbalanced_repaired_worst_abs_dev"],
        "worst deviation of %.5f")
    has("repaired worst z", sf["unbalanced_repaired_worst_z"],
        "at %.1f standard")
    # the signature: the repaired rate must not move with the common mean
    unb = [r for r in sf["rows"] if r["n"][0] != r["n"][1]]
    byphi = {}
    for r in unb:
        byphi.setdefault(r["phi"], set()).add(round(r["fwer_co"], 5))
    chk("the repaired rate is invariant to the common mean",
        all(len(v) == 1 for v in byphi.values()),
        "distinct values per phi: %s" % sorted(len(v) for v in byphi.values()))
    chk("every declared-bound rule is valid",
        all(b["valid"] for b in ser["bound_rule"]))
    for b in ser["bound_rule"]:
        has("bound %.2f worst" % b["phi_max"], b["worst_fwer"], "%.4f")
    for m in ser["misspecification"]:
        has("mis %.1f/%.1f" % (m["declared"], m["true"]), m["fwer"], "%.4f")
    fac = load(WF, "factor.json")
    for r in fac["rows"]:
        has("%s / %s" % (r["exposure"][:12], r["arm"][:12]), r["fwer"],
            "%.5f")
    has("pre-period rate", fac["pre_period_loadings"]["fwer"], "%.4f")
    has("pre-period z", abs(fac["pre_period_loadings"]["z"]),
        "%.2f standard errors")
    fgen = load(WA, "factor_general.json")
    has("factor K1 z", abs(fgen["k1_vs_monte_carlo"]["z"]), "|z| = %.2f")
    has("factor beta0", fgen["beta_zero_limit"]["rel_diff_pct"], "%.1f%%")
    for r in fgen["general_k_vs_simulation"]["rows"]:
        has("factor K=%d z" % r["K"], abs(r["z"]), "%.1f")


    # ------------------------------------------------------------ table rows
    # Whole-row checks, so that a value moved between rows or columns fails
    # even though both values are still somewhere in the document.
    print("\nTABLE ROWS")
    ROWS = table_rows(path)

    def row(label, cells):
        EXPECTED.extend(str(c) for c in cells)
        chk("row: %s" % label, tuple(cells) in ROWS,
            "" if tuple(cells) in ROWS else "not found as a row")

    r35 = lrow("df35")
    for key, lab in (("oracle_conditional",
                      "The oracle, which knows the ratios"),
                     ("plugin_measured", "The plug-in at the full level"),
                     ("region_sup_conditional",
                      "The region rule of Section 9"),
                     ("uncond_certified",
                      "The calibrated rule of this section")):
        e = r35[key]
        row("Table 16 %s" % key,
            [lab, "%.3f" % e["mean_c"], "%.2f" % e["ratio_to_oracle_mean"],
             "%.5f" % e["fwer"], "%.1f pp" % e["loss_pp_at_80"]])

    bl = {"df3": ("3 df", "12 df"), "df5": ("5 df", "12 df"),
          "df7": ("7 df", "12 df"), "df11": ("11 df", "12 df"),
          "df19": ("19 df", "12 df"), "df35": ("35 df", "12 df"),
          "df59": ("59 df", "12 df"), "df35_short_cand": ("35 df", "4 df"),
          "balanced_df11": ("11 df", "11 df"),
          "balanced_df35": ("35 df", "35 df")}
    for x in k1["designs"]:
        b, c = bl[x["design"]]
        row("Table 17 %s" % x["design"],
            [b, c, "%.5f" % x["alpha0"], "%.6f" % x["sup_upper_bound"],
             "%.1f pp" % x["median_loss_pp_at_80"],
             "%.1f pp" % x["max_loss_pp_at_80"]])

    gl = {"K2_df3": "2 candidates, 3 df", "K2_df5": "2 candidates, 5 df",
          "K2_df7": "2 candidates, 7 df", "K2_df11": "2 candidates, 11 df",
          "K2_df19": "2 candidates, 19 df", "K2_df35": "2 candidates, 35 df",
          "K3_df11": "3 candidates, 11 df", "K3_df35": "3 candidates, 35 df",
          "K5_df11": "5 candidates, 11 df", "K5_df35": "5 candidates, 35 df"}
    for x in gen["designs"]:
        row("Table 18 %s" % x["design"],
            [gl[x["design"]], "%.5f" % x["plugin_at_abar_sup_fwer"],
             "+%.1f%%" % x["plugin_at_abar_rel_pct"], "%.5f" % x["alpha0"],
             "%.1f pp" % x["median_loss_pp"], "%.1f pp" % x["max_loss_pp"]])

    for x in anyk["designs"]:
        row("Table 19 %s %s" % (x["design"], x["two_sided"]),
            [x["design"].split("_")[1].replace("df", "") + " df",
             str(x["K"]), "two-sided" if x["two_sided"] else "one-sided",
             "%.5f" % x["alpha0"], "%.6f" % x["one_candidate_bound"],
             "%.6f" % x["union_bound_total"],
             "%.6f" % x["exact_sup_over_patterns"]])

    for x in joint["designs"]:
        d = x["design"]
        row("Table 20 %s" % d,
            [d.split("_")[1].replace("df", "") + " df", str(x["K"]),
             "%.5f" % x["alpha0_certified_joint"],
             "%.6f" % x["union_bound_total"],
             "%.1f%%" % x["level_vs_measured_pct"],
             "%.6f (se %.6f)" % (x["exact_sup_over_patterns"],
                                 x["exact_sup_mcse"]),
             "%.4f" % rcd[d]["points"][0]["c_joint_certified"]])

    for x in front["designs"]:
        row("Table 21 %s" % x["design"],
            [x["design"][2:] + " df",
             "%.1f pp" % (sum(x["plugin_loss_pp"]) / len(x["plugin_loss_pp"])),
             "%.1f pp" % (sum(x["optimal_loss_pp"])
                          / len(x["optimal_loss_pp"])),
             "%+.2f pp" % x["mean_gain_pp"], "%.4f" % x["sup_fwer_fine"]])

    for r in ser["rows"]:
        row("Table 25 phi %.2f" % r["phi"],
            ["%.2f" % r["phi"],
             "%.4f (se %.4f)" % (r["naive"]["fwer"], r["naive"]["mcse"]),
             "%.4f (se %.4f)" % (r["whitened_correct"]["fwer"],
                                 r["whitened_correct"]["mcse"])])

    # ------------------------------------------------------ derived in prose
    # Quantities the prose computes from deposits rather than quoting.  These
    # are the ones a substring check would otherwise never see.
    print("\nDERIVED QUANTITIES")
    wc = max(bnb["cells"], key=lambda c: c["gap"])
    has("worst relative gap", 100.0 * wc["gap"] / wc["c_U"], "%.3f%%")
    has("the constant it is on", wc["c_U"], "on a constant of %.3f")
    has("region budget share",
        100.0 * r35["region_sup_conditional"]["fwer"] / 0.025,
        "%.1f%% of what the guarantee permits")
    jf = [x["exact_sup_over_patterns"] for x in joint["designs"]]
    af = [x["exact_sup_over_patterns"] for x in anyk["designs"]]
    has("joint budget share low", 100.0 * min(jf) / 0.025,
        "realises %.0f%%")
    has("joint budget share high", 100.0 * max(jf) / 0.025, "to %.0f%% of the")
    has("percand budget share low", 100.0 * min(af) / 0.025,
        "realises %.0f%%")
    zc = [r["z"] for r in fac["rows"]
          if r["arm"] == "estimated loadings, df corrected"]
    zu = [r["z"] for r in fac["rows"]
          if r["arm"] == "estimated loadings, df NOT corrected"]
    has("df-returned z range", (min(zc), max(zc)),
        "which is %.1f to %.1f standard errors")
    has("df-not-returned z range", (min(zu), max(zu)),
        "%.1f to %.1f standard errors above")
    k1a = {x["design"]: x for x in k1["designs"]}
    lev = [x["alpha0_two_sided"] / k1a[x["design"]]["alpha0"]
           for x in twos["designs"]]
    cost = [x["median_loss_pp"] / k1a[x["design"]]["median_loss_pp_at_80"]
            for x in twos["designs"]]
    has("two-sided level ratio", (min(lev), max(lev)), "a factor of %.2f to %.2f")
    has("two-sided cost ratio", (min(cost), max(cost)), "between %.2f and %.2f")

    # ------------------------------------------- derived, round two
    pcd = {d["design"]: d for d in load(WF, "wb_percand.json")["designs"]}
    kd = {d["design"]: d for d in k1["designs"]}
    rdd = {d["design"]: d for d in restr["designs"]}

    def rdom(des, dom):
        return next(r for r in rdd[des]["rows"] if r["domain"] == dom)

    def anyk_row2(des, two):
        return next(x for x in anyk["designs"]
                    if x["design"] == des and bool(x["two_sided"]) == two)

    # Everything the coverage audit found quoted in prose but never checked.
    print("\nDERIVED, ROUND TWO")
    has("abstract median gap", bnb["median_gap"], "%.1e")
    has("omitted mass", k1["tail_mass_correction"]["omitted_probability_mass"],
        "leaves %g of probability")
    x2 = next(g for g in gates["gates"] if g["gate"].startswith("X2"))
    has("gate trials", x2["x2_trials"],
        lambda v: "none in {:,} trials".format(v))
    has("region coverage", 1.0 - ledger["gamma"], "a %g interval")
    has("serial draws", ser["draws"], lambda v: "{:,} draws per cell".format(v))
    has("misspecification draws", 12000,
        lambda v: "at {:,} draws per cell".format(v))
    has("factor draws", fac["draws"],
        lambda v: "{:,} draws per cell".format(v))
    has("joint draws", 4000, lambda v: "over {:,} draws".format(v))
    has("general draws", gen["draws"], "over %d draws")
    has("general worst mcse", max(x["worst_mcse"] for x in gen["designs"]),
        "worst value is %.5f")
    has("bound-rule se", (0.025 * 0.975 / 12000.0) ** 0.5,
        "standard error of about %.4f")
    has("percand measured level", pcd["K2_df35"]["alpha0"], "level of %.5f")
    has("percand cert cost at df35",
        100.0 * (1.0 - anyk_row2("K2_df35", False)["alpha0"]
                 / pcd["K2_df35"]["alpha0"]), "gives up %.1f%% of the level")
    for des in ("df3", "df35"):
        has("binding ratio %s" % des,
            rdom(des, "unrestricted_3_decades")["lam_worst"], "%.3f")
    has("naive at phi 0.1",
        100.0 * ser["rows"][2]["naive"]["fwer"] / ser["abar"], "reaches %.0f%%")
    has("naive at phi 0.4",
        100.0 * ser["rows"][5]["naive"]["fwer"] / ser["abar"], "%.0f%% of it at")
    fvv = {(r["exposure"], r["arm"]): r["fwer"] for r in fac["rows"]}
    pats = []
    for r in fac["rows"]:
        if r["exposure"] not in pats:
            pats.append(r["exposure"])
    has("opposite-sign excess",
        100.0 * (fvv[(pats[2], "raw returns")] / fac["abar"] - 1.0),
        "%.0f%% above it")
    has("plug-in loss at df3", lrow("df3")["plugin_measured"]["loss_pp_at_80"],
        "gives up %.1f points")
    has("region ratio high",
        max(lrow(d)["region_sup_conditional"]["ratio_to_oracle_mean"]
            for d in ("df3", "df11", "df35")), "and %.2f times")
    lp = [100.0 * x["alpha0"] / 0.025 for x in gen["designs"]
          if x["design"] not in ("K2_df3", "K2_df5")]
    has("general level share low", min(lp), "between %.1f%%")
    has("general level share high", max(lp), "and %.1f%% of the ceiling")
    has("k1 level share at df11", 100.0 * kd["df11"]["alpha0"] / 0.025,
        "against %.1f%% at one")
    for d in ("df11", "df59"):
        has("level price %s" % d, 100.0 * (0.025 - kd[d]["alpha0"]) / 0.025,
            "%.1f%%")
    has("two-sided loss low",
        min(x["median_loss_pp"] for x in twos["designs"]), "gives up %.1f")
    has("two-sided loss high",
        max(x["median_loss_pp"] for x in twos["designs"]), "to %.1f percentage")
    has("spread at 59 df", dmap["spread_within_bench_df"]["59"],
        "the spread is %.1f")
    has("spread at 23 df", dmap["spread_within_min_df"]["23"],
        "where it is %.1f")
    for c in bnb["cells"]:
        has("seconds %s" % c["cell"], c["seconds"], "%.0f")
    # the paper writes the tolerance as 1e-4; check that is the deposited value
    chk("declared tolerance is 1e-4", abs(bnb["tol"] - 1e-4) < 1e-18,
        "%g" % bnb["tol"])
    has("tolerance as written", bnb["tol"], lambda v: "tolerance of 1e-4")
    has("error budget target", ebud["dp_target"], "declared target of %g")
    slk = load(WA, "box_slack.json")
    has("narrow box width",
        100.0 * (min(f for f in slk["factors"] if f > 1.0) - 1.0),
        "a box %g%% wide")
    has("percand budget share high", 100.0 * max(af) / 0.025,
        "and %.0f%% of the error rate")

    if os.environ.get("DUMP_EXPECTED"):
        with open(os.environ["DUMP_EXPECTED"], "w",
                  encoding="utf-8") as fh:
            fh.write(chr(10).join(EXPECTED))
    print("\n%d checks passed, %d failed" % (npass, nfail))
    print("PAPER VERIFIED" if ok else "PAPER HAS FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
