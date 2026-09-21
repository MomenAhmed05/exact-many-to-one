"""Check every number quoted in REPORT-FRONTIER.md against its deposit.

Same discipline as the preceding programmes, and the same reason: a check is
only worth running if it can fail for the right reason.  Each expected string is
FORMATTED FROM A DEPOSIT and then looked for in the report, so editing the
report cannot make a check pass and editing a deposit cannot make one pass
either - only agreement between them does.
"""

from __future__ import annotations

import io
import os
import re
import sys

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))
TXT = io.open(os.path.join(HERE, "REPORT-FRONTIER.md"), encoding="utf-8").read()
TXT = re.sub(r"\s+", " ", TXT)

ok = True
npass = nfail = 0


def chk(label, cond, detail=""):
    global ok, npass, nfail
    cond = bool(cond)
    ok &= cond
    npass += cond
    nfail += not cond
    print(("  OK   " if cond else "  FAIL ") + label +
          ("   " + detail if detail else ""))


def has(label, value, fmt="%s"):
    v = fmt % value if not callable(fmt) else fmt(value)
    chk("%s -> %s" % (label, v), v in TXT)


def main():
    print("CERTIFIED K=1 RULE")
    d = F.read_deposit("wb_certified_rule.json")
    chk("all certified", d["all_certified"])
    for x in d["designs"]:
        if x["design"] in ("df3", "df5", "df7", "df11", "df19", "df35",
                           "df59", "balanced_df11"):
            has("%s alpha0" % x["design"], x["alpha0"], "%.5f")
            has("%s median loss" % x["design"], x["median_loss_pp_at_80"],
                "%.1f")

    print("\nGENERAL-K MEASURED RULE")
    g = F.read_deposit("wb_general.json")
    chk("all within ceiling", g["all_within_ceiling"])
    for x in g["designs"]:
        has("%s alpha0" % x["design"], x["alpha0"], "%.5f")
        has("%s max loss" % x["design"], x["max_loss_pp"], "%.1f")

    print("\nFRONTIER / LOWER BOUND")
    w = F.read_deposit("wc_frontier.json")
    chk("all valid", w["all_valid"])
    for x in w["designs"]:
        has("%s gain" % x["design"], x["mean_gain_pp"], "%+.2f")

    print("\nDESIGN MAP")
    m = F.read_deposit("wd_designmap.json")
    chk("binds on min df", m["binds_on_min_df"])
    for k, v in m["planning"].items():
        chk("planning %s = %s" % (k, v), str(v) in TXT)

    print("\nLEDGER")
    ld = F.read_deposit("wa_ledger.json")
    for des in ld["designs"]:
        if des["design"] not in ("df3", "df11", "df35"):
            continue
        for row in des["rows"]:
            if abs(row["lam"] - 1.0) > 1e-9:
                continue
            for lab in ("oracle_conditional", "plugin_measured",
                        "region_sup_conditional", "uncond_certified"):
                e = row[lab]
                has("%s %s mean c" % (des["design"], lab[:9]),
                    e["mean_c"], "%.3f")
                has("%s %s loss" % (des["design"], lab[:9]),
                    e["loss_pp_at_80"], "%.1f")

    print("\nRESTRICTED DOMAINS")
    rs = F.read_deposit("wd_restricted.json")
    for des in rs["designs"]:
        has("%s best recovery" % des["design"], des["best_recovery_pp"],
            "%.1f")
        for r in des["rows"]:
            if r["domain"] in ("unrestricted_3_decades", "factor_of_4",
                               "factor_of_2", "vol_targeted_1.5"):
                has("%s %s loss" % (des["design"], r["domain"][:12]),
                    r["median_loss_pp"], "%.1f")

    print("\nPER-CANDIDATE CLASS")
    try:
        pc = F.read_deposit("wb_percand.json")
    except Exception:  # noqa: BLE001
        pc = None
    if pc:
        chk("per-candidate within ceiling", pc["all_within_ceiling"])
        for x in pc["designs"]:
            has("%s percand alpha0" % x["design"], x["alpha0"], "%.5f")
            has("%s percand max loss" % x["design"], x["max_loss_pp"], "%.1f")
    else:
        chk("per-candidate deposit present", False, "not written yet")

    print("\nTWO-SIDED (M6)")
    ts = F.read_deposit("twosided.json")
    chk("two-sided reduction exact", ts["reduction"]["pass"],
        "cond %.1e unc %.1e" % (ts["reduction"]["worst_conditional_abs_diff"],
                                ts["reduction"]["worst_unconditional_abs_diff"]))
    chk("reduction FAILS above K=1 (the falsifier can fire)",
        ts["no_reduction_above_k1"]["pass"],
        "ratio %.4f" % ts["no_reduction_above_k1"]["ratio"])
    chk("all two-sided within ceiling", ts["all_within_ceiling"])
    for x in ts["designs"]:
        has("%s 2s alpha0" % x["design"], x["alpha0_two_sided"], "%.5f")
        has("%s 2s median loss" % x["design"], x["median_loss_pp"], "%.1f")

    print("\nSERIAL DEPENDENCE (M6)")
    sr = F.read_deposit("serial.json")
    for r in sr["rows"]:
        has("phi %.2f naive" % r["phi"], r["naive"]["fwer"], "%.4f")
        has("phi %.2f whitened" % r["phi"], r["whitened_correct"]["fwer"],
            "%.4f")
    chk("a declared-bound rule exists and is valid",
        all(b["valid"] for b in sr["bound_rule"]))

    print("\nFACTOR (M6)")
    fc = F.read_deposit("factor.json")
    for r in fc["rows"]:
        has("%s / %s" % (r["exposure"][:14], r["arm"][:14]), r["fwer"],
            "%.5f")
    has("pre-period loadings", fc["pre_period_loadings"]["fwer"], "%.5f")

    print("\nCERTIFIED AT ANY K (M3/M6)")
    ak = F.read_deposit("certify_anyk.json")
    chk("all any-K certificates pass", ak["all_certified"],
        "%d certificates" % len(ak["designs"]))
    for x in ak["designs"]:
        has("%s %s alpha0" % (x["design"], "2s" if x["two_sided"] else "1s"),
            x["alpha0"], "%.5f")
        chk("%s %s union total <= abar"
            % (x["design"], "2s" if x["two_sided"] else "1s"),
            x["union_bound_total"] <= ak["abar"] + 1e-12,
            "%.6f" % x["union_bound_total"])
    has("worst union slack", ak["worst_union_slack"], "%.2e")

    print("\nJOINT RULE CERTIFIED (K = 2, 3)")
    jc = F.read_deposit("certify_joint.json")
    chk("all joint certificates pass", jc["all_certified"],
        "%d designs" % len(jc["designs"]))
    chk("every joint bound dominates the deployed rule's own FWER",
        jc["all_bounds_dominate"])
    for x in jc["designs"]:
        has("%s joint alpha0" % x["design"], x["alpha0_certified_joint"],
            "%.5f")
        has("%s joint gain over per-candidate" % x["design"],
            x["level_gain_vs_percand_x"], "%.2fx")
        has("%s joint realised FWER" % x["design"],
            x["exact_sup_over_patterns"], "%.6f")

    print("\nSPEED LAYER")
    fj = F.read_deposit("fast.json")
    chk("parallel map is bit-identical to serial",
        fj["selftest"]["worst_parallel_vs_serial"] == 0.0)
    chk("tabulated inversion within tolerance",
        fj["selftest"]["worst_inversion_vs_rootfind"] < 1e-3,
        "%.1e" % fj["selftest"]["worst_inversion_vs_rootfind"])
    has("sweep speedup", fj["bench"]["speedup_x"], "%.0fx")

    print("\nCROSS-CHECK (M7)")
    cc = F.read_deposit("crosscheck.json")
    chk("all cross-checks pass", cc["all_pass"])
    chk("certified bound never breached by a sample", cc["xd_all_dominated"],
        "%d points" % len(cc["xd"]))
    has("cross-check worst |z| on FWER", cc["xb_worst_abs_z"], "%.2f")
    has("cross-check worst |z| on power", cc["xc_worst_abs_z"], "%.2f")

    print("\nGATES")
    gt = F.read_deposit("gates.json")
    chk("gates all pass", gt["n_pass"] == gt["n_gates"],
        "%d/%d" % (gt["n_pass"], gt["n_gates"]))

    # ------------------------------------------------ the desk-facing document
    print("\nDESK PROCEDURE")
    DESK = re.sub(r"\\s+", " ", io.open(
        os.path.join(HERE, "DESK-PROCEDURE.md"), encoding="utf-8").read())

    def dhas(label, value, fmt="%s"):
        v = fmt % value
        chk("desk: %s -> %s" % (label, v), v in DESK)

    kd1 = F.read_deposit("wb_certified_rule.json")
    for x in kd1["designs"][:7]:
        dhas("%s alpha0" % x["design"], x["alpha0"], "%.5f")
        dhas("%s bound" % x["design"], x["sup_upper_bound"], "%.6f")
        dhas("%s loss" % x["design"], x["median_loss_pp_at_80"], "%.1f points")
    ak2 = F.read_deposit("certify_anyk.json")
    for x in ak2["designs"]:
        if not x["two_sided"]:
            dhas("%s level" % x["design"], x["alpha0"], "%.5f")
            dhas("%s own rate" % x["design"], x["exact_sup_over_patterns"],
                 "%.5f")
    dm2 = F.read_deposit("wd_designmap.json")
    for k, v in dm2["planning"].items():
        chk("desk: planning %s = %s" % (k, v), "| %d |" % v in DESK)
    ld2 = F.read_deposit("wa_ledger.json")
    row = next(q for d in ld2["designs"] if d["design"] == "df35"
               for q in d["rows"] if q["lam"] == 1.0)
    dhas("region loss", row["region_sup_conditional"]["loss_pp_at_80"],
         "%.1f points")
    dhas("calibrated loss", row["uncond_certified"]["loss_pp_at_80"],
         "%.1f points")
    sf2 = F.read_deposit("serial_falsify.json")
    dhas("serial breach", max(r["fwer"] for r in sf2["rows"]), "%.5f")
    ex2 = [x["exact_sup_over_patterns"] for x in ak2["designs"]]
    dhas("spend low", 100.0 * min(ex2) / 0.025, "%.0f%%")
    dhas("spend high", 100.0 * max(ex2) / 0.025, "%.0f%%")

    print("\n%d checks passed, %d failed" % (npass, nfail))
    print("REPORT VERIFIED" if ok else "REPORT HAS FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
