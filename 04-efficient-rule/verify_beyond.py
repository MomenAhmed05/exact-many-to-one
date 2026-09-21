"""Check every number in REPORT-BEYOND.md against its deposit.

Same discipline as `verify_frontier.py`: each expected string is FORMATTED FROM
A DEPOSIT and then looked for in the report, so editing the report cannot make a
check pass and editing a deposit cannot either.  Only agreement between them
does.

Two things learned from the preceding programme are built in.

  * Rows that matter are matched WHOLE.  Substring containment cannot see a
    transposition inside a table row, because both values are still somewhere
    in the document.
  * A missing or stale deposit FAILS its section and the verifier carries on.
    A crash stops the remaining checks and so hides them, which is the opposite
    of what a verifier is for.

    python verify_beyond.py             check
    python verify_beyond.py probes      check, then run the corruption probes
"""

from __future__ import annotations

import io
import os
import re
import sys

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(HERE, "REPORT-BEYOND.md")
PROTOCOL = os.path.join(HERE, "PROMOTION-PROTOCOL.md")

ok = True
npass = nfail = 0
FAILED = []


def load(path):
    return re.sub(r"\s+", " ", io.open(path, encoding="utf-8").read())


TXT = load(REPORT) if os.path.exists(REPORT) else ""
PROT = load(PROTOCOL) if os.path.exists(PROTOCOL) else ""


class Stale(Exception):
    """A deposit is missing or was written by an older version of its module."""


def chk(label, cond, detail=""):
    global ok, npass, nfail
    cond = bool(cond)
    ok &= cond
    npass += cond
    nfail += not cond
    if not cond:
        FAILED.append(label)
    print(("  OK   " if cond else "  FAIL ") + label +
          ("   " + detail if detail else ""))


EXPECTED = []


def has(label, value, fmt="%s", txt=None):
    v = fmt % value
    EXPECTED.append(v)
    chk("%s -> %s" % (label, v), v in (TXT if txt is None else txt))


def row(label, cells):
    """Match a whole table row, so a transposition inside it cannot pass."""
    s = "| " + " | ".join(cells) + " |"
    EXPECTED.append(s)
    chk("row: %s" % label, s in TXT, "" if s in TXT else s[:90])


def dget(name, *keys):
    try:
        d = F.read_deposit(name)
    except Exception as exc:                                   # noqa: BLE001
        chk("deposit present: %s" % name, False, str(exc)[:60])
        raise Stale(name)
    missing = [k for k in keys if k not in d]
    if missing:
        chk("deposit current: %s" % name, False,
            "missing " + ", ".join(missing))
        raise Stale(name)
    return d


# ------------------------------------------------------------------ sections
def s_invariant():
    d = dget("w1_invariant.json", "max_hidden_pp", "bin_resolution")
    has("sufficient dim", d["sufficient_dim"], "%d")
    has("group dim", d["group_dim"], "%d")
    has("invariant dim", d["maximal_invariant_dim"], "%d")
    has("invariance drift", d["worst_invariance_drift"], "%.2e")
    has("recovery error", d["worst_recovery_error"], "%.2e")
    has("recovery formula", d["recovery_formula"], "%s")
    has("max hidden", d["max_hidden_pp"], "%.3f points")
    chk("hunt-stein flagged as derivation",
        "derivation, not a computation" in TXT
        or "derivation, not a proof" in TXT)
    for x in d["bin_resolution"]:
        r = dict((q["nbins"], q["mean_loss_pp"]) for q in x["rows"])
        row("bins %s" % x["design"],
            [x["design"], "%.2f pp" % r[36], "%.2f pp" % r[72],
             "%.2f pp" % r[144], "%.3f pp" % x["hidden_pp"]])


def s_region():
    d = dget("w1_region.json", "worst_overshoot",
             "min_gain_pp_from_dropping_monotonicity")
    has("min gain", d["min_gain_pp_from_dropping_monotonicity"], "%+.2f to")
    has("max gain", d["max_gain_pp_from_dropping_monotonicity"], "%+.2f points")
    has("t-grid sensitivity", d["worst_t_grid_sensitivity_pp"], "%.3f points")
    dfl = [x["loss_pp_deterministic_floor"] - x["loss_pp_monotone"]
           for x in d["designs"]]
    has("det-floor min", min(dfl), "%.2f to")
    has("det-floor max", max(dfl), "%.2f points")
    has("overshoot", d["worst_overshoot"], "%.2e")
    has("overshoot pct", d["worst_overshoot_pct_of_abar"], "%.2f%% of the")
    chk("spurious first result disclosed", "spurious 1.47-point gain" in TXT)
    for x in d["designs"]:
        row("region %s" % x["design"],
            [x["design"], "%.2f pp" % x["loss_pp_monotone"],
             "%.2f pp" % x["loss_pp_randomised"],
             "%.2f pp" % x["loss_pp_deterministic_floor"],
             "%+.2f pp" % x["gain_pp_from_dropping_monotonicity"]])


def s_lowdf():
    d = dget("w1_lowdf.json", "mean_gain_by_min_df", "operational_rule")
    has("material threshold", d["material_threshold_pp"], "%.2f points")
    has("material cells", d["n_material_cells"], "%d of ")
    has("max gain", d["max_gain_pp"], "%.2f points")
    for k, v in d["mean_gain_by_min_df"].items():
        row("gain min nu %s" % k, [k, "%.2f pp" % v])


def s_factor():
    d = dget("w2_factor.json", "falsifiers", "critical_table")
    chk("reduction holds", d["reduction_holds"])
    has("outer dim", d["outer_dimension"], "%d-dimensional outer rule")
    has("worst z", d["worst_abs_z"], "%.2f")
    has("brute draws", "{:,}".format(d["brute_draws"]), "%s draws")
    has("refinement", d["refinement_gap"], "%.1e")
    has("departure", 100 * d["max_relative_departure_from_independent_law"],
        "%.1f%%")
    for r in d["falsifiers"]:
        row("falsifier %s" % r["cell"],
            [r["cell"], "%.5f" % r["brute"], "%.5f" % r["exact"],
             "%+.2f" % r["z"], "%.5f" % r["independent_law"]])
    for t in d["critical_table"]:
        row("critical %s" % t["design"], [t["design"], "%.4f" % t["c_alpha"]])


def s_cert():
    d = dget("w2_cert.json", "mechanism", "price_rows",
             "n_opposite_sign_liberal")
    has("liberal cells", d["n_liberal_cells"], "%d of ")
    has("worst iid fwer", d["worst_fwer_iid_rule_on_aligned_data"], "%.5f")
    has("worst ratio", d["worst_ratio_to_abar"], "%.2f times")
    has("worst same sign", d["worst_same_sign_fwer"], "%.5f")
    has("opposite-sign liberal", d["n_opposite_sign_liberal"], "%d of ")
    has("box points", d["box"]["points"], "%d points")
    has("certified constant", d["certified_constant_over_box"], "%.4f")
    has("oracle min", d["oracle_constant_min"], "%.4f")
    has("oracle max", d["oracle_constant_max"], "%.4f")
    has("mechanism", d["mechanism"], "%s")
    has("price vs worst oracle", d["certificate_price_vs_worst_oracle_x"],
        "%.2f times")
    has("price vs best oracle", d["certificate_price_vs_best_oracle_x"],
        "%.2f times")
    chk("grid maximum labelled", "grid maximum" in TXT)
    chk("three-dimensional nuisance stated",
        "three-dimensional per candidate" in TXT)
    for r in d["price_rows"]:
        row("price lam %.2f b0 %.1f b1 %.1f" % (r["lam"], r["b0"], r["b1"]),
            ["%.2f" % r["lam"], "%.1f" % r["b0"], "%.1f" % r["b1"],
             "%.5f" % r["fwer"], "%+.1f" % r["z"]])


def s_unequal():
    d = dget("w2_unequal.json", "saturates_at_prehistory",
             "full_history_power_gain_by_pre")
    has("overlap", d["n_overlap"], "overlap of %d observations")
    has("two-window gain", d["two_window_max_power_gain_pp"], "%.2f points")
    has("full-history gain", d["full_history_max_power_gain_pp"], "%.2f points")
    has("broken fwer", d["broken_stationarity_worst_fwer"], "%.5f")
    has("broken ratio", d["broken_ratio_to_abar"], "%.2f times")
    chk("both break directions reported",
        "makes the test powerless" in TXT and "makes the test liberal" in TXT)
    for k, v in d["full_history_power_gain_by_pre"].items():
        row("prehistory %s" % k, [k, "%+.2f pp" % v])


def s_tiers():
    d = dget("w3_tiers.json", "serial_tierA_worst_fwer")
    for key in ("serial", "factor"):
        for cl in (36, 120, 360):
            has("%s tierB %d" % (key, cl),
                d["%s_tierB_cal%d_worst_fwer" % (key, cl)], "%.5f")
        has("%s tierA" % key, d["%s_tierA_worst_fwer" % key], "%.5f")
        has("%s tierC" % key, d["%s_tierC_worst_fwer" % key], "%.5f")
    for cl in (36, 120, 360):
        row("tierB %d" % cl,
            ["%d" % cl, "%.5f" % d["serial_tierB_cal%d_worst_fwer" % cl],
             "%.5f" % d["factor_tierB_cal%d_worst_fwer" % cl]])


def s_repair():
    d = dget("w3_serial_repair.json", "pooled_spread_over_mu_tierC",
             "separate_worst_at_mu5")
    has("separate worst", d["separate_worst_fwer"], "%.5f")
    has("pooled worst", d["pooled_worst_fwer"], "%.5f")
    has("separate ratio", d["separate_worst_ratio"], "%.1f")
    has("pooled ratio", d["pooled_worst_ratio"], "%.1f")
    has("spread B", d["pooled_spread_over_mu_tierB"], "%.5f")
    has("spread C", d["pooled_spread_over_mu_tierC"], "%.5f")
    chk("pooling not called a repair",
        "Pooling is an improvement, not a repair" in TXT)
    chk("tier-A-only conclusion stated", "is a tier-A result" in TXT)
    row("repair worst",
        ["worst realised rate", "%.5f" % d["separate_worst_fwer"],
         "%.5f" % d["pooled_worst_fwer"]])
    row("repair mu0",
        ["worst at a zero mean", "%.5f" % d["separate_worst_at_mu0"],
         "%.5f" % d["pooled_worst_at_mu0"]])
    row("repair mu5",
        ["worst at a mean of 5", "%.5f" % d["separate_worst_at_mu5"],
         "%.5f" % d["pooled_worst_at_mu5"]])
    # the mean-invariance demonstration, matched whole
    wc = max(d["rows"], key=lambda r: r["fwer"])
    trip = [r for r in d["rows"] if r["design"] == wc["design"]
            and r["phi"] == wc["phi"] and r["tier"] == wc["tier"]]
    for m in sorted(set(r["mu"] for r in trip)):
        sep = [r["fwer"] for r in trip if r["mu"] == m and not r["pooled"]]
        poo = [r["fwer"] for r in trip if r["mu"] == m and r["pooled"]]
        row("invariance mu %.0f" % m,
            ["%.0f" % m, "%.5f" % sep[0], "%.5f" % poo[0]])


def s_asymptotic():
    d = dget("w3_asymptotic.json", "mu", "methods")
    has("cells per method", d["cells_per_method"], "%d cells per method")
    has("common mean", d["mu"], "common mean %.1f")
    for mth in d["methods"]:
        row("asym %s" % mth,
            [mth, "%.5f" % d["worst_size"][mth],
             "%.2f" % d["worst_ratio_to_abar"][mth],
             "%d of %d" % (d["breaches"][mth], d["cells_per_method"])])
    chk("declared-phi reference holds",
        d["worst_ratio_to_abar"]["the exact rule, coefficient declared"] < 1.05)
    chk("reference labelled as not a competitor",
        "REFERENCE, not a competitor" in TXT)
    est = [m for m in d["methods"] if m != "the exact rule, coefficient declared"]
    chk("no estimated-coefficient method holds",
        all(d["worst_ratio_to_abar"][m] > 1.5 for m in est))


def s_stress():
    d = dget("w3_stress.json", "breaching_cells")
    has("breach se", d["breach_se"], "%.0f standard errors")
    has("breaches", d["n_breaches"], "%d of ")
    has("cells", d["n_cells"], "%d cells breach")
    has("worst", d["worst_fwer"], "%.5f")
    has("worst ratio", d["worst_ratio_to_abar"], "%.2f times")
    has("worst cell", d["worst_cell"], "`%s`")
    for c in d["breaching_cells"]:
        has("breaching cell", c, "%s")
    for c in d["surviving_cells"]:
        has("surviving cell", c, "%s")
    chk("short variance window called out",
        "short variance window" in d["breaching_cells"]
        and "rolling risk model" in TXT)


def s_protocol():
    d = dget("w4_manifest.json", "probes", "verdict_labels")
    chk("all probes refused", d["all_probes_refused"])
    has("probes refused", d["n_refused"], "%d of ")
    has("probe count", d["n_probes"], "%d refusal probes")
    has("example verdict", d["example_verdict"], "**%s**")
    has("tier c verdict", d["tier_c_verdict"], "**%s**")
    for p in d["probes"]:
        chk("probe listed: %s" % p["probe"], p["probe"] in TXT)
    for v in d["verdict_labels"]:
        chk("protocol defines %s" % v, v in PROT)
    for key in ("Exact in model", "Certified bound", "Measured",
                "Asymptotic", "Diagnostic"):
        chk("taxonomy: %s" % key, key in PROT)


def s_probes():
    d = dget("probe_beyond.json", "n_transposition_caught", "probes")
    chk("all live probes caught", d["all_caught"])
    has("probes caught", d["n_caught"], "%d of ")
    has("live probes", d["n_live"], "%d probes caught")
    has("transpositions caught", d["n_transposition_caught"], "%d of ")
    for r in d["probes"]:
        if not r["skipped"]:
            row("probe %s" % r["probe"],
                [r["probe"], "yes" if r["caught"] else "no"])
    chk("the vacuous probe is disclosed",
        "A vacuous probe is worse than no probe" in TXT)


def s_coverage():
    d = dget("coverage_beyond.json", "n_uncovered", "excused")
    chk("no uncovered numeric token", d["n_uncovered"] == 0)
    has("tokens covered", d["n_covered"], "%d of ")
    has("tokens total", d["n_tokens"], "%d distinct tokens")
    has("tokens excused", d["n_excused"], "%d excused")
    for t in d["excused"]:
        chk("excusal named: %s" % t, ("`%s`, " % t) in TXT)


def s_labels():
    for lab in ("exact in model", "measured", "diagnostic"):
        chk("label used: %s" % lab, "*(%s)*" % lab in TXT)
    # "certified bound" must NOT appear as a label, and the report has to say
    # why: nothing in this programme earns it.  A future run that produces a
    # real certificate flips both of these together.
    chk("no result claims a certified bound",
        "*(certified bound)*" not in TXT)
    chk("the absence is stated, not merely true",
        "No result in this programme is a *certified bound*" in TXT)
    chk("grid versus continuum distinction drawn", "holds over a GRID" in TXT)


SECTIONS = [
    ("W1 INVARIANCE", s_invariant),
    ("W1 REGION LP", s_region),
    ("W1 LOW DF", s_lowdf),
    ("W2 ALIGNED FACTOR", s_factor),
    ("W2 CERTIFICATE", s_cert),
    ("W2 UNEQUAL HISTORIES", s_unequal),
    ("W3 TIERS", s_tiers),
    ("W3 SERIAL REPAIR", s_repair),
    ("W3 ASYMPTOTIC ALTERNATIVES", s_asymptotic),
    ("W3 STRESS", s_stress),
    ("W4 PROTOCOL", s_protocol),
    ("CORRUPTION PROBES", s_probes),
    ("COVERAGE", s_coverage),
    ("TAXONOMY LABELS IN THE REPORT", s_labels),
]


def main():
    if not TXT:
        print("REPORT-BEYOND.md not found; run w_report.py first")
        return 1
    for title, fn in SECTIONS:
        print("\n" + title)
        try:
            fn()
        except Stale:
            pass                       # the failure is already recorded
    # the coverage hook sits AFTER every section, or it measures a stale
    # set and reports coverage the checks do not actually have
    dump = os.environ.get("DUMP_EXPECTED")
    if dump:
        with io.open(os.path.join(HERE, dump), "w",
                     encoding="utf-8") as fh:
            for e in EXPECTED:
                fh.write(e + chr(10))
        print("wrote %d expected strings to %s" % (len(EXPECTED), dump))

    print("\n%d checks passed, %d failed" % (npass, nfail))
    if FAILED:
        print("failures:")
        for f in FAILED[:30]:
            print("   " + f)
    print("REPORT VERIFIED" if ok else "REPORT HAS FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
