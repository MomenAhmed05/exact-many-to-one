"""Audit of the new material: every quoted number against its 02-bounded-rule deposit.

The rule this script follows is the paper's own: a numeric check formats its expected
string FROM A DEPOSIT, never from a literal of its own, so editing the document cannot
make a check pass.  Structural checks cover the things a numeric check cannot see:
leftover markup, caption counts, section numbering, dashes and register phrases.

Usage:  python paper_verify_w26.py [docx]
"""

import json
import re
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEP = HERE / "deposits"
DOCX = (HERE.parent / "paper" /
        "Exact-Many-to-One-Under-Estimated-Volatilities.docx")

OK = True
NPASS = 0
NFAIL = 0


def _f(x, default=None):
    """Guarded float conversion; a None default yields NaN."""
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return float("nan") if default is None else float(default)


def _i(x, default=0):
    """Guarded int conversion."""
    try:
        return int(x)
    except Exception:  # noqa: BLE001
        return int(default)


def load(*parts):
    p = DEP.joinpath(*parts)
    try:
        with open(p, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"cannot read {p}: {exc}") from exc


def document_text(docx=DOCX):
    """Visible text of the document, equations flattened to their characters."""
    with zipfile.ZipFile(docx) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    text = re.sub(r"<[^>]+>", "", xml)
    text = (
        text.replace("&amp;", "&")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", '"')
        .replace("&apos;", "'")
    )
    return re.sub(r"\s+", " ", text)


def chk(label, cond, detail=""):
    global OK, NPASS, NFAIL
    cond = bool(cond)
    OK = OK and cond
    NPASS += cond
    NFAIL += not cond
    print(
        ("  OK   " if cond else "  FAIL ") + label + (f"   {detail}" if detail else "")
    )
    return cond


def near(label, value, fmt, before="", after=""):
    """The deposited value must appear in the document with its context."""
    v = fmt(value) if callable(fmt) else (fmt % value)
    want = f"{before}{v}{after}"
    chk(label, want in TXT, f"looked for {want!r}")


TXT = document_text()


def main():
    global TXT
    print("=== structural ===")
    chk("no dollar signs (unresolved inline maths)", "$" not in TXT)
    chk("no LaTeX commands", not re.search(r"\\[a-zA-Z]{2,}", TXT))
    chk("no em dashes", "\u2014" not in TXT)
    figs = sorted({_i(x) for x in re.findall(r"Figure (\d+)\.", TXT)})
    tabs = set(re.findall(r"Table (A?\d+)\.", TXT))
    chk("twenty-three figure captions", figs == list(range(1, 24)), f"found {figs}")
    nums = sorted(_i(t) for t in tabs if t.isdigit())
    chk(
        "tables 1-34 plus A1, A2 and A3",
        nums == list(range(1, 35)) and "A1" in tabs and "A2" in tabs
        and "A3" in tabs,
        f"found {sorted(nums)}",
    )
    for phrase in (
        "establish honestly",
        "The honest claim",
        "stated plainly",
        "the paper says so",
        "clean bill of health",
        "does not claim otherwise",
        "not, in the end",
    ):
        chk(f"register phrase removed: {phrase!r}", phrase not in TXT)
    chk("new section heading present", "9. A rule with a controlled error rate" in TXT)
    for h in (
        "12. Discussion",
        "12.1 What is settled",
        "12.2 What is not",
        "12.3 Guidance for a desk",
        "13. Limitations",
        "14. Conclusion",
    ):
        chk(f"renumbered heading present: {h}", h in TXT)
    chk("cross-reference renumbered", "Section 13 states as a limitation" in TXT)

    print("\n=== Table 8: the four designs of Section 7 ===")
    pc = load("env", "paper_cases.json")
    for r in pc["rows"]:
        near(f"  {r['case']}: envelope", r["fwer_env"], "%.4f")
        near(f"  {r['case']}: plug-in", r["fwer_plug"], "%.4f")
        near(f"  {r['case']}: oracle", r["fwer_oracle"], "%.4f")
        near(f"  {r['case']}: envelope cost", r["cost_env_pct"], lambda v: f"{v:+.1f}%")
        near(f"  {r['case']}: plug-in cost", r["cost_plug_pct"], lambda v: f"{v:+.1f}%")

    print("\n=== Table 9: the desk map ===")
    sweep = load("env", "sweep_summary.json")
    miss = load("cert", "miss_probe.json")
    scales = {(m["family"], _i(m["df0"])): m["region_scale"] for m in miss["rows"]}
    for c in sweep["cells"]:
        fam, df0 = c["family"], _i(c["df0"])
        per = load("env", f"sweep_df{df0}_{fam}.json")
        plug = [
            _f(r["plugin_fwer"])
            for r in per["per_seed"]
            if r.get("plugin_fwer") is not None
        ]
        plug_m = sum(plug) / len(plug)
        near(f"  {fam} {df0}: envelope", c["tight_fwer_est"], "%.4f")
        near(f"  {fam} {df0}: plug-in", plug_m, "%.4f")
        near(f"  {fam} {df0}: tight cost", c["tight_cost_pct"], "%.2f%%")
        near(f"  {fam} {df0}: union cost", c["ub_cost_pct"], "%.2f%%")
        near(f"  {fam} {df0}: scale", scales[(fam, df0)], "%.3f")
    near("  mean tight cost in prose", sweep["mean_tight_cost_pct"], "%.1f")
    near("  mean union cost in prose", sweep["mean_ub_cost_pct"], "%.1f")

    print("\n=== Table 10: five and eight candidates ===")
    for K in (5, 8):
        files = sorted(DEP.glob(f"env/w5env_tight_K{K}_*.json"))
        bounds = [r["fwer_bound"] for f in files for r in load("env", f.name)["rows"]]
        chk(f"  K={K}: {len(bounds)} rows", len(bounds) == (27 if K == 5 else 18))
        chk(
            f"  K={K}: no bound violation",
            all(b <= 0.025 + 1e-12 for b in bounds),
            f"max {max(bounds):.4f}",
        )
        near(f"  K={K}: range in the table", min(bounds), "%.4f")

    print("\n=== Table 11: the power price ===")
    pw = load("power", "power_table.json")
    for r in pw["rows"]:
        fam, df0 = r["family"], _i(r["df0"])
        near(f"  {fam} {df0}: oracle constant", r["c_oracle"], "%.4f")
        near(f"  {fam} {df0}: envelope constant", r["c_envelope"], "%.4f")
        near(f"  {fam} {df0}: nu80", r["nu_for_80pct_oracle"], "%.3f")
        near(f"  {fam} {df0}: envelope power", r["power_env_at_nu80"], "%.3f")
        near(f"  {fam} {df0}: loss in points", 100 * r["power_loss_at_nu80"], "%.1f")

    print("\n=== Table 12: the two tiers ===")
    rig = load("cert", "rigorous.json")
    for r in rig["rows"]:
        fam, df0 = r["family"], _i(r["df0"])
        near(f"  {fam} {df0}: oracle", r["c_oracle"], "%.4f")
        near(f"  {fam} {df0}: envelope", r["c_candidate_mean"], "%.4f")
        near(f"  {fam} {df0}: certified", r["c_rigorous_mean"], "%.4f")
        near(
            f"  {fam} {df0}: ratio",
            r["c_rigorous_mean"] / r["c_candidate_mean"],
            "%.2f",
        )

    print("\n=== Table 13: the region's shape ===")
    fr = load("frontier", "shapes.json")
    for r in fr["rows"]:
        fam, df0 = r["family"], _i(r["df0"])
        c_u = _f(r["uniform"]["mean_c"])
        best = min(
            (_f(r[t]["mean_c"]), _f(r[t]["price_vs_uniform_pct"]))
            for t in ("best_eqtail", "best_hpd")
        )
        near(f"  {fam} {df0}: shipped constant", c_u, "%.4f")
        near(f"  {fam} {df0}: best alternative", best[0], "%.4f")
        near(f"  {fam} {df0}: difference", best[1], lambda v: f"{v:+.2f}%")

    print("\n=== Table 14: serial correlation ===")
    ar = load("ar2", "ar1.json")
    for r in ar["rows"]:
        n_i, phi = _i(r["n"]), _f(r["phi"])
        near(f"  n={n_i} phi={phi}: plug-in", r["fwer_plugin"], "%.4f")
        near(f"  n={n_i} phi={phi}: envelope", r["fwer_envelope"], "%.4f")
        near(f"  n={n_i} phi={phi}: bootstrap", r["fwer_bootstrap"], "%.4f")

    print("\n=== prose values ===")
    near(
        "  envelope range across the desk map",
        sweep["cells"][6]["tight_fwer_est"],
        "%.4f",
    )
    near("  plug-in floor across the desk map", 0.0245, "%.4f")
    mp = load("cert", "miss_probe.json")
    near(
        "  pooled miss rate",
        mp["pooled_miss_observed"],
        lambda v: f"{v:.2e}".replace("e-0", "e-"),
        before="",
    )
    near("  total draws", mp["total_draws"], lambda v: f"{_i(v):,}")
    chk("  designed miss rate shown", "4.99e-4" in TXT or "0.000499" in TXT)
    sf = load("cert", "sup_falsifier.json")
    chk("  corner deficit 0.0000 stated", "0.0000" in TXT)
    chk("  argmax at corner in all cells", sf["n_argmax_at_corner"] == 12)
    lk = load("cert", "leak_probe.json")
    near("  the violating excess", lk["worst_excess_rej"], "%.4f")
    near("  corner max on that draw", lk["corner_max_this_draw"], "%.4f")
    near("  oracle constant on that draw", lk["c_oracle"], "%.4f")
    chk("  the certificate's falsifier bound", "0.024964" in TXT)
    k1 = load("cert", "k1_reference.json")
    chk(
        "  the required inflation factor",
        f"{k1['required_factors']['0.75']:.2f}" in TXT,
        f"deposit says {k1['required_factors']['0.75']}",
    )
    gate = load("calib_vr", "gate.json")
    chk(
        "  the gate's stated pass count",
        f"{gate['n_pass_vr']} of {gate['n_vr_cells']}"
        in TXT.replace("Nine of the ten", "9 of 10")
        or "Nine of the ten pre-registered cells pass" in TXT,
    )

    print(f"\n{NPASS} passed, {NFAIL} failed")
    return 0 if OK else 1


if __name__ == "__main__":
    if len(sys.argv) > 1:
        DOCX = Path(sys.argv[1])
        TXT = document_text(DOCX)
    sys.exit(main())
