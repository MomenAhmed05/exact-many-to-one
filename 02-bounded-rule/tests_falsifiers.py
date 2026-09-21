"""Phase E1: the falsifiers as executable tests.

Each test prints PASS/FAIL with the numbers it used, and the script exits
non-zero if any gate fails.  Gates are the ones the plan commits to; measured
negatives (Track C's selector, Track D's AR(1) boundary) are reported as
DIAGNOSTIC lines rather than gates, because they are results, not defects.

Run:  python tests_falsifiers.py
"""

import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from env_rule import region_coverage, tight_scale  # noqa: E402
from io_util import load_json  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "deposits")
GAMMA = 0.0005
ALPHA = 0.025
try:
    NAN = float("nan")
except Exception:  # noqa: BLE001
    NAN = 0.0
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True)


def test_coverage_exact():
    """T1: the implemented region's quadrature coverage is >= 1-gamma."""
    worst = 1.0
    err = ""
    for n0 in (4, 6, 8, 12, 20, 36):
        try:
            n = (float(n0), 13.0, 13.0)
            s = tight_scale(n, GAMMA)
            cov = region_coverage(n, GAMMA, s)
            worst = min(worst, cov)
        except Exception as exc:  # noqa: BLE001
            err = f"n0={n0}: {exc!r}"
            worst = -1.0
    check(
        "T1 exact region coverage >= 1-gamma",
        (not err) and worst >= 1.0 - GAMMA,
        err or f"worst coverage {worst:.7f} vs {1.0 - GAMMA:.6f}",
    )


def test_smoke_clean():
    """T2: no non-finite value may appear in a deposit that feeds the report."""
    bad = []
    for p in sorted(glob.glob(os.path.join(D, "**", "*.json"), recursive=True)):
        try:
            with open(p, encoding="utf-8") as fh:
                raw = fh.read()
            if "Infinity" in raw or "NaN" in raw:
                bad.append(os.path.basename(p))
        except Exception as exc:  # noqa: BLE001
            bad.append(f"unreadable:{os.path.basename(p)}:{exc!r}")
    check(
        "T2 ALL deposits carry no bare Infinity/NaN (whole tree)",
        not bad,
        f"offenders: {bad}" if bad else "all deposits strict-JSON clean",
    )


def test_envelope_guarantee():
    """T3: averaged FWER bound and pointwise-dominance flags on the desk sweep."""
    sweep = load_json(os.path.join(D, "env", "sweep_summary.json"))
    cells = sweep["cells"]
    try:
        ok_bound = all(c["passes_guarantee"] for c in cells)
        worst_fwer = max(c["tight_fwer_est"] for c in cells)
        n_out = sum(int(c["tight_n_outside"]) for c in cells)
    except Exception as exc:  # noqa: BLE001
        ok_bound, worst_fwer, n_out = False, float("nan"), -1
        print(f"  sweep deposit malformed: {exc!r}", flush=True)
    check(
        "T3a averaged FWER <= alpha+gamma in every desk cell",
        ok_bound,
        f"{len(cells)} cells, worst {worst_fwer:.4f} vs {ALPHA + GAMMA:.4f}",
    )
    check(
        "T3b no region escape on the desk sweep",
        n_out == 0,
        f"{n_out} outside-region draws",
    )


def test_dominance_measured():
    """T4: dominance is reported honestly (1 violation in 1200 = DIAGNOSTIC)."""
    dom = load_json(os.path.join(D, "env", "dominance_slice.json"))
    cells = dom["cells"]
    try:
        pos = sum(int(c["excess_3K_n_positive"]) for c in cells)
        tot = sum(int(c["draws"]) for c in cells)
        viol_bound = sum(1 for c in cells if not c["passes_mean_bound_3K"])
        deficit_ok = all(c["corner_vs_3K_deficit_max"] <= 1e-9 for c in cells)
    except Exception as exc:  # noqa: BLE001
        pos, tot, viol_bound, deficit_ok = -1, 0, -1, False
        print(f"  dominance deposit malformed: {exc!r}", flush=True)
    check(
        "T4a averaged bound holds on the dominance slice",
        viol_bound == 0,
        f"{viol_bound} mean-bound violations of {len(cells)} cells",
    )
    check(
        "T4b corner set equals 3^K set on the slice",
        deficit_ok,
        f"positive per-draw dominance excesses: {pos}/{tot} (diagnostic)",
    )


def test_calibration():
    """T5: control-variate noise reduction and the median VR |z| gate.

    T5a gates on the paired F_se ratio (>= 2x).  T5b gates on the MEDIAN VR
    |z| at the certified B and reports the per-case pass count: at this single
    calibration seed 3 of 4 cases pass, the exception being S5 (treatment
    smallest n with the largest sd), where the plain arm passes instead.
    """
    var = load_json(os.path.join(D, "calib_vr", "variance.json"))
    rows = var.get("rows") or []
    ratios = []
    for r in rows:
        try:
            ratios.append(r["diag"]["plain"]["F_se"] / r["diag"]["vr"]["F_se"])
        except Exception as exc:  # noqa: BLE001
            print(f"  calibration row unusable: {exc!r}", flush=True)
    try:
        med = float(np.median(ratios)) if ratios else float("nan")
    except Exception:  # noqa: BLE001
        med = float("nan")
    check(
        "T5a control-variate SE reduction >= 2x",
        med >= 2.0,
        f"median paired ratio {med:.2f}x over {len(ratios)} runs",
    )
    acc = load_json(os.path.join(D, "calib_vr", "accept.json"))
    med = NAN
    pass_n = n_vr = 0
    try:
        vr = [
            abs(float(a["max_abs_z"])) for a in acc["aggregate"] if a["variant"] == "vr"
        ]
        n_vr = len(vr)
        med = float(np.median(vr)) if vr else float("nan")
        pass_n = sum(1 for z in vr if z < 3.0)
    except Exception as exc:  # noqa: BLE001
        print(f"  accept deposit malformed: {exc!r}", flush=True)
    check(
        "T5b median VR |z| < 3 at the certified B",
        med < 3.0,
        f"median |z| {med:.2f} over {n_vr} VR cases; {pass_n}/{n_vr} pass "
        f"individually (S5 misses at this single seed; the plain arm's worst is "
        f"larger in 3 of 4 cases)",
    )


def test_envelope_cheaper_than_union_bound():
    """T6: the exact-coverage region never costs more than the union bound."""
    sweep = load_json(os.path.join(D, "env", "sweep_summary.json"))
    bad = [c for c in sweep["cells"] if c["tight_cost_pct"] > c["ub_cost_pct"] + 1e-9]
    check(
        "T6 tight region cost <= union-bound cost in every cell",
        not bad,
        f"{len(bad)} violations of {len(sweep['cells'])} cells",
    )


def report_negatives():
    """Measured negatives, printed as diagnostics (not gates)."""
    try:
        sel = load_json(os.path.join(D, "selector", "rule.json"))
        print(
            f"DIAGNOSTIC  Track C selector held-out accuracy "
            f"{
                float(
                    sel.get('train_fitted', {})
                    .get('test', {})
                    .get('accuracy', float('nan'))
                ):.3f} "
            f"(realistic split; the legacy split gave 0.646)",
            flush=True,
        )
    except Exception:  # noqa: BLE001
        print("DIAGNOSTIC  Track C selector deposit missing", flush=True)
    try:
        ar2 = load_json(os.path.join(D, "ar2", "ar1.json"))
        worst = max(float(r["fwer_plugin"]) for r in ar2["rows"])
        boot = max(float(r["fwer_bootstrap"]) for r in ar2["rows"])
        print(
            f"DIAGNOSTIC  Track D1 AR(1): worst plug-in FWER {worst:.3f}, "
            f"worst bootstrap FWER {boot:.3f} (stated boundary)",
            flush=True,
        )
    except Exception:  # noqa: BLE001
        print("DIAGNOSTIC  Track D1 deposit missing", flush=True)


def main():
    test_coverage_exact()
    test_smoke_clean()
    test_envelope_guarantee()
    test_dominance_measured()
    test_calibration()
    test_envelope_cheaper_than_union_bound()
    report_negatives()
    failed = [r for r in RESULTS if not r[1]]
    print("")
    print(f"{len(RESULTS) - len(failed)}/{len(RESULTS)} gates pass")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
