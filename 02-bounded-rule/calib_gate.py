"""S3.1: the pre-registered multi-seed calibration gate.

Criterion frozen in `plans/desk-guarantee-program.md` BEFORE this run:

  seeds 101,202,303,404,505,606,707,808 ; B = 3000 ; D = 1200 ; alpha = 0.025
  statistic: mean over seeds of the end-to-end residual (FWER - alpha), with
  between-seed SE = sd/sqrt(8)
  acceptance: |mean residual| <= max(0.0002, 2*SE)
              AND mean + 2*SE <= 0.0005
  B-sweep {600, 1500, 3000, 5000}: smallest B satisfying both

The sweep runs on the two informative cases (S1 balanced; S5 treatment-smallest-n
with the largest sd, which was the single-seed exception); S4 and S6 are checked at
B = 3000 only.  The plain arm is measured on the same seeds for contrast and carries
no criterion.

Outputs: deposits/calib_vr/gate.json
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calib_vr import CASES, calibrate, fwer_case  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "calib_vr")

# --- frozen pre-registration (do not edit after the run starts) --------------
SEEDS = (101, 202, 303, 404, 505, 606, 707, 808)
B_MAIN = 3000
B_SWEEP = (600, 1500, 3000, 5000)
D = 1200
ALPHA = 0.025
GAMMA = 0.0005
BAR_ABS = 0.0002
BAR_UPPER = 0.0005
MAIN_CASES = ("k2 S1 balanced, equal sd", "k2 S5 trt smallest n, largest sd")
EXTRA_CASES = (
    "k2 S4 control smallest n and sd",
    "k2 S6 control smallest n, largest sd",
)
# ---------------------------------------------------------------------------


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


POOL_N = _env_int("W26_POOL", 3)


def run_case(name, n, sig, b_list, pool):
    rows = []
    for B in b_list:
        for seed in SEEDS:
            rec = {"B": i64(B), "seed": i64(seed)}
            for arm in ("vr", "plain"):
                try:
                    a_star, _diag = calibrate(n, sig, ALPHA, B, seed, arm, pool=pool)
                    r = fwer_case(n, sig, a_star, D, seed, pool=pool)
                    rec[arm] = {
                        "a_star": f64(a_star),
                        "residual": r["residual"],
                        "z": r["z"],
                        "mcse": r["mcse"],
                    }
                except Exception as exc:  # noqa: BLE001
                    rec[arm] = {
                        "a_star": None,
                        "residual": None,
                        "z": None,
                        "mcse": None,
                        "error": repr(exc),
                    }
            rows.append(rec)
    agg = []
    for B in b_list:
        for arm in ("vr", "plain"):
            vals = []
            for r in rows:
                if r["B"] == B and r[arm]["residual"] is not None:
                    vals.append(float(r[arm]["residual"]))
            try:
                mean = float(np.mean(vals)) if vals else float("nan")
                se = (
                    float(np.std(vals, ddof=1) / np.sqrt(len(vals)))
                    if len(vals) > 1
                    else float("nan")
                )
                upper = mean + 2.0 * se if np.isfinite(se) else float("nan")
                lower = mean - 2.0 * se if np.isfinite(se) else float("nan")
                thr = max(BAR_ABS, 2.0 * se) if np.isfinite(se) else BAR_ABS
                passes = bool(
                    np.isfinite(mean) and abs(mean) <= thr and upper <= BAR_UPPER
                )
            except Exception:  # noqa: BLE001
                mean = se = upper = lower = float("nan")
                passes = False
            agg.append(
                {
                    "B": i64(B),
                    "arm": arm,
                    "n_seeds": i64(len(vals)),
                    "mean_residual": f64(mean),
                    "se": f64(se),
                    "lower2se": f64(lower),
                    "upper2se": f64(upper),
                    "passes": passes,
                }
            )
    return {
        "case": name,
        "n": list(n),
        "sig": list(sig),
        "rows": rows,
        "aggregate": agg,
    }


def main():
    t0 = time.time()
    out = {
        "seeds": [i64(s) for s in SEEDS],
        "B_main": i64(B_MAIN),
        "B_sweep": [i64(b) for b in B_SWEEP],
        "D": i64(D),
        "alpha": ALPHA,
        "gamma": GAMMA,
        "bar_abs": f64(BAR_ABS),
        "bar_upper": f64(BAR_UPPER),
        "preregistered_in": "plans/desk-guarantee-program.md",
        "cases": [],
    }
    lookup = {c[0]: c for c in CASES}
    with mp.Pool(POOL_N) as pool:
        for name in MAIN_CASES:
            if name not in lookup:
                raise SystemExit(f"case not found: {name}")
            _n, n, sig = lookup[name]
            res = run_case(name, n, sig, B_SWEEP, pool)
            out["cases"].append(res)
            for a in res["aggregate"]:
                print(
                    f"{name[:24]:<24s} B={a['B']:<5d} {a['arm']:<5s} mean "
                    f"{a['mean_residual']:+.6f} +/-2SE [{a['lower2se']:+.6f}, "
                    f"{a['upper2se']:+.6f}] pass {a['passes']}",
                    flush=True,
                )
        for name in EXTRA_CASES:
            if name not in lookup:
                continue
            _n, n, sig = lookup[name]
            res = run_case(name, n, sig, (B_MAIN,), pool)
            out["cases"].append(res)
            for a in res["aggregate"]:
                print(
                    f"{name[:24]:<24s} B={a['B']:<5d} {a['arm']:<5s} mean "
                    f"{a['mean_residual']:+.6f} +/-2SE [{a['lower2se']:+.6f}, "
                    f"{a['upper2se']:+.6f}] pass {a['passes']}",
                    flush=True,
                )
    passing = [
        a["B"]
        for c in out["cases"]
        for a in c["aggregate"]
        if a["arm"] == "vr" and a["passes"]
    ]
    try:
        out["smallest_passing_B"] = i64(min(passing)) if passing else i64(-1)
    except Exception:  # noqa: BLE001
        out["smallest_passing_B"] = i64(-1)
    out["n_pass_vr"] = i64(
        sum(
            1
            for c in out["cases"]
            for a in c["aggregate"]
            if a["arm"] == "vr" and a["passes"]
        )
    )
    out["n_vr_cells"] = i64(
        sum(1 for c in out["cases"] for a in c["aggregate"] if a["arm"] == "vr")
    )
    out["seconds"] = f64(time.time() - t0)
    dump_json(out, os.path.join(OUT, "gate.json"))
    print(
        f"vr cells passing: {out['n_pass_vr']}/{out['n_vr_cells']} | smallest "
        f"passing B: {out['smallest_passing_B']} | {time.time() - t0:.0f}s"
    )


if __name__ == "__main__":
    main()
