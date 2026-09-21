"""Track B acceptance: end-to-end |z| at B = 3000, plain vs control-variate alpha*.

The W4 falsifier failed on single-estimate z (S1 B600 -5.86, S4 B600 +2.70,
B2400 -3.58).  With a 2.1-2.4x smaller estimator standard error at the same B,
the same test should pass: this script measures it at B = 3000, D = 1200.
"""

import multiprocessing as mp
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calib_vr import CASES, POOL_N, calibrate, fwer_case  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "deposits", "calib_vr")


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


B_ACC = _env_int("W26_ACC_B", 3000)
D_ACC = _env_int("W26_ACC_D", 1200)
OUTER_SEEDS = (7, 8)


def _finite_or_none(x):
    """JSON-safe value: non-finite becomes null (bare NaN is invalid strict JSON)."""
    try:
        v = float(x)
        return v if v == v and abs(v) != float("inf") else None
    except Exception:  # noqa: BLE001
        return None


def main():
    rows = []
    with mp.Pool(POOL_N) as pool:
        for name, n, sig in CASES:
            for variant in ("plain", "vr"):
                a_star, diag = calibrate(n, sig, 0.025, B_ACC, 777, variant, pool=pool)
                for s in OUTER_SEEDS:
                    r = fwer_case(n, sig, a_star, D_ACC, s, pool=pool)
                    payload = {
                        "case": name,
                        "variant": variant,
                        "a_star": f64(a_star),
                        "B": i64(B_ACC),
                        "outer_seed": i64(s),
                        "D": i64(D_ACC),
                        "F_se": diag["F_se"],
                        "control_r2": _finite_or_none(diag["control_r2"]),
                        "fwer": r["fwer"],
                        "mcse": r["mcse"],
                        "residual": r["residual"],
                        "z": r["z"],
                    }
                    rows.append(payload)
                    print(
                        f"{name[:22]:<22s} {variant:<6s} B={B_ACC} "
                        f"seed={s} a*={a_star:.5f} "
                        f"resid={payload['residual']:+.5f} "
                        f"z={payload['z']:+.2f}",
                        flush=True,
                    )
    agg = []
    for name, _n, _s in CASES:
        for variant in ("plain", "vr"):
            sub = [r for r in rows if r["case"] == name and r["variant"] == variant]
            try:
                zs = [float(r["z"]) for r in sub]
                res = [float(r["residual"]) for r in sub]
                worst = max(abs(z) for z in zs)
                agg.append(
                    {
                        "case": name,
                        "variant": variant,
                        "mean_residual": f64(float(np.mean(res))),
                        "max_abs_z": f64(worst),
                        "n_outer": i64(len(sub)),
                    }
                )
            except Exception:  # noqa: BLE001
                agg.append(
                    {
                        "case": name,
                        "variant": variant,
                        "mean_residual": f64(float("nan")),
                        "max_abs_z": f64(float("nan")),
                        "n_outer": i64(0),
                    }
                )
    passes = [
        a
        for a in agg
        if a["variant"] == "vr" and np.isfinite(a["max_abs_z"]) and a["max_abs_z"] < 3.0
    ]
    out = {
        "B": i64(B_ACC),
        "D": i64(D_ACC),
        "outer_seeds": [i64(s) for s in OUTER_SEEDS],
        "rows": rows,
        "aggregate": agg,
        "vr_cases_passing_z3": i64(len(passes)),
        "n_cases": i64(len(CASES)),
        "all_vr_pass": bool(len(passes) == len(CASES)),
        "note": (
            "plain and vr use the same B and the same outer seeds; the "
            "comparison is paired in B, so a smaller |z| at the same B "
            "is the variance-reduction effect"
        ),
    }
    dump_json(out, os.path.join(OUT, "accept.json"))
    for a in agg:
        print(
            f"{a['case'][:22]:<22s} {a['variant']:<6s} mean_resid "
            f"{a['mean_residual']:+.5f} max|z| {a['max_abs_z']:.2f}",
            flush=True,
        )
    print(f"vr cases with max|z| < 3: {len(passes)}/{len(CASES)}")


if __name__ == "__main__":
    main()
