"""Track B probe: is the S5 VR regression bias, or one unlucky calibration seed?

The corrected calib run left S5 (treatment smallest n, largest sd) with a VR
residual of +0.00164 (|z| = 8.38) while the plain arm sat at +0.00026.  A single
seed cannot distinguish (a) an unlucky realisation from (b) a systematic bias in
the control-variate arm for that configuration.  This script calibrates S5 at
several independent seeds, at the certified B, and reports the distribution of
the VR and plain residuals with each arm's outer-noise z.

Outputs: deposits/calib_vr/s5_probe.json
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
B = 3000
D = 1200
OUTER_SEED = 7
SEEDS = (777, 1313, 2024, 4242, 5555, 6789, 8008, 9009)
NAME = "k2 S5 trt smallest n, largest sd"


def _finite_or_none(x):
    """JSON-safe value: non-finite becomes null.

    Bare NaN/Infinity is invalid strict JSON, and the earlier smoke deposit was
    polluted exactly this way; the plain arm's control_r2 is not applicable by
    construction, so it must serialise as null rather than NaN.
    """
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except Exception:  # noqa: BLE001
        return None


def _run_variant(n, sig, seed, variant, pool):
    try:
        a_star, diag = calibrate(n, sig, 0.025, B, seed, variant, pool=pool)
        r = fwer_case(n, sig, a_star, D, OUTER_SEED, pool=pool)
        return {
            "a_star": _finite_or_none(a_star),
            "residual": _finite_or_none(r["residual"]),
            "z": _finite_or_none(r["z"]),
            "F_se": _finite_or_none(diag["F_se"]),
            "control_r2": _finite_or_none(diag["control_r2"]),
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "a_star": None,
            "residual": None,
            "z": None,
            "F_se": None,
            "control_r2": None,
            "error": repr(exc),
        }


def _col(rows, variant, key):
    vals = []
    for r in rows:
        try:
            v = r[variant][key]
            if v is not None:
                vals.append(float(v))
        except Exception:  # noqa: BLE001
            print(f"  skipping unreadable {variant}/{key} row", flush=True)
    return vals


def _stats(vals):
    try:
        if not vals:
            return {"mean": f64(float("nan")), "sd": f64(float("nan")), "n": i64(0)}
        sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else float("nan")
        return {"mean": f64(float(np.mean(vals))), "sd": f64(sd), "n": i64(len(vals))}
    except Exception:  # noqa: BLE001
        return {"mean": f64(float("nan")), "sd": f64(float("nan")), "n": i64(0)}


def _verdict(vals):
    try:
        if len(vals) < 3:
            return "inconclusive (too few seeds)"
        mean = float(np.mean(vals))
        se = float(np.std(vals, ddof=1) / np.sqrt(len(vals)))
        if se <= 0:
            return "inconclusive (zero spread)"
        t = mean / se
        if abs(t) > 2.5:
            return (
                f"systematic bias: mean {mean:+.5f} = {t:+.1f} x its own SE, "
                f"so the VR arm is off for this configuration"
            )
        return (
            f"consistent with zero mean (mean {mean:+.5f}, t {t:+.1f}): "
            f"the single B=3000 seed was unlucky rather than biased"
        )
    except Exception:  # noqa: BLE001
        return "inconclusive"


def main():
    case = [c for c in CASES if c[0] == NAME]
    if not case:
        raise SystemExit(f"case not found: {NAME}")
    _name, n, sig = case[0]
    rows = []
    with mp.Pool(POOL_N) as pool:
        for seed in SEEDS:
            rec = {
                "seed": i64(seed),
                "plain": _run_variant(n, sig, seed, "plain", pool),
                "vr": _run_variant(n, sig, seed, "vr", pool),
            }
            rows.append(rec)
            print(
                f"seed {seed}: vr resid {rec['vr']['residual']} "
                f"z {rec['vr']['z']} | plain resid {rec['plain']['residual']} "
                f"z {rec['plain']['z']}",
                flush=True,
            )
    vr_res = _col(rows, "vr", "residual")
    pl_res = _col(rows, "plain", "residual")
    summary = {
        "vr_residual": _stats(vr_res),
        "plain_residual": _stats(pl_res),
        "vr_z": _stats(_col(rows, "vr", "z")),
        "plain_z": _stats(_col(rows, "plain", "z")),
    }
    try:
        delta = (
            f64(float(np.mean(pl_res) - np.mean(vr_res)))
            if (pl_res and vr_res)
            else f64(float("nan"))
        )
    except Exception:  # noqa: BLE001
        delta = f64(float("nan"))
    out = {
        "case": NAME,
        "B": i64(B),
        "D": i64(D),
        "outer_seed": i64(OUTER_SEED),
        "seeds": [i64(s) for s in SEEDS],
        "rows": rows,
        "summary": summary,
        "verdict": _verdict(vr_res),
        "mean_plain_minus_vr_residual": delta,
    }
    dump_json(out, os.path.join(OUT, "s5_probe.json"))
    print("verdict:", out["verdict"])
    print(
        "vr residual mean/sd:",
        summary["vr_residual"]["mean"],
        summary["vr_residual"]["sd"],
    )
    print(
        "plain residual mean/sd:",
        summary["plain_residual"]["mean"],
        summary["plain_residual"]["sd"],
    )


if __name__ == "__main__":
    main()
