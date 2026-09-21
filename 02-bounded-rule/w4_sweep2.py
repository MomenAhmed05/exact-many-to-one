"""W4 reduced sweep: single-alpha* end-to-end falsifier vs B, cases S1 and S4.

Re-runs the paper's F1b arm with the corrected parallel calibrator at
B in {600, 2400} and D=1200, and reports the residual and its z.  The
extension (mean over independent alpha*) is the definitional test; this sweep
is the paper's single-estimate test.
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA  # noqa: E402
from io_util import dump_json, f64  # noqa: E402
from w4_calibration import calibrate_par, oracle_cal_fwer  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w4")
CASES = [
    ("k2 S1 balanced, equal sd", [10, 10, 10], [30.0, 30.0, 30.0]),
    ("k2 S4 control smallest n and sd", [4, 13, 13], [10.0, 10.0, 50.0]),
]
BLIST = (600, 2400)
D = 1200


def _main(pool):
    out = {"alpha": ALPHA, "B_list": list(BLIST), "draws": D, "cases": []}
    for label, n_l, sig_l in CASES:
        n = np.asarray(n_l, float)
        sig = np.asarray(sig_l, float)
        rec = {"case": label, "n": n_l, "sigma": sig_l, "by_B": {}}
        for B in BLIST:
            t0 = time.time()
            a_star, F = calibrate_par(
                n, sig, ALPHA, B, seed=31337, workers=3, pool=pool
            )
            arm = oracle_cal_fwer(n, sig, a_star, D, 9911, workers=3, pool=pool)
            rec["by_B"][str(B)] = {
                "alpha_star": f64(a_star),
                "F_at_alpha": f64(F(ALPHA)),
                "fwer": arm["fwer"],
                "residual": arm["residual"],
                "mcse": arm["mcse"],
                "z": arm["z"],
                "seconds": time.time() - t0,
            }
            print(
                f"{label}: B={B} alpha*={a_star:.5f} "
                f"FWER={arm['fwer']:.5f} resid={arm['residual']:+.5f} "
                f"z={arm['z']:+.2f}",
                flush=True,
            )
        out["cases"].append(rec)
    dump_json(out, os.path.join(OUTDIR, "calibration_B2.json"))
    print("wrote", os.path.join(OUTDIR, "calibration_B2.json"))


if __name__ == "__main__":
    with mp.Pool(3) as pool:
        _main(pool)
