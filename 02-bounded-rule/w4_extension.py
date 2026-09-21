"""W4 extension: separate the calibration falsifier's noise from any bias.

deviation_k = mean(nr_o - nr) = FWER(rule) - alpha, the paired control-variate
definition used in core.rb_fwer.

The end-to-end falsifier (calibrate at the truth, apply that level to the
plug-in rule, FWER must be alpha) uses ONE alpha* estimate per scenario, so at
finite B its deviation is dominated by the estimation noise in alpha*.  Here
K independent alpha* estimates are drawn on the SAME data stream, and the mean
deviation is reported with its own standard error.  The mean is an estimate of
the systematic part of the falsifier's failure; the spread is the finite-B
noise the paper attributes it to.

Outputs: deposits/w4/extension.json
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, batch_c, cond_noreject, draw_var, oracle_c  # noqa
from io_util import dump_json, f64  # noqa
from w4_calibration import calibrate_par  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w4")
CASES = [
    ("k2 S1 balanced, equal sd", [10, 10, 10], [30.0, 30.0, 30.0]),
    ("k2 S4 control smallest n and sd", [4, 13, 13], [10.0, 10.0, 50.0]),
]
K_REPS = 8
D = 1200


def main(B=600):
    import multiprocessing as mp

    with mp.Pool(4) as pool:
        _main(B, pool)


def _main(B, pool):
    out = {"alpha": ALPHA, "B": B, "K_reps": K_REPS, "draws": D, "cases": []}
    for label, n_l, sig_l in CASES:
        n = np.asarray(n_l, float)
        sig = np.asarray(sig_l, float)
        lam = list((sig / sig[0]) ** 2)[1:]
        s0, si, lh = draw_var(n, lam, D, 9911)
        c_or = oracle_c(n, lam, ALPHA)
        nr_o = cond_noreject(np.full(D, c_or), s0, si, n, lam)
        t0 = time.time()
        stars, devs = [], []
        for k in range(K_REPS):
            a_star, _F = calibrate_par(
                n, sig, ALPHA, B, seed=5000 + 37 * k, workers=4, pool=pool
            )
            stars.append(a_star)
            c_rule = batch_c(n, lh, alpha=a_star, workers=4, pool=pool)
            nr = cond_noreject(c_rule, s0, si, n, lam)
            devs.append(f64((nr_o - nr).mean()))
        devs = np.asarray(devs)
        stars = np.asarray(stars)
        rec = {
            "case": label,
            "n": n_l,
            "sigma": sig_l,
            "alpha_star_mean": f64(stars.mean()),
            "alpha_star_sd": f64(stars.std(ddof=1)),
            "deviation_mean": f64(devs.mean()),
            "deviation_sd_across_reps": f64(devs.std(ddof=1)),
            "deviation_se_of_mean": f64(devs.std(ddof=1) / np.sqrt(len(devs))),
            "z_mean": f64(devs.mean() / (devs.std(ddof=1) / np.sqrt(len(devs)))),
            "deviations": [f64(x) for x in devs],
            "fwer_mean": f64(ALPHA + devs.mean()),
            "seconds": time.time() - t0,
        }
        out["cases"].append(rec)
        print(
            f"{label}: B={B} alpha* sd {rec['alpha_star_sd']:.5f} "
            f"mean dev {rec['deviation_mean']:+.6f} "
            f"(SE {rec['deviation_se_of_mean']:.6f}, z {rec['z_mean']:+.2f})",
            flush=True,
        )
    dump_json(out, os.path.join(OUTDIR, "extension.json"))
    print("wrote", os.path.join(OUTDIR, "extension.json"))


if __name__ == "__main__":

    def _cli_int(argv):
        try:
            return int(argv[0]) if argv else 600
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(f"bad B: {argv}") from exc

    main(_cli_int(sys.argv[1:]))
