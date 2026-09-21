"""W4: rebuild the calibration and drive its failed end-to-end falsifier.

The paper's F1b falsifier: calibrate at the TRUE ratios, apply that level to the
plug-in rule, and the resulting FWER must be alpha.  The paper reports it failing
at |z| = 3.8 with B = 600 bootstrap replicates and a residual of ~0.0006, and
attributes the residual to estimation noise in alpha*.  Here the bootstrap is
parallelised, B is raised through 600 -> 9600, and the same falsifier is re-run
on paired draws.  Pass bar is the paper's own: worst |z| < 3.

Outputs: deposits/w4/calibration_B.json (this full sweep was superseded by
the reduced w4_sweep2.py, which writes calibration_B2.json)
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, cond_noreject  # noqa: E402
from h12_calibrated import BOOT, NGRID, SPAN, plugin_sigma  # noqa: E402
from h12_hetero_exact import HeteroMaxT  # noqa: E402
from io_util import dump_json, f64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w4")

CASES = [
    ("k2 S1 balanced, equal sd", [10, 10, 10], [30.0, 30.0, 30.0]),
    ("k2 S4 control smallest n and sd", [4, 13, 13], [10.0, 10.0, 50.0]),
    ("k2 S5 trt smallest n, largest sd", [13, 13, 4], [10.0, 10.0, 50.0]),
    ("k2 S6 control smallest n, largest sd", [4, 13, 13], [50.0, 10.0, 10.0]),
]


def _boot_task(args):
    """One bootstrap replicate of the calibrator: (G, H) on the c-grid."""
    n, sw, s2, cgrid = args
    mod = HeteroMaxT(
        np.asarray(n, float), plugin_sigma(np.asarray(s2, float)), grid=BOOT
    )
    g = np.array([1.0 - mod.cdf(f64(c)) for c in cgrid])
    h = 1.0 - cond_noreject(
        np.asarray(cgrid, float),
        np.array([s2[0]]),
        np.asarray(s2, float)[None, 1:],
        np.asarray(n, float),
        np.asarray(sw, float)[1:] ** 2,
    )
    return np.asarray(g, float), np.asarray(h, float)


def calibrate_par(
    n, sigma_world, alpha, B, seed, workers=10, ngrid=NGRID, span=SPAN, pool=None
):
    """Parallel version of h12_calibrated.calibrate; same definition."""
    import multiprocessing as mp

    n = np.asarray(n, float)
    sw = np.asarray(sigma_world, float) / np.asarray(sigma_world, float)[0]
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    c_ref = HeteroMaxT(n, sw, grid=BOOT).critical_value(alpha)
    cgrid = c_ref + np.linspace(span[0], span[1], ngrid)
    s2s = [(rng.chisquare(nu, size=n.size) / nu) * sw**2 for _ in range(B)]
    args = [(tuple(n), tuple(sw), tuple(s2), tuple(cgrid)) for s2 in s2s]
    if pool is None:
        with mp.Pool(workers) as p:
            out = p.map(_boot_task, args)
    else:
        out = pool.map(_boot_task, args)
    G = np.array([o[0] for o in out])
    H = np.array([o[1] for o in out])

    def F(a):
        tot = 0.0
        for b in range(B):
            c = np.interp(a, G[b][::-1], cgrid[::-1])
            tot += np.interp(c, cgrid, H[b])
        return tot / B

    lo, hi = 1e-5, 0.5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if F(mid) < alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), F


def oracle_cal_fwer(n, sig, a_star, D, seed, alpha=ALPHA, workers=10, pool=None):
    """FWER of: calibrate at the truth, apply that level to the plug-in rule."""
    from core import batch_c, draw_var, oracle_c

    n = np.asarray(n, float)
    sig = np.asarray(sig, float)
    lam = list((sig / sig[0]) ** 2)[1:]
    s0, si, lh = draw_var(n, lam, D, seed)
    c_or = oracle_c(n, lam, alpha)
    c_rule = batch_c(n, lh, alpha=a_star, workers=workers, pool=pool)
    d = cond_noreject(np.full(D, c_or), s0, si, n, lam) - cond_noreject(
        c_rule, s0, si, n, lam
    )
    f = alpha + f64(d.mean())
    return {
        "fwer": f,
        "mcse": f64(d.std(ddof=1) / np.sqrt(D)),
        "residual": f64(f - alpha),
        "z": f64((f - alpha) / (d.std(ddof=1) / np.sqrt(D))),
    }


def main(B_list=(600, 1200, 2400, 4800), draws=1200, seeds=(9911,)):
    import multiprocessing as mp

    with mp.Pool(6) as pool:
        return _main(B_list, draws, seeds, pool)


def _main(B_list, draws, seeds, pool):
    res = {
        "alpha": ALPHA,
        "B_list": list(B_list),
        "draws": draws,
        "seeds": list(seeds),
        "cases": {},
    }
    for label, n_l, sig_l in CASES:
        n = np.asarray(n_l, float)
        sig = np.asarray(sig_l, float)
        rec = {"n": n_l, "sigma": sig_l, "by_B": {}}
        for B in B_list:
            t0 = time.time()
            a_star, F = calibrate_par(
                n, sig, ALPHA, B, seed=31337, workers=6, pool=pool
            )
            arms = []
            for sd in seeds:
                arms.append(
                    oracle_cal_fwer(n, sig, a_star, draws, sd, workers=6, pool=pool)
                )
            z = max(abs(a["z"]) for a in arms)
            rec["by_B"][str(B)] = {
                "alpha_star": f64(a_star),
                "F_at_alpha_star": f64(F(a_star)),
                "F_at_alpha": f64(F(ALPHA)),
                "arms": arms,
                "worst_abs_z": f64(z),
                "seconds": time.time() - t0,
            }
            print(
                f"{label}: B={B} alpha*={a_star:.5f} "
                f"worst|z|={z:.2f} "
                f"(resid {arms[0]['residual']:+.5f}) "
                f"{time.time() - t0:.0f}s",
                flush=True,
            )
        res["cases"][label] = rec
    # predicted sd of alpha* at each B, measured from independent calibrations
    pred = {}
    _label, n_l, sig_l = CASES[0]
    n = np.asarray(n_l, float)
    sig = np.asarray(sig_l, float)
    for B in B_list:
        vals = [
            calibrate_par(n, sig, ALPHA, B, seed=1000 + 17 * r, workers=6, pool=pool)[0]
            for r in range(6)
        ]
        pred[str(B)] = {
            "sd_alpha_star": f64(np.std(vals, ddof=1)),
            "mean": f64(np.mean(vals)),
        }
    res["alpha_star_sd_vs_B"] = pred
    dump_json(res, os.path.join(OUTDIR, "calibration_B.json"))
    print("wrote", os.path.join(OUTDIR, "calibration_B.json"))


if __name__ == "__main__":
    main()
