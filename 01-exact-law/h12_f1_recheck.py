"""F1 re-run: the calibration's own falsifier, at a replicate count that can pass it.

WHY THIS EXISTS. h12_calibrated.py carries an arm called `oracle_cal`: calibrate
at the TRUE variance ratios, then apply that level to the plug-in rule. By the
definition of alpha* its familywise error rate must be exactly alpha, which makes
it the falsifier for the calibration definition itself. In the main run it came
back at 0.01863 and 0.03507 on two scenarios, i.e. -25 per cent and +40 per cent.

That is a flaw in the check, not evidence against the calibration, and the cause
is identifiable: alpha* for that arm is estimated ONCE per scenario at B = 200
bootstrap replicates, and the same single estimate is then used for every data
set. Its error therefore passes straight through to the reported rate instead of
averaging out, unlike the deployable arm where alpha* is redrawn per data set. At
B = 200 the standard deviation of alpha* is about 0.0017 and dFWER/dalpha* is
close to 1, so an error of that size lands directly on a quantity whose target is
0.025. The arm was measuring its own estimator.

So it is rerun here two ways, both of which can fail:
  F1a  DEFINITION. Using the calibrator's own high-replicate F(a) - the true size
       of the level-a plug-in rule - check that F(alpha*) = alpha and, separately,
       that F(alpha) reproduces the uncalibrated size measured by the entirely
       independent head-to-head code.
  F1b  END TO END. Recompute alpha* at B = 600 and measure the oracle_cal arm on a
       fresh data stream. With the estimate an order of magnitude tighter, the arm
       must now land on alpha.

Also reported: dFWER/dalpha*, so the replicate count needed for any future run can
be read off rather than guessed.
"""

import json
import os
from math import sqrt

import numpy as np

from h12_calibrated import FAST, calibrate, plugin_sigma
from h12_hetero_exact import HeteroMaxT
from h12_quad import normal_nodes
from h12_tamhane_headtohead import ALPHA, cond_reject

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "series", "h12_f1_recheck.json")

# uncalibrated sizes measured by h12_tamhane_headtohead.py, for the cross-check
H2H_UNCAL = {"k2 S4": 0.02842, "k2 S5": 0.02884, "k2 S6": 0.02752,
             "k2 S1": 0.02407}
CASES = [("k2 S1", [10, 10, 10], [30., 30., 30.]),
         ("k2 S4", [4, 13, 13], [10., 10., 50.]),
         ("k2 S5", [13, 13, 4], [10., 10., 50.]),
         ("k2 S6", [4, 13, 13], [50., 10., 10.])]


def oracle_cal_fwer(n, sig, a_star, n_draws, seed, alpha=ALPHA):
    """FWER of: calibrate at the truth, apply that level to the plug-in rule."""
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    un = normal_nodes(20, 12)
    u_nodes = (un[0] * sig[0] / sqrt(n[0]), un[1])
    c_or = HeteroMaxT(n, sig, grid=FAST).critical_value(alpha)
    oc, orc = [], []
    for _ in range(n_draws):
        s2 = (rng.chisquare(nu, size=n.size) / nu) * sig ** 2
        mod = HeteroMaxT(n, plugin_sigma(s2), grid=FAST)
        oc.append(cond_reject(n, sig, s2, mod.critical_value(a_star), u_nodes))
        orc.append(cond_reject(n, sig, s2, c_or, u_nodes))
    d = np.asarray(oc) - np.asarray(orc)
    return alpha + float(d.mean()), float(d.std(ddof=1) / sqrt(len(d)))


def main(n_draws=500, b_big=600):
    res = {"alpha": ALPHA, "b_big": b_big, "draws": n_draws, "cases": []}
    print("F1a  DEFINITION: the calibrator's own F(a), at B = %d" % b_big)
    print("     %-8s %10s %10s %10s %12s" % ("case", "alpha*", "F(alpha*)",
                                             "F(alpha)", "head-to-head"))
    for label, n_l, sig_l in CASES:
        n = np.asarray(n_l, float)
        sig = np.asarray(sig_l, float)
        a_star, F = calibrate(n, sig, ALPHA, np.random.default_rng(31337),
                              B=b_big, ngrid=11)
        rec = {"case": label, "n": n_l, "sigma": sig_l, "alpha_star": a_star,
               "F_at_alpha_star": F(a_star), "F_at_alpha": F(ALPHA),
               "h2h_uncalibrated": H2H_UNCAL[label]}
        rec["F_alpha_minus_h2h"] = rec["F_at_alpha"] - rec["h2h_uncalibrated"]
        # local slope, for sizing any future calibration
        h = 0.002
        rec["dFWER_dalphastar"] = (F(a_star + h) - F(a_star - h)) / (2 * h)
        res["cases"].append(rec)
        print("     %-8s %10.5f %10.5f %10.5f %12.5f"
              % (label, a_star, rec["F_at_alpha_star"], rec["F_at_alpha"],
                 rec["h2h_uncalibrated"]))
    w_def = max(abs(r["F_at_alpha_star"] - ALPHA) for r in res["cases"])
    w_x = max(abs(r["F_alpha_minus_h2h"]) for r in res["cases"])
    print("     worst |F(alpha*) - alpha|            = %.2e" % w_def)
    print("     worst |F(alpha) - head-to-head size| = %.5f  (independent code)"
          % w_x)
    print("     dFWER/dalpha* = %s"
          % " ".join("%.2f" % r["dFWER_dalphastar"] for r in res["cases"]))

    print("\nF1b  END TO END: alpha* at B = %d, fresh data stream, %d draws"
          % (b_big, n_draws))
    print("     %-8s %10s %10s %8s %10s" % ("case", "alpha*", "FWER", "MCSE", "z"))
    for rec in res["cases"]:
        n = np.asarray(rec["n"], float)
        sig = np.asarray(rec["sigma"], float)
        f, se = oracle_cal_fwer(n, sig, rec["alpha_star"], n_draws, seed=9911)
        rec["oracle_cal_fwer"] = f
        rec["oracle_cal_mcse"] = se
        rec["oracle_cal_z"] = (f - ALPHA) / se
        print("     %-8s %10.5f %10.5f %8.5f %+10.2f"
              % (rec["case"], rec["alpha_star"], f, se, rec["oracle_cal_z"]))
    w_z = max(abs(r["oracle_cal_z"]) for r in res["cases"])
    res["summary"] = {"worst_abs_F_minus_alpha": w_def,
                      "worst_abs_F_vs_headtohead": w_x,
                      "worst_abs_z_end_to_end": w_z,
                      "f1a_passes": bool(w_def < 1e-4),
                      "f1b_passes": bool(w_z < 3.0)}
    print("     worst |z| = %.2f" % w_z)
    print("\nF1a passes: %s     F1b passes: %s"
          % (res["summary"]["f1a_passes"], res["summary"]["f1b_passes"]))
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
