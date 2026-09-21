"""Four diagnostics the head-to-head left open, each of which could change its reading.

The head-to-head concluded that Tamhane and Xi's Method 3 controls the error rate
more accurately than the plug-in exact rule. That conclusion rests on their
scenario grid, at their alpha, with their method reimplemented. Four things about
it were not established, and all four are cheap to settle.

D1  WHERE IS THE CROSSOVER? Their designs put four observations on a group, i.e.
    three degrees of freedom, where a variance estimate is nearly worthless. The
    earlier plug-in study used 5 to 19 d.f. and found the exact rule holding size.
    Both can be true, with the winner flipping somewhere between. Nobody has
    measured where. A trading desk with 36 monthly returns has 35 d.f. and would
    be far from the flip point, so this is the single most decision-relevant
    number in the whole file. Swept here by varying the control's sample size with
    everything else held fixed.

D2  IS THE CANCELLATION STORY CLEAN? The head-to-head compared Method 3's critical
    value, computed at PLUG-IN estimates, against the oracle value. That number
    therefore contains Method 3's approximation error AND its estimation error
    together, so "two errors cancel" mixes three ingredients, not two. Evaluating
    Method 3's critical value at the TRUE variances isolates its approximation
    error, and the difference from the plug-in version isolates its estimation
    error. Only then is the decomposition honest.

D3  HOW MUCH OF METHOD 3 IS ITS OWN MONTE CARLO NOISE? Its critical value is a
    simulated quantile, so it is a random variable even given the data. If that
    noise is comparable to the biases being discussed, the comparison is partly
    measuring seed choice. Checked by re-running the arm at ten times the
    replicate count.

D4  DOES THE GAP SURVIVE AT alpha = 0.05? The head-to-head ran at one-sided 0.025
    to match their tables, where relative errors are amplified. Re-checked at 0.05.
"""

import json
import os
from math import sqrt

import numpy as np

from h12_calibrated import FAST, calibrate, plugin_sigma
from h12_hetero_exact import HeteroMaxT
from h12_quad import normal_nodes
from h12_tamhane_headtohead import ALPHA, c_method3, cond_reject

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "series", "h12_crossover.json")


def arms_fwer(n_l, sig_l, n_draws, seed, alpha=ALPHA, m3_sims=20000,
              with_calib=False):
    """Paired FWER of oracle / uncalibrated / Method 3 (and optionally the
    calibrated rule), with the oracle as an exactly known control variate."""
    n = np.asarray(n_l, float)
    sig = np.asarray(sig_l, float)
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    un = normal_nodes(20, 12)
    u_nodes = (un[0] * sig[0] / sqrt(n[0]), un[1])
    c_or = HeteroMaxT(n, sig, grid=FAST).critical_value(alpha)

    keys = ["oracle", "uncalibrated", "method3"]
    if with_calib:
        keys.append("calibrated")
    acc = {k: [] for k in keys}
    for _ in range(n_draws):
        s2 = (rng.chisquare(nu, size=n.size) / nu) * sig ** 2
        mod = HeteroMaxT(n, plugin_sigma(s2), grid=FAST)
        cs = {"oracle": c_or,
              "uncalibrated": mod.critical_value(alpha),
              "method3": c_method3(n, s2, rng, alpha, sims=m3_sims)}
        if with_calib:
            a_hat, _ = calibrate(n, plugin_sigma(s2), alpha, rng)
            cs["calibrated"] = mod.critical_value(a_hat)
        for k, c in cs.items():
            acc[k].append(cond_reject(n, sig, s2, c, u_nodes))

    orc = np.asarray(acc["oracle"])
    out = {}
    for k in keys:
        v = np.asarray(acc[k])
        if k == "oracle":
            out[k] = {"fwer": alpha, "direct": float(v.mean()),
                      "direct_mcse": float(v.std(ddof=1) / sqrt(len(v)))}
        else:
            d = v - orc
            out[k] = {"fwer": alpha + float(d.mean()),
                      "mcse": float(d.std(ddof=1) / sqrt(len(d)))}
    return out


# --------------------------------------------------------------------------
def d1_crossover(n_draws=300):
    """Vary the control's sample size; everything else fixed."""
    rows = []
    for n0 in (4, 5, 6, 8, 12, 20, 36, 60):
        r = arms_fwer([n0, 13, 13], [10., 10., 50.], n_draws, seed=5100 + n0)
        e_unc = 100 * (r["uncalibrated"]["fwer"] - ALPHA) / ALPHA
        e_m3 = 100 * (r["method3"]["fwer"] - ALPHA) / ALPHA
        rows.append({"n0": n0, "df0": n0 - 1,
                     "uncalibrated": r["uncalibrated"], "method3": r["method3"],
                     "uncal_rel_err_pct": e_unc, "m3_rel_err_pct": e_m3,
                     "exact_closer": abs(e_unc) < abs(e_m3)})
    return rows


def d2_decompose(n_draws=400):
    """Split Method 3's critical-value bias into approximation and estimation."""
    rows = []
    cases = [("S4 control smallest n and sd", [4, 13, 13], [10., 10., 50.]),
             ("S5 trt smallest n, largest sd", [13, 13, 4], [10., 10., 50.]),
             ("S6 control smallest n, largest sd", [4, 13, 13], [50., 10., 10.]),
             ("S1 balanced, equal sd", [10, 10, 10], [30., 30., 30.])]
    for label, n_l, sig_l in cases:
        n = np.asarray(n_l, float)
        sig = np.asarray(sig_l, float)
        nu = n - 1.0
        rng = np.random.default_rng(6100 + len(label))
        c_or = HeteroMaxT(n, sig, grid=FAST).critical_value(ALPHA)
        # Method 3 at the TRUE variances: pure approximation error, no estimation
        s2_true = sig ** 2
        c3_true = float(np.mean([c_method3(n, s2_true, rng, ALPHA, sims=60000)
                                 for _ in range(12)]))
        ce, c3 = [], []
        for _ in range(n_draws):
            s2 = (rng.chisquare(nu, size=n.size) / nu) * sig ** 2
            ce.append(HeteroMaxT(n, plugin_sigma(s2), grid=FAST).critical_value(ALPHA))
            c3.append(c_method3(n, s2, rng, ALPHA, sims=20000))
        ce, c3 = np.asarray(ce), np.asarray(c3)
        rows.append({
            "case": label, "n": n_l, "sigma": sig_l, "c_oracle": c_or,
            "exact_plugin_bias": float(ce.mean() - c_or),
            "exact_plugin_sd": float(ce.std(ddof=1)),
            "m3_total_bias": float(c3.mean() - c_or),
            "m3_approximation_bias": c3_true - c_or,
            "m3_estimation_bias": float(c3.mean()) - c3_true,
            "m3_sd": float(c3.std(ddof=1))})
    return rows


def d3_m3_noise(n_draws=150):
    """Is Method 3's own Monte Carlo noise material at the size it is run?"""
    rows = []
    for label, n_l, sig_l in (("S4", [4, 13, 13], [10., 10., 50.]),
                              ("S6", [4, 13, 13], [50., 10., 10.])):
        rec = {"case": label}
        for sims in (20000, 100000):
            r = arms_fwer(n_l, sig_l, n_draws, seed=7100, m3_sims=sims)
            rec["sims_%d" % sims] = r["method3"]
        rec["shift"] = rec["sims_100000"]["fwer"] - rec["sims_20000"]["fwer"]
        rec["shift_in_mcse"] = abs(rec["shift"]) / rec["sims_20000"]["mcse"]
        rows.append(rec)
    return rows


def d4_alpha05(n_draws=300):
    """Does the head-to-head gap survive at a two-sided-style 5 per cent level?"""
    rows = []
    for label, n_l, sig_l in (("S4", [4, 13, 13], [10., 10., 50.]),
                              ("S5", [13, 13, 4], [10., 10., 50.]),
                              ("S6", [4, 13, 13], [50., 10., 10.]),
                              ("S1", [10, 10, 10], [30., 30., 30.])):
        r = arms_fwer(n_l, sig_l, n_draws, seed=8100, alpha=0.05)
        rows.append({"case": label,
                     "uncal_rel_err_pct": 100 * (r["uncalibrated"]["fwer"] - 0.05) / 0.05,
                     "m3_rel_err_pct": 100 * (r["method3"]["fwer"] - 0.05) / 0.05,
                     "uncalibrated": r["uncalibrated"], "method3": r["method3"]})
    return rows


def main():
    res = {"alpha": ALPHA}

    print("D1  CROSSOVER: vary the control's sample size, n=(n0,13,13), "
          "sd=(10,10,50)")
    print("    %-5s %-5s | %12s %12s | %s" % ("n0", "df0", "exact plug-in",
                                              "Method 3", "closer to nominal"))
    res["d1_crossover"] = d1_crossover()
    for r in res["d1_crossover"]:
        print("    %-5d %-5d | %+11.1f%% %+11.1f%% | %s"
              % (r["n0"], r["df0"], r["uncal_rel_err_pct"], r["m3_rel_err_pct"],
                 "EXACT" if r["exact_closer"] else "Method 3"))
    flips = [r["n0"] for r in res["d1_crossover"] if r["exact_closer"]]
    res["d1_exact_closer_from_n0"] = min(flips) if flips else None
    print("    exact is closer from n0 = %s upward" % res["d1_exact_closer_from_n0"])

    print("\nD2  DECOMPOSE Method 3's critical-value bias against the oracle value")
    print("    %-34s %8s | %9s | %9s %9s %9s"
          % ("case", "oracle c", "exact", "M3 total", "M3 appx", "M3 est"))
    res["d2_decompose"] = d2_decompose()
    for r in res["d2_decompose"]:
        print("    %-34s %8.4f | %+9.4f | %+9.4f %+9.4f %+9.4f"
              % (r["case"], r["c_oracle"], r["exact_plugin_bias"],
                 r["m3_total_bias"], r["m3_approximation_bias"],
                 r["m3_estimation_bias"]))

    print("\nD3  Method 3's own Monte Carlo noise: 20,000 vs 100,000 replicates")
    res["d3_m3_noise"] = d3_m3_noise()
    for r in res["d3_m3_noise"]:
        print("    %-4s FWER %.5f -> %.5f   shift %+.5f = %.2f MCSE"
              % (r["case"], r["sims_20000"]["fwer"], r["sims_100000"]["fwer"],
                 r["shift"], r["shift_in_mcse"]))

    print("\nD4  Does the gap survive at alpha = 0.05?")
    res["d4_alpha05"] = d4_alpha05()
    for r in res["d4_alpha05"]:
        print("    %-4s exact plug-in %+7.1f%%   Method 3 %+7.1f%%   closer: %s"
              % (r["case"], r["uncal_rel_err_pct"], r["m3_rel_err_pct"],
                 "EXACT" if abs(r["uncal_rel_err_pct"]) < abs(r["m3_rel_err_pct"])
                 else "Method 3"))

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    main()
