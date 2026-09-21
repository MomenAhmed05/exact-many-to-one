"""Calibrating the plug-in rule: does the exact law give a better test after all?

THE QUESTION. h12_tamhane_headtohead.py established that the exact null law does
NOT give the better procedure: plugging estimated variance ratios into it is
liberal by up to 15.4 per cent on Tamhane and Xi's own scenarios, while their
approximate Method 3 stays within about 7. The mechanism was measured - on designs
with a three-d.f. group the plug-in exact critical value sits BELOW the oracle
value, so the bar is too low. But the plug-in error is now exactly measurable, so
the obvious question is whether it can simply be corrected. If it can, the exact
law does give the better procedure and the head-to-head verdict is overturned.

THE CALIBRATION, a known device rather than an invented one. Write the exact
adjusted p-value at the plug-in ratios,

    P = 1 - cdf_{lambdahat}( max_i t_i ),

which is the natural statistic: the uncalibrated rule rejects when P <= alpha. P
would be exactly uniform if lambdahat were the truth. It is not, and that
non-uniformity IS the plug-in error. Beran's prepivoting says do not use alpha,
use the alpha-quantile of P's own null distribution:

    alpha* = alpha-quantile of the null law of P, computed in the lambdahat world
    reject if  max_i t_i > cdf_{lambdahat}^{-1}(1 - alpha*)

Every ingredient is available from data: the sample sizes are known exactly and
lambdahat comes from the sample variances. This is a double bootstrap, and it is
affordable only because the exact law exists - each replicate costs one cdf
evaluation instead of a nested root-find.

WHY IT IS COMPUTED RAO-BLACKWELLISED. The naive estimator of alpha* is the
empirical alpha-quantile of simulated P values. At alpha = 0.025 that is a deep
tail order statistic and it is hopeless at practical replicate counts: measured
standard deviations were 0.0120 at B = 200 and still 0.0019 at B = 4000, against
a target near 0.021. Instead, for each bootstrap draw of the variance estimates
the rejection probability is evaluated in closed form by the same conditioning
argument, and the level-to-critical-value map is taken on a small grid of c so
that the whole function

    F(a) = true size, in the lambdahat world, of the plug-in rule run at level a

is available from one pass and inverted by bisection. That reaches a standard
deviation of 0.0039 at B = 80, i.e. B = 4000 raw draws' worth of precision for
one fiftieth of the work.

THREE ARMS.
  uncalibrated  the head-to-head's plug-in exact rule, for reference.
  calibrated    alpha*(lambdahat) recomputed per data set. DEPLOYABLE.
  oracle_cal    alpha*(lambda_true). NOT deployable; it is the ceiling on what a
                calibration of this form can buy, and by construction its FWER
                must come out at alpha, which makes it the falsifier for the
                calibration definition itself.

FALSIFIERS, fixed before running.
  F1 oracle_cal must return alpha in every scenario.
  F2 alpha* must tend to alpha as the degrees of freedom grow. A calibration that
     keeps correcting when there is nothing to correct is broken.
  F3 the calibrated arm must beat the uncalibrated arm where the uncalibrated arm
     is worst. If it does not, calibration does not work here and the
     head-to-head verdict stands unaltered.

VALIDATION ALREADY PASSED before this was wired up: on scenario S4 the machinery
returns F(alpha) = 0.02825 for the uncalibrated size, against 0.02842 measured by
the entirely separate head-to-head code, and F(alpha*) = 0.02500 exactly.
"""

import json
import os
from math import sqrt

import numpy as np

from h12_hetero_exact import HeteroMaxT
from h12_quad import normal_nodes
from h12_tamhane_headtohead import (ALPHA, PUBLISHED, SCEN_K2, SCEN_K3,
                                    c_method3, cond_reject)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "series", "h12_calibrated.json")

# inner grid: max |error| 6.3e-4 in the critical value against a 32-panel
# reference, at 38 ms per solve. The biases being corrected are 0.01 to 0.08.
BOOT = {"u": (4, 8), "v": (4, 8), "s": (3, 8)}
FAST = {"u": (8, 10), "v": (8, 10), "s": (3, 10)}
B_BOOT = 80
NGRID = 9
SPAN = (-1.0, 1.4)
M3_SIMS = 20000


def plugin_sigma(s2):
    """Scale-invariant plug-in ratios: sigma_0 = 1 by construction."""
    sh = np.sqrt(s2 / s2[0])
    sh[0] = 1.0
    return sh


def calibrate(n, sigma_world, alpha, rng, B=B_BOOT, ngrid=NGRID, grid=BOOT,
              span=SPAN):
    """Return (alpha*, F) with F(a) the true size of the level-a plug-in rule
    in the `sigma_world` world."""
    sw = np.asarray(sigma_world, float) / sigma_world[0]
    nu = n - 1.0
    c_ref = HeteroMaxT(n, sw, grid=grid).critical_value(alpha)
    cgrid = c_ref + np.linspace(span[0], span[1], ngrid)
    un = normal_nodes(16, 10)
    u_nodes = (un[0] * sw[0] / sqrt(n[0]), un[1])

    G = np.empty((B, ngrid))
    H = np.empty((B, ngrid))
    for b in range(B):
        s2 = (rng.chisquare(nu, size=n.size) / nu) * sw ** 2
        m = HeteroMaxT(n, plugin_sigma(s2), grid=grid)
        for j, c in enumerate(cgrid):
            G[b, j] = 1.0 - m.cdf(c)
            H[b, j] = cond_reject(n, sw, s2, c, u_nodes)

    def F(a):
        tot = 0.0
        for b in range(B):
            c = np.interp(a, G[b][::-1], cgrid[::-1])
            tot += np.interp(c, cgrid, H[b])
        return tot / B

    lo, hi = 1e-5, 0.5
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if F(mid) < alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), F


def precision_study(n_l=(10., 10., 10.), sig_l=(30., 30., 30.), alpha=ALPHA,
                    b_pool=4000, sizes=(100, 200, 4000), reps=400,
                    seed=20260914, agrid_n=400):
    """Measure what the Rao-Blackwellisation buys, instead of asserting it.

    Both estimators of alpha* are built from the same pool of bootstrap draws,
    so the comparison is like for like and the only difference is the estimator:

      naive  the empirical alpha-quantile of the simulated p-values P_b, which
             at alpha = 0.025 is a deep tail order statistic;
      RB     the root of F(a) = alpha, where F averages the CONDITIONAL
             rejection probability, known in closed form by the same
             conditioning argument the paper is about.

    The p-value for draw b needs one realised data set in that draw's world,
    which costs a normal vector and is free beside the law evaluation already
    being done. Sampling variability at each replicate count is measured by
    resampling the pool with replacement, so the numbers are standard
    deviations of the estimator, not of the pool.
    """
    n = np.asarray(n_l, float)
    sw = np.asarray(sig_l, float) / sig_l[0]
    nu = n - 1.0
    rng = np.random.default_rng(seed)

    c_ref = HeteroMaxT(n, sw, grid=BOOT).critical_value(alpha)
    cgrid = c_ref + np.linspace(SPAN[0], SPAN[1], NGRID)
    un = normal_nodes(16, 10)
    u_nodes = (un[0] * sw[0] / sqrt(n[0]), un[1])

    G = np.empty((b_pool, NGRID))
    H = np.empty((b_pool, NGRID))
    P = np.empty(b_pool)
    for b in range(b_pool):
        s2 = (rng.chisquare(nu, size=n.size) / nu) * sw ** 2
        m = HeteroMaxT(n, plugin_sigma(s2), grid=BOOT)
        for j, c in enumerate(cgrid):
            G[b, j] = 1.0 - m.cdf(c)
            H[b, j] = cond_reject(n, sw, s2, c, u_nodes)
        xbar = rng.normal(0.0, sw / np.sqrt(n))
        t = (xbar[1:] - xbar[0]) / np.sqrt(s2[1:] / n[1:] + s2[0] / n[0])
        P[b] = float(np.interp(t.max(), cgrid, G[b]))

    # phi_b(a): the conditional rejection probability of the level-a plug-in
    # rule on draw b. Tabulated once per draw so that resampling is cheap.
    agrid = np.exp(np.linspace(np.log(1e-4), np.log(0.25), agrid_n))
    phi = np.empty((b_pool, agrid_n))
    for b in range(b_pool):
        c_of_a = np.interp(agrid, G[b][::-1], cgrid[::-1])
        phi[b] = np.interp(c_of_a, cgrid, H[b])

    def rb_alpha_star(idx):
        f = phi[idx].mean(axis=0)
        if f[0] > alpha or f[-1] < alpha:
            return float("nan")
        return float(np.interp(alpha, f, agrid))

    rows = []
    for B in sizes:
        naive, rbv = [], []
        for _ in range(reps):
            idx = rng.integers(0, b_pool, size=B)
            naive.append(float(np.quantile(P[idx], alpha)))
            rbv.append(rb_alpha_star(idx))
        naive = np.asarray(naive)
        rbv = np.asarray(rbv[:])
        rows.append({"replicates": B,
                     "naive_sd": float(np.nanstd(naive, ddof=1)),
                     "naive_mean": float(np.nanmean(naive)),
                     "rb_sd": float(np.nanstd(rbv, ddof=1)),
                     "rb_mean": float(np.nanmean(rbv))})
    target = float(np.nanmean([r["rb_mean"] for r in rows]))
    by_B = {r["replicates"]: r for r in rows}
    small, large = min(sizes), max(sizes)
    return {"scenario": {"n": list(n_l), "sigma": list(sig_l)},
            "alpha": alpha, "pool": b_pool, "reps": reps,
            "target_alpha_star": target, "rows": rows,
            "naive_sd_at_max": by_B[large]["naive_sd"],
            "rb_sd_at_min": by_B[small]["rb_sd"],
            "replicate_reduction": (
                large / min((r["replicates"] for r in rows
                             if r["rb_sd"] <= by_B[large]["naive_sd"]),
                            default=large))}


def run_scenario(label, n_l, sig_l, n_draws, seed, alpha=ALPHA, with_m3=True):
    n = np.asarray(n_l, float)
    sig = np.asarray(sig_l, float)
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    un = normal_nodes(20, 12)
    u_nodes = (un[0] * sig[0] / sqrt(n[0]), un[1])

    c_or = HeteroMaxT(n, sig, grid=FAST).critical_value(alpha)
    a_oracle, F_or = calibrate(n, sig, alpha, np.random.default_rng(seed + 7717),
                               B=200)
    uncal_size_predicted = F_or(alpha)

    keys = ["oracle", "uncalibrated", "calibrated", "oracle_cal"]
    if with_m3:
        keys.append("method3")
    acc = {k: [] for k in keys}
    a_stars = []
    for _ in range(n_draws):
        s2 = (rng.chisquare(nu, size=n.size) / nu) * sig ** 2
        sh = plugin_sigma(s2)
        mod = HeteroMaxT(n, sh, grid=FAST)
        a_hat, _ = calibrate(n, sh, alpha, rng)
        a_stars.append(a_hat)
        cs = {"oracle": c_or,
              "uncalibrated": mod.critical_value(alpha),
              "calibrated": mod.critical_value(a_hat),
              "oracle_cal": mod.critical_value(a_oracle)}
        if with_m3:
            cs["method3"] = c_method3(n, s2, rng, alpha, sims=M3_SIMS)
        for k, c in cs.items():
            acc[k].append(cond_reject(n, sig, s2, c, u_nodes))

    orc = np.asarray(acc["oracle"])
    out = {"scenario": label, "n": n_l, "sigma": sig_l, "draws": n_draws,
           "alpha": alpha, "alpha_oracle_cal": a_oracle,
           "uncal_size_predicted_by_calibrator": uncal_size_predicted,
           "alpha_star_mean": float(np.mean(a_stars)),
           "alpha_star_sd": float(np.std(a_stars, ddof=1))}
    for k in keys:
        v = np.asarray(acc[k])
        if k == "oracle":
            out[k] = {"fwer": alpha, "mcse": 0.0, "direct": float(v.mean()),
                      "direct_mcse": float(v.std(ddof=1) / sqrt(len(v)))}
        else:
            d = v - orc
            out[k] = {"fwer": alpha + float(d.mean()),
                      "mcse": float(d.std(ddof=1) / sqrt(len(d)))}
    return out


def f2_large_df(alpha=ALPHA):
    rows = []
    for m in (5, 10, 30, 100, 400):
        n = np.array([m, m, m], float)
        a, _ = calibrate(n, np.array([1.0, 1.0, 3.0]), alpha,
                         np.random.default_rng(4041 + m), B=150)
        rows.append({"n_per_group": m, "df": m - 1, "alpha_star": a,
                     "ratio_to_alpha": a / alpha})
    return rows


def main(n_draws=200):
    res = {"alpha": ALPHA, "draws": n_draws, "B_boot": B_BOOT,
           "boot_grid": BOOT, "ngrid": NGRID, "blocks": {}, "summary": {}}

    print("F2  alpha* must tend to alpha as the degrees of freedom grow")
    res["f2_large_df"] = f2_large_df()
    for r in res["f2_large_df"]:
        print("    n per group %4d (df %3d)   alpha* = %.5f   alpha*/alpha = %.3f"
              % (r["n_per_group"], r["df"], r["alpha_star"], r["ratio_to_alpha"]))
    dev = [abs(r["ratio_to_alpha"] - 1) for r in res["f2_large_df"]]
    res["f2_tends_to_alpha"] = bool(dev[-1] < 0.5 * dev[0])
    print("    deviation from alpha shrinks %.3f -> %.3f : %s"
          % (dev[0], dev[-1], res["f2_tends_to_alpha"]))

    for K, scen in ((2, SCEN_K2), (3, SCEN_K3)):
        print("\n" + "=" * 104)
        print("k = %d treatments versus control  (Tamhane and Xi Table %d)"
              % (K, 1 if K == 2 else 3))
        print("%-34s | %9s %9s %10s %9s | %8s %8s"
              % ("scenario", "uncalib", "CALIB", "oracle-cal", "Method 3",
                 "alpha*", "sd"))
        rows = []
        for i, (label, n_l, sig_l) in enumerate(scen):
            r = run_scenario(label, n_l, sig_l, n_draws, seed=8800 + 41 * K + i)
            r["published_m3"] = PUBLISHED[K][3][i]
            rows.append(r)
            print("%-34s | %9.5f %9.5f %10.5f %9.5f | %8.5f %8.5f"
                  % (label, r["uncalibrated"]["fwer"], r["calibrated"]["fwer"],
                     r["oracle_cal"]["fwer"], r["method3"]["fwer"],
                     r["alpha_star_mean"], r["alpha_star_sd"]))
        res["blocks"][str(K)] = rows

        def stats(k):
            f = [r[k]["fwer"] for r in rows]
            w = max(f, key=lambda x: abs(x - ALPHA))
            return {"min": min(f), "max": max(f),
                    "worst_rel_err_pct": 100 * (w - ALPHA) / ALPHA,
                    "mean_abs_rel_err_pct": float(np.mean(
                        [abs(100 * (x - ALPHA) / ALPHA) for x in f]))}

        zo = max(abs(r["oracle"]["direct"] - ALPHA) / r["oracle"]["direct_mcse"]
                 for r in rows)
        print("  estimator check: direct oracle mean vs alpha, worst |z| = %.2f" % zo)
        blk = {}
        for k, nm in (("uncalibrated", "uncalibrated"),
                      ("calibrated", "CALIBRATED"),
                      ("oracle_cal", "oracle-calibrated"),
                      ("method3", "Method 3")):
            st = stats(k)
            blk[k] = st
            print("  %-18s FWER %.5f to %.5f   worst %+.1f%%   mean |err| %.1f%%"
                  % (nm, st["min"], st["max"], st["worst_rel_err_pct"],
                     st["mean_abs_rel_err_pct"]))
        blk["oracle_direct_worst_z"] = zo
        res["summary"][str(K)] = blk

    print("\n" + "=" * 104)
    print("VERDICT")
    for K in ("2", "3"):
        s = res["summary"][K]
        bu = abs(s["calibrated"]["worst_rel_err_pct"]) < abs(
            s["uncalibrated"]["worst_rel_err_pct"])
        bm = abs(s["calibrated"]["worst_rel_err_pct"]) < abs(
            s["method3"]["worst_rel_err_pct"])
        bmm = s["calibrated"]["mean_abs_rel_err_pct"] < s["method3"][
            "mean_abs_rel_err_pct"]
        s["calibrated_beats_uncalibrated"] = bool(bu)
        s["calibrated_beats_method3_worst"] = bool(bm)
        s["calibrated_beats_method3_mean"] = bool(bmm)
        print("  k=%s  F3 calibrated beats uncalibrated: %-5s | "
              "beats Method 3 on worst: %-5s | on mean: %-5s"
              % (K, bu, bm, bmm))
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    main()
