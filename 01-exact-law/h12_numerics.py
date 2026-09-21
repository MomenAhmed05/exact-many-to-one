"""The numerical-realisation facts the paper states, computed and deposited.

Section 4 of the paper quotes the failure of the textbook rules, the geometry
of the hidden step and the cost of an evaluation; Section 6.3 quotes a worked
desk example. Those numbers previously lived only in module docstrings and in
the prose, so they could not be checked against anything. This script computes
each of them and writes them to series/h12_numerics.json, which the paper's
verification script then reads.

Nothing here is used by the method. It is diagnostic: it reproduces the two
ways the classical rules fail, so that the claim "they fail" is reproducible
rather than asserted.
"""

import json
import os
import time

import numpy as np
from scipy.special import gamma as gamma_fn
from scipy.special import roots_genlaguerre
from scipy.stats import chi2

import h12_hetero_exact as HX
from h12_hetero_compare import c_onefactor_t, c_smm, factor_loadings, welch_df
from h12_calibrated import precision_study
from h12_hetero_exact import GRID, ZMAX, HeteroMaxT

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "series", "h12_numerics.json")


# ----------------------------------------------------------------- Laguerre
def chi2_over_df_laguerre(nu, n_nodes=24, order=None, zmax=None):
    """Generalized Gauss-Laguerre nodes for chi2_nu/nu, the textbook choice.

    X = chi2_nu/nu has density proportional to x^(nu/2-1) exp(-nu x/2), so with
    t = nu x / 2 the expectation is a generalized Laguerre integral of index
    alpha = nu/2 - 1. The weights of that rule sum to Gamma(alpha+1) = Gamma(nu/2),
    which is where the overflow comes from: it is a property of the rule, not of
    the implementation. `order` and `zmax` are accepted and ignored so that this
    can be substituted for the composite rule of the same name.
    """
    a = 0.5 * nu - 1.0
    t, w = roots_genlaguerre(n_nodes, a)
    tot = gamma_fn(0.5 * nu)
    return 2.0 * t / nu, w / tot


def laguerre_overflow():
    """The degrees of freedom at which the Laguerre weights stop summing."""
    rows = []
    for nu in (300, 340, 342, 343, 344, 350, 400):
        g = float(gamma_fn(0.5 * nu))
        rows.append({"nu": nu, "gamma_nu_over_2": g if np.isfinite(g) else None,
                     "finite": bool(np.isfinite(g))})
    last_finite = max(r["nu"] for r in rows if r["finite"])
    first_infinite = min(r["nu"] for r in rows if not r["finite"])
    return {"rows": rows, "last_finite_nu": last_finite,
            "first_infinite_nu": first_infinite}


def laguerre_refinement(nu0=5, alpha=0.05):
    """Refine the chi-square rule and watch the answer fail to settle.

    The composite rule is swapped out for the Laguerre one at the module level,
    so everything else about the evaluation is held fixed and the only thing
    changing is the chi-square quadrature and its node count.
    """
    n = [nu0 + 1, 21, nu0 + 1, nu0 + 1]
    sigma = [1.0, 0.25, 4.0, 1.0]
    saved = HX.chi2_over_df_nodes
    rows = []
    try:
        for k in (24, 48, 96, 192, 384):
            HX.chi2_over_df_nodes = (
                lambda nu, *a, _k=k, **kw: chi2_over_df_laguerre(nu, _k))
            try:
                c = float(HeteroMaxT(n, sigma).critical_value(alpha))
                ok = bool(np.isfinite(c))
            except Exception:
                c, ok = None, False
            rows.append({"nodes": k, "critical_value": c if ok else None,
                         "finite": ok})
    finally:
        HX.chi2_over_df_nodes = saved
    vals = [r["critical_value"] for r in rows if r["finite"]]
    return {"design": "n=%s sigma=%s" % (n, sigma), "rows": rows,
            "finite_values": vals,
            "spread": (max(vals) - min(vals)) if len(vals) > 1 else None,
            "monotone": bool(vals == sorted(vals) or vals == sorted(vals,
                                                                   reverse=True)),
            "first_nonfinite_nodes": next((r["nodes"] for r in rows
                                           if not r["finite"]), None)}


# ------------------------------------------------------------ step geometry
def step_geometry(n=(10., 4000., 4000., 4000.), sigma=(1., 1e-3, 1e-3, 1e-3)):
    """How much narrower the transition is than a base quadrature panel.

    Both are expressed in u, which is measured in the benchmark's own standard
    deviations by the z-scale the outer rule runs on.
    """
    n = np.asarray(n, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    sd0 = sigma[0] / np.sqrt(n[0])
    panel_z = 2.0 * ZMAX / GRID["u"][0]
    panel_u = panel_z * sd0
    step_u = sigma[1] / np.sqrt(n[1])
    return {"design": {"n": n.tolist(), "sigma": sigma.tolist()},
            "u_panels": GRID["u"][0], "zmax": ZMAX,
            "benchmark_sd": float(sd0),
            "panel_width_u": float(panel_u),
            "step_width_u": float(step_u),
            "ratio": float(panel_u / step_u),
            "ratio_rounded_to_hundreds": int(round(panel_u / step_u, -2))}


# -------------------------------------------------------------------- cost
def timings(alpha=0.05, repeats=1):
    designs = [("balanced n=10, homoscedastic", [10.] * 4, [1.] * 4),
               ("balanced n=20, sd 1:3:3:3", [20.] * 4, [1., 3., 3., 3.]),
               ("small benchmark n0=6, sd 1:.25:4:1", [6., 21., 6., 6.],
                [1., 0.25, 4., 1.]),
               ("large df n0=400, sd 1:2:.5:3", [400., 60., 200., 60.],
                [1., 2., 0.5, 3.])]
    rows = []
    for name, n, s in designs:
        best = None
        for _ in range(repeats):
            t0 = time.perf_counter()
            c = float(HeteroMaxT(n, s).critical_value(alpha))
            dt = time.perf_counter() - t0
            best = dt if best is None else min(best, dt)
        rows.append({"design": name, "critical_value": c, "seconds": best})
    secs = [r["seconds"] for r in rows]
    return {"rows": rows, "min_seconds": min(secs), "max_seconds": max(secs)}


# ----------------------------------------------------------- worked example
def worked_example(alpha=0.05):
    """One benchmark at 12 per cent and five candidates at 8 to 30 per cent."""
    n = np.array([36.] * 6)
    sigma = np.array([0.12, 0.08, 0.10, 0.15, 0.20, 0.30])
    K = 5
    law = HeteroMaxT(n, sigma)
    c_exact = float(law.critical_value(alpha))

    # The two rules are priced exactly as the arms in h12_hetero_compare.py
    # price them, so this example and Table 2 are the same two procedures:
    # T3 takes its own constant per comparison at that comparison's Welch
    # degrees of freedom, and Dunnett's route keeps the loadings but carries
    # the pooled degrees of freedom across.
    dfs = np.array([welch_df(n, sigma, i) for i in range(1, K + 1)])
    mean_df = float(dfs.mean())
    h = factor_loadings(n, sigma)

    c_t3 = np.array([c_smm(K, d, alpha) for d in dfs])
    pooled_df = float((n - 1).sum())
    c_dun = float(c_onefactor_t(h, pooled_df, alpha))

    fwer_t3 = float(1.0 - law.cdf(c_t3))
    fwer_dun = float(1.0 - law.cdf(c_dun))
    fwer_exact = float(1.0 - law.cdf(c_exact))
    return {"n": n.tolist(), "sigma_annualised": sigma.tolist(),
            "alpha": alpha, "loadings": [float(x) for x in h],
            "homoscedastic_loading": float(np.sqrt(0.5)),
            "welch_dfs": [float(x) for x in dfs], "mean_welch_df": mean_df,
            "pooled_df": pooled_df,
            "c_exact": c_exact, "c_t3_smm": [float(x) for x in c_t3],
            "c_dunnett_carried": c_dun,
            "fwer_exact": fwer_exact, "fwer_t3_smm": fwer_t3,
            "fwer_dunnett_carried": fwer_dun,
            "t3_budget_shortfall_pct": 100.0 * (alpha - fwer_t3) / alpha,
            "dunnett_budget_overshoot_pct": 100.0 * (fwer_dun - alpha) / alpha}


def cross_run_consistency():
    """Tables 3 and 7 are two independent runs. Measure their agreement.

    The head to head runs 800 data sets per scenario and the calibration study
    runs 200, so the two tables report the same quantity for the uncalibrated
    exact rule and for Method 3 from different samples. Rather than leave the
    reader to notice, compare them arm by arm against the combined standard
    error the two runs themselves report.
    """
    here = os.path.dirname(OUT)
    h2h = json.load(open(os.path.join(here, "h12_tamhane_headtohead.json")))
    cal = json.load(open(os.path.join(here, "h12_calibrated.json")))
    pairs = (("exact", "uncalibrated"), ("method3", "method3"))
    rows = []
    for K in ("2", "3"):
        by_lab = {r["scenario"]: r for r in cal["blocks"][K]}
        for r3 in h2h["blocks"][K]:
            r7 = by_lab.get(r3["scenario"])
            if r7 is None:
                continue
            for a3, a7 in pairs:
                d = r3[a3]["fwer"] - r7[a7]["fwer"]
                se = float(np.hypot(r3[a3]["mcse"], r7[a7]["mcse"]))
                rows.append({"K": int(K), "scenario": r3["scenario"],
                             "arm": a7, "headtohead_fwer": r3[a3]["fwer"],
                             "calibration_fwer": r7[a7]["fwer"],
                             "difference": d, "combined_mcse": se,
                             "z": d / se if se > 0 else None})
    zs = [abs(r["z"]) for r in rows if r["z"] is not None]
    return {"headtohead_draws": h2h["draws"], "calibration_draws": cal["draws"],
            "rows": rows, "n_comparisons": len(rows),
            "worst_abs_z": max(zs), "mean_abs_z": float(np.mean(zs)),
            "fraction_within_3z": float(np.mean([z < 3.0 for z in zs]))}


def reimplementation_faithfulness():
    """How far the reimplementation sits from the published table, in its own
    standard errors.

    Both sides carry sampling error: the published figures come from 100,000
    data sets per cell, and the recomputation from 800 with a Rao-Blackwellised
    estimator whose standard error is deposited per cell. The right yardstick
    is the two combined, and both halves of it are read from deposits rather
    than assumed.
    """
    here = os.path.dirname(OUT)
    h2h = json.load(open(os.path.join(here, "h12_tamhane_headtohead.json")))
    arms = ("method1", "method2", "method3")
    diffs, mcse = [], []
    for K in ("2", "3"):
        for r in h2h["blocks"][K]:
            pub = r["published"]
            for j, arm in enumerate(arms):
                key = str(j + 1)          # published column 0 is not a method
                if key not in pub:
                    continue
                diffs.append(r[arm]["fwer"] - pub[key])
                mcse.append(r[arm]["mcse"])
    diffs, mcse = np.asarray(diffs), np.asarray(mcse)
    pub_se = float(np.sqrt(0.025 * 0.975 / h2h["published_sims"]))
    rec_se = float(np.sqrt((mcse ** 2).mean()))
    comb = float(np.hypot(pub_se, rec_se))
    worst = float(np.abs(diffs).max())
    return {"cells": int(diffs.size),
            "published_sims": h2h["published_sims"],
            "published_se": pub_se,
            "recomputation_se_rms": rec_se,
            "recomputation_se_mean": float(mcse.mean()),
            "recomputation_se_max": float(mcse.max()),
            "combined_se": comb,
            "worst_abs_difference": worst,
            "worst_in_combined_se": worst / comb,
            "n_positive": int((diffs > 0).sum()),
            "n_negative": int((diffs < 0).sum())}


def main():
    res = {}

    print("Laguerre weights: Gamma(nu/2) and where it stops being finite")
    res["laguerre_overflow"] = laguerre_overflow()
    for r in res["laguerre_overflow"]["rows"]:
        print("    nu=%-4d finite=%s" % (r["nu"], r["finite"]))
    print("    last finite %d, first infinite %d"
          % (res["laguerre_overflow"]["last_finite_nu"],
             res["laguerre_overflow"]["first_infinite_nu"]))

    print("\nLaguerre refinement: the answer must settle, and does not")
    res["laguerre_refinement"] = laguerre_refinement()
    for r in res["laguerre_refinement"]["rows"]:
        print("    %3d nodes -> %s"
              % (r["nodes"], ("%.6f" % r["critical_value"]) if r["finite"]
                 else "non-finite"))
    print("    spread over the finite answers: %.2e"
          % res["laguerre_refinement"]["spread"])

    print("\nStep geometry")
    res["step_geometry"] = step_geometry()
    g = res["step_geometry"]
    print("    panel %.4f wide in u, step %.3e wide, ratio %.0f"
          % (g["panel_width_u"], g["step_width_u"], g["ratio"]))

    print("\nCost of one critical value")
    res["timings"] = timings()
    for r in res["timings"]["rows"]:
        print("    %-40s %.2f s" % (r["design"], r["seconds"]))
    print("    %.2f to %.2f seconds"
          % (res["timings"]["min_seconds"], res["timings"]["max_seconds"]))

    print("\nWorked desk example")
    res["worked_example"] = worked_example()
    w = res["worked_example"]
    print("    exact c = %.4f, delivers %.3f per cent"
          % (w["c_exact"], 100 * w["fwer_exact"]))
    print("    studentised max modulus delivers %.3f per cent (%.0f per cent "
          "of the budget thrown away)"
          % (100 * w["fwer_t3_smm"], w["t3_budget_shortfall_pct"]))
    print("    Dunnett carried across delivers %.3f per cent (%.0f per cent over)"
          % (100 * w["fwer_dunnett_carried"], w["dunnett_budget_overshoot_pct"]))
    print("    loadings %s vs %.3f"
          % (" ".join("%.3f" % x for x in w["loadings"]),
             w["homoscedastic_loading"]))

    print("\nWhat the Rao-Blackwellisation buys, measured")
    res["precision_study"] = precision_study(
        sizes=(50, 100, 200, 1000, 4000), reps=400)
    ps = res["precision_study"]
    print("    target alpha* %.5f on n=%s"
          % (ps["target_alpha_star"], ps["scenario"]["n"]))
    for r in ps["rows"]:
        print("    B=%-5d naive sd %.4f   conditional sd %.4f"
              % (r["replicates"], r["naive_sd"], r["rb_sd"]))
    print("    %.0f-fold reduction in the replicate count"
          % ps["replicate_reduction"])

    print("\nReimplementation against the published table")
    res["reimplementation_faithfulness"] = reimplementation_faithfulness()
    f = res["reimplementation_faithfulness"]
    print("    %d cells, worst |difference| %.5f (%d up, %d down)"
          % (f["cells"], f["worst_abs_difference"], f["n_positive"],
             f["n_negative"]))
    print("    published se %.5f, recomputation se %.5f (rms), combined %.5f"
          % (f["published_se"], f["recomputation_se_rms"], f["combined_se"]))
    print("    worst cell sits at %.2f combined standard errors"
          % f["worst_in_combined_se"])

    print("\nAgreement between the two runs that share arms")
    res["cross_run_consistency"] = cross_run_consistency()
    x = res["cross_run_consistency"]
    print("    %d comparisons, %d draws against %d"
          % (x["n_comparisons"], x["headtohead_draws"],
             x["calibration_draws"]))
    print("    worst |z| %.2f, mean |z| %.2f, within 3 se: %.0f per cent"
          % (x["worst_abs_z"], x["mean_abs_z"], 100 * x["fraction_within_3z"]))

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    main()
