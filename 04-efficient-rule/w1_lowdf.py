"""W1, Question 3: the low-degrees-of-freedom regime as a first-class result.

Section 10's headline is that above eleven degrees of freedom on the benchmark
the residual loss is not removable by any threshold function of the estimated
ratio.  The same computation says something else that was reported as an aside:
at three, five and seven degrees of freedom there IS recoverable power, 4.24,
1.30 and 0.49 points respectively.  A desk that reads only the headline will
deploy the calibrated table at a sample size where a different table is
measurably better.

This module maps that regime properly.  Both group sizes move, not just the
benchmark, because the cost is set by the smaller of the two and a design map
over one axis cannot see that.  For every cell it reports:

  * the calibrated plug-in's loss against the known-ratio oracle;
  * the optimised threshold function's loss;
  * the gain, which is what the desk is leaving on the table;
  * where the gain comes from, measured as the shape change in the threshold.

and from the grid it derives the operational statement: the minimum degrees of
freedom below which the shipped table should not be used unaltered.

Deposit: deposits/w1_lowdf.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L
import wc_frontier as WC

ABAR = F.ABAR
LAM_DOMAIN = (1e-3, 1e3)
ALT_LAMS = (0.04, 0.25, 1.0, 4.0, 25.0)
GAIN_MATERIAL_PP = 0.5      # declared in advance: what counts as worth acting on


def cell(nu0, nu1, nbins=36):
    """One (nu0, nu1) design: calibrated rule against the optimised one."""
    t0 = time.time()
    n = np.array([nu0 + 1.0, nu1 + 1.0], float)
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], 0.8)) for l in ALT_LAMS]

    coarse = F.log_grid(*LAM_DOMAIN, 25)
    cal = L.calibrate_level(n, coarse, ABAR, npanel=12, order=8, iters=10)
    tab = cal["table"]
    plug = [F.power_exact_k1(tab, n, l, d, 24, 12) for l, d in alts]

    r = WC.solve_frontier(n, alts, nbins=nbins, c_init=tab)
    rule = WC.step_rule(r["c"], r["edges"])
    sup, lw, _ = F.sup_fwer_k1(rule, n, F.log_grid(*LAM_DOMAIN, 97), 24, 12)
    scale, tries = 1.0, 0
    while sup > ABAR and tries < 25:
        scale *= 1.004
        rule = WC.step_rule(r["c"] * scale, r["edges"])
        sup, lw, _ = F.sup_fwer_k1(rule, n, F.log_grid(*LAM_DOMAIN, 97), 24, 12)
        tries += 1
    opt = [F.power_exact_k1(rule, n, l, d, 24, 12) for l, d in alts]

    # where the gain comes from: how much the optimised threshold tilts across
    # the ratio range relative to the calibrated one, and in which direction
    mid = np.concatenate([[r["edges"][0] / 3.0],
                          np.sqrt(r["edges"][:-1] * r["edges"][1:]),
                          [r["edges"][-1] * 3.0]])
    c_opt = np.asarray(rule(mid), float)
    c_cal = np.asarray(tab(mid), float)
    d = c_opt - c_cal
    lo = float(np.mean(d[:len(d) // 3]))
    hi = float(np.mean(d[-len(d) // 3:]))

    return {
        "nu0": int(nu0), "nu1": int(nu1), "n": n.tolist(),
        "alpha0_plugin": float(cal["alpha0"]),
        "plugin_loss_pp": float(100 * (0.8 - np.mean(plug))),
        "optimal_loss_pp": float(100 * (0.8 - np.mean(opt))),
        "gain_pp": float(100 * (np.mean(opt) - np.mean(plug))),
        "material": bool(100 * (np.mean(opt) - np.mean(plug))
                         >= GAIN_MATERIAL_PP),
        "sup_fwer": float(sup), "lam_worst": float(lw),
        "shrink_scale": float(scale),
        "valid": bool(sup <= ABAR),
        "tilt_low_ratio": lo, "tilt_high_ratio": hi,
        "tilt_asymmetry": float(hi - lo),
        "max_abs_threshold_change": float(np.max(np.abs(d))),
        "seconds": float(time.time() - t0),
    }


def main():
    t0 = time.time()
    dfs = (3, 4, 5, 7, 9, 11)
    rows = []
    print("%4s %4s %10s %10s %8s %9s %8s" %
          ("nu0", "nu1", "plugLoss", "optLoss", "gain", "tiltAsym", "sec"))
    for nu0 in dfs:
        for nu1 in dfs:
            r = cell(nu0, nu1)
            rows.append(r)
            print("%4d %4d %10.2f %10.2f %8.2f %9.3f %8.0f"
                  % (nu0, nu1, r["plugin_loss_pp"], r["optimal_loss_pp"],
                     r["gain_pp"], r["tilt_asymmetry"], r["seconds"]))

    mat = [r for r in rows if r["material"]]
    # the operational rule: the largest min(nu0, nu1) at which the gain is
    # still material, so a desk above it can deploy the shipped table
    if mat:
        worst_min_df = max(min(r["nu0"], r["nu1"]) for r in mat)
    else:
        worst_min_df = 0
    by_min = {}
    for r in rows:
        k = min(r["nu0"], r["nu1"])
        by_min.setdefault(k, []).append(r["gain_pp"])

    out = {
        "abar": ABAR, "df_grid": list(dfs), "alt_lams": list(ALT_LAMS),
        "material_threshold_pp": GAIN_MATERIAL_PP,
        "cells": rows,
        "max_gain_pp": float(max(r["gain_pp"] for r in rows)),
        "n_material_cells": len(mat),
        "n_cells": len(rows),
        "max_min_df_with_material_gain": int(worst_min_df),
        "operational_rule": (
            "the shipped table is deployable unaltered when the SMALLER group "
            "has more than %d degrees of freedom; at or below that a "
            "design-specific threshold function is worth %0.2f points or more"
            % (worst_min_df, GAIN_MATERIAL_PP)),
        "mean_gain_by_min_df": dict(
            (str(k), float(np.mean(v))) for k, v in sorted(by_min.items())),
        "all_valid": bool(all(r["valid"] for r in rows)),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_lowdf.json", out, "frontier_measured")
    print("\nmaterial in %d of %d cells; shipped table is safe above "
          "min(nu) = %d" % (len(mat), len(rows), worst_min_df))
    print("mean gain by min(nu):", out["mean_gain_by_min_df"])
    print("all valid: %s ; %.0f s" % (out["all_valid"], out["seconds"]))


if __name__ == "__main__":
    main()
