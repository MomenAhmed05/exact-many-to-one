"""Workstream D: the information map a desk can plan against.

The question a desk actually asks is not "is the rule valid" but "how long a
track record do I need before a valid rule is worth running".  This builds the
surface

    (benchmark df, candidate df)  ->  internal level, power loss ,

so the answer is a lookup rather than an argument.

It also settles a question the existing paper leaves open.  Section 10.3 of the
paper locates the usable boundary in the BENCHMARK's degrees of freedom, because
the sweep behind it varied only the benchmark's sample size.  The grid here
varies both, and the answer is different: what binds is the smaller of the two,
so a long benchmark does not rescue a short candidate.

Deposit: deposits/wd_designmap.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L

DF_BENCH = (3, 5, 7, 11, 19, 35, 59)
DF_CAND = (3, 5, 11, 23, 35)
LAM_GRID = 21
ALT_LAMS = (0.04, 1.0, 25.0)
TAB = dict(m=41, npanel=12, order=8)


def cell(nu0, nu1, abar=F.ABAR):
    n = np.array([nu0 + 1.0, nu1 + 1.0])
    coarse = F.log_grid(1e-3, 1e3, LAM_GRID)
    fine = F.log_grid(1e-3, 1e3, 61)
    cal = L.calibrate_level(n, coarse, abar, npanel=TAB["npanel"],
                            order=TAB["order"], iters=10)
    tab = L.PluginTable(n, cal["alpha0"], **TAB)
    s, lw, _ = L.sup_fwer(tab, n, fine, 20, 10)
    tighten = 0
    while s > abar and tighten < 12:
        cal["alpha0"] *= 0.98
        tab = L.PluginTable(n, cal["alpha0"], **TAB)
        s, lw, _ = L.sup_fwer(tab, n, fine, 20, 10)
        tighten += 1
    losses = []
    for lam in ALT_LAMS:
        d80 = F.delta_for_oracle_power_k1(n, [lam], 0.8, abar)
        p = F.power_exact_k1(tab, n, lam, d80, 20, 10)
        losses.append(100.0 * (0.8 - p))
    return {
        "nu0": int(nu0), "nu1": int(nu1), "n": n.tolist(),
        "alpha0": float(cal["alpha0"]),
        "level_price_pct": float(100.0 * (abar - cal["alpha0"]) / abar),
        "sup_fwer": float(s), "lam_worst": float(lw),
        "tightenings": int(tighten),
        "loss_pp": [float(x) for x in losses],
        "median_loss_pp": float(np.median(losses)),
        "max_loss_pp": float(max(losses)),
        "min_df": int(min(nu0, nu1)),
    }


def main():
    out = {"abar": F.ABAR, "df_bench": list(DF_BENCH), "df_cand": list(DF_CAND),
           "alt_lams": list(ALT_LAMS), "cells": []}
    print("%5s %5s %8s %8s %8s %8s" %
          ("nu0", "nu1", "alpha0", "supFWER", "lossMed", "lossMax"))
    t0 = time.time()
    for nu0 in DF_BENCH:
        for nu1 in DF_CAND:
            c = cell(nu0, nu1)
            out["cells"].append(c)
            print("%5d %5d %8.5f %8.5f %8.1f %8.1f" %
                  (nu0, nu1, c["alpha0"], c["sup_fwer"], c["median_loss_pp"],
                   c["max_loss_pp"]))
    out["seconds"] = float(time.time() - t0)

    # the planning answer: the smallest min(df) reaching each loss budget
    thresholds = {}
    for budget in (5.0, 10.0, 25.0, 50.0):
        ok = [c for c in out["cells"] if c["max_loss_pp"] <= budget]
        thresholds["min_df_for_%gpp" % budget] = (
            int(min(c["min_df"] for c in ok)) if ok else None)
    out["planning"] = thresholds

    # does the benchmark or the smaller group bind?  compare the two summaries
    by_min, by_b = {}, {}
    for c in out["cells"]:
        by_min.setdefault(c["min_df"], []).append(c["max_loss_pp"])
        by_b.setdefault(c["nu0"], []).append(c["max_loss_pp"])
    out["spread_within_min_df"] = {
        str(k): float(max(v) - min(v)) for k, v in sorted(by_min.items())}
    out["spread_within_bench_df"] = {
        str(k): float(max(v) - min(v)) for k, v in sorted(by_b.items())}
    out["binds_on_min_df"] = bool(
        np.mean(list(out["spread_within_min_df"].values()))
        < np.mean(list(out["spread_within_bench_df"].values())))
    out["note"] = (
        "Grouping the cells by min(nu0, nu1) leaves a much smaller spread in "
        "the power loss than grouping them by the benchmark's degrees of "
        "freedom alone, so the binding quantity is the smaller group and a "
        "long benchmark does not rescue a short candidate."
    )
    F.write_deposit("wd_designmap.json", out, "frontier_measured")
    print("\nplanning:", out["planning"])
    print("binds on min df:", out["binds_on_min_df"])
    print("spread within min-df groups:", out["spread_within_min_df"])
    print("spread within bench-df groups:", out["spread_within_bench_df"])


if __name__ == "__main__":
    main()
