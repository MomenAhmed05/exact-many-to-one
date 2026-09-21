"""THE DESK PROCEDURE: a certified, tabulated, uniformly valid threshold rule.

This is the deliverable the whole programme was for.  What a desk runs:

    1. compute the sample variances and lambda_hat_i = S_i^2 / S_0^2;
    2. look up c from the deposited table for its (n, K, abar);
    3. promote candidate i iff T_i > c.

No region, no branch-and-bound at decision time, no simulation, no seed.  The
table is computed once per design and carries a certificate:

    sup_{lambda in [1e-3, 1e3]}  FWER(lambda)  <=  abar ,

proven by the rigorous enclosure of `wb_certify.py` over the CONTINUUM of ratio
values, not over a grid.

WHY THIS IS CHEAPER THAN THE REGION RULE.  The region rule protects, on every
draw, against every ratio in an admission region; it therefore pays for a
supremum it does not need.  The guarantee is an average over draws, so the level
a draw does not spend is available to another draw.  Calibrating one internal
level against the exact unconditional constraint spends the budget where it is
worth something.

Deposit: deposits/wb_certified_rule.json  (tier unconditional_valid_certified)
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L
import wb_certify as C

LAM_LO, LAM_HI = 1e-3, 1e3
NBIN = 60
CELLS = 800
TABLE = dict(m=41, npanel=12, order=8)

DESIGNS = [
    ("df3", [4.0, 13.0]),
    ("df5", [6.0, 13.0]),
    ("df7", [8.0, 13.0]),
    ("df11", [12.0, 13.0]),
    ("df19", [20.0, 13.0]),
    ("df35", [36.0, 13.0]),
    ("df59", [60.0, 13.0]),
    ("df35_short_cand", [36.0, 5.0]),
    ("balanced_df11", [12.0, 12.0]),
    ("balanced_df35", [36.0, 36.0]),
]
ALT_LAMS = (0.04, 0.25, 1.0, 4.0, 25.0)


def step_edges(nbin=NBIN):
    return F.log_grid(LAM_LO * 3.0, LAM_HI / 3.0, nbin - 1)


def rule_from_level(n, alpha0, edges):
    tab = L.PluginTable(n, alpha0, **TABLE)
    mids = np.concatenate([[edges[0] / 3.0],
                           np.sqrt(edges[:-1] * edges[1:]),
                           [edges[-1] * 3.0]])
    return C.StepRule(tab(mids), edges)


def certified_level(n, abar=F.ABAR, iters=9, cells=CELLS):
    """The largest internal level whose step rule certifies over the continuum."""
    edges = step_edges()
    cache = {}

    def ok(a):
        if a not in cache:
            r = rule_from_level(n, a, edges)
            cert = C.certify(r, n, LAM_LO, LAM_HI, abar, m=cells,
                             max_boxes=640)
            cache[a] = (cert, r)
        return cache[a]

    hi = abar
    cert, rule = ok(hi)
    if cert["certified"]:
        return hi, cert, rule, edges
    lo = 0.004
    for _ in range(iters):
        m = 0.5 * (lo + hi)
        if ok(m)[0]["certified"]:
            lo = m
        else:
            hi = m
    cert, rule = ok(lo)
    return lo, cert, rule, edges


def price(rule, n, abar=F.ABAR):
    rows = []
    for lam in ALT_LAMS:
        for t in (0.5, 0.8, 0.95):
            dl = F.delta_for_oracle_power_k1(n, [lam], t, abar)
            p = F.power_exact_k1(rule, n, lam, dl, 24, 12)
            rows.append({"lam": float(lam), "target": float(t),
                         "delta": float(dl), "rule_power": float(p),
                         "loss_pp": float(100.0 * (t - p))})
    return rows


def main():
    out = {"abar": F.ABAR, "lam_domain": [LAM_LO, LAM_HI], "nbins": NBIN,
           "cells_per_axis": CELLS, "table": TABLE,
           "alt_lams": list(ALT_LAMS), "designs": []}
    print("%-16s %9s %9s %11s %8s %8s %6s" %
          ("design", "alpha0", "supUpper", "certified", "loss80m", "loss80x",
           "sec"))
    for name, n in DESIGNS:
        t0 = time.time()
        nn = np.asarray(n, float)
        a0, cert, rule, edges = certified_level(nn)
        rows = price(rule, nn)
        l80 = [r["loss_pp"] for r in rows if r["target"] == 0.8]
        rec = {
            "design": name, "n": nn.tolist(), "alpha0": float(a0),
            "level_price_pct": float(100.0 * (F.ABAR - a0) / F.ABAR),
            "certificate": {k: v for k, v in cert.items()
                            if k != "undecided"},
            "certified": bool(cert["certified"]),
            "sup_upper_bound": cert["sup_upper_bound"],
            "edges": [float(x) for x in edges],
            "c": [float(x) for x in rule.c],
            "rows": rows,
            "median_loss_pp_at_80": float(np.median(l80)),
            "max_loss_pp_at_80": float(max(l80)),
            "seconds": float(time.time() - t0),
        }
        out["designs"].append(rec)
        print("%-16s %9.5f %9.5f %11s %8.1f %8.1f %6.0f" %
              (name, a0, cert["sup_upper_bound"] or float("nan"),
               cert["certified"], rec["median_loss_pp_at_80"],
               rec["max_loss_pp_at_80"], rec["seconds"]))
    out["all_certified"] = all(d["certified"] for d in out["designs"])
    out["worst_loss_pp_at_80"] = max(d["max_loss_pp_at_80"]
                                     for d in out["designs"])
    out["falsifiers"] = C.selftest()
    out["note"] = (
        "The deployed object is the step function c(lambda_hat) deposited in "
        "each design's `c`/`edges`.  The certificate bounds its unconditional "
        "FWER above by a monotone corner substitution on a probability-scale "
        "partition of the state space, over a branch-and-bound partition of "
        "the ratio interval; every accepted box carries its own upper bound."
    )
    F.write_deposit("wb_certified_rule.json", out,
                    "unconditional_valid_certified")
    print("\nall certified:", out["all_certified"],
          "| worst loss at 80%%: %.1f pp" % out["worst_loss_pp_at_80"])


if __name__ == "__main__":
    main()
