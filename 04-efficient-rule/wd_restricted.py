"""Workstream D section 6.2: what an ex ante ratio restriction buys.

The uniform constraint is taken over a declared ratio domain.  The default is
three decades either side, which is deliberately generous.  A desk that runs
volatility targeting, or has contractual risk limits, can often justify a much
narrower domain BEFORE seeing any performance data - and the narrower the
domain, the less the rule has to protect against.

This measures the trade.  For each restriction the internal level is
recalibrated with the supremum taken over that domain only, and the power is
re-priced.  Nothing here is free: the guarantee becomes conditional on the
restriction holding, so the restriction has to be auditable and fixed in advance.
A restriction chosen after seeing returns is worth nothing, and the deposit
records that requirement alongside the numbers.

Deposit: deposits/wd_restricted.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L

DOMAINS = [
    ("unrestricted_3_decades", 1e-3, 1e3),
    ("two_decades", 1e-2, 1e2),
    ("one_decade", 0.1, 10.0),
    ("factor_of_4", 0.25, 4.0),
    ("factor_of_2", 0.5, 2.0),
    ("vol_targeted_1.5", 1.0 / 1.5, 1.5),
]
DESIGNS = [("df3", [4.0, 13.0]), ("df5", [6.0, 13.0]),
           ("df11", [12.0, 13.0]), ("df35", [36.0, 13.0])]
TAB = dict(m=41, npanel=12, order=8)


def run(n, lo, hi, abar=F.ABAR):
    n = np.asarray(n, float)
    coarse = F.log_grid(lo, hi, 21)
    fine = F.log_grid(lo, hi, 81)
    cal = L.calibrate_level(n, coarse, abar, npanel=TAB["npanel"],
                            order=TAB["order"], iters=11)
    tab = L.PluginTable(n, cal["alpha0"], **TAB)
    s, lw, _ = L.sup_fwer(tab, n, fine, 20, 10)
    t = 0
    while s > abar and t < 12:
        cal["alpha0"] *= 0.98
        tab = L.PluginTable(n, cal["alpha0"], **TAB)
        s, lw, _ = L.sup_fwer(tab, n, fine, 20, 10)
        t += 1
    # price only at ratios inside the restriction: outside it the rule makes no
    # claim, so quoting its power there would be meaningless
    alts = [x for x in (0.04, 0.25, 1.0, 4.0, 25.0) if lo <= x <= hi] or [1.0]
    losses = []
    for lam in alts:
        d80 = F.delta_for_oracle_power_k1(n, [lam], 0.8, abar)
        p = F.power_exact_k1(tab, n, lam, d80, 20, 10)
        losses.append(100.0 * (0.8 - p))
    return {"lo": lo, "hi": hi, "alpha0": float(cal["alpha0"]),
            "sup_fwer": float(s), "lam_worst": float(lw),
            "alts_priced": alts,
            "median_loss_pp": float(np.median(losses)),
            "max_loss_pp": float(max(losses)),
            "tightenings": int(t)}


def main():
    out = {"abar": F.ABAR, "domains": [d[0] for d in DOMAINS],
           "designs": [],
           "governance": (
               "A restriction is worth its power only if it is declared before "
               "the performance data are seen and cannot be revised after. "
               "Admissible sources: a published volatility target, a "
               "contractual risk limit, a strategy-construction constraint, or "
               "a pre-evaluation volatility history not reused for selection. "
               "A ratio bound inferred from the same returns being tested is "
               "circular and voids the guarantee.")}
    print("%-8s %-24s %9s %9s %8s %8s" %
          ("design", "domain", "alpha0", "supFWER", "lossMed", "lossMax"))
    for name, n in DESIGNS:
        rec = {"design": name, "n": list(n), "rows": []}
        base = None
        for dname, lo, hi in DOMAINS:
            t0 = time.time()
            r = run(n, lo, hi)
            r["domain"] = dname
            if base is None:
                base = r["median_loss_pp"]
            r["loss_recovered_pp"] = float(base - r["median_loss_pp"])
            r["seconds"] = float(time.time() - t0)
            rec["rows"].append(r)
            print("%-8s %-24s %9.5f %9.5f %8.1f %8.1f" %
                  (name, dname, r["alpha0"], r["sup_fwer"],
                   r["median_loss_pp"], r["max_loss_pp"]))
        rec["best_recovery_pp"] = float(max(x["loss_recovered_pp"]
                                            for x in rec["rows"]))
        out["designs"].append(rec)
    out["note"] = (
        "The level and the loss are both monotone in the width of the declared "
        "domain, which is the whole point: the uniform constraint binds at the "
        "least favourable ratio, and removing far-away ratios from the domain "
        "removes the state the rule was protecting against."
    )
    F.write_deposit("wd_restricted.json", out, "frontier_measured")
    print("\nbest recovery per design:",
          {d["design"]: round(d["best_recovery_pp"], 1) for d in out["designs"]})


if __name__ == "__main__":
    main()
