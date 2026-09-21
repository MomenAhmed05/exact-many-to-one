"""Workstream B at K = 1: a uniformly valid rule that is not the region rule.

THE CONSTRUCTION.  The certified region rule protects, on every draw, against
every ratio in an admission region.  That is far more than the guarantee asks
for.  The guarantee is UNCONDITIONAL:

    sup_lambda  E_S[ P_lambda( max_i T_i > c(S) | S ) ]  <=  abar .

So take the plug-in threshold the paper already studies - the exact fixed-ratio
constant evaluated at the estimated ratios - and run it at an internal level
alpha0 chosen so that the exact uniform constraint above holds.  One scalar, no
region, no supremum over nuisance values per draw.

At K = 1 the constraint is an exact 2-D quadrature and the known-lambda oracle
power is closed form, so the whole calibration and its price are deterministic.

Deposit: deposits/wb_k1.json  (tier unconditional_candidate_measured; the
certification of the uniform supremum is wb_certify.py, and only that promotes
the tier).
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L

# the declared nuisance domain: five decades of ratio, which covers every desk
# and stress design in the programme (02-bounded-rule's widest is a factor of 25)
LAM_DOMAIN = (1e-3, 1e3)
LAM_COARSE = 25       # search grid
LAM_FINE = 97         # verification grid
TAB = dict(m=41, npanel=12, order=8)
VER = dict(npanel=24, order=12)

DESIGNS = [
    ("df3", [4.0, 13.0]),
    ("df5", [6.0, 13.0]),
    ("df7", [8.0, 13.0]),
    ("df11", [12.0, 13.0]),
    ("df19", [20.0, 13.0]),
    ("df35", [36.0, 13.0]),
    ("df59", [60.0, 13.0]),
    ("df35_short_cand", [36.0, 5.0]),
    ("df11_long_cand", [12.0, 61.0]),
    ("balanced_df11", [12.0, 12.0]),
]
ALT_LAMS = (0.04, 0.25, 1.0, 4.0, 25.0)
TARGETS = (0.5, 0.8, 0.95)


def run_design(name, n, abar=F.ABAR):
    n = np.asarray(n, float)
    t0 = time.time()
    coarse = F.log_grid(*LAM_DOMAIN, LAM_COARSE)
    fine = F.log_grid(*LAM_DOMAIN, LAM_FINE)

    # the uncalibrated plug-in, for the price comparison
    tab0 = L.PluginTable(n, abar, **TAB)
    s0, lw0, _ = L.sup_fwer(tab0, n, fine, **VER)

    cal = L.calibrate_level(n, coarse, abar, npanel=TAB["npanel"],
                            order=TAB["order"], iters=12)
    # re-verify the accepted level on the FINE grid at the fine quadrature: the
    # search may sit on a grid artefact, and a rule that fails here is not valid
    tab = L.PluginTable(n, cal["alpha0"], **TAB)
    s1, lw1, vals = L.sup_fwer(tab, n, fine, **VER)
    tighten = 0
    while s1 > abar and tighten < 14:        # shrink until the fine check holds
        cal["alpha0"] *= 0.97
        tab = L.PluginTable(n, cal["alpha0"], **TAB)
        s1, lw1, vals = L.sup_fwer(tab, n, fine, **VER)
        tighten += 1

    rows = []
    for lam in ALT_LAMS:
        for t in TARGETS:
            dl = F.delta_for_oracle_power_k1(n, [lam], t, abar)
            p = F.power_exact_k1(tab, n, lam, dl, **VER)
            p0 = F.power_exact_k1(tab0, n, lam, dl, **VER)
            rows.append({
                "lam": float(lam), "target": float(t), "delta": float(dl),
                "oracle_power": float(t),
                "calibrated_power": float(p),
                "uncalibrated_plugin_power": float(p0),
                "loss_pp": float(100.0 * (t - p)),
                "loss_pp_uncalibrated": float(100.0 * (t - p0)),
            })
    losses = [r["loss_pp"] for r in rows]
    l80 = [r["loss_pp"] for r in rows if r["target"] == 0.8]
    return {
        "design": name, "n": n.tolist(),
        "abar": abar,
        "plugin_at_abar": {"sup_fwer": s0, "lam_worst": lw0,
                           "rel_excess_pct": 100.0 * (s0 - abar) / abar},
        "alpha0": float(cal["alpha0"]),
        "level_price_pct": float(100.0 * (abar - cal["alpha0"]) / abar),
        "sup_fwer_fine": float(s1), "lam_worst_fine": float(lw1),
        "fine_grid_tightenings": int(tighten),
        "valid_on_fine_grid": bool(s1 <= abar),
        "rows": rows,
        "max_loss_pp": float(max(losses)),
        "median_loss_pp": float(np.median(losses)),
        "max_loss_pp_at_80": float(max(l80)),
        "median_loss_pp_at_80": float(np.median(l80)),
        "seconds": float(time.time() - t0),
    }


def main():
    out = {"lam_domain": list(LAM_DOMAIN), "lam_coarse": LAM_COARSE,
           "lam_fine": LAM_FINE, "table": TAB, "verify": VER,
           "alt_lams": list(ALT_LAMS), "targets": list(TARGETS),
           "abar": F.ABAR, "designs": []}
    print("%-16s %8s %9s %9s %9s %9s" %
          ("design", "supFWER0", "alpha0", "supFWER", "loss80med", "loss80max"))
    for name, n in DESIGNS:
        r = run_design(name, n)
        out["designs"].append(r)
        print("%-16s %8.5f %9.5f %9.5f %9.1f %9.1f   [%.0fs]" %
              (name, r["plugin_at_abar"]["sup_fwer"], r["alpha0"],
               r["sup_fwer_fine"], r["median_loss_pp_at_80"],
               r["max_loss_pp_at_80"], r["seconds"]))
    out["all_valid_on_fine_grid"] = all(d["valid_on_fine_grid"]
                                        for d in out["designs"])
    out["worst_loss_pp_at_80"] = max(d["max_loss_pp_at_80"]
                                     for d in out["designs"])
    out["note"] = (
        "alpha0 is the internal level of the plug-in exact-law rule; the "
        "constraint sup_lambda FWER <= abar is evaluated exactly by quadrature "
        "on the fine grid over the declared five-decade ratio domain.  The "
        "tier is measured until wb_certify.py encloses the supremum over the "
        "CONTINUUM rather than a grid."
    )
    F.write_deposit("wb_k1.json", out, "unconditional_candidate_measured")
    print("\nall valid on fine grid:", out["all_valid_on_fine_grid"])
    print("worst power loss at the 80%% point: %.1f pp"
          % out["worst_loss_pp_at_80"])


if __name__ == "__main__":
    main()
