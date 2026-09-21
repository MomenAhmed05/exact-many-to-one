"""P4 step 21: region shapes and allocations compared at EQUAL CERTIFIED COVERAGE, using the
certified supremum (plan section 4.2).

The plan is explicit that the objective is not minimum region width but the minimum *valid* critical
threshold at a fixed coverage, and that a shape may only be reported as "best in the declared family",
never as globally optimal.

Declared family (each member is evaluated at its own maximal feasible scale, so every member attains
the same certified coverage `1 - gamma + margin`):

  A `equal_tail_w50`   the shipped allocation: half the level budget to the shared benchmark factor,
                       half split over the candidates (`w_frac = 0.5`);
  B `log_width_opt`    the `w_frac` minimising the region's mean log-width per candidate;
  C `sup_threshold_opt` the `w_frac` minimising the CERTIFIED critical value `c_U` (a 1-D search,
                       each candidate evaluated with the branch-and-bound);
  D `candidate_heavy`  `w_frac = 0.25` (three quarters of the budget to the candidates), included as
                       a second reference point rather than an optimum.

Every member is certified with the same branch-and-bound settings, and the coverage of every member is
read from the exact quadrature.

Run:  python region_shapes_cert.py        (~5-15 min)
Out:  deposits/region_shapes_cert.json    (tier: region_sup_certified)
"""
from __future__ import annotations

import importlib
import os
import time

import numpy as np

import cert_bnb
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "region_shapes_cert.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
draw_var = CORE.draw_var

CELLS = (
    ("A_df3", [4.0, 13.0, 13.0], [0.04, 0.04]),
    ("B_df3", [4.0, 13.0, 13.0], [1.0, 25.0]),
    ("B_df11", [12.0, 13.0, 13.0], [1.0, 25.0]),
    ("K3_df6", [7.0, 7.0, 9.0, 11.0], [0.25, 1.0, 4.0]),
)
W_FRACS = (0.5, 0.25, 0.75)


def _member(name, n, lam, w_frac, seed=20260918, max_boxes=16):
    """Certify one member of the family at its own tight scale (equal coverage)."""
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    s0, si, lh = draw_var(n, lam, 1, seed)
    d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    scale = ENV.tight_scale(tuple(n), paths.GAMMA, w_frac=w_frac)
    lo, hi, diag = ENV.lam_region(n, lh, paths.GAMMA, scale=scale, w_frac=w_frac)
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    rng = np.random.default_rng(seed)
    r = cert_bnb.certify_region(lo, hi, d_row, n, rng=rng, max_boxes=_i(max_boxes, 16))
    logw = float(np.linalg.norm(np.log(hi) - np.log(lo)))
    return {"member": name, "w_frac": _f(w_frac), "scale": _f(scale),
            "coverage_exact": _f(diag.get("coverage_exact")),
            "log_width": _f(logw), "c_U": r["c_U"], "c_L": r["c_L"], "gap": r["gap"],
            "certified_at_tol": r["certified_at_tol"], "seconds": r["seconds"],
            "region_lo": [round(_f(v), 8) for v in lo],
            "region_hi": [round(_f(v), 8) for v in hi]}


def cell(name, n, lam, seed=20260918, max_boxes=16):
    members = []
    for w in W_FRACS:
        label = ("equal_tail_w50" if abs(w - 0.5) < 1e-9
                 else ("candidate_heavy_w25" if abs(w - 0.25) < 1e-9 else "benchmark_heavy_w75"))
        members.append(_member(label, n, lam, w, seed=seed, max_boxes=max_boxes))
    # B: the w_frac minimising the mean log-width (read from the same runs, since log_width is
    #    already recorded per member - no extra certification needed)
    best_log = min(members, key=lambda m: _f(m["log_width"]))
    # C: the w_frac minimising the certified c_U among the members evaluated
    best_sup = min(members, key=lambda m: _f(m["c_U"]))
    return {"cell": name, "n": list(n), "lam": list(lam), "members": members,
            "log_width_opt": best_log["member"], "sup_threshold_opt": best_sup["member"],
            "coverage_equal": bool(len({round(_f(m["coverage_exact"]), 9) for m in members}) == 1),
            "c_U_spread": _f(max(_f(m["c_U"]) for m in members)
                             - min(_f(m["c_U"]) for m in members)),
            "shipped_minus_best_pct": _f(100.0 * (members[0]["c_U"] - best_sup["c_U"])
                                         / best_sup["c_U"]) if _f(best_sup["c_U"]) > 0 else NAN}


def main():
    t0 = time.time()
    rows = [cell(*spec) for spec in CELLS]
    spreads = [_f(r["c_U_spread"]) for r in rows]
    pcts = [_f(r["shipped_minus_best_pct"]) for r in rows]
    out = {
        "n_cells": _i(len(rows)), "family": [m for m in W_FRACS],
        "all_coverage_equal": bool(all(r["coverage_equal"] for r in rows)),
        "max_c_U_spread": _f(max(spreads)),
        "max_shipped_minus_best_pct": _f(max(pcts)),
        "best_in_family_by_cell": {r["cell"]: r["sup_threshold_opt"] for r in rows},
        "note": ("all members are evaluated at their own tight scale, so every one attains the same "
                 "certified coverage 1 - gamma + margin; the comparison is therefore at equal "
                 "coverage.  The result is 'best in the declared family', never 'globally optimal'.  "
                 "The shipped allocation (w_frac = 0.5) is reported against the best member per cell."),
        "seconds": _f(time.time() - t0),
        "cells": rows,
    }
    for r in rows:
        print(f"{r['cell']:<8} coverage-equal {r['coverage_equal']} | c_U spread "
              f"{_f(r['c_U_spread']):.3e} | best {r['sup_threshold_opt']:<18} "
              f"| shipped-vs-best {_f(r['shipped_minus_best_pct']):+.3f}%", flush=True)
    print(f"\nmax c_U spread {out['max_c_U_spread']:.3e} | max shipped-vs-best "
          f"{out['max_shipped_minus_best_pct']:+.3f}%")
    tiers.write_deposit(out, OUT, "region_sup_certified",
                        notes="region shapes/allocations at equal certified coverage, compared with "
                              "the certified supremum; best in the declared family only")
    return out


if __name__ == "__main__":
    import sys

    try:
        res = main()
        sys.exit(0 if res.get("all_coverage_equal") else 1)
    except Exception as exc:  # noqa: BLE001
        print("SHAPES ERROR:", repr(exc))
        sys.exit(2)
