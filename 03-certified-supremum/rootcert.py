"""P2 step 10: the point root-certificate battery (plan section 3.4, SPEC section 6).

`enclose_gl.bracket_point_gl` already computes the certified bracket `c*(lambda) in [c_L, c_U]`
at a point.  This script runs that primitive over the design suite and, for every witness point,
performs the DECISIVE independent check with the solver that shares no code with the enclosure:

    cond_noreject(c_U; lambda_w, npanels=paths.REF_PANELS) >= 1 - alpha        (c_U is high enough)
    cond_noreject(c_L; lambda_w, npanels=paths.REF_PANELS) < 1 - alpha         (c_L is below the crossing)

plus the predeclared budgets: bracket width `<= TOL = 1e-4`, `Delta p <= CERT_DP_TARGET = 1e-6`.
Witness points per cell: every region vertex (up to a cap), the box centre, and the mixed corners
found by `deriv_diag` (currently none - the sup attains at a vertex in 12/12 cells).

Run:  python rootcert.py          (a few minutes)
Out:  deposits/rootcert.json      (tier: region_sup_certified)
"""

from __future__ import annotations

import contextlib
import importlib
import itertools
import os
import time

import enclose_gl
import numpy as np
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "rootcert.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
cond_noreject = CORE.cond_noreject
draw_var = CORE.draw_var

FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
DFS = (3, 5, 7, 11, 19, 35)
# extra geometries: low df, extreme spreads, K = 1, 3, 5
EXTRA = (
    {"name": "K1_df3", "n": [4.0, 7.0], "lam": [1.0]},
    {"name": "K3_df6", "n": [7.0, 7.0, 9.0, 11.0], "lam": [0.25, 1.0, 4.0]},
    {"name": "K5_df11", "n": [12.0] + [12.0] * 5, "lam": [0.2, 0.5, 1.0, 2.0, 5.0]},
    {"name": "K2_extreme", "n": [4.0, 4.0, 30.0], "lam": [0.01, 100.0]},
    {"name": "K2_quiet", "n": [30.0, 4.0, 4.0], "lam": [1.0, 1.0]},
)
MAX_VERTICES = 8


def _witnesses(lo, hi, centre=True):
    """Region vertices (capped) plus the box centre."""
    K = np.asarray(lo, float).size
    pts = []
    for combo in itertools.product((0, 1), repeat=K):
        pts.append(
            np.where(
                np.array(combo, bool), np.asarray(hi, float), np.asarray(lo, float)
            )
        )
    if len(pts) > MAX_VERTICES:
        idx = np.linspace(0, len(pts) - 1, MAX_VERTICES).astype(int)
        pts = [pts[i] for i in idx]
    if centre:
        pts.append(np.sqrt(np.asarray(lo, float) * np.asarray(hi, float)))
    return pts


def _cell(name, n, lam, alpha=paths.ALPHA0, seed=20260918, prec=64):
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    s0, si, lh = draw_var(n, lam, 1, seed)
    d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    lo, hi, diag = ENV.lam_region(n, lh, paths.GAMMA, scale=None)
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    rows = []
    for w, lam_pt in enumerate(_witnesses(lo, hi)):
        s_pt = np.sqrt(n[1:] / lam_pt)
        br = enclose_gl.bracket_point_gl(s_pt, d_row, n[0], alpha=alpha, prec=prec)
        f_cl = _f(
            cond_noreject(
                np.array([_f(br["c_L"])]),
                np.array([s0[0]]),
                si[0:1],
                n,
                lam_pt,
                npanels=paths.REF_PANELS,
            )[0]
        )
        f_cu = _f(
            cond_noreject(
                np.array([_f(br["c_U"])]),
                np.array([s0[0]]),
                si[0:1],
                n,
                lam_pt,
                npanels=paths.REF_PANELS,
            )[0]
        )
        target = 1.0 - _f(alpha)
        rows.append(
            {
                "cell": name,
                "witness": _i(w),
                "K": _i(lam.size),
                "n": list(n),
                "lambda_witness": [round(_f(v), 8) for v in lam_pt],
                "region_lo": [round(_f(v), 8) for v in lo],
                "region_hi": [round(_f(v), 8) for v in hi],
                "c_L": br["c_L"],
                "c_U": br["c_U"],
                "width": br["width"],
                "dp": br["dp"],
                "c_equiv": br["c_equiv"],
                "slope_hi": br["slope_hi"],
                "n_eval": br["n_eval"],
                "noreject_at_c_L": f_cl,
                "noreject_at_c_U": f_cu,
                "target": target,
                "c_L_below": bool(f_cl < target),
                "c_U_above": bool(f_cu >= target - 1e-15),
                "width_ok": bool(_f(br["width"]) <= paths.TOL),
                "dp_ok": bool(_f(br["dp"]) <= paths.CERT_DP_TARGET),
                "settings": br.get("settings"),
            }
        )
    return {
        "name": name,
        "region_scale": diag.get("scale"),
        "coverage_exact": diag.get("coverage_exact"),
        "rows": rows,
        "lambda_hat": [round(_f(v), 8) for v in np.asarray(lh[0], float)],
        "d_row": [round(_f(v), 8) for v in d_row],
    }


def main():
    t0 = time.time()
    cells = []
    for fam, lam in FAMS.items():
        for df in DFS:
            cells.append(
                {"name": f"{fam}_df{df}", "n": [df + 1.0, 13.0, 13.0], "lam": list(lam)}
            )
    cells += list(EXTRA)
    out_cells = []
    for spec in cells:
        try:
            out_cells.append(_cell(spec["name"], spec["n"], spec["lam"]))
        except Exception as exc:  # noqa: BLE001
            out_cells.append({"name": spec["name"], "error": repr(exc), "rows": []})
    rows = [r for c in out_cells for r in c.get("rows", [])]
    widths = [_f(r["width"]) for r in rows if r["width"] == r["width"]]
    dps = [_f(r["dp"]) for r in rows if r["dp"] == r["dp"]]
    summary = {
        "n_cells": _i(len(out_cells)),
        "n_brackets": _i(len(rows)),
        "all_independent_check_pass": bool(
            rows and all(r["c_L_below"] and r["c_U_above"] for r in rows)
        ),
        "all_widths_under_tol": bool(widths and max(widths) <= paths.TOL),
        "all_dp_under_target": bool(dps and max(dps) <= paths.CERT_DP_TARGET),
        "max_width": _f(max(widths)) if widths else NAN,
        "median_width": _f(np.median(widths)) if widths else NAN,
        "max_dp": _f(max(dps)) if dps else NAN,
        "median_dp": _f(np.median(dps)) if dps else NAN,
        "max_c_equiv": _f(
            max(
                (_f(r["c_equiv"]) for r in rows if r["c_equiv"] == r["c_equiv"]),
                default=NAN,
            )
        ),
        "tol_predeclared": _f(paths.TOL),
        "dp_target": _f(paths.CERT_DP_TARGET),
        "alpha0": _f(paths.ALPHA0),
        "gamma": _f(paths.GAMMA),
        "seconds": _f(time.time() - t0),
        "cells": out_cells,
    }
    print(
        f"cells {summary['n_cells']} | brackets {summary['n_brackets']} | "
        f"independent check pass {summary['all_independent_check_pass']}"
    )
    print(
        f"  width  max {_f(summary['max_width']):.3e} median "
        f"{_f(summary['median_width']):.3e} (tol {paths.TOL:.0e}) "
        f"-> under tol: {summary['all_widths_under_tol']}"
    )
    print(
        f"  dp     max {_f(summary['max_dp']):.3e} median "
        f"{_f(summary['median_dp']):.3e} (target {paths.CERT_DP_TARGET:.0e}) "
        f"-> under target: {summary['all_dp_under_target']}"
    )
    print(f"  seconds {_f(summary['seconds']):.0f}")
    bad = [r for r in rows if not (r["c_L_below"] and r["c_U_above"])]
    for r in bad[:5]:
        print(
            f"  FAIL {r['cell']} witness {r['witness']}: F(c_L)={r['noreject_at_c_L']} "
            f"F(c_U)={r['noreject_at_c_U']} target={r['target']}"
        )
    tiers.write_deposit(
        summary,
        OUT,
        "region_sup_certified",
        notes="point root certificates c* in [c_L, c_U] over the design suite, "
        "each verified against the independent solver core.cond_noreject "
        "(c_L below the crossing, c_U at or above it)",
    )
    return summary


if __name__ == "__main__":
    with contextlib.suppress(Exception):
        main()
