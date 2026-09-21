"""P3 step 16: the per-certificate error budget (SPEC section 9, plan section 3.2).

For each certificate point the enclosure width `Delta p = F_hi - F_lo` is decomposed into its
actual components, each measured rather than asserted:

  * **analyticity remainder** - the summed Bernstein-ellipse panel bounds (`diag["err_total"]`);
  * **node/weight radius** - the contribution of the float64 node/weight perturbation, measured by
    re-running with `NODE_RAD = WEIGHT_RAD = 0` and differencing the widths (this replaced the old,
    unfounded `2*2*eps*M` term: the perturbation is now carried by the interval arithmetic itself);
  * **tail mass** - the explicit tail terms the enclosure adds (`Phi(-z_max)` and `1 - Phi(z_max)`);
  * **Arb ball radii** - measured by re-running at a lower precision (`prec = 32` vs `64`) and
    differencing;
  * **summation slack** - the explicit `GL_SUM_SLACK` added to both endpoints;
  * **bisection width** - the search's own bracket width (`c_hi - c_lo` at the returned c);
  * **Delta p and its c-equivalent** - `F_hi - F_lo` at the certificate point and `Delta p / slope`
    with a local slope bound on `dF/dc`.

The budget is reported per point and summarised; the predeclared targets are `Delta p <= 1e-6` and
`c_U - c_L <= 1e-4`.

Run:  python error_budget.py        (~1 min)
Out:  deposits/error_budget.json    (tier: region_sup_certified)
"""

from __future__ import annotations

import importlib
import os
import time

import enclose_gl
import numpy as np
import paths
import tiers
from paths import _f, _i

OUT = os.path.join(paths.ensure_dirs(), "error_budget.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
draw_var = CORE.draw_var

POINTS = (
    ("desk_B_df35", [36.0, 13.0, 13.0], [1.0, 25.0], 2.6),
    ("desk_A_df3", [4.0, 13.0, 13.0], [0.04, 0.04], 2.4),
    ("sharp_A_df3_lo", [4.0, 13.0, 13.0], [0.04, 0.04], 5.0),
    ("K5_balanced", [13.0] * 6, [1.0, 0.5, 0.25, 2.0, 4.0], 3.0),
)


def _width_at(c, d_row, n0, s_pt, prec=64, node_rad=None, weight_rad=None):
    """Enclosure width at c, optionally overriding the node/weight radii."""
    keep = (enclose_gl.NODE_RAD, enclose_gl.WEIGHT_RAD)
    try:
        if node_rad is not None:
            enclose_gl.NODE_RAD = _f(node_rad, 0.0)
        if weight_rad is not None:
            enclose_gl.WEIGHT_RAD = _f(weight_rad, 0.0)
        lo, hi, diag = enclose_gl.bounds_gl_arb(
            c, d_row, n0, s_pt, s_pt, prec=prec, return_diag=True
        )
        return _f(hi) - _f(lo), diag, lo, hi
    finally:
        enclose_gl.NODE_RAD, enclose_gl.WEIGHT_RAD = keep


def point(name, n, lam, c, seed=20260918, prec=64):
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    s0, si, _lh = draw_var(n, lam, 1, seed)
    d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    s_pt = np.sqrt(n[1:] / lam)
    w_full, diag, lo, hi = _width_at(c, d_row, n[0], s_pt, prec=prec)
    w_nor, _, _, _ = _width_at(
        c, d_row, n[0], s_pt, prec=prec, node_rad=0.0, weight_rad=0.0
    )
    w_low, _, _, _ = _width_at(c, d_row, n[0], s_pt, prec=32)
    br = enclose_gl.bracket_point_gl(s_pt, d_row, n[0], prec=prec)
    tails = diag.get("tail_terms", {}) or {}
    return {
        "point": name,
        "n": list(n),
        "lam": list(lam),
        "c": _f(c),
        "K": _i(lam.size),
        "Delta_p": _f(w_full),
        "analyticity_remainder": _f(diag.get("err_total")),
        "node_weight_radius": _f(w_full - w_nor),
        # upper bound on the Arb ball-radius contribution: the width GROWS when precision is lowered,
        # so (w_low - w_full) bounds what the radii contribute at prec (measured by differencing).
        "arb_ball_radii_upper": _f(w_low - w_full),
        "summation_slack": _f(2.0 * enclose_gl.GL_SUM_SLACK),
        "tail_mass_low": _f(tails.get("low_mass")),
        "tail_mass_up": _f(tails.get("up_mass")),
        "n_panels": _i(diag.get("n_panels")),
        "bisection_width": br.get("bracket_width")
        if br.get("bracket_width") is not None
        else br.get("width"),
        "slope_hi": br.get("slope_hi"),
        "c_equivalent": br.get("c_equiv"),
        "c_L": br.get("c_L"),
        "c_U": br.get("c_U"),
        "settings": diag.get("settings"),
    }


def main():
    t0 = time.time()
    rows = [point(*spec) for spec in POINTS]
    dps = [_f(r["Delta_p"]) for r in rows]
    comps = {
        k: _f(max(_f(r[k]) for r in rows))
        for k in (
            "analyticity_remainder",
            "node_weight_radius",
            "arb_ball_radii_upper",
            "summation_slack",
            "tail_mass_low",
            "tail_mass_up",
        )
    }
    out = {
        "n_points": _i(len(rows)),
        "max_Delta_p": _f(max(dps)),
        "median_Delta_p": _f(np.median(dps)),
        "components_max": comps,
        "dp_target": _f(paths.CERT_DP_TARGET),
        "all_dp_under_target": bool(all(d <= paths.CERT_DP_TARGET for d in dps)),
        "tol": _f(paths.TOL),
        "note": (
            "Delta p is decomposed into components that are MEASURED (by differencing runs with "
            "the radius/precision overridden) rather than asserted.  The node/weight radius "
            "replaced the earlier unfounded 2*2*eps*M term (pass-2 defect #20); it is now "
            "carried inside the interval arithmetic."
        ),
        "seconds": _f(time.time() - t0),
        "rows": rows,
    }
    for r in rows:
        print(
            f"{r['point']:<16} dp {_f(r['Delta_p']):.3e} | analyticity "
            f"{_f(r['analyticity_remainder']):.3e} | radius {_f(r['node_weight_radius']):.3e} "
            f"| arb<= {_f(r['arb_ball_radii_upper']):.3e} | panels {r['n_panels']} "
            f"| c-equiv {_f(r['c_equivalent']):.3e}"
        )
    print(
        f"\nmax Delta p {out['max_Delta_p']:.3e} (target {paths.CERT_DP_TARGET:.0e}) "
        f"-> under target: {out['all_dp_under_target']}"
    )
    tiers.write_deposit(
        out,
        OUT,
        "region_sup_certified",
        notes="per-certificate error budget with measured components",
    )
    return out


if __name__ == "__main__":
    import sys

    try:
        res = main()
        sys.exit(0 if res.get("all_dp_under_target") else 1)
    except Exception as exc:  # noqa: BLE001
        print("BUDGET ERROR:", repr(exc))
        sys.exit(2)
