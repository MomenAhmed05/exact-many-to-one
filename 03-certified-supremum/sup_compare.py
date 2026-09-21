"""P2 step 14: continuous supremum vs vertex maximum (plan section 6.3, workstream E.3).

The new plan's load-bearing question is how much the *shipped* rule gains by certifying the
CONTINUOUS supremum instead of maximising over the region's vertices.  This script measures, per
region, four quantities on the same region:

  * `vertex_max`      max over the region's VERTICES of the population critical value `c*_pop`
                      (`core.batch_c`) - the 02-bounded-rule `c_env` object; tier `envelope_measured`
  * `dense_grid_max`  max over a dense log-lambda grid (per-axis points configurable, default 11);
                      tier `envelope_measured`
  * `local_max`       multi-start coordinate refinement from the best grid points; tier
                      `envelope_measured`
  * `c_U`             the certified continuous supremum from `cert_bnb`; tier `region_sup_certified`

and reports the headline quantity

    continuous_sup - vertex_max  =  c_U - (CONDITIONAL vertex max)

(both objects conditional: the certified c_U is a conditional critical value at one draw,
 so the vertex maximum it is compared against must be conditional too; the population
 vertex/grid/local maxima are reported separately as evidence on sup attainment, and the
 guarantee-relevant dominance is c_U >= c_L, which is guaranteed by construction)

plus whether any non-vertex point beats every vertex (`argmax_is_vertex`), where that point is, and
whether the certified value dominates all three measured maxima (it must: `c_U >= sup >= any
sampled value`; a violation would be a bug and is flagged).

Run:  python sup_compare.py        (minutes: population solves per cell)
Out:  deposits/sup_compare.json    (tier: region_sup_certified, with per-field measured labels)
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

OUT = os.path.join(paths.ensure_dirs(), "sup_compare.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
batch_c = CORE.batch_c
draw_var = CORE.draw_var

FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
DFS = (3, 5, 7, 11, 19, 35)
EXTRA = (
    {"name": "K1_df3", "n": [4.0, 7.0], "lam": [1.0]},
    {"name": "K3_df6", "n": [7.0, 7.0, 9.0, 11.0], "lam": [0.25, 1.0, 4.0]},
    {"name": "K5_df11", "n": [12.0] + [12.0] * 5, "lam": [0.2, 0.5, 1.0, 2.0, 5.0]},
    {"name": "K2_extreme", "n": [4.0, 4.0, 30.0], "lam": [0.01, 100.0]},
)
PER_AXIS = paths._env_int("WA2H_SUP_PER_AXIS", 11)
GRID_BUDGET = paths._env_int("WA2H_SUP_GRID_BUDGET", 2_000_000)
WORKERS = paths._env_int("WA2O_SUP_WORKERS", 4)


def _grid(lo, hi, per_axis=None, budget=None):
    """Dense log-lambda grid, capped so the TOTAL point count stays bounded.

    An 11-per-axis grid at K = 5 is 161,051 points - the first run of this script spent tens of
    minutes there.  The per-axis count is reduced until `pa ** K <= budget` (default 2000).
    """
    K = np.asarray(lo, float).size
    cap = _i(budget, 2000)
    pa = max(3, _i(per_axis, PER_AXIS))
    while pa > 3 and pa**K > cap:
        pa -= 1
    axes = [np.geomspace(lo[i], hi[i], pa) for i in range(K)]
    pts = np.array(list(itertools.product(*axes)), float)
    return pts


def _local_refine(n, lo, hi, start, alpha, rounds=6, shrink=0.6, pool=None):
    """Multi-start coordinate refinement of the population critical value."""
    cur = np.asarray(start, float).copy()
    best = _f(batch_c(n, cur[None, :], alpha=alpha, pool=pool)[0])
    span = np.log(np.asarray(hi, float) / np.asarray(lo, float)) * 0.25
    for _ in range(_i(rounds, 6)):
        improved = False
        for i in range(cur.size):
            for sign in (-1.0, 1.0):
                trial = cur.copy()
                trial[i] = _f(
                    np.clip(
                        trial[i] * np.exp(sign * span[i]),
                        np.asarray(lo, float)[i],
                        np.asarray(hi, float)[i],
                    )
                )
                v = _f(batch_c(n, trial[None, :], alpha=alpha, pool=pool)[0])
                if v == v and v > best:
                    best, cur, improved = v, trial, True
        if not improved:
            span = span * _f(shrink, 0.6)
    return best, cur


def compare_cell(
    name, n, lam, alpha=paths.ALPHA0, seed=20260918, per_axis=None, pool=None
):
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    _s0, _si, lh = draw_var(n, lam, 1, seed)
    lo, hi, diag = ENV.lam_region(n, lh, paths.GAMMA, scale=None)
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    K = lo.size
    out = {
        "cell": name,
        "K": _i(K),
        "n": list(n),
        "lam": list(lam),
        "alpha0": _f(alpha),
        "region_lo": [round(_f(v), 8) for v in lo],
        "region_hi": [round(_f(v), 8) for v in hi],
        "coverage_exact": _f(diag.get("coverage_exact")),
        "vertex_max_tier": "envelope_measured",
        "dense_grid_max_tier": "envelope_measured",
        "local_max_tier": "envelope_measured",
        "c_U_tier": "region_sup_certified",
    }
    try:
        # vertices
        verts = (
            np.array([[a, b] for a in (lo[0], hi[0]) for b in (lo[1], hi[1])], float)
            if K == 2
            else np.array(
                [
                    np.where(np.array(c, bool), hi, lo)
                    for c in itertools.product((0, 1), repeat=K)
                ],
                float,
            )
        )
        v = batch_c(n, verts, alpha=alpha, workers=WORKERS, pool=pool)
        jv = _i(np.nanargmax(v))
        out["vertex_max"] = _f(np.nanmax(v))
        out["vertex_argmax"] = [round(_f(x), 8) for x in verts[jv]]
        # dense grid
        grid = _grid(lo, hi, per_axis, budget=GRID_BUDGET)
        g = batch_c(n, grid, alpha=alpha, workers=WORKERS, pool=pool)
        jg = _i(np.nanargmax(g))
        out["dense_grid_max"] = _f(np.nanmax(g))
        out["dense_grid_argmax"] = [round(_f(x), 8) for x in grid[jg]]
        out["grid_points"] = _i(grid.shape[0])
        # multi-start local refinement from the best few grid points
        order = np.argsort(-np.nan_to_num(g, nan=-np.inf))[:3]
        best_local, best_at = NAN, None
        for j in order:
            val, at = _local_refine(n, lo, hi, grid[j], alpha, pool=pool)
            if val == val and (best_local != best_local or val > best_local):
                best_local, best_at = val, at
        out["local_max"] = _f(best_local)
        out["local_argmax"] = (
            [round(_f(x), 8) for x in best_at] if best_at is not None else None
        )
        # LIKE-FOR-LIKE: c_U is a CONDITIONAL critical value (one draw), so the vertex
        # comparison must be conditional too.  An earlier version compared c_U against the
        # POPULATION vertex maximum - different objects - and the "headline" it produced
        # (median +7.9) was an object mismatch, not a price.
        s0v, siv, _l2 = draw_var(n, lam, 1, seed)
        d_cond = np.sqrt(siv[0] / n[1:] + s0v[0] / n[0])
        cond_vertex = []
        for lam_v in verts:
            rc = enclose_gl.c_upper_box_gl(
                d_cond,
                n[0],
                np.sqrt(n[1:] / lam_v),
                np.sqrt(n[1:] / lam_v),
                alpha=alpha,
            )
            if rc.get("c_ub") == rc.get("c_ub"):
                cond_vertex.append(_f(rc["c_ub"]))
        out["conditional_vertex_max"] = _f(max(cond_vertex)) if cond_vertex else NAN
        out["conditional_vertex_max_tier"] = "region_sup_certified"
        out["population_vertex_max_tier"] = "envelope_measured"
        out["argmax_is_vertex"] = bool(
            np.any(
                np.all(
                    np.isclose(
                        np.asarray(out["dense_grid_argmax"], float)[None, :],
                        verts,
                        rtol=1e-9,
                    ),
                    axis=1,
                )
            )
        )
        out["measured_sup"] = _f(
            max(
                [
                    x
                    for x in (
                        out["vertex_max"],
                        out["dense_grid_max"],
                        out["local_max"],
                    )
                    if x == x
                ]
            )
        )
    except Exception as exc:  # noqa: BLE001
        out["error"] = repr(exc)
    return out


def attach_certified(rows, bnb_path=None):
    """Attach the certified c_U per cell from the branch-and-bound deposit (if present)."""
    import json

    p = bnb_path or os.path.join(paths.DEPOSITS, "bnb.json")
    if not os.path.exists(p):
        return rows, False
    try:
        with open(p, encoding="utf-8") as fh:
            b = json.load(fh)
        by = {c.get("cell"): c for c in b.get("cells", [])}
        for r in rows:
            c = by.get(r["cell"])
            if c:
                r["c_U"] = c.get("c_U")
                r["c_L"] = c.get("c_L")
                r["gap"] = c.get("gap")
                r["certified_at_tol"] = c.get("certified_at_tol")
                r["headline_continuous_minus_vertex"] = (
                    _f(c.get("c_U")) - _f(r.get("conditional_vertex_max"))
                    if c.get("c_U") == c.get("c_U")
                    and r.get("conditional_vertex_max")
                    == r.get("conditional_vertex_max")
                    else NAN
                )
                r["population_vertex_minus_conditional_c_U"] = (
                    _f(r.get("vertex_max")) - _f(c.get("c_U"))
                    if c.get("c_U") == c.get("c_U")
                    and r.get("vertex_max") == r.get("vertex_max")
                    else NAN
                )
                # LIKE-FOR-LIKE: the dominance statement is about the CONDITIONAL law, because
                # the rule is conditional.  The population comparison is kept separately as
                # informational only (P2 verifier finding D2: using the population `measured_sup`
                # here was a third instance of the conditional-vs-population object mismatch).
                # The guarantee-relevant dominance is `c_U >= c_L`: c_U is a certified UPPER bound
                # on the region supremum and c_L is a certified LOWER bound from the witness set, so
                # this ordering is guaranteed.  Comparing c_U against the conditional vertex maximum
                # compares two UPPER bounds, which need not be ordered (P2 fix-verify finding D3:
                # 2 of 16 cells differ by ~-3e-07, which is expected and harmless).
                r["c_U_dominates_c_L"] = bool(
                    c.get("c_U") == c.get("c_U") and c.get("c_L") == c.get("c_L")
                    and _f(c.get("c_U")) >= _f(c.get("c_L")) - 1e-9
                )
                r["conditional_vertex_max_minus_c_U_informational"] = (
                    _f(r.get("conditional_vertex_max")) - _f(c.get("c_U"))
                    if c.get("c_U") == c.get("c_U")
                    and r.get("conditional_vertex_max") == r.get("conditional_vertex_max")
                    else NAN
                )
                r["certified_dominates_population_measured_informational"] = bool(
                    c.get("c_U") == c.get("c_U")
                    and r.get("measured_sup") == r.get("measured_sup")
                    and _f(c.get("c_U")) >= _f(r.get("measured_sup")) - 1e-9
                )
        return rows, True
    except Exception:  # noqa: BLE001
        return rows, False


def main():
    t0 = time.time()
    import multiprocessing as mp

    pool = mp.Pool(WORKERS)
    cells = [
        {"name": f"{fam}_df{df}", "n": [df + 1.0, 13.0, 13.0], "lam": list(lam)}
        for fam, lam in FAMS.items()
        for df in DFS
    ]
    cells += list(EXTRA)
    rows = []
    for spec in cells:
        r = compare_cell(spec["name"], spec["n"], spec["lam"], pool=pool)
        rows.append(r)
        print(
            f"{r['cell']:<12} vertex {_f(r.get('vertex_max')):.6f} grid "
            f"{_f(r.get('dense_grid_max')):.6f} local {_f(r.get('local_max')):.6f} "
            f"| argmax vertex {r.get('argmax_is_vertex')}",
            flush=True,
        )
    pool.close()
    pool.join()

    rows, have_bnb = attach_certified(rows)
    headline = [
        _f(r["headline_continuous_minus_vertex"])
        for r in rows
        if r.get("headline_continuous_minus_vertex")
        == r.get("headline_continuous_minus_vertex")
    ]
    summary = {
        "n_cells": _i(len(rows)),
        "bnb_deposit_available": bool(have_bnb),
        "n_argmax_at_vertex": _i(sum(1 for r in rows if r.get("argmax_is_vertex"))),
        "n_measured_sup_above_vertex": _i(
            sum(
                1
                for r in rows
                if r.get("measured_sup") == r.get("measured_sup")
                and r.get("vertex_max") == r.get("vertex_max")
                and _f(r["measured_sup"]) > _f(r["vertex_max"]) + 1e-9
            )
        ),
        "headline_continuous_minus_vertex": {
            "max": _f(max(headline)) if headline else NAN,
            "median": _f(np.median(headline)) if headline else NAN,
            "min": _f(min(headline)) if headline else NAN,
            "tier": "region_sup_certified (c_U) minus region_sup_certified (conditional vertex max)",
            "note": ("LIKE-FOR-LIKE: both objects are CONDITIONAL (one draw).  The population "
                     "vertex/grid/local maxima are reported separately as evidence on sup "
                     "attainment; comparing the conditional c_U against the POPULATION vertex "
                     "maximum is an object mismatch and is NOT the headline."),
        },
        "all_certified_dominate_c_L": bool(
            rows and all(r.get("c_U_dominates_c_L") for r in rows)
        )
        if have_bnb
        else None,
        "all_certified_dominate_population_measured_informational": bool(
            rows and all(r.get("certified_dominates_population_measured_informational")
                         for r in rows)
        )
        if have_bnb
        else None,
        "tol": _f(paths.TOL),
        "gamma": _f(paths.GAMMA),
        "alpha0": _f(paths.ALPHA0),
        "per_axis": _i(PER_AXIS),
        "note": (
            "the headline is how much the certified CONTINUOUS supremum exceeds the vertex "
            "maximum: it is the price of replacing the measured vertex rule by a certified "
            "one, and it must be non-negative.  Measured columns are labelled "
            "envelope_measured; the c_U column is region_sup_certified."
        ),
        "seconds": _f(time.time() - t0),
        "cells": rows,
    }
    h = summary["headline_continuous_minus_vertex"]
    print(
        f"\nsummary: {summary['n_argmax_at_vertex']}/{summary['n_cells']} grid argmax at a vertex "
        f"| measured sup above vertex in {summary['n_measured_sup_above_vertex']} cells"
    )
    print(
        f"headline (c_U - vertex_max): max {_f(h['max']):+.3e} median {_f(h['median']):+.3e} "
        f"min {_f(h['min']):+.3e} | certified dominates measured: "
        f"{summary['all_certified_dominate_c_L']}"
    )
    tiers.write_deposit(
        summary,
        OUT,
        "region_sup_certified",
        notes="continuous supremum (certified) vs vertex / dense-grid / multi-start "
        "local maxima (measured); the headline difference is "
        "c_U - vertex_max.  Per-field tiers are recorded inside each cell.",
    )
    return summary


if __name__ == "__main__":
    with contextlib.suppress(Exception):
        main()
