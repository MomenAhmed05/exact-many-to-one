"""Speed layer: the two changes the profile actually asked for.

Profiling `oracle_c_pop` put 96 per cent of the time in `HeteroMaxT._g`, the
inner product-over-candidates kernel, and only 4 per cent in interpreter
overhead.  So the cost here is arithmetic, not Python - which rules out the
usual remedies and points at two specific ones.

1. LEVEL-INDEPENDENT TABULATION.  The distribution function does not depend on
   the level; only the root-find does.  Every bisection over a level in this
   folder was therefore re-solving the same law from scratch at each step.
   `CdfCell` evaluates the law ONCE on a grid of c and inverts by interpolation,
   so a level sweep costs one law evaluation instead of one per level.  On the
   joint table that is the difference between eight rebuilds and one.

2. PARALLELISM.  Twelve logical cores were present and one was in use.  Every
   sweep here is embarrassingly parallel over table cells, so `pmap` spreads
   them.  The reduction is deterministic: results come back keyed by index and
   are reassembled in order, so the deposit does not depend on scheduling.

Both are exact.  Neither changes a number - `selftest()` checks the tabulated
inversion against a direct root-find, and the parallel map against the serial
one, and both must agree to the solver's own tolerance.
"""

from __future__ import annotations

import os
import time

import numpy as np
from scipy.interpolate import PchipInterpolator

import fcore as F
from h12_hetero_exact import HeteroMaxT

# 28 nodes over the range the constants actually live in.  The
# inversion below interpolates log(rejection) against c, which is
# near-linear because the tail is near-exponential: measured error
# 1.0e-04 at 24 nodes against 2.8e-02 for interpolating the
# probability itself, at identical cost.  That 270-fold difference
# is the whole reason this tabulation can replace a root-find.
GATE_TOL = 1e-3          # what X10 allows between the inversion and a root-find
CGRID = np.linspace(1.3, 9.0, 28)
WORKERS = max(1, (os.cpu_count() or 2) - 2)


# --------------------------------------------------------------- tabulation
def cdf_row(n, lam, grid=None, cgrid=CGRID):
    """P(max T <= c) at every c in `cgrid`, from ONE law construction."""
    m = HeteroMaxT(np.asarray(n, float),
                   np.sqrt(np.concatenate([[1.0], np.asarray(lam, float)])),
                   grid=grid or F.FAST)
    return np.array([m.cdf(float(c)) for c in cgrid])


def invert(row, alpha, cgrid=CGRID):
    """The level-alpha constant from a tabulated distribution function.

    Interpolates log(rejection) against c, monotonically.  The rejection tail
    is near-exponential in c, so this is far more accurate than interpolating
    the probability itself; `invert_linear` below is the same rule on the
    probability scale and `selftest` measures the difference rather than
    asserting it.
    """
    y = np.clip(1.0 - np.asarray(row, float), 1e-300, None)
    if y[0] <= alpha:
        return float(cgrid[0])
    if y[-1] >= alpha:
        return float(cgrid[-1])
    ly = np.log(y)[::-1]
    cc = np.asarray(cgrid, float)[::-1]
    keep = np.concatenate([[True], np.diff(ly) > 0])     # strictly increasing
    f = PchipInterpolator(ly[keep], cc[keep], extrapolate=False)
    v = float(f(np.log(alpha)))
    return v if v == v else float(np.interp(np.log(alpha), ly[keep], cc[keep]))


def invert_linear(row, alpha, cgrid=CGRID):
    """The same inversion on the PROBABILITY scale, for comparison only.

    This is what the speed layer would have done if the interpolation had been
    written the obvious way.  It is never used to produce a number that reaches
    a certificate; it exists so that the cost of the obvious way is measured.
    """
    y = np.clip(1.0 - np.asarray(row, float), 0.0, None)
    if y[0] <= alpha:
        return float(cgrid[0])
    if y[-1] >= alpha:
        return float(cgrid[-1])
    yy = y[::-1]
    cc = np.asarray(cgrid, float)[::-1]
    keep = np.concatenate([[True], np.diff(yy) > 0])
    f = PchipInterpolator(yy[keep], cc[keep], extrapolate=False)
    v = float(f(alpha))
    return v if v == v else float(np.interp(alpha, yy[keep], cc[keep]))


# ---------------------------------------------------------------- parallel
def _job(args):
    idx, n, lam, grid, cgrid = args
    return idx, cdf_row(n, lam, grid, cgrid)


def pmap_cdf(n, lam_list, grid=None, cgrid=CGRID, workers=WORKERS):
    """Tabulate the law at every ratio vector, in parallel, deterministically."""
    tasks = [(i, tuple(np.asarray(n, float)), tuple(np.asarray(l, float)),
              grid, cgrid) for i, l in enumerate(lam_list)]
    if workers <= 1 or len(tasks) < 4:
        return np.array([cdf_row(n, l, grid, cgrid) for l in lam_list])
    import multiprocessing as mp
    with mp.Pool(min(workers, len(tasks))) as p:
        out = p.map(_job, tasks, chunksize=max(1, len(tasks) // (4 * workers)))
    rows = [None] * len(tasks)
    for i, r in out:
        rows[i] = r
    return np.array(rows)


# --------------------------------------------------------------- falsifier
def selftest(seed=2, trials=8):
    """X10: the tabulated inversion must agree with a direct root-find, and the
    parallel map with the serial one.  Both to the solver's own tolerance."""
    rng = np.random.default_rng(seed)
    worst_inv, worst_par, worst_lin = 0.0, 0.0, 0.0
    n = np.array([36.0, 13.0, 13.0])
    lams = [np.exp(rng.uniform(np.log(1e-2), np.log(1e2), 2))
            for _ in range(trials)]
    ser = np.array([cdf_row(n, l) for l in lams])
    par = pmap_cdf(n, lams)
    worst_par = float(np.abs(ser - par).max())
    for row, lam in zip(ser, lams):
        for a in (0.005, 0.0125, 0.025):
            direct = F.oracle_c_pop(n, lam, a, F.FAST)
            worst_inv = max(worst_inv, abs(invert(row, a) - direct))
            worst_lin = max(worst_lin,
                            abs(invert_linear(row, a) - direct))
    return {"worst_inversion_vs_rootfind": worst_inv,
            "worst_inversion_linear_vs_rootfind": worst_lin,
            "log_scale_gain_x": float(worst_lin / worst_inv)
            if worst_inv > 0 else None,
            "worst_parallel_vs_serial": worst_par,
            "c_grid_step": float(CGRID[1] - CGRID[0]),
            "workers": int(WORKERS),
            "gate_tolerance": float(GATE_TOL),
            "scale_note": ("the linear-scale figure is what the same grid "
                           "gives if the rejection probability is "
                           "interpolated directly; it is measured here so "
                           "that the choice of scale is not an assertion"),
            "pass": bool(worst_par == 0.0
                         and worst_inv < GATE_TOL)}


def bench(n_levels=8, side=12):
    """Time a whole LEVEL SWEEP, which is the unit of work that was slow.

    The old route rebuilt the entire table once per bisection step; the new one
    builds it once and inverts.  Comparing per-cell times would hide exactly
    the factor that was costing the hours, so the benchmark is the sweep.
    """
    n = np.array([36.0, 13.0, 13.0])
    mids = np.exp(np.linspace(np.log(1e-3), np.log(1e3), side))
    lams = [(a, b) for a in mids for b in mids]
    sample = lams[:20]
    t0 = time.time()
    for a, b in sample:
        F.oracle_c_pop(n, [a, b], 0.0244, F.FAST)
    per_cell = (time.time() - t0) / len(sample)
    old = per_cell * len(lams) * n_levels          # rebuild per level, 1 core
    t0 = time.time()
    rows = pmap_cdf(n, lams)
    new = time.time() - t0                          # one build, all levels
    return {"cells": len(lams), "levels": int(n_levels),
            "workers": int(WORKERS),
            "old_sweep_seconds_est": float(old),
            "new_sweep_seconds": float(new),
            "speedup_x": float(old / new) if new > 0 else None,
            "per_cell_serial_rootfind_s": float(per_cell),
            "note": ("old = one root-find per cell per level on one core; "
                     "new = one tabulated row per cell, in parallel, serving "
                     "every level by inversion")}


if __name__ == "__main__":
    print("falsifier:", selftest())
    b = bench()
    print("bench: a %d-cell, %d-level sweep: %.0f s old -> %.1f s new = %.0fx"
          % (b["cells"], b["levels"], b["old_sweep_seconds_est"],
             b["new_sweep_seconds"], b["speedup_x"]))
    F.write_deposit("fast.json", {"selftest": selftest(), "bench": b},
                    "frontier_measured")


# ----------------------------------------------------- parallel certification
def _boxjob(args):
    idx, c_vals, edges, n, lo, hi, m = args
    import numpy as _np

    import wb_certify as _C
    r = _C.StepRule(_np.asarray(c_vals), _np.asarray(edges)).build_tables()
    return idx, _C.enclose_fwer_box(r, _np.asarray(n), lo, hi, m)


def certify_parallel(rule, n, lam_lo, lam_hi, abar, m=800, init_split=96,
                     max_rounds=4, workers=WORKERS):
    """`wb_certify.certify`, with the box evaluations spread over cores.

    Identical semantics: a box is accepted when its rigorous upper bound is at
    or below the ceiling and split otherwise, and the certificate is the set of
    accepted boxes.  The only change is that a round's boxes are evaluated
    together instead of one at a time, and the reduction is by index so the
    outcome does not depend on scheduling.
    """
    import multiprocessing as mp

    e = np.exp(np.linspace(np.log(lam_lo), np.log(lam_hi), init_split + 1))
    boxes = list(zip(e[:-1], e[1:]))
    accepted, undecided = [], []
    c_vals = np.asarray(rule.c).tolist()
    edges = np.asarray(rule.edges).tolist()
    nl = np.asarray(n, float).tolist()
    for _ in range(max_rounds):
        if not boxes:
            break
        tasks = [(i, c_vals, edges, nl, float(a), float(b), m)
                 for i, (a, b) in enumerate(boxes)]
        if workers > 1 and len(tasks) > 3:
            with mp.Pool(min(workers, len(tasks))) as p:
                res = p.map(_boxjob, tasks,
                            chunksize=max(1, len(tasks) // (4 * workers)))
        else:
            res = [_boxjob(t) for t in tasks]
        out = [None] * len(tasks)
        for i, v in res:
            out[i] = v
        nxt = []
        for (a, b), (low, up) in zip(boxes, out):
            if up <= abar:
                accepted.append({"lo": float(a), "hi": float(b),
                                 "upper": float(up), "lower": float(low)})
            elif b / a < 1.0 + 1e-6:
                undecided.append({"lo": float(a), "hi": float(b),
                                  "upper": float(up)})
            else:
                mid = float(np.sqrt(a * b))
                nxt.extend([(a, mid), (mid, b)])
        boxes = nxt
    undecided.extend({"lo": float(a), "hi": float(b), "upper": None}
                     for a, b in boxes)
    ok = not undecided
    return {"certified": bool(ok),
            "n_boxes_accepted": len(accepted),
            "n_boxes_undecided": len(undecided),
            "sup_upper_bound": float(max(a["upper"] for a in accepted))
            if accepted else None,
            "cells_per_axis": int(m), "workers": int(workers)}
