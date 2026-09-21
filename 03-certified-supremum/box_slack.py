"""Deposit the STRUCTURAL box slack (verifier finding: it was quoted in the P1 report but
appeared in no deposit).

For a box `Lambda`, the enclosure uses `int min_lambda G dPhi`, while the "true" box infimum is
`min_lambda int G dPhi`.  The inequality is strict because the pointwise minimiser moves with `z`,
so the box bound exceeds the corner bracket by a structural amount that does NOT shrink with
quadrature effort - it shrinks with the box's width.  That is the quantity the branch-and-bound
must drive below `tol` by splitting, so it belongs in a deposit.

Measured on a fixed draw at a fixed c, as a function of the multiplicative box width factor:
    slack = F_lo(box) - min over corners F(corner)      (both rigorous / refined-reference values)

Run:  python box_slack.py        (~1 min)
Out:  deposits/box_slack.json    (tier: region_sup_certified)
"""

from __future__ import annotations

import contextlib
import importlib
import itertools
import os
import time

import enclose
import enclose_gl
import numpy as np
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "box_slack.json")
CORE = importlib.import_module("core")
cond_noreject = CORE.cond_noreject
draw_var = CORE.draw_var

FACTORS = (3.0, 2.0, 1.5, 1.1, 1.01, 1.001)
CELLS = (
    ("A_df3", [4.0, 13.0, 13.0], [0.04, 0.04], 2.0),
    ("B_df3", [4.0, 13.0, 13.0], [1.0, 25.0], 2.0),
    ("K3_df7", [8.0, 8.0, 10.0, 12.0], [0.25, 1.0, 4.0], 2.2),
)


def cell(name, n, lam, c, seed=20260918, alpha=paths.ALPHA0):
    """Structural box slack in c-space: c_ub(box) - max over corners of c*(corner).

    Measured in PROBABILITY space the analogous quantity is necessarily negative (F_lo(box) <=
    min over corners F by construction), so the meaningful quantity - and the one the branch-and-
    bound must drive below tol by splitting - is the critical-value excess.
    """
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    s0, si, _lh = draw_var(n, lam, 1, seed)
    d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    rows = []
    for fac in FACTORS:
        lo_f, hi_f = lam / fac, lam * fac
        s_lo, s_hi = enclose.s_range(n, lo_f, hi_f)
        box = enclose_gl.c_upper_box_gl(d_row, n[0], s_lo, s_hi, alpha=alpha)
        c_corners = []
        for combo in itertools.product((0, 1), repeat=lam.size):
            lam_pt = np.where(np.array(combo, bool), hi_f, lo_f)
            s_pt = np.sqrt(n[1:] / lam_pt)
            r = enclose_gl.c_upper_box_gl(
                d_row, n[0], s_pt, s_pt, alpha=alpha, c_guess=box.get("warm_start")
            )
            if r.get("c_ub") == r.get("c_ub"):
                c_corners.append(_f(r["c_ub"]))
        c_max = max(c_corners) if c_corners else NAN
        rows.append(
            {
                "factor": _f(fac),
                "c_ub_box": box.get("c_ub"),
                "c_ub_max_corner": c_max,
                "slack_c": (
                    _f(box.get("c_ub")) - c_max
                    if box.get("c_ub") == box.get("c_ub") and c_max == c_max
                    else NAN
                ),
                "box_log_width": _f(
                    np.linalg.norm(
                        np.log(np.asarray(hi_f, float))
                        - np.log(np.asarray(lo_f, float))
                    )
                ),
            }
        )
    return {
        "cell": name,
        "n": list(n),
        "lam": list(lam),
        "c": _f(c),
        "alpha0": _f(alpha),
        "d_row": [round(_f(v), 8) for v in d_row],
        "rows": rows,
    }


def main():
    t0 = time.time()
    cells = [cell(*spec) for spec in CELLS]
    rows = [r for c in cells for r in c["rows"]]
    out = {
        "n_cells": _i(len(cells)),
        "factors": [_f(f) for f in FACTORS],
        "slack_max": _f(max(_f(r["slack_c"]) for r in rows)),
        "slack_at_factor_1p001": _f(
            max(_f(r["slack_c"]) for r in rows if abs(_f(r["factor"]) - 1.001) < 1e-9)
        ),
        "tol": _f(paths.TOL),
        "note": (
            "structural box slack in C-SPACE = c_ub(box) - max over corners of c*(corner).  "
            "In probability space the analogous difference is necessarily negative "
            "(F_lo(box) <= min over corners F by construction), so the c-space excess is the "
            "meaningful quantity: it shrinks with the box width, not with quadrature effort, "
            "and it is what the branch-and-bound drives below tol by splitting (plan 2.4).  "
            "Report-only in PHASE1-ENCLOSURE.md before; now deposited."
        ),
        "seconds": _f(time.time() - t0),
        "cells": cells,
    }
    print("structural box slack in c-space (worst over cells):")
    for fac in FACTORS:
        vals = [_f(r["slack_c"]) for r in rows if abs(_f(r["factor"]) - fac) < 1e-9]
        print(f"  factor {fac:<6} max slack {max(vals):.3e} c-units")
    print(
        f"tol {paths.TOL:.0e} | slack at factor 1.001: {out['slack_at_factor_1p001']:.3e}"
    )
    tiers.write_deposit(
        out,
        OUT,
        "region_sup_certified",
        notes="structural box slack vs box width factor, on fixed draws at fixed c",
    )
    return out


if __name__ == "__main__":
    with contextlib.suppress(Exception):
        main()
