"""P1 step 8: derivative diagnostics, the predeclared Lipschitz materiality test, and the
mixed-corner geometry that seeds the P2 branch-and-bound.

Three deliverables (plan sections 2.1-2.3 and 2.5):

1. **Derivative check.** `cert_sup.dcdloglam` evaluates the implicit-function derivative

       d c* / d log lambda_i = - (dF/d log lambda_i) / (dF/dc)

   at a fixed draw.  It is compared here against central finite differences of the
   *conditional* critical value, solved by bisection on the independent solver
   `core.cond_noreject` - the same object, so the comparison is apples to apples.  (The
   first draft of the P1 falsifiers compared a conditional bracket with the *population*
   oracle; that was a test-design error, recorded in PROGRESS as finding 10.)

2. **Predeclared Lipschitz materiality test** (plan section 2.3.5).  For a uniform grid of
   `N^K` log-ratio cells the Lipschitz-grid certificate would add

       correction = max_C L_C * diam(C),      L_C = max |d c* / d log lambda| over the cell,

   to the grid maximum.  The test asks whether that correction can be driven below the
   predeclared `TOL = 1e-4`, at what grid size, and whether that grid is feasible for
   `K = 1..5`.  The answer is recorded either way; the Lipschitz route is a diagnostic and
   a potential accelerator, never the foundation of the theorem (the certificate path uses
   the branch-and-bound of plan section 2.4 instead).  Note the honest limitation: `L_C` is
   sampled at cell corners here, so this measures whether the route is *worth pursuing*
   with rigorous enclosures, not whether a certificate exists.

3. **Mixed-corner geometry.**  For each desk cell the population critical value is maximised
   over a dense log-lambda grid and the argmax is classified as a vertex (every coordinate
   at a bound) or mixed/interior.  This is the seeding information for P2: the B&B witness
   set starts from vertices, face midpoints, the box centre and the mixed corners found
   here.

Run:  python deriv_diag.py        (a few minutes)
Out:  deposits/deriv_diag.json
"""

from __future__ import annotations

import contextlib
import importlib
import os
import time

import derivs
import numpy as np
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "deriv_diag.json")
CORE = importlib.import_module("core")
CERT = importlib.import_module("cert_sup")
ENV = importlib.import_module("env_rule")
cond_noreject = CORE.cond_noreject
draw_var = CORE.draw_var

FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
DFS = (3, 5, 7, 11, 19, 35)


def c_cond(s0, si, n, lam, alpha, lo=0.3, hi=12.0, tol=1e-10):
        """Conditional critical value at a fixed draw: the smallest c whose conditional
        NO-REJECTION probability reaches 1 - alpha (bisection on the independent solver).

        The first version compared the *rejection* tail against `1 - alpha`, which is never true,
        so every call returned NaN and the derivative columns came back empty (the FD check
        reported FAIL with a zero deviation - a green-by-silence failure worth remembering).
        """
        try:

                def rej(c):
                        return 1.0 - _f(
                                cond_noreject(
                                        np.array([c]),
                                        np.array([s0]),
                                        np.asarray(si, float)[None, :],
                                        np.asarray(n, float),
                                        np.asarray(lam, float),
                                )[0]
                        )

                if rej(hi) > alpha:  # even c = hi rejects too often: no bracket
                        return NAN
                if rej(lo) <= alpha:  # lo already attains the level
                        return lo
                a, b = lo, hi
                while b - a > tol:
                        m = 0.5 * (a + b)
                        if rej(m) <= alpha:
                                b = m
                        else:
                                a = m
                return b
        except Exception:  # noqa: BLE001
                return NAN


def fd_check(alpha=paths.ALPHA0, h=2e-3, seed=20260918):
        """`dcdloglam` vs central finite differences of the conditional critical value."""
        rows = []
        for name, lam0 in FAMS.items():
                for df in DFS:
                        n = np.array([df + 1.0, 13.0, 13.0])
                        lam = np.asarray(lam0, float)
                        s0, si, _ = draw_var(n, lam, 1, seed + df)
                        d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
                        c0 = c_cond(s0[0], si[0], n, lam, alpha)
                        an = derivs.dcdloglam(n, d_row, lam, c0) if c0 == c0 else []
                        if not an:
                                rows.append(
                                        {
                                                "family": name,
                                                "df0": _i(df),
                                                "n": list(n),
                                                "lam": list(lam),
                                                "c_conditional": c0,
                                                "analytic": [],
                                                "finite_difference": [],
                                                "abs_dev": [],
                                                "max_rel_dev": None,
                                                "agree": False,
                                                "error": "derivative unavailable (c_cond NaN or empty result)",
                                        }
                                )
                                continue
                        fd = []
                        for i in range(lam.size):
                                up = lam.copy()
                                dn = lam.copy()
                                up[i] = lam[i] * np.exp(h)
                                dn[i] = lam[i] * np.exp(-h)
                                cu = c_cond(s0[0], si[0], n, up, alpha)
                                cd = c_cond(s0[0], si[0], n, dn, alpha)
                                fd.append(_f((cu - cd) / (2.0 * h)))
                        dev = (
                                [abs(_f(an[i]) - fd[i]) for i in range(len(fd))]
                                if an
                                else []
                        )
                        rel = [dev[i] / max(1e-12, abs(fd[i])) for i in range(len(dev))]
                        rows.append(
                                {
                                        "family": name,
                                        "df0": _i(df),
                                        "n": list(n),
                                        "lam": list(lam),
                                        "c_conditional": c0,
                                        "analytic": [round(_f(v), 6) for v in an],
                                        "finite_difference": [
                                                round(_f(v), 6) for v in fd
                                        ],
                                        "abs_dev": [round(v, 8) for v in dev],
                                        "max_rel_dev": round(max(rel), 6)
                                        if rel
                                        else None,
                                        "agree": bool(rel and max(rel) < 5e-3),
                                }
                        )
        return {
                "rows": rows,
                "h": _f(h),
                "alpha": _f(alpha),
                "all_agree": bool(all(r["agree"] for r in rows)),
                "note": (
                        "dcdloglam is the implicit-function derivative at a fixed draw; the "
                        "reference is a central difference of the conditional critical value "
                        "solved by bisection on core.cond_noreject, so both objects are "
                        "conditional on the same draw"
                ),
        }


def _max_abs_deriv(cell_lo, cell_hi, n, d_row, c_mid, alpha):
        """Sampled max |dc*/dlog lam| over the cell's corners (an estimate, not a bound)."""
        best = 0.0
        for combo in np.ndindex(*([2] * len(cell_lo))):
                pt = np.where(
                        np.array(combo, bool),
                        np.asarray(cell_hi, float),
                        np.asarray(cell_lo, float),
                )
                try:
                        g = derivs.dcdloglam(n, d_row, pt, c_mid)
                except Exception:  # noqa: BLE001
                        g = []
                if g:
                        best = max(best, max(abs(_f(v)) for v in g))
        return best


def lipschitz_materiality(alpha=paths.ALPHA0, seed=20260918, Ns=(2, 4, 8, 16)):
        """Would a Lipschitz-grid certificate reach the predeclared tolerance, and at what cost?"""
        rows = []
        for name, lam0 in FAMS.items():
                for df in (3, 7, 19):
                        n = np.array([df + 1.0, 13.0, 13.0])
                        lam = np.asarray(lam0, float)
                        s0, si, lh = draw_var(n, lam, 1, seed + df)
                        d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
                        lo, hi, _ = ENV.lam_region(n, lh, paths.GAMMA, scale=None)
                        lo = np.asarray(lo[0], float)
                        hi = np.asarray(hi[0], float)
                        K = lo.size
                        c_mid = c_cond(s0[0], si[0], n, np.sqrt(lo * hi), alpha)
                        for N in Ns:
                                edges = [
                                        np.geomspace(lo[i], hi[i], N + 1)
                                        for i in range(K)
                                ]
                                worst = 0.0
                                for combo in np.ndindex(*([N] * K)):
                                        cl = np.array(
                                                [edges[i][combo[i]] for i in range(K)]
                                        )
                                        ch = np.array(
                                                [
                                                        edges[i][combo[i] + 1]
                                                        for i in range(K)
                                                ]
                                        )
                                        diam = _f(
                                                np.linalg.norm(np.log(ch) - np.log(cl))
                                        )
                                        L = _max_abs_deriv(
                                                cl, ch, n, d_row, c_mid, alpha
                                        )
                                        worst = max(worst, L * diam)
                                cells = N**K
                                rows.append(
                                        {
                                                "family": name,
                                                "df0": _i(df),
                                                "K": _i(K),
                                                "N_per_axis": _i(N),
                                                "cells": _i(cells),
                                                "correction": _f(worst),
                                                "reaches_tol": bool(worst <= paths.TOL),
                                                "evals_needed": _i(cells * 1),
                                                "feasible": bool(cells <= 100000),
                                        }
                                )
        ok = [r for r in rows if r["reaches_tol"] and r["feasible"]]
        return {
                "rows": rows,
                "tol": _f(paths.TOL),
                "min_correction": _f(min(r["correction"] for r in rows)),
                "feasible_and_within_tol": len(ok),
                "verdict": (
                        "the Lipschitz route can reach the predeclared tolerance on a feasible "
                        "grid for these cells"
                        if ok
                        else "the Lipschitz route does NOT reach the predeclared tolerance on a "
                        "feasible grid for these cells: it is a diagnostic, not the certificate"
                ),
                "caveat": (
                        "L_C is sampled at cell corners, so this measures whether the route is "
                        "worth pursuing with rigorous enclosures, not whether a certificate "
                        "exists; the certificate path is the branch-and-bound of plan 2.4"
                ),
        }


def argmax_geometry(alpha=paths.ALPHA0, seed=20260918, per_axis=9):
        """Where does the population critical value attain its grid maximum, and at what cost?"""
        rows = []
        for name, lam0 in FAMS.items():
                for df in DFS:
                        n = np.array([df + 1.0, 13.0, 13.0])
                        lam = np.asarray(lam0, float)
                        _s0, _si, lh = draw_var(n, lam, 1, seed + df)
                        lo, hi, _ = ENV.lam_region(n, lh, paths.GAMMA, scale=None)
                        lo = np.asarray(lo[0], float)
                        hi = np.asarray(hi[0], float)
                        K = lo.size
                        axes = [np.geomspace(lo[i], hi[i], per_axis) for i in range(K)]
                        grid = np.array(list(np.ndindex(*([per_axis] * K))), float)
                        pts = np.empty_like(grid)
                        for i in range(K):
                                pts[:, i] = axes[i][grid[:, i].astype(int)]
                        c = CORE.batch_c(n, pts, alpha=alpha, workers=2)
                        j = _i(np.nanargmax(c))
                        arg = pts[j]
                        at_lo = np.isclose(arg, lo, rtol=1e-9)
                        at_hi = np.isclose(arg, hi, rtol=1e-9)
                        rows.append(
                                {
                                        "family": name,
                                        "df0": _i(df),
                                        "K": _i(K),
                                        "grid_points": _i(pts.shape[0]),
                                        "sup_grid": _f(np.nanmax(c)),
                                        "argmax": [round(_f(v), 8) for v in arg],
                                        "argmax_at_lo": [bool(v) for v in at_lo],
                                        "argmax_at_hi": [bool(v) for v in at_hi],
                                        "is_vertex": bool(np.all(at_lo | at_hi)),
                                        "n_coords_at_bound": _i(
                                                np.count_nonzero(at_lo | at_hi)
                                        ),
                                        "region_lo": [round(_f(v), 8) for v in lo],
                                        "region_hi": [round(_f(v), 8) for v in hi],
                                }
                        )
                        del c
        return {
                "rows": rows,
                "per_axis": _i(per_axis),
                "n_vertex": _i(sum(1 for r in rows if r["is_vertex"])),
                "n_mixed": _i(sum(1 for r in rows if not r["is_vertex"])),
                "note": (
                        "P2 seeds its witness set with vertices, face midpoints, the box centre and "
                        "any mixed corners found here"
                ),
        }


def main():
        t0 = time.time()
        res: dict = {"alpha0": paths.ALPHA0, "gamma": paths.GAMMA, "tol": paths.TOL}
        res["fd_check"] = fd_check()
        print(
                "FD check:",
                "PASS" if res["fd_check"]["all_agree"] else "FAIL",
                "max rel dev",
                max((r["max_rel_dev"] or 0) for r in res["fd_check"]["rows"]),
        )
        res["w26_derivative_defect"] = derivs.w26_version_defect()
        res["lipschitz_materiality"] = lipschitz_materiality()
        lm = res["lipschitz_materiality"]
        print(
                f"Lipschitz materiality: min correction {lm['min_correction']:.3e} vs tol "
                f"{lm['tol']:.0e} | feasible+within tol: {lm['feasible_and_within_tol']}"
        )
        print("   ", lm["verdict"])
        res["argmax_geometry"] = argmax_geometry()
        ag = res["argmax_geometry"]
        print(f"argmax geometry: vertices {ag['n_vertex']}, mixed {ag['n_mixed']}")
        res["seconds"] = _f(time.time() - t0)
        tiers.write_deposit(
                res,
                OUT,
                "envelope_measured",
                notes="diagnostics for the measured-tier geometry: finite-difference "
                "checks of the implicit-function derivative, a SAMPLED Lipschitz "
                "materiality estimate (not a bound), and the argmax geometry used "
                "to seed the P2 branch-and-bound",
        )
        print(f"wrote {OUT} ({_f(res['seconds']):.0f}s)")
        return res


if __name__ == "__main__":
        with contextlib.suppress(Exception):
                main()
