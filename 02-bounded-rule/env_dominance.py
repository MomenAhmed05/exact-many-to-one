"""Honest dominance falsifier for the envelope's sup-over-region step.

The guarantee needs  rej(c_env, s^2) <= rej(c_star_pop, s^2)  on draws whose
lambda_true lies inside the region.  Any finite candidate set can miss the
supremum of c_star_pop over the region, so this cannot be certified by
candidate construction alone.  This script measures, on the sweep's own per-cell
streams for the small-df cells:

  * the c-deficit of the corner set relative to the 3^K candidate set;
  * the distribution of the per-draw rejection excess over the population oracle;
  * how often that excess is positive (a dominance failure) and by how much.

It does NOT establish dominance; it bounds how often it fails on these streams
and is deposited so the report can state the caveat with numbers.
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, cond_noreject, draw_var, oracle_c  # noqa: E402
from env_rule import c_env  # noqa: E402
from env_sweep import cell_seed  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "env")
FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
CELLS = (("A", 3), ("B", 3), ("A", 5), ("B", 5), ("A", 7), ("B", 7))
DRAWS = 200


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


def cell(fam, df, pool):
    lam = np.asarray(FAMS[fam], float)
    n = np.array([df + 1, 13, 13], float)
    s0, si, lh = draw_var(n, lam, DRAWS, cell_seed(fam, df) + 101)
    c3 = c_env(n, s0, si, lh, tight=True, pool=pool)
    cc = c_env(n, s0, si, lh, tight=True, pool=pool, max_vertices=4)
    c_or = oracle_c(n, list(lam), ALPHA)
    try:
        rej_or = 1.0 - cond_noreject(np.full(DRAWS, c_or), s0, si, n, lam)
        rej_3 = 1.0 - cond_noreject(c3, s0, si, n, lam)
        rej_c = 1.0 - cond_noreject(cc, s0, si, n, lam)
        exc_c = rej_c - rej_or
        exc_3 = rej_3 - rej_or
        deficit = cc - c3
        payload = {
            "family": fam,
            "df0": df,
            "draws": DRAWS,
            "corner_vs_3K_deficit_max": f64(float(deficit.max())),
            "corner_vs_3K_deficit_n": i64(int(np.count_nonzero(deficit > 1e-9))),
            "excess_corners_max": f64(float(exc_c.max())),
            "excess_3K_max": f64(float(exc_3.max())),
            "excess_corners_n_positive": i64(int(np.count_nonzero(exc_c > 1e-9))),
            "excess_3K_n_positive": i64(int(np.count_nonzero(exc_3 > 1e-9))),
            "excess_corners_mean": f64(float(exc_c.mean())),
            "excess_3K_mean": f64(float(exc_3.mean())),
            "fwer_est_corners": f64(ALPHA + float(exc_c.mean())),
            "fwer_est_3K": f64(ALPHA + float(exc_3.mean())),
            "passes_mean_bound_corners": bool(
                ALPHA + float(exc_c.mean()) <= ALPHA + 0.0005
            ),
            "passes_mean_bound_3K": bool(ALPHA + float(exc_3.mean()) <= ALPHA + 0.0005),
        }
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"dominance aggregation failed for {fam} df{df}") from exc
    return payload


def main():
    t0 = time.time()
    rows = []
    with mp.Pool(_env_int("W26_POOL", 4)) as pool:
        for fam, df in CELLS:
            r = cell(fam, df, pool)
            rows.append(r)
            print(
                f"{fam} df{df:>2}: corners-3K deficit "
                f"{r['corner_vs_3K_deficit_max']:+.2e} "
                f"(n={r['corner_vs_3K_deficit_n']}) | excess max "
                f"{r['excess_3K_max']:+.4f} "
                f"(n+={r['excess_3K_n_positive']}) | mean excess "
                f"{r['excess_3K_mean']:+.5f} | fwer_est "
                f"{r['fwer_est_3K']:.4f}",
                flush=True,
            )
    fail = [r for r in rows if not r["passes_mean_bound_3K"]]
    dump_json(
        {
            "cells": rows,
            "alpha": ALPHA,
            "gamma": 0.0005,
            "draws_per_cell": DRAWS,
            "n_mean_bound_violations": i64(len(fail)),
            "note": (
                "dominance over the population oracle is NOT certified by the "
                "candidate set; this deposit measures how often and by how much "
                "the per-draw rejection exceeds the oracle, and whether the "
                "AVERAGED FWER estimate still respects alpha+gamma."
            ),
            "seconds": f64(time.time() - t0),
        },
        os.path.join(OUT, "dominance_slice.json"),
    )
    print(f"mean-bound violations {len(fail)}/{len(rows)}")


if __name__ == "__main__":
    main()
