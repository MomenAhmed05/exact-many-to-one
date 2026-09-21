"""S1.2: does the candidate set attain the POPULATION supremum over the region?

The guarantee needs c_env >= c_star_pop(lam) for every lam in the region, where
c_star_pop is the population critical value (draw-free, `core.batch_c`).  The
candidate envelope takes the max over the region's corners.  This script measures,
per desk cell, the population supremum on a dense log-lam grid (11 points per
coordinate, corners included) and reports:

  * sup_pop and the argmax lam, and whether the argmax is a corner;
  * the candidate envelope on the same region (and on draws);
  * the certified envelope (cert_sup.c_rigorous);
  * deficit_candidate = cand - sup_pop   (negative => the candidate set misses the sup)
  * margin_certified  = certified - sup_pop (must be >= 0 for the certificate)

An interior-grid sup that beats every corner IS the leak class, so this is the
direct falsifier for the sup-at-candidate step.

Outputs: deposits/cert/sup_falsifier.json
"""
import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cert_sup import DFS, FAMS, GAMMA, certify_draws  # noqa: E402
from core import ALPHA, batch_c, draw_var  # noqa: E402
from env_rule import c_env, lam_region  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "cert")
def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


P_GRID = _env_int("W26_FALS_P", 11)
DRAWS = _env_int("W26_FALS_D", 20)
POOL_N = _env_int("W26_POOL", 2)


def cell(fam, df, pool):
    lam = np.asarray(FAMS[fam], float)
    n = np.array([df + 1, 13, 13], float)
    s0, si, lh = draw_var(n, lam, DRAWS, 4242 + 97 * df)
    lo, hi, _ = lam_region(n, lh, GAMMA, scale=None)
    # dense log-grid over the region (per-coordinate geometric spacing)
    ax = [np.geomspace(lo[0][i], hi[0][i], P_GRID) for i in range(len(lam))]
    grid = np.array([[a, b] for a in ax[0] for b in ax[1]], float)
    c_pop = batch_c(n, grid, alpha=ALPHA, pool=pool)
    try:
        j_max = int(np.argmax(c_pop))
        sup_pop = float(c_pop[j_max])
        argmax = grid[j_max]
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("sup grid argmax failed") from exc
    corners = np.array([[a, b] for a in (lo[0][0], hi[0][0])
                        for b in (lo[0][1], hi[0][1])], float)
    c_corners = batch_c(n, corners, alpha=ALPHA, pool=pool)
    is_corner = bool(np.any(np.all(np.isclose(argmax[None, :], corners,
                                              rtol=1e-9), axis=1)))
    c_cand = c_env(n, s0, si, lh, gamma=GAMMA, alpha=ALPHA, pool=pool)
    c_rig = certify_draws(n, s0, si, lo, hi, alpha=ALPHA, pool=pool)
    try:
        cand_mean = float(np.mean(c_cand))
        corners_max = float(np.max(c_corners))
        rig_mean = float(np.mean(c_rig))
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("cell aggregation failed") from exc
    return {
        "family": fam, "df0": df, "n": list(n), "lam_true": list(lam),
        "grid_points": i64(len(grid)), "draws": i64(DRAWS),
        "sup_pop_grid": f64(sup_pop), "argmax_lam": [f64(v) for v in argmax],
        "argmax_is_corner": is_corner,
        "corners_max": f64(corners_max),
        "c_candidate_mean": f64(cand_mean),
        "c_rigorous_mean": f64(rig_mean),
        "deficit_candidate": f64(cand_mean - sup_pop),
        "deficit_corners_only": f64(corners_max - sup_pop),
        "margin_certified": f64(rig_mean - sup_pop),
        "candidate_covers_sup": bool(cand_mean >= sup_pop - 1e-9),
        "corners_cover_sup": bool(corners_max >= sup_pop - 1e-9),
        "certified_covers_sup": bool(rig_mean >= sup_pop - 1e-9),
    }


def main():
    t0 = time.time()
    rows = []
    with mp.Pool(POOL_N) as pool:
        for df in DFS:
            for fam in FAMS:
                r = cell(fam, df, pool)
                rows.append(r)
                print(f"{fam} df{df:>2}: sup_pop {r['sup_pop_grid']:.4f} "
                      f"(corner {r['argmax_is_corner']}) | corners {r['corners_max']:.4f} "
                      f"| cand {r['c_candidate_mean']:.4f} "
                      f"(deficit {r['deficit_corners_only']:+.4f}) | certified "
                      f"{r['c_rigorous_mean']:.4f} (margin {r['margin_certified']:+.4f})",
                      flush=True)
    summary = {
        "rows": rows, "alpha": ALPHA, "gamma": GAMMA, "grid_per_axis": i64(P_GRID),
        "draws_per_cell": i64(DRAWS),
        "all_corners_cover_sup": bool(all(r["corners_cover_sup"] for r in rows)),
        "all_certified_cover_sup": bool(all(r["certified_covers_sup"] for r in rows)),
        "worst_corner_deficit": f64(min(r["deficit_corners_only"] for r in rows)),
        "n_argmax_at_corner": i64(sum(1 for r in rows if r["argmax_is_corner"])),
        "seconds": f64(time.time() - t0),
        "note": ("A negative deficit means an interior grid point beat every corner, "
                 "which is the leak class the certificate exists to remove; the "
                 "certified column must never be negative."),
    }
    dump_json(summary, os.path.join(OUT, "sup_falsifier.json"))
    print(f"corners cover sup: {summary['all_corners_cover_sup']} | certified cover sup: "
          f"{summary['all_certified_cover_sup']} | worst corner deficit "
          f"{summary['worst_corner_deficit']:+.4f} | argmax at corner "
          f"{summary['n_argmax_at_corner']}/{len(rows)} | {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
