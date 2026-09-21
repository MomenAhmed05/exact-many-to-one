"""S1.4 leak probe: reproduce the family-B df-3 dominance violation, explain it.

The Stage-A dominance slice found one draw (of 1200) on which the envelope's
conditional rejection at the TRUE lambda exceeded the population oracle's by
+0.0118.  A certificate that cannot explain that draw is not a certificate, so
this probe

  1. re-derives the slice (family B, df 3, the slice's own seed) and locates the
     violating draw(s) by paired conditional rejection (cheap, no grid);
  2. for the worst violator, evaluates the POPULATION critical value on a fine
     log-lambda grid inside that draw's own region and compares it with the corner
     maximum (the candidate set), naming the argmax lambda;
  3. reports the deficit in c and the implied rejection excess, and whether the
     argmax sits at a corner (attained) or strictly inside (the leak class).

Outputs: deposits/cert/leak_probe.json
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, batch_c, cond_noreject, draw_var, oracle_c  # noqa: E402
from env_rule import c_env, lam_region  # noqa: E402
from env_sweep import cell_seed  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "cert")
GAMMA = 0.0005
DRAWS = 200
def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


P_FINE = _env_int("W26_LEAK_P", 25)
POOL_N = _env_int("W26_POOL", 3)


def main():
    lam = np.array([1.0, 25.0])
    n = np.array([4.0, 13.0, 13.0])
    seed = cell_seed("B", 3) + 101
    s0, si, lh = draw_var(n, lam, DRAWS, seed)
    lo, hi, reg = lam_region(n, lh, GAMMA, scale=None)
    c_cand = c_env(n, s0, si, lh, gamma=GAMMA, alpha=ALPHA)
    c_or = oracle_c(n, lam, ALPHA)
    rej_cand = 1.0 - cond_noreject(c_cand, s0, si, n, lam)
    rej_or = 1.0 - cond_noreject(np.full(DRAWS, c_or), s0, si, n, lam)
    excess = rej_cand - rej_or
    try:
        j = int(np.argmax(excess))
        n_pos = int(np.count_nonzero(excess > 1e-9))
        worst_excess = float(excess[j])
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("violation scan failed") from exc
    inside = bool(np.all((lam >= lo[j]) & (lam <= hi[j])))

    # fine population supremum inside THIS draw's region
    ax = [np.geomspace(lo[j][i], hi[j][i], P_FINE) for i in range(2)]
    grid = np.array([[a, b] for a in ax[0] for b in ax[1]], float)
    c_pop = batch_c(n, grid, alpha=ALPHA, workers=POOL_N)
    try:
        k = int(np.argmax(c_pop))
        sup_fine = float(c_pop[k])
        argmax_lam = grid[k]
        corners = np.array([[a, b] for a in (lo[j][0], hi[j][0])
                            for b in (lo[j][1], hi[j][1])], float)
        c_corners = batch_c(n, corners, alpha=ALPHA, workers=POOL_N)
        corner_max = float(np.max(c_corners))
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("fine-grid sup failed") from exc
    at_corner = bool(np.any(np.all(np.isclose(argmax_lam[None, :], corners,
                                              rtol=1e-9), axis=1)))

    # implied rejection excess at the argmax lambda, holding the draw fixed
    def rej_at(cvec, lam_pt):
        try:
            return 1.0 - float(cond_noreject(np.array([cvec]), np.array([s0[j]]),
                                             si[j:j + 1], n, lam_pt)[0])
        except Exception:  # noqa: BLE001
            return float("nan")

    implied = rej_at(c_cand[j], argmax_lam) - rej_at(c_or, argmax_lam)
    out = {
        "family": "B", "df0": 3, "n": list(n), "lam_true": list(lam),
        "slice_seed": i64(seed), "draws": i64(DRAWS),
        "c_oracle": f64(c_or),
        "n_dominance_violations": i64(n_pos),
        "n_possible": i64(DRAWS),
        "worst_excess_rej": f64(worst_excess),
        "worst_index": i64(j),
        "worst_inside_region": inside,
        "worst_lambda_hat": [f64(v) for v in lh[j]],
        "region_lo": [f64(v) for v in lo[j]],
        "region_hi": [f64(v) for v in hi[j]],
        "region_scale": f64(reg["scale"]),
        "corner_max_this_draw": f64(corner_max),
        "fine_grid_points": i64(len(grid)),
        "sup_pop_fine": f64(sup_fine),
        "argmax_lam": [f64(v) for v in argmax_lam],
        "argmax_is_corner": at_corner,
        "corner_deficit_c": f64(corner_max - sup_fine),
        "implied_rej_excess_at_argmax": f64(implied),
        "explanation": ("corner set attains the sup on this draw" if at_corner else
                        "an interior lambda attains a higher population c* than every "
                        "corner on this draw: the leak class the certificate removes"),
    }
    dump_json(out, os.path.join(OUT, "leak_probe.json"))
    print(f"slice: {n_pos}/{DRAWS} dominance violations | worst excess "
          f"{worst_excess:+.5f} at index {j}")
    print(f"that draw: corner max {corner_max:.4f} | fine sup {sup_fine:.4f} "
          f"(argmax corner={at_corner}) | deficit {corner_max - sup_fine:+.5f} | "
          f"implied rej excess {implied:+.5f}")
    print("explanation:", out["explanation"])


if __name__ == "__main__":
    main()
