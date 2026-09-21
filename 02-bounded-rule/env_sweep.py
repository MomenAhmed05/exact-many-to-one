"""12-cell A/B envelope sweep: dual region arms, per-cell seeds, deterministic check.

Primary evidence is deliberately NOT the 200-draw FWER estimate.  Its Monte Carlo
standard error (1.3e-4 .. 9.4e-4, from the deposits) is comparable to or larger
than the 5e-4 slack at df 3, so it cannot certify the guarantee.  The primary
check is deterministic per draw: evaluate the CONDITIONAL rejection probability
at the TRUE lambda under c_env, restricted to draws where lambda_true lies inside
the region.  On that event the envelope argument requires it to be <= alpha
pointwise; `max_rej_excess_inside` is the direct falsifier, and the empirical
region-miss rate enters FWER <= alpha + P(lambda not in Lambda).

Two arms per cell: `tight` (exact-coverage region, scale > 1) and `ub` (the
union-bound allocation, scale = 1), so the power cost of the tightening is
measured on the same draws.
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, batch_c, cond_noreject, draw_var, oracle_c  # noqa: E402
from env_rule import c_env, fwer_of, lam_region  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "env")
FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
FAM_IDX = {"A": 1, "B": 2}
DFS = (3, 5, 7, 11, 19, 35)
GAMMA = 0.0005
SEEDS = (101, 202, 303, 404)
D_PER_SEED = 50


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


DRAWS = D_PER_SEED * len(SEEDS)
POOL_N = _env_int("W26_POOL", 4)


def cell_seed(fam, df):
    """Distinct, documented seed per cell (the audit flagged a shared stream)."""
    return 97 * df + 131 * FAM_IDX[fam] + 4242


def _arm(n, lam, s0, si, lh, alpha, gamma, tight, pool, c_or):
    """One region arm, with the CORRECT falsifier.

    The guarantee's pointwise step is `rej(c_env, s^2) <= rej(c_star_pop, s^2)` on
    every draw whose lambda_true lies inside the region -- NOT `rej <= alpha`
    (the population oracle is calibrated only on average over s^2, so on a single
    draw it can sit far above alpha).  `excess_vs_oracle_max` is therefore the
    invariant, and `alpha + mean(excess)` estimates the realised FWER.
    """
    c = c_env(n, s0, si, lh, alpha=alpha, gamma=gamma, pool=pool, tight=tight)
    lo, hi, reg = lam_region(n, lh, gamma, scale=(None if tight else 1.0))
    try:
        rej = 1.0 - cond_noreject(c, s0, si, n, lam)
        c_or_vec = np.full(int(np.size(rej)), float(c_or))
        rej_or = 1.0 - cond_noreject(c_or_vec, s0, si, n, lam)
        inside = np.all((lam[None, :] >= lo) & (lam[None, :] <= hi), axis=1)
        n_out = int(np.count_nonzero(~inside))
        ndraw = int(np.size(rej))
        d = rej - rej_or
        worst = float(np.max(d[inside])) if bool(np.any(inside)) else float("nan")
        payload = {
            "c_mean": f64(float(np.mean(c))),
            "rej_mean": f64(float(np.mean(rej))),
            "excess_vs_oracle_max": f64(worst),
            "excess_vs_oracle_mean": f64(float(np.mean(d))),
            "fwer_est": f64(alpha + float(np.mean(d))),
            "n_outside": i64(n_out),
            "frac_outside": f64(n_out / ndraw),
            "region_scale": f64(reg["scale"]),
            "region_coverage": f64(reg["coverage_exact"]),
        }
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("arm aggregation failed") from exc
    return payload


def cell(fam, df, pool):
    lam = np.asarray(FAMS[fam], float)
    n = np.array([df + 1, 13, 13], float)
    c_or = oracle_c(n, lam, ALPHA)
    per_seed = []
    for s in SEEDS:
        seed = cell_seed(fam, df) + s
        s0, si, lh = draw_var(n, lam, D_PER_SEED, seed)
        t = _arm(n, lam, s0, si, lh, ALPHA, GAMMA, True, pool, c_or)
        u = _arm(n, lam, s0, si, lh, ALPHA, GAMMA, False, pool, c_or)
        c_p = batch_c(n, lh, alpha=ALPHA, pool=pool)
        r_p = fwer_of(n, lam, c_p, s0, si)
        t["power_cost_pct"] = f64(100 * (t["c_mean"] - c_or) / c_or)
        u["power_cost_pct"] = f64(100 * (u["c_mean"] - c_or) / c_or)
        t["seed"] = i64(seed)
        t["plugin_fwer"] = r_p["fwer"]
        t["plugin_cost_pct"] = r_p["power_cost_pct"]
        per_seed.append({**t, "ub": u})
    return {
        "family": fam,
        "df0": df,
        "lam": list(lam),
        "n": list(n),
        "alpha": ALPHA,
        "gamma": GAMMA,
        "draws_per_seed": D_PER_SEED,
        "seeds": [i64(s) for s in SEEDS],
        "c_oracle": f64(c_or),
        "per_seed": per_seed,
    }


def _aggregate(fam, df, res):
    """Cell aggregates; all conversions happen here under one guard."""
    rows = res["per_seed"]
    try:
        tight_cost = float(np.mean([r["power_cost_pct"] for r in rows]))
        ub_cost = float(np.mean([r["ub"]["power_cost_pct"] for r in rows]))
        rej_mean = float(np.mean([r["rej_mean"] for r in rows]))
        fwer_tight = float(np.mean([r["fwer_est"] for r in rows]))
        exc = [float(r["excess_vs_oracle_max"]) for r in rows]
        exc = [e for e in exc if np.isfinite(e)]
        max_exc = max(exc) if exc else float("nan")
        n_out = sum(int(r["n_outside"]) for r in rows)
        ndraw = D_PER_SEED * len(SEEDS)
        bound = fwer_tight + n_out / ndraw
        cond_ok = all(
            (not np.isfinite(float(r["excess_vs_oracle_max"])))
            or float(r["excess_vs_oracle_max"]) <= 1e-9
            for r in rows
        )
        return {
            "family": fam,
            "df0": df,
            "c_oracle": res["c_oracle"],
            "tight_cost_pct": f64(tight_cost),
            "ub_cost_pct": f64(ub_cost),
            "tight_rej_mean": f64(rej_mean),
            "tight_fwer_est": f64(fwer_tight),
            "tight_excess_vs_oracle_max": f64(max_exc),
            "tight_n_outside": i64(n_out),
            "window_draws": i64(ndraw),
            "guarantee": f64(ALPHA + GAMMA),
            "passes_guarantee": bool(bound <= ALPHA + GAMMA),
            "passes_pointwise": bool(cond_ok),
        }
    except Exception as exc2:  # noqa: BLE001
        raise RuntimeError(f"aggregate failed for {fam} df{df}") from exc2


def main():
    t0 = time.time()
    summary = []
    with mp.Pool(POOL_N) as pool:
        for df in DFS:
            for fam in FAMS:
                res = cell(fam, df, pool)
                dump_json(res, os.path.join(OUT, f"sweep_df{df}_{fam}.json"))
                agg = _aggregate(fam, df, res)
                summary.append(agg)
                print(
                    f"{fam} df{df:>2}: tight {agg['tight_cost_pct']:+6.2f}% "
                    f"(ub {agg['ub_cost_pct']:+6.2f}%) | fwer "
                    f"{agg['tight_fwer_est']:.4f} | exc_vs_oracle "
                    f"{agg['tight_excess_vs_oracle_max']:+.2e} | out "
                    f"{agg['tight_n_outside']} | pointwise "
                    f"{agg['passes_pointwise']}",
                    flush=True,
                )
    viol = [s for s in summary if not s["passes_guarantee"]]
    cviol = [s for s in summary if not s["passes_pointwise"]]
    try:
        mt = float(np.mean([s["tight_cost_pct"] for s in summary]))
        mu = float(np.mean([s["ub_cost_pct"] for s in summary]))
    except Exception:  # noqa: BLE001
        mt = mu = float("nan")
    out = {
        "alpha": ALPHA,
        "gamma": GAMMA,
        "seeds": list(SEEDS),
        "draws_per_seed": D_PER_SEED,
        "cells": summary,
        "mean_tight_cost_pct": f64(mt),
        "mean_ub_cost_pct": f64(mu),
        "violations": viol,
        "conditional_violations": cviol,
        "all_pass": bool(not viol and not cviol),
        "note": (
            "primary check is the per-draw conditional rejection at the "
            "true lambda (deterministic); the FWER MC cannot certify "
            "the 5e-4 slack at df 3"
        ),
    }
    dump_json(out, os.path.join(OUT, "sweep_summary.json"))
    print(
        f"mean cost: tight {mt:+.2f}% vs ub {mu:+.2f}% | guarantee viol "
        f"{len(viol)} | conditional viol {len(cviol)} | {time.time() - t0:.0f}s"
    )


if __name__ == "__main__":
    main()
