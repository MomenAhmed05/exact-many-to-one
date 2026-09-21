"""W2: the self-calibrating, direction-aware confidence-bound rule (CCB).

Idea: the CB family is c(lam_hat / q) with q > 1 the conservative direction when
c decreases in lam, and c(lam_hat * q) when it increases.  Rather than fixing
gamma (and therefore q) by assumption, choose it per data set by matching the
rule's estimated size to alpha in the lam_hat world:

  candidates   q in {1/F_gamma} U {F_gamma}, gamma in a grid, plus q = 1
  c_j          the exact critical value at the shrunk ratios
  FWER_est_j   mean_b [ 1 - P(no reject | inner draw b, world = lam_hat) ]
  pick the c_j whose FWER_est is the largest value <= alpha; if none, the
  smallest FWER_est (most conservative).

The inner draws are cheap (variance draws only) and the conditional rejection
probability is the exact 1-D integral, so the whole rule costs |candidates|
solves per data set - about fifteen times cheaper than the full calibrator.

Usage: python w2_adaptive.py [df ...]
Outputs: deposits/w2/adapt_df<df>_<fam>.json
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import (
    ALPHA,
    batch_c,
    cond_noreject,
    draw_var,
    f_quantile,  # noqa
    oracle_c,
)
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w2")
FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
GAMMAS = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95)


def candidate_matrix(n, lh, alpha=ALPHA, pool=None):
    """Per-draw conservative direction, then the gamma ladder on that side.

    Direction is detected at gamma = 0.8: if c(lh/q) > c(lh*q) the lower side is
    conservative for that draw.  Candidate rows are then lh/q_g (lower side) or
    lh*q_g (upper side) for the gamma ladder, plus the plugin row.
    """
    nu = np.asarray(n, float) - 1.0
    D, K = lh.shape
    q80 = np.array([f_quantile(nu[i + 1], nu[0], 0.8) for i in range(K)])
    c_lo = batch_c(n, lh / q80[None, :], alpha=alpha, pool=pool)
    c_hi = batch_c(n, lh * q80[None, :], alpha=alpha, pool=pool)
    down = c_lo > c_hi
    rows = []
    labels = []
    for g in GAMMAS:
        q = np.array([f_quantile(nu[i + 1], nu[0], g) for i in range(K)])
        rows.append(np.where(down[:, None], lh / q[None, :], lh * q[None, :]))
        labels.append(f"g{g}" + ("_down" if down.all() else ""))
    rows.append(lh)
    labels.append("plugin")
    return np.vstack(rows), labels, down


def candidates(n, lh, alpha, pool=None):
    """Return (lam_stack, cvecs, direction) for the CB candidate set."""
    lam_stack, _labels, down = candidate_matrix(n, lh, alpha=alpha, pool=pool)
    cs = batch_c(n, lam_stack, alpha=alpha, workers=10, pool=pool)
    D = lh.shape[0]
    m = lam_stack.shape[0] // D
    return lam_stack.reshape(m, D, lh.shape[1]), cs.reshape(m, D), down


def inner_world(n, lh, B, seed):
    """B variance draws in the world whose true ratios are lam_hat."""
    n = np.asarray(n, float)
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    D, K = lh.shape
    s0 = np.empty((D, B))
    si = np.empty((D, B, K))
    for d in range(D):
        s0[d] = rng.chisquare(nu[0], B) / nu[0]
        for i in range(K):
            si[d, :, i] = lh[d, i] * rng.chisquare(nu[i + 1], B) / nu[i + 1]
    return s0, si


def choose_ccb(n, lh, B, alpha, seed, pool=None, taus=None):
    """Return the chosen critical value per outer draw plus diagnostics."""
    n = np.asarray(n, float)
    D, K = lh.shape
    lam_stack, cs, down = candidates(n, lh, alpha, pool=pool)
    s0, si = inner_world(n, lh, B, seed)
    m = lam_stack.shape[0]
    est = np.empty((m, D))
    for d in range(D):
        cvec = np.repeat(cs[:, d], B)
        s0_t = np.tile(s0[d], m)
        si_t = np.tile(si[d], (m, 1))
        no = cond_noreject(cvec, s0_t, si_t, n, lh[d], npanels=24, order=8)
        est[:, d] = 1.0 - no.reshape(m, B).mean(axis=1)
    pick = np.empty(D, int)
    for d in range(D):
        below = np.flatnonzero(est[:, d] <= alpha)
        if below.size:
            j = below[np.argmax(est[below, d])]
        else:
            j = i64(np.argmin(est[:, d]))
        pick[d] = j
    c_out = cs[pick, np.arange(D)]
    return c_out, est, pick, down


def run(
    n, lam_true, alpha=ALPHA, D=100, seeds=(11, 22, 33), B=200, workers=10, pool=None
):
    n = np.asarray(n, float)
    lam_true = np.asarray(lam_true, float)
    c_or = oracle_c(n, lam_true, alpha)
    per = {
        "plugin": [],
        "cb_adapt": [],
        "cb0.9": [],
        "cb0.8": [],
        "cb0.7": [],
        "cbdir0.9": [],
        "cap0.02": [],
        "cap0.05": [],
        "cap0.10": [],
    }
    t0 = time.time()
    picks_hist = np.zeros(len(GAMMAS) + 1)
    for seed in seeds:
        s0, si, lh = draw_var(n, lam_true, D, seed)
        cplug = batch_c(n, lh, alpha=alpha, workers=workers, pool=pool)
        c_adapt, _est, pick, _down = choose_ccb(n, lh, B, alpha, seed + 913, pool=pool)
        for p in pick:
            picks_hist[p] += 1
        nu = n - 1.0
        q09 = np.array([f_quantile(nu[i + 1], nu[0], 0.9) for i in range(len(nu) - 1)])
        q08 = np.array([f_quantile(nu[i + 1], nu[0], 0.8) for i in range(len(nu) - 1)])
        c_cb9 = batch_c(n, lh / q09[None, :], alpha=alpha, workers=workers, pool=pool)
        c_cb8 = batch_c(n, lh / q08[None, :], alpha=alpha, workers=workers, pool=pool)
        c_cb7 = batch_c(
            n,
            lh
            / np.array([f_quantile(nu[i + 1], nu[0], 0.7) for i in range(len(nu) - 1)])[
                None, :
            ],
            alpha=alpha,
            workers=workers,
            pool=pool,
        )
        c_hi9 = batch_c(n, lh * q09[None, :], alpha=alpha, workers=workers, pool=pool)
        c_cons = np.maximum(c_cb9, c_hi9)
        capvecs = {
            name: np.minimum(c_cons, cplug + cap)
            for name, cap in (("cap0.02", 0.02), ("cap0.05", 0.05), ("cap0.10", 0.10))
        }
        nr_o = cond_noreject(np.full(D, c_or), s0, si, n, lam_true)
        for rule, cv in (
            ("plugin", cplug),
            ("cb_adapt", c_adapt),
            ("cb0.9", c_cb9),
            ("cb0.8", c_cb8),
            ("cb0.7", c_cb7),
            ("cbdir0.9", c_cons),
            *capvecs.items(),
        ):
            nr = cond_noreject(cv, s0, si, n, lam_true)
            d = nr_o - nr
            per[rule].append(
                (alpha + f64(d.mean()), f64(d.std(ddof=1) / np.sqrt(D)), f64(cv.mean()))
            )
    out = {
        "n": list(n),
        "lam_true": list(lam_true),
        "alpha": alpha,
        "D": D,
        "seeds": list(seeds),
        "B_inner": B,
        "c_oracle": c_or,
        "gamma_grid": list(GAMMAS),
        "rules": {},
        "pick_counts": [f64(x) for x in picks_hist],
    }
    for rule in per:
        rows = np.array(per[rule])
        f = rows[:, 0]
        out["rules"][rule] = {
            "fwer": f64(f.mean()),
            "seed_se_pct": f64(100 * f.std(ddof=1) / np.sqrt(len(f)) / alpha),
            "rel_err_pct": 100 * (f64(f.mean()) - alpha) / alpha,
            "mean_c": f64(rows[:, 2].mean()),
            "power_cost_pct": 100 * (f64(rows[:, 2].mean()) - c_or) / c_or,
        }
    out["seconds"] = time.time() - t0
    return out


def main(dfs=(3, 5, 7, 11, 19, 35)):
    import multiprocessing as mp

    with mp.Pool(i64(os.environ.get("W26_POOL", "3"))) as pool:
        for df in dfs:
            for fam, lam in FAMS.items():
                p = os.path.join(OUTDIR, f"adapt_df{df}_{fam}.json")
                if os.path.exists(p):
                    continue
                n = [df + 1, 13, 13]
                res = run(n, lam, pool=pool)
                res["family"] = fam
                res["benchmark_df"] = df
                dump_json(res, p)
                line = " ".join(
                    f"{k}:{v['rel_err_pct']:+.1f}%/{v['power_cost_pct']:+.2f}%"
                    for k, v in res["rules"].items()
                )
                print(f"df{df} {fam} {line} picks={res['pick_counts']}", flush=True)


if __name__ == "__main__":
    main()
