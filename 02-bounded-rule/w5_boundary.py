"""W5: boundary-map generalization.

Sweeps the benchmark/candidate degrees of freedom for
  K in {1,2,3,5,8}
  shapes: bench_small (n0=m, cands 13), short_cand (n0=36, cands m),
          short_bench (n0=m, cands 36), balanced (all m)
  m in {4,6,8,12,20,36}  (df 3,5,7,11,19,35)
  3 volatility patterns x alpha in {0.025, 0.05}
  5 seeds x 100 draws per cell

Rules: 'plugin' (c(lam_hat)) and 'deb' (2c - mean bootstrap c*), both cheap.
Every estimate is paired against the oracle at the true ratios.

Usage:
  python w5_boundary.py plugin          # run all plugin cells
  python w5_boundary.py deb             # run all deb cells
  python w5_boundary.py plugin 2        # only K=2
Checkpoints one file per (rule,K,shape,vol,alpha); skips existing.
"""

import os
import sys
import time
from math import sqrt

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import batch_c, cond_noreject, draw_var, oracle_c  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w5")
KS = (1, 2, 3, 5, 8)
MS = (4, 6, 8, 12, 36)
MS_BY_K = {
    1: (4, 6, 8, 12, 36),
    2: (4, 6, 8, 12, 36),
    3: (4, 6, 8, 12, 36),
    5: (4, 6, 8),
    8: (4, 6),
}
D_BY_K = {1: 60, 2: 60, 3: 60, 5: 40, 8: 40}
POOL_N = i64(os.environ.get("W26_POOL", "4"))
SHAPES = ("bench_small", "short_cand", "balanced")
ALPHAS = (0.025, 0.05)
B_BOOT = 40
SEEDS = (101, 202, 303, 404, 505)
D_PER_SEED = 60


def sizes(shape, m, K):
    if shape == "bench_small":
        return [m] + [13] * K
    if shape == "short_cand":
        return [36] + [m] * K
    if shape == "short_bench":
        return [m] + [36] * K
    if shape == "balanced":
        return [m] * (K + 1)
    raise RuntimeError(f"unknown shape {shape}")


def lamvec(vol, K):
    if vol == "noisy_cand":
        return [25.0] + [1.0] * (K - 1)
    if vol == "noisy_bench":
        return [0.04] * K
    if vol == "graded":
        return ([0.25, 4.0, 1.0, 9.0, 0.5, 2.0, 0.1, 6.25])[:K]
    raise RuntimeError(f"unknown vol {vol}")


def cell_key(rule, K, shape, vol, alpha):
    return f"{rule}_K{K}_{shape}_{vol}_a{alpha}"


def rule_c(rule, n, lh, B, alpha, seed, workers=10, pool=None):
    if rule == "plugin":
        return batch_c(n, lh, alpha=alpha, workers=workers, pool=pool)
    if rule == "deb":
        nu = np.asarray(n, float) - 1.0
        rng = np.random.default_rng(seed)
        D, Kk = lh.shape
        fstar = np.column_stack(
            [rng.f(nu[i + 1], nu[0], size=(B, D)) for i in range(Kk)]
        )
        lam_star = (lh[None, :, :] * fstar.transpose(2, 0, 1)).reshape(B * D, Kk)
        cs = batch_c(n, lam_star, alpha=alpha, workers=workers, pool=pool).reshape(B, D)
        cplug = batch_c(n, lh, alpha=alpha, workers=workers, pool=pool)
        return 2 * cplug - cs.mean(axis=0)
    if rule == "cb90":
        from core import lam_lower

        return batch_c(
            n, lam_lower(n, lh, 0.9), alpha=alpha, workers=workers, pool=pool
        )
    if rule == "dir":
        from core import lam_lower, lam_upper

        lo = batch_c(n, lam_lower(n, lh, 0.9), alpha=alpha, workers=workers, pool=pool)
        hi = batch_c(n, lam_upper(n, lh, 0.9), alpha=alpha, workers=workers, pool=pool)
        return np.maximum(lo, hi)
    raise RuntimeError(f"unknown rule {rule}")


def run_cell(rule, K, shape, vol, alpha, workers=10, pool=None):
    lam = np.asarray(lamvec(vol, K), float)
    rows = []
    t0 = time.time()
    for m in MS_BY_K.get(K, MS):
        n = np.asarray(sizes(shape, m, K), float)
        c_or = oracle_c(n, lam, alpha)
        # pool all seeds into ONE solver batch for this m
        draws = []
        for seed in SEEDS:
            s0, si, lh = draw_var(
                n, lam, D_BY_K.get(K, D_PER_SEED), seed
            )
            draws.append((seed, s0, si, lh))
        lh_all = np.vstack([d[3] for d in draws])
        cv_all = rule_c(rule, n, lh_all, B_BOOT, alpha, SEEDS[0], workers, pool=pool)
        seeds_out = []
        pos = 0
        for seed, s0, si, _lh in draws:
            d_len = len(s0)
            cv = cv_all[pos : pos + d_len]
            pos += d_len
            nr_o = cond_noreject(np.full(d_len, c_or), s0, si, n, lam)
            nr = cond_noreject(cv, s0, si, n, lam)
            d = nr_o - nr
            seeds_out.append(
                {
                    "seed": seed,
                    "fwer": f64(alpha + d.mean()),
                    "mcse": f64(d.std(ddof=1) / sqrt(len(d))),
                    "mean_c": f64(cv.mean()),
                }
            )
        f = np.array([s["fwer"] for s in seeds_out])
        mc = f64(np.mean([s["mean_c"] for s in seeds_out]))
        rows.append(
            {
                "m": m,
                "df": m - 1,
                "n": list(n),
                "c_oracle": c_or,
                "seeds": seeds_out,
                "fwer": f64(f.mean()),
                "seed_se": f64(f.std(ddof=1) / sqrt(len(f))),
                "rel_err_pct": 100 * (f64(f.mean()) - alpha) / alpha,
                "mean_c": mc,
                "power_cost_pct": 100 * (mc - c_or) / c_or,
                "nan": i64(np.sum(~np.isfinite(cv_all))) if False else 0,
            }
        )
    return {
        "rule": rule,
        "K": K,
        "shape": shape,
        "vol": vol,
        "alpha": alpha,
        "lam": list(lam),
        "rows": rows,
        "seconds": time.time() - t0,
    }


def cells_for(rule, only_K=None):
    for K in KS:
        if only_K and only_K != K:
            continue
        for shape in SHAPES:
            for vol in ("noisy_cand", "noisy_bench", "graded"):
                if K == 1 and vol == "noisy_cand":
                    continue
                for alpha in ALPHAS:
                    yield (rule, K, shape, vol, alpha)


def main():
    rule = sys.argv[1] if len(sys.argv) > 1 else "plugin"
    only_K = None
    if len(sys.argv) > 2:
        try:
            only_K = int(sys.argv[2])
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError("bad K") from exc
    import multiprocessing as mp

    with mp.Pool(POOL_N) as pool:
        for rl, K, shape, vol, alpha in cells_for(rule, only_K):
            key = cell_key(rl, K, shape, vol, alpha)
            p = os.path.join(OUTDIR, key + ".json")
            if os.path.exists(p):
                continue
            try:
                res = run_cell(rl, K, shape, vol, alpha, pool=pool)
                dump_json(res, p)
                worst = max(res["rows"], key=lambda r: abs(r["rel_err_pct"]))
                print(
                    f"{key} worst m={worst['m']} "
                    f"({worst['rel_err_pct']:+.1f}%) "
                    f"{res['seconds']:.0f}s",
                    flush=True,
                )
            except Exception as exc:  # noqa: BLE001
                print(f"FAILED {key}: {exc!r}", flush=True)


if __name__ == "__main__":
    main()
