"""W5 best-rule slice: the W2 CB arms on the crossing region.

Computes plugin + cb0.9 + cb0.7 on a reduced grid (K=2, shapes
bench_small/short_cand/balanced, vols noisy_cand/noisy_bench/graded,
m in {4,6,8}, both alphas) with 5 seeds x 60 draws, common draws across arms.
The full-grid plugin map is in the other W5 deposits; this slice supplies the
"same map for the best rule" comparison where the crossing happens.

Outputs: deposits/w5best_<K>_<shape>_<vol>_a<alpha>.json
"""
import multiprocessing as mp
import os
import sys
import time
from math import sqrt

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import (batch_c, cond_noreject, draw_var, f_quantile,  # noqa
                  oracle_c)
from io_util import dump_json, f64  # noqa: E402
from w5_boundary import lamvec, sizes  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'deposits')
MS = (4, 6, 8)
SEEDS = (101, 202, 303, 404, 505)
D = 60
SHAPES = ('bench_small', 'short_cand')
VOLS = ('noisy_cand', 'noisy_bench', 'graded')
ALPHAS = (0.025, 0.05)


def run_cell(K, shape, vol, alpha, pool):
    lam = np.asarray(lamvec(vol, K), float)
    rows = []
    t0 = time.time()
    for m in MS:
        n = np.asarray(sizes(shape, m, K), float)
        nu = n - 1.0
        q09 = np.array([f_quantile(nu[i + 1], nu[0], 0.9)
                        for i in range(K)])
        q07 = np.array([f_quantile(nu[i + 1], nu[0], 0.7)
                        for i in range(K)])
        c_or = oracle_c(n, lam, alpha)
        draws = [draw_var(n, lam, D, s) for s in SEEDS]
        lh_all = np.vstack([d[2] for d in draws])
        solved = batch_c(n, np.vstack([lh_all, lh_all / q09[None, :],
                                       lh_all / q07[None, :]]),
                         alpha=alpha, workers=3, pool=pool).reshape(3, -1)
        seeds_out = []
        pos = 0
        for (s0, si, _lh) in draws:
            dd = len(s0)
            cplug = solved[0, pos:pos + dd]
            c09 = solved[1, pos:pos + dd]
            c07 = solved[2, pos:pos + dd]
            pos += dd
            nr_o = cond_noreject(np.full(dd, c_or), s0, si, n, lam)
            rec = {}
            for rule, cv in (('plugin', cplug), ('cb0.9', c09),
                             ('cb0.7', c07)):
                nr = cond_noreject(cv, s0, si, n, lam)
                d = nr_o - nr
                rec[rule] = {
                    'fwer': f64(alpha + d.mean()),
                    'mcse': f64(d.std(ddof=1) / sqrt(dd)),
                    'mean_c': f64(cv.mean())}
            seeds_out.append(rec)
        row = {'m': m, 'df': m - 1, 'n': list(n), 'c_oracle': c_or}
        for rule in ('plugin', 'cb0.9', 'cb0.7'):
            f = np.array([s[rule]['fwer'] for s in seeds_out])
            row[rule] = {
                'fwer': f64(f.mean()),
                'rel_err_pct': 100 * (f64(f.mean()) - alpha) / alpha,
                'seed_se': f64(f.std(ddof=1) / sqrt(len(f))),
                'mean_c': f64(np.mean([s[rule]['mean_c'] for s in seeds_out])),
                'power_cost_pct': 100 * (f64(np.mean(
                    [s[rule]['mean_c'] for s in seeds_out])) - c_or) / c_or}
        rows.append(row)
    return {'K': K, 'shape': shape, 'vol': vol, 'alpha': alpha,
            'lam': list(lam), 'rows': rows, 'seconds': time.time() - t0}


def _main(pool):
    for K in (2,):
        for shape in SHAPES:
            for vol in VOLS:
                for alpha in ALPHAS:
                    p = os.path.join(OUTDIR,
                                     f'w5best_K{K}_{shape}_{vol}_a{alpha}.json')
                    if os.path.exists(p):
                        continue
                    res = run_cell(K, shape, vol, alpha, pool)
                    dump_json(res, p)
                    line = ' '.join(
                        f"{r['df']}:{r['plugin']['rel_err_pct']:+.1f}/"
                        f"{r['cb0.9']['rel_err_pct']:+.1f}/"
                        f"{r['cb0.7']['rel_err_pct']:+.1f}"
                        for r in res['rows'])
                    print(f"K{K} {shape} {vol} a{alpha} {line}", flush=True)


if __name__ == '__main__':
    with mp.Pool(2) as pool:
        _main(pool)
