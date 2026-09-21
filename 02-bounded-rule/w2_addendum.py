"""W2 addendum: per-component MIXED-direction CB,
evaluated on Tamhane and Xi's own fourteen scenarios.

The global direction map shrinks every ratio the same way.  Here the sign of
dc/dlam_i is estimated per component by finite differences (lam_i * q versus
lam_i / q, others fixed) and each component is bounded in its own conservative
direction, then all components are shrunk by the F-quantile of gamma.

Rules per scenario: plugin, cb0.9, mixcb0.8, mixcb0.9, cap0.10, m3.
Paired against the oracle at the true ratios; 3 seeds x 100 draws.

Outputs: deposits/w2/addendum_<i>.json
"""
import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import (ALPHA, batch_c, cond_noreject, draw_var, f_quantile,  # noqa
                  oracle_c)
from h12_tamhane_headtohead import SCEN_K2, SCEN_K3, c_method3  # noqa: E402
from io_util import dump_json, f64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'deposits', 'w2')
SEEDS = (11, 22, 33)
D = 100
RULES = ['plugin', 'cb0.9', 'mixcb0.8', 'mixcb0.9', 'cap0.10', 'm3']


def _mix_rule(n, lh, gamma, alpha, pool):
    """Per-component mixed-direction CB at confidence gamma."""
    n = np.asarray(n, float)
    nu = n - 1.0
    K = lh.shape[1]
    q = np.array([f_quantile(nu[i + 1], nu[0], gamma) for i in range(K)])
    Dn = lh.shape[0]
    rows = [lh]
    for i in range(K):
        up = lh.copy()
        up[:, i] = lh[:, i] * q[i]
        dn = lh.copy()
        dn[:, i] = lh[:, i] / q[i]
        rows += [up, dn]
    solved = batch_c(n, np.vstack(rows), alpha=alpha, workers=3,
                     pool=pool).reshape(len(rows), Dn)
    base = solved[0]
    use_up = np.zeros((Dn, K), dtype=bool)
    for i in range(K):
        c_up = solved[1 + 2 * i]
        c_dn = solved[2 + 2 * i]
        use_up[:, i] = c_up > c_dn
    mix = lh.copy()
    for i in range(K):
        mix[:, i] = np.where(use_up[:, i], lh[:, i] * q[i], lh[:, i] / q[i])
    c_mix = batch_c(n, mix, alpha=alpha, workers=3, pool=pool)
    return base, c_mix, mix


def run_scenario(idx, label, n_l, sig_l, pool):
    n = np.asarray(n_l, float)
    sig = np.asarray(sig_l, float)
    lam = np.asarray((sig / sig[0]) ** 2, float)[1:]
    c_or = oracle_c(n, lam, ALPHA)
    t0 = time.time()
    per = {r: [] for r in RULES}
    for seed in SEEDS:
        s0, si, lh = draw_var(n, lam, D, seed)
        cplug, c_mix8, _m8 = _mix_rule(n, lh, 0.8, ALPHA, pool)
        _, c_mix9, _m9 = _mix_rule(n, lh, 0.9, ALPHA, pool)
        q09 = np.array([f_quantile(n.size - 1 if False else n[i + 1] - 1,
                                   n[0] - 1, 0.9)
                        for i in range(len(n) - 1)])
        c09 = batch_c(n, lh / q09[None, :], alpha=ALPHA, workers=3, pool=pool)
        chead = np.vstack([lh, lh / q09[None, :], lh * q09[None, :]])
        c_h = batch_c(n, chead, alpha=ALPHA, workers=3, pool=pool).reshape(3, D)
        ccons = np.maximum(c_h[1], c_h[2])
        ccap = np.minimum(ccons, cplug + 0.10)
        c_m3 = _m3_vec(n, s0, si, seed)
        cvecs = {'plugin': cplug, 'cb0.9': c09, 'mixcb0.8': c_mix8,
                 'mixcb0.9': c_mix9, 'cap0.10': ccap, 'm3': c_m3}
        nr_o = cond_noreject(np.full(D, c_or), s0, si, n, lam)
        for rule, cv in cvecs.items():
            ok = np.isfinite(cv)
            nr = cond_noreject(cv[ok], s0[ok], si[ok], n, lam)
            d = nr_o[ok] - nr
            per[rule].append((ALPHA + f64(d.mean()),
                              f64(d.std(ddof=1) / np.sqrt(len(d))),
                              f64(cv[ok].mean())))
    out = {'idx': idx, 'scenario': label, 'n': n_l, 'sigma': sig_l,
           'lam': [f64(x) for x in lam], 'alpha': ALPHA, 'D': D,
           'seeds': list(SEEDS), 'c_oracle': c_or, 'rules': {},
           'seconds': 0.0}
    for rule in RULES:
        rows = np.array(per[rule])
        f = rows[:, 0]
        out['rules'][rule] = {
            'fwer': f64(f.mean()),
            'seed_se_pct': f64(100 * f.std(ddof=1) / np.sqrt(len(f)) / ALPHA),
            'rel_err_pct': 100 * (f64(f.mean()) - ALPHA) / ALPHA,
            'mean_c': f64(rows[:, 2].mean()),
            'power_cost_pct': 100 * (f64(rows[:, 2].mean()) - c_or) / c_or}
    out['seconds'] = time.time() - t0
    return out


def _m3_vec(n, s0, si, seed):
    rng = np.random.default_rng(seed * 104729)
    vals = []
    for d in range(len(s0)):
        s2 = np.concatenate([[s0[d]], si[d]])
        vals.append(c_method3(n, s2, rng, ALPHA, sims=20000))
    return np.asarray(vals, float)


def _main(pool):
    scen = [('k2 ' + a, b, c) for a, b, c in SCEN_K2] + \
           [('k3 ' + a, b, c) for a, b, c in SCEN_K3]
    for idx, (label, n_l, sig_l) in enumerate(scen):
        p = os.path.join(OUTDIR, f'addendum_{idx:02d}.json')
        if os.path.exists(p):
            continue
        res = run_scenario(idx, label, n_l, sig_l, pool)
        dump_json(res, p)
        line = ' '.join(f"{k}:{v['rel_err_pct']:+.1f}%"
                        for k, v in res['rules'].items())
        print(f'{idx:02d} {label} {line}', flush=True)


if __name__ == '__main__':
    with mp.Pool(3) as pool:
        _main(pool)
