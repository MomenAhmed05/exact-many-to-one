"""H12: exact max-t with estimated scales and correlation (equicorrelated class).

Model: K strategies, T obs, returns X ~ N(0, R), R equicorrelated rho, unit var.
Statistic (the field's practice): t_i = sqrt(T)*mean_i / s_i with s_i the sample sd.
EXACT FORM (derived; Basu: sample means independent of the sample scales):
  given the K sample scales s=(s_1..s_K), t | s has the same correlation
  structure as X/s: Corr(t_i,t_j) = rho for ALL i!=j (scales cancel), with
  marginal scale 1/s_i. Hence
      P(max t <= c) = E_s [ E_z [ PROD_i Phi((c*s_i - sqrt(rho)*z)/sqrt(1-rho)) ] ]
  (z~N(0,1); the expectation over s uses the ACTUAL joint law of the sample
  scales, obtained by simulating returns -> no independence assumption).
  This is the "correlated chi-square mixtures" object the 2023 literature
  simulates. Dunnett's known-scale/known-rho critical value is the same with
  s_i = 1 and known rho; the gap is the estimation penalty.
Plug-in with estimated rho (Basu: r_hat independent of the means): c_plug(r_hat)
  solves the known-scale formula at r_hat; realized size measured under the truth.
Cells: K in {5,10} x T in {60,250,1000} x rho in {0,0.5}; 3 seeds.
Falsifier: MC-certified size above nominal (beyond 3 MCSE) for c_exact.
"""
import json
import math
import numpy as np
from math import erf, sqrt

OUT = 'tools/audit/series/h12_dunnett.json'
SEEDS = (67110, 67111, 67112)
ALPHA = 0.05


def Phi(x):
    return 0.5 * (1 + np.vectorize(lambda z: erf(z))(np.asarray(x, float) / sqrt(2)))


ZG = np.linspace(-9, 9, 4001)
WZ = np.exp(-0.5 * ZG ** 2)
WZ = WZ / WZ.sum() / (ZG[1] - ZG[0]) * (ZG[1] - ZG[0])  # normalized discrete weights


def cond_tail_vec(c, scales, rho):
    """E_s[P(max<=c|s)] for a PRECOMPUTED (R_s,K) scale array, vectorized."""
    rho = min(max(rho, 0.0), 0.9999)
    a = sqrt(1 - rho)
    cs = c * scales[:, :, None]                      # (Rs,K,1)
    arg = (cs - sqrt(rho) * ZG[None, None, :]) / a   # (Rs,K,G)
    prod = np.prod(Phi(arg), axis=1)                 # (Rs,G)
    w = np.exp(-0.5 * ZG ** 2)
    w = w / w.sum()
    return float((prod * w).sum(axis=1).mean())


def solve_c_exact(rho, K, T, rng, R_s=300, alpha=ALPHA):
    """Exact CV: invert E_s[cond_tail(c|s)] = 1-alpha on precomputed scales."""
    scales = np.array([X.std(axis=1, ddof=1) for X in (sim(K, T, rho, rng) for _ in range(R_s))])
    lo, hi = 0.0, 12.0
    for _ in range(24):
        mid = (lo + hi) / 2
        if cond_tail_vec(mid, scales, rho) < 1 - alpha:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def dunnett_c(rho, K, alpha=ALPHA, R=400000):
    """Known-scale, known-rho critical value."""
    rng = np.random.default_rng(4242)
    z = rng.standard_normal(R)
    a = sqrt(1 - min(max(rho, 0.0), 0.9999))
    def tail(c):
        return float(np.mean(Phi((c - sqrt(rho) * z) / a) ** K))
    lo, hi = 0.0, 12.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if tail(mid) < 1 - alpha:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def sim(K, T, rho, rng):
    z = rng.standard_normal((K, T))
    c = rng.standard_normal((1, T))
    return sqrt(rho) * c + sqrt(1 - rho) * z


def rhat(X):
    C = np.corrcoef(X)
    K = C.shape[0]
    return float(np.mean(C[np.triu_indices(K, 1)]))


_CP_TAB = {}


def _cp_table(K, alpha=ALPHA, R=200000):
    if K in _CP_TAB:
        return _CP_TAB[K]
    rng = np.random.default_rng(8484)
    z = rng.standard_normal(R)
    rhos = np.linspace(0.0, 0.9, 121)
    cs = np.empty_like(rhos)
    for j, rh in enumerate(rhos):
        a = sqrt(1 - rh)
        def tail(c):
            return float(np.mean(Phi((c - sqrt(rh) * z) / a) ** K))
        lo, hi = 0.0, 12.0
        for _ in range(34):
            mid = (lo + hi) / 2
            if tail(mid) < 1 - alpha:
                lo = mid
            else:
                hi = mid
        cs[j] = (lo + hi) / 2
    _CP_TAB[K] = (rhos, cs)
    return _CP_TAB[K]


def c_plug(rho_hat, K, alpha=ALPHA):
    rhos, cs = _cp_table(K, alpha)
    return float(np.interp(max(rho_hat, 0.0), rhos, cs))


def run():
    res = {'alpha': ALPHA, 'cells': {}, 'falsified': []}
    R_cert = 4000
    for K in (5, 10):
        for T in (60, 250, 1000):
            for rho in (0.0, 0.5):
                key = f'K{K}/T{T}/rho{rho}'
                c_dun = dunnett_c(rho, K)
                c_ex = np.mean([solve_c_exact(rho, K, T, np.random.default_rng(s * 31 + K + T)) for s in SEEDS])
                # MC certificate: simulate returns, max-t at both CVs
                size_dun = size_ex = size_plug = 0.0
                plugs = []
                for s in SEEDS:
                    rng = np.random.default_rng(s * 777 + K * 13 + T + int(rho * 10))
                    n = R_cert // len(SEEDS)
                    rej_d = rej_e = rej_p = 0
                    for _ in range(n):
                        X = sim(K, T, rho, rng)
                        t = sqrt(T) * X.mean(axis=1) / X.std(axis=1, ddof=1)
                        m = float(np.max(t))
                        rej_d += m > c_dun
                        rej_e += m > c_ex
                        cp = c_plug(rhat(X), K)
                        plugs.append(cp)
                        rej_p += m > cp
                    size_dun += rej_d / n
                    size_ex += rej_e / n
                    size_plug += rej_p / n
                size_dun /= len(SEEDS)
                size_ex /= len(SEEDS)
                size_plug /= len(SEEDS)
                mcse = math.sqrt(ALPHA * (1 - ALPHA) / R_cert)
                res['cells'][key] = {
                    'c_dunnett': round(c_dun, 4), 'c_exact': round(float(c_ex), 4),
                    'penalty_pct': round((c_ex - c_dun) / c_dun * 100, 3),
                    'c_plug_mean': round(float(np.mean(plugs)), 4), 'c_plug_sd': round(float(np.std(plugs, ddof=1)), 4),
                    'size_dunnett': round(size_dun, 4), 'size_exact': round(size_ex, 4),
                    'size_plug': round(size_plug, 4), 'mcse': round(mcse, 4)}
                row = res['cells'][key]
                if row['size_exact'] > ALPHA + 3 * mcse:
                    res['falsified'].append(f'{key} size_exact={row["size_exact"]}')
                print(f"{key}: c_dun={row['c_dunnett']} c_exact={row['c_exact']} "
                      f"(penalty {row['penalty_pct']}%) c_plug={row['c_plug_mean']}±{row['c_plug_sd']} | "
                      f"size dun={row['size_dunnett']} exact={row['size_exact']} plug={row['size_plug']}", flush=True)
    print('falsified:', res['falsified'], flush=True)
    with open(OUT, 'w') as f:
        json.dump(res, f, indent=1)
    print('wrote', OUT, flush=True)


if __name__ == '__main__':
    run()

