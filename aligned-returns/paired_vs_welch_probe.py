# -*- coding: utf-8 -*-
"""paired_vs_welch_probe.py — the paired test against the Welch rule on aligned dates (Section 7).

On aligned one-factor data r0_t = b0 f_t + z0_t, r1_t = b1 f_t + sqrt(lam) z1_t + delta,
the daily differences d_t = r1_t - r0_t are iid N(delta, (b1-b0)^2 + 1 + lam), so the
paired t statistic sqrt(n) dbar / s_d is EXACTLY t_{n-1} under delta = 0, at every theta.

Declared before running:
  size: paired t at the t_35 upper 0.025 quantile must have rate 0.025 within 4 SE at
        every theta (sanity anchor for the exactness claim);
  power: at the delta where the paired test has about 80% power, report the power of
         2B's certified Welch rule (constant Uc = 2.7174, the box-wide constant) and of
         the Welch statistic at its own known-theta constant c*(theta) (the oracle, from
         the 2B engine's certified bracket), on the SAME draws.
N = 2,000,000 per theta (SE about 1.1e-4 on a size of 0.025).
"""
import numpy as np
from scipy.stats import t as tdist

n, N = 36, 2_000_000
UC = 2.7174
cpaired = tdist.isf(0.025, n - 1)
# certified c*(theta) brackets from the 2B engine (_cstar_spread_probe.py)
thetas = {'same sign (1,1,1)': ((1.0, 1.0, 1.0), 1.416272),
          'no loadings (1,0,0)': ((1.0, 0.0, 0.0), 1.994437),
          'opposite (1,1,-1)': ((1.0, 1.0, -1.0), 2.453056),
          'corner (0.04,-1.5,1.5)': ((0.04, -1.5, 1.5), 2.717340)}
rng = np.random.default_rng(20261006)
print('paired critical value t_35 upper 0.025 = %.6f' % cpaired)
for name, ((lam, b0, b1), cstar) in thetas.items():
    vd = (b1 - b0) ** 2 + 1 + lam
    delta80 = (cpaired + 0.8416) * np.sqrt(vd / n)      # approx 80% paired power
    res = {}
    for d_lab, delta in (('null', 0.0), ('alt', delta80)):
        hits = np.zeros(3, dtype=np.int64)
        done = 0
        while done < N:
            m = min(200_000, N - done)
            f = rng.standard_normal((m, n))
            r0 = b0 * f + rng.standard_normal((m, n))
            r1 = b1 * f + np.sqrt(lam) * rng.standard_normal((m, n)) + delta
            d = r1 - r0
            tp = np.sqrt(n) * d.mean(1) / d.std(1, ddof=1)
            tw = (r1.mean(1) - r0.mean(1)) / np.sqrt(r1.var(1, ddof=1) / n + r0.var(1, ddof=1) / n)
            hits += [(tp > cpaired).sum(), (tw > UC).sum(), (tw > cstar).sum()]
            done += m
        res[d_lab] = hits / N
    se = np.sqrt(0.025 * 0.975 / N)
    print('CLAIM PVW-%s: size paired / Welch@Uc / Welch@c*(theta) = %.5f (z %+.1f) / %.5f / %.5f ;'
          ' power at paired-80%% effect = %.3f / %.3f / %.3f'
          % (name, res['null'][0], (res['null'][0] - 0.025) / se, res['null'][1], res['null'][2],
             res['alt'][0], res['alt'][1], res['alt'][2]))
