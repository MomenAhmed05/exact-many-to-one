# -*- coding: utf-8 -*-
"""Independent check of the aligned-factor corner, written from the model alone (Section 6.2).

r_{g,t} = beta_g f_t + eps_{g,t},  f_t ~ N(0,1),  Var eps_0 = 1,  Var eps_1 = lam,
n = 36 per group, Welch statistic T = (xbar_1 - xbar_0)/sqrt(s1^2/n + s0^2/n),
one-sided rejection T > c. Shares no code with the other scripts in this folder.

Writes deposits/corner_check_20m_claims.txt.  About 35 s at 20,000,000 draws.
"""
import pathlib
import sys
import time

import numpy as np

LAM, B0, B1 = 0.04, -1.5, 1.5
N = 36
CS = {'published 11.6 constant': 2.7366485595703125,
      'new certified U': 13587 / 5000,
      'new certified L': 271731 / 100000}
DRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 20_000_000
CHUNK = 250_000

rng = np.random.default_rng(987654321)
hits = {k: 0 for k in CS}
keep = []                       # the upper tail, to read off the empirical quantile
t0 = time.time()
done = 0
while done < DRAWS:
    m = min(CHUNK, DRAWS - done)
    f = rng.standard_normal((m, N))
    e0 = rng.standard_normal((m, N))
    e1 = rng.standard_normal((m, N)) * np.sqrt(LAM)
    r0 = B0 * f + e0
    r1 = B1 * f + e1
    t = (r1.mean(1) - r0.mean(1)) / np.sqrt(r1.var(1, ddof=1) / N + r0.var(1, ddof=1) / N)
    for k, c in CS.items():
        hits[k] += int(np.count_nonzero(t > c))
    keep.append(t[t > 2.5])
    done += m

tail = np.sort(np.concatenate(keep))[::-1]
k = int(round(0.025 * DRAWS))
cstar = tail[k - 1]
rows = []
for name, c in CS.items():
    p = hits[name] / DRAWS
    se = (p * (1 - p) / DRAWS) ** 0.5
    line = '%-26s c=%.7f  P(T>c)=%.6f  se=%.1e  (alpha - p)/se = %+.1f' % (name, c, p, se, (0.025 - p) / se)
    print(line)
    rows.append('CLAIM CC20-%s: c, P(T > c), its standard error and (alpha - P)/se = %.7f, %.6f, %.1e, %+.1f'
                % (name.replace(' ', '-'), c, p, se, (0.025 - p) / se))
# a rough CI for the empirical 97.5% quantile from order statistics
sd = (0.025 * 0.975 * DRAWS) ** 0.5
lo, hi = tail[int(k + 2 * sd) - 1], tail[int(k - 2 * sd) - 1]
print('empirical critical constant c* = %.5f   (approx 95%% CI %.5f to %.5f)' % (cstar, lo, hi))
print('draws %d, %.0f s' % (DRAWS, time.time() - t0))
rows.append('CLAIM CC20-CSTAR: empirical critical constant and approximate 95%% interval = %.5f [%.5f, %.5f]'
            % (cstar, lo, hi))
rows.append('CLAIM CC20-DRAWS: draws = %d' % DRAWS)
rows.append('CLAIM CC20-SEED: seed = 987654321')
out = pathlib.Path(__file__).resolve().parent / 'deposits' / 'corner_check_20m_claims.txt'
out.write_text('\n'.join(rows) + '\n', encoding='utf-8')
