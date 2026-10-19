# -*- coding: utf-8 -*-
"""p_verify_t4k5.py — Track A, T4 at K = 5 (continuation of p_verify.py A).

The first run of this section (p_verify.py A) used a raw-data inner draw, ~2.4e10 normals per
cell, and was stopped for cost after the EQ0 structure; its EQ0 claims stand in
deposits/p_verify_A_claims.txt. This continuation uses the Bartlett inner draw, valid whenever
nu = n - 1 >= K (all n >= 6 here), with new seeds, and adds the replication triage of the one
anomalous first-run cell (EQ0, n = 12: 0.02960, z +4.2) at three times the outer sample.
Declared before the run: outer 20,000 (SE 1.1e-3), inner 20,000; replication outer 60,000
(SE 6.4e-4). T4 criterion (declared): the smallest n at which the rule is within 1.05 alpha.
"""
import math
import time

import numpy as np

import p_killround_k5 as K5
from p_verify import MIX

ALPHA = 0.025
rows = []


def claim(cid, label, value):
    rows.append("CLAIM %s: %s = %s" % (cid, label, value))
    print(rows[-1], flush=True)


def plugin(R, n, NO, M, seed):
    rng = np.random.default_rng(seed)
    rng2 = np.random.default_rng(seed + 1)
    L = np.linalg.cholesky(R)
    hits = 0
    for _ in range(NO):
        D = rng2.standard_normal((n, R.shape[0])) @ L.T
        Sh = np.cov(D, rowvar=False)
        sd = np.sqrt(np.diag(Sh))
        Rh = Sh / np.outer(sd, sd)
        tt = math.sqrt(n) * D.mean(0) / sd
        Zi, si = K5.draw_stats(Rh + 1e-9 * np.eye(R.shape[0]), n, M, rng)
        hits += int(tt.max() > np.quantile((Zi / si).max(1), 1 - ALPHA))
    return hits / NO


import sys
t0 = time.time()
S = dict(K5.structures(), MIX=MIX())
name = sys.argv[1]                       # one structure per process (run in parallel)
si = list(S).index(name)
R = S[name]
se = math.sqrt(ALPHA * (1 - ALPHA) / 20_000)
cells = []
for n in (6, 8, 12, 24, 36):
    p = plugin(R, n, 20_000, 20_000, 410000 + 100 * si + n)
    cells.append("n%d:%.5f(z%+.1f)" % (n, p, (p - ALPHA) / se))
    print(cells[-1], flush=True)
claim("PVA2-T4-K5-%s" % name, "plug-in realised FWER by n, Bartlett inner (SE 1.1e-3)", "; ".join(cells))
claim("PVA2-WALL-%s" % name, "seconds", "%.0f" % (time.time() - t0))
open("deposits/p_verify_A2_%s_claims.txt" % name, "w", encoding="utf-8").write(chr(10).join(rows) + chr(10))
