"""Micro-timing: per-row cost of c_env at K=2 vs K=5 vs K=8 (for run-scope decisions)."""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import draw_var  # noqa: E402
from env_rule import c_env  # noqa: E402

for K, maxv in ((2, 32), (5, 32), (8, 19)):
    n = np.array([4.0] + [13.0] * K, float)
    lam = np.array([1.0, 25.0] + [4.0] * (K - 2), float)  # candidates only (sigma_0=1)
    s0, si, lh = draw_var(n, lam, 1, 7)
    t0 = time.time()
    try:
        c = c_env(n, s0, si, lh, max_vertices=maxv)
        dt = time.time() - t0
        print(
            f"K={K} MAXV={maxv}: 1 row, {dt:.2f}s -> {dt:.2f}s/row "
            f"(c={float(np.asarray(c).ravel()[0]):.4f})",
            flush=True,
        )
    except Exception as exc:  # noqa: BLE001
        print(
            f"K={K} MAXV={maxv}: FAILED {exc!r} after {time.time() - t0:.1f}s",
            flush=True,
        )
