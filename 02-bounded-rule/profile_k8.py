"""Diagnostic: WHERE does the K=5/K=8 critical-value cost actually go?

Distinguishes two very different diagnoses:

  interpreter-bound  the inner work is small per call and Python overhead
                     dominates -> a compiled language (or batching) buys 10-100x
  compute-bound      the inner work is already numpy-vectorised and large ->
                     a compiled language buys ~1x, and the only levers are the
                     quadrature size, batching across rows, or fewer root steps

Achieved element-throughput is the tell: numpy/erf-based kernels reach ~1e8-1e9
element-ops/s, an interpreter-bound loop sits near 1e6.

Outputs: deposits/cert/profile_k8.json
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import GRIDS  # noqa: E402
from h12_hetero_exact import HeteroMaxT  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "cert")
N_REP = 20


def _f(x, default=0.0):
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return float(default)


def probe(K, grid_name, alpha=0.025):
    grid = GRIDS[grid_name]
    n = np.array([4.0] + [13.0] * K, float)
    lam = np.array([1.0, 25.0] + [4.0] * (K - 2), float)
    sig = np.sqrt([1.0] + list(lam))
    t0 = time.perf_counter()
    obj = HeteroMaxT(n, sig, grid=grid)
    t_init = time.perf_counter() - t0
    # cdf cost
    c_probe = 3.0
    val = float("nan")
    t0 = time.perf_counter()
    for _ in range(N_REP):
        val = obj.cdf(c_probe)
    t_cdf = (time.perf_counter() - t0) / N_REP
    # quadrature sizes
    m_v = _f(obj.v.size, 0)
    m_u = _f(obj._u_nodes(c_probe)[0].shape[1], 0)
    m_s = _f(np.asarray(obj.ws[0]).size, 0)
    elems = m_v * m_u * m_s * K  # elementwise work per cdf: (v,u,s) per treatment
    t0 = time.perf_counter()
    c = obj.critical_value(alpha)
    t_root = time.perf_counter() - t0
    return {
        "K": i64(K),
        "grid": grid_name,
        "t_init_s": f64(t_init),
        "t_cdf_s": f64(t_cdf),
        "t_critical_value_s": f64(t_root),
        "implied_cdf_calls": f64(t_root / t_cdf if t_cdf > 0 else float("nan")),
        "m_v": i64(m_v),
        "m_u": i64(m_u),
        "m_s": i64(m_s),
        "elem_per_cdf": f64(elems),
        "elem_per_s": f64(elems / t_cdf if t_cdf > 0 else float("nan")),
        "c_value": f64(c),
        "cdf_at_probe": f64(val),
    }


def main():
    rows = []
    for K in (2, 5, 8):
        for grid_name in ("fast", "boot"):
            try:
                r = probe(K, grid_name)
                rows.append(r)
                print(
                    f"K={r['K']} grid={r['grid']}: init {r['t_init_s']:.3f}s | "
                    f"cdf {r['t_cdf_s'] * 1000:.2f}ms | root "
                    f"{r['t_critical_value_s']:.2f}s ({r['implied_cdf_calls']:.0f} "
                    f"cdf calls) | grid {r['m_v']}x{r['m_s']}x{r['m_u']} = "
                    f"{r['elem_per_cdf']:.2e} elem/cdf | achieved "
                    f"{r['elem_per_s']:.2e} elem/s | c={r['c_value']:.4f}",
                    flush=True,
                )
            except Exception as exc:  # noqa: BLE001
                print(f"K={K} grid={grid_name}: FAILED {exc!r}", flush=True)
    dump_json(
        {
            "rows": rows,
            "n_rep": i64(N_REP),
            "note": "achieved element-throughput ~1e8-1e9 => compute-bound "
            "(numpy); ~1e6 => interpreter-bound (batching/compiled wins)",
        },
        os.path.join(OUT, "profile_k8.json"),
    )
    print("wrote deposits/cert/profile_k8.json")


if __name__ == "__main__":
    main()
