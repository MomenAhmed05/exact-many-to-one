"""Vectorised exact-law solver for the bulk W5 sweep.

Same reduction as core.HeteroMaxT (probability-scale composite Gauss-Legendre,
planted edges unused here because W5 designs stay away from the F5 corner), but
solves a whole batch of lambda vectors at once with one shared bracket.  The
grid is validated against core.HeteroMaxT(FAST) per design family in
w5_validate_fast.py; the max |dc| is deposited with the W5 results.
"""

import os
import sys

import numpy as np
from scipy.special import ndtr
from scipy.stats import chi2, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_util import f64  # noqa: E402


def _gl(a, b, npanels, order):
    x, w = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(a, b, npanels + 1)
    half = 0.5 * (edges[1:] - edges[:-1])
    mid = 0.5 * (edges[1:] + edges[:-1])
    z = (mid[:, None] + half[:, None] * x[None, :]).ravel()
    return z, (w[None, :] * half[:, None]).ravel()


def _grid(npanels, order, zmax=8.0):
    z, w = _gl(-zmax, zmax, npanels, order)
    return z, w * norm.pdf(z)


def cdf_vec(n, lam_mat, cvec, g=(24, 6, 6), orders=(5, 5, 8), zmax=8.0):
    """P(max_i t_i <= c_i) for a batch of draws.  cvec length D or scalar."""
    n = np.asarray(n, float)
    lam = np.asarray(lam_mat, float)
    D, K = lam.shape
    n0 = n[0]
    nu0 = n0 - 1.0
    c = np.asarray(cvec, float)
    if c.ndim == 0:
        c = np.full(D, f64(c))
    zu, wu = _grid(g[0], orders[0], zmax)
    zv, wv = _grid(g[1], orders[1], zmax)
    v = chi2.isf(norm.sf(zv), nu0) / (nu0 * n0)
    acc = np.ones((D, len(zu), len(zv)))
    for i in range(K):
        zi, wi = _grid(g[2], orders[2], zmax)
        yb = chi2.isf(norm.sf(zi), n[i + 1] - 1) / (n[i + 1] - 1)
        pref = np.sqrt(n[i + 1] / lam[:, i])
        inner = np.zeros((D, len(zu), len(zv)))
        for m in range(len(yb)):
            s2n = lam[:, i] * yb[m] / n[i + 1]
            arg = pref[:, None, None] * (
                zu[None, :, None] / np.sqrt(n0)
                + c[:, None, None] * np.sqrt(s2n[:, None, None] + v[None, None, :])
            )
            inner += wi[m] * ndtr(arg)
        acc = acc * inner
    return (acc * wu[None, :, None] * wv[None, None, :]).sum(axis=(1, 2))


def batch_c_fast(
    n,
    lam_mat,
    alpha=0.025,
    c_ref=None,
    width=0.7,
    iters=28,
    g=(24, 6, 6),
    orders=(5, 5, 8),
    zmax=8.0,
):
    """Vectorised bisection.  Returns (c, flags) with flags for edge hits."""
    n = np.asarray(n, float)
    lam = np.asarray(lam_mat, float)
    D = lam.shape[0]
    if c_ref is None:
        c_ref = 2.5
    lo = np.full(D, c_ref - width)
    hi = np.full(D, c_ref + width)
    target = 1.0 - alpha
    for _ in range(6):  # widen until every draw brackets the root
        f_lo = cdf_vec(n, lam, lo, g, orders, zmax) - target
        f_hi = cdf_vec(n, lam, hi, g, orders, zmax) - target
        bad_lo = f_lo > 0
        bad_hi = f_hi < 0
        if not (bad_lo.any() or bad_hi.any()):
            break
        lo = np.where(bad_lo, lo - width, lo)
        hi = np.where(bad_hi, hi + width, hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        f = cdf_vec(n, lam, mid, g, orders, zmax) - target
        up = f < 0
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    c = 0.5 * (lo + hi)
    flags = {
        "left_hit": _int(np.sum(lo <= c_ref - width + 1e-9)),
        "right_hit": _int(np.sum(hi >= c_ref + width - 1e-9)),
    }
    return c, flags


def _int(x):
    try:
        return int(x)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("not an int") from exc


def clip_lam(lam_mat, lo=1e-4, hi=1e4):
    return np.clip(np.asarray(lam_mat, float), lo, hi)


if __name__ == "__main__":
    import time

    from core import batch_c, draw_var

    rows = []
    for n_l, lam_l in (
        ([4, 13, 13], [1.0, 25.0]),
        ([4, 13, 13], [0.04, 0.04]),
        ([13, 13, 4], [1.0, 25.0]),
        ([36] * 6, [0.25, 4.0, 1.0, 9.0, 0.5]),
        ([8] * 9, [0.25, 4.0, 1.0, 9.0, 0.5, 2.0, 0.1, 6.25]),
    ):
        n = np.asarray(n_l, float)
        lam = np.asarray(lam_l, float)
        s0, si, lh = draw_var(n, lam, 24, 4242)
        ref = batch_c(n, lh, alpha=0.025, grid="fast", workers=10)
        t0 = time.time()
        got, flags = batch_c_fast(n, lh, 0.025, c_ref=f64(ref.mean()))
        dt = time.time() - t0
        rows.append(
            {
                "n": n_l,
                "lam": lam_l,
                "draws": 24,
                "max_abs_dc": f64(np.abs(got - ref).max()),
                "mean_abs_dc": f64(np.abs(got - ref).mean()),
                "mean_dc": f64((got - ref).mean()),
                "flags": flags,
                "seconds_per_draw": f64(dt / 24),
            }
        )
        print(rows[-1])
