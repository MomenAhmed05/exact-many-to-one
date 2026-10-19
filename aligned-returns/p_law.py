# -*- coding: utf-8 -*-
"""p_law.py — Unit P: the exact null law of the paired many-to-one max-t at K = 2.

D_t iid N(0, R) in R^2, t_i = sqrt(n) Dbar_i / s_i (own denominators), R = [[1, rho], [rho, 1]].
With A = sum (D_t - Dbar)(D_t - Dbar)' ~ Wishart(nu = n - 1, R), Z = sqrt(n) Dbar ~ N(0, R)
independent of A, and t_i <= c  <=>  Z_i <= c sqrt(A_ii / nu):

    A_11 = u ~ chi2_nu,  A_22 = (rho sqrt(u) + sqrt(1 - rho^2) z)^2 + (1 - rho^2) w,
    z ~ N(0,1), w ~ chi2_{nu-1}, independent;

    F(c) = P(t_1 <= c, t_2 <= c) = E[ Phi2( c sqrt(u/nu), c sqrt(A_22/nu); rho ) ].

The bivariate normal CDF is evaluated exactly through Owen's T. The expectation is a tensor
rule: composite Gauss-Legendre panels on the standard-normal score of each variable
(u = chi2.ppf(Phi(x)), etc.), truncated at |x| <= 8 (mass outside < 1.3e-15 per axis).

Gates E-1..E-3 were declared before the run.
    python p_law.py gates
"""
import math
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import owens_t, ndtr
from scipy.stats import chi2, norm, t as tdist

X_TRUNC = 8.0


def phi2(h, k, rho):
    """Bivariate standard normal CDF P(X <= h, Y <= k), corr rho, via Owen's T (vectorised)."""
    h, k = np.broadcast_arrays(np.asarray(h, float), np.asarray(k, float))
    r = math.sqrt(1.0 - rho * rho)
    out = np.empty(h.shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        ah = (k - rho * h) / (h * r)
        ak = (h - rho * k) / (k * r)
    nz = (h != 0) & (k != 0)
    beta = ((h * k < 0) | ((h * k == 0) & (h + k < 0))).astype(float) * 0.5
    out[nz] = (0.5 * ndtr(h[nz]) + 0.5 * ndtr(k[nz]) - owens_t(h[nz], ah[nz])
               - owens_t(k[nz], ak[nz]) - beta[nz])
    if (~nz).any():                       # an argument exactly 0: closed form via ndtr of the other
        hh, kk = h[~nz], k[~nz]
        # P(X <= 0, Y <= k) = int_{-inf}^0 phi(x) Phi((k - rho x)/r) dx, done by the general
        # formula with a tiny shift (the CDF is continuous); exact zeros do not occur for c > 0
        out[~nz] = phi2(np.where(hh == 0, 1e-300, hh), np.where(kk == 0, 1e-300, kk), rho)
    return out


def score_rule(panels, order):
    """Composite Gauss-Legendre nodes/weights in x on [-8, 8], weights include phi(x)."""
    g, wg = leggauss(order)
    edges = np.linspace(-X_TRUNC, X_TRUNC, panels + 1)
    xs, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        xs.append(0.5 * (b - a) * g + 0.5 * (a + b))
        ws.append(0.5 * (b - a) * wg)
    x = np.concatenate(xs)
    w = np.concatenate(ws) * norm.pdf(x)
    return x, w


class Law2:
    """F(c; rho, n) for K = 2 on a fixed tensor rule."""

    def __init__(self, n, rho, panels=16, order=8):
        self.n, self.rho, self.nu = n, rho, n - 1
        x, w = score_rule(panels, order)
        nu = self.nu
        u = chi2.ppf(ndtr(x), nu)
        z = x
        ww = chi2.ppf(ndtr(x), nu - 1) if nu > 1 else np.zeros_like(x)
        if nu == 1:                       # chi2_0 is the point mass at 0
            wwts = np.zeros_like(w)
            wwts[len(w) // 2] = 1.0
            ww = np.zeros_like(x)
        else:
            wwts = w
        U, Zz, W = np.meshgrid(u, z, ww, indexing="ij")
        self.wt = (w[:, None, None] * w[None, :, None] * wwts[None, None, :]).ravel()
        self.su = np.sqrt(U.ravel() / nu)
        a22 = (rho * np.sqrt(U) + math.sqrt(1 - rho * rho) * Zz) ** 2 + (1 - rho * rho) * W
        self.s2 = np.sqrt(a22.ravel() / nu)
        self.norm_ = self.wt.sum()

    def F(self, c):
        """P(t1 <= c, t2 <= c)."""
        return float(phi2(c * self.su, c * self.s2, self.rho) @ self.wt / self.norm_)

    def rect(self, c):
        """P(|t1| <= c, |t2| <= c)."""
        h, k, r = c * self.su, c * self.s2, self.rho
        v = phi2(h, k, r) - phi2(-h, k, r) - phi2(h, -k, r) + phi2(-h, -k, r)
        return float(v @ self.wt / self.norm_)

    def cstar(self, alpha=0.025, lo=0.5, hi=20.0, tol=1e-10):
        """c with P(max t_i > c) = alpha (monotone bisection)."""
        target = 1.0 - alpha
        for _ in range(200):
            m = 0.5 * (lo + hi)
            if self.F(m) < target:
                lo = m
            else:
                hi = m
            if hi - lo < tol:
                break
        return 0.5 * (lo + hi)


def brute(n, rho, c, N, seed, batch=200_000):
    rng = np.random.default_rng(seed)
    L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]))
    hits = done = 0
    while done < N:
        m = min(batch, N - done)
        D = rng.standard_normal((m, n, 2)) @ L.T
        t = math.sqrt(n) * D.mean(1) / D.std(1, ddof=1)
        hits += int((t.max(1) <= c).sum())
        done += m
    p = hits / N
    return p, math.sqrt(p * (1 - p) / N)


def gates():
    rows = []

    def claim(cid, label, value):
        rows.append("CLAIM %s: %s = %s" % (cid, label, value))
        print(rows[-1], flush=True)

    # E-0: Owen's-T bivariate CDF against scipy's Genz routine at random points
    from scipy.stats import multivariate_normal
    rng = np.random.default_rng(5)
    worst = 0.0
    for rho in (-0.9, -0.3, 0.0, 0.5, 0.95):
        mvn = multivariate_normal([0, 0], [[1, rho], [rho, 1]])
        pts = rng.uniform(-3, 3, size=(40, 2))
        a = phi2(pts[:, 0], pts[:, 1], rho)
        b = np.array([mvn.cdf(p) for p in pts])
        worst = max(worst, float(np.abs(a - b).max()))
    claim("PLAW-E0", "max |Owen-T Phi2 - scipy Genz| over 200 points", "%.2e" % worst)
    # E-1: independence anchor
    c = 2.3
    for n in (12, 36):
        L = Law2(n, 0.0)
        d = abs(L.F(c) - tdist.cdf(c, n - 1) ** 2)
        claim("PLAW-E1-n%d" % n, "|F(rho=0) - P(t<=c)^2| (pass <= 1e-9)", "%.2e pass %s" % (d, d <= 1e-9))
    # E-2 refinement and E-3 brute force
    k = 0
    for n in (12, 36):
        for rho in (-0.5, 0.0, 0.5, 0.95):
            a = Law2(n, rho).F(c)
            b = Law2(n, rho, panels=32, order=8).F(c)
            pb, se = brute(n, rho, c, 4_000_000, 91000 + k)
            k += 1
            z = (pb - a) / se
            claim("PLAW-E23-n%d-r%+.2f" % (n, rho),
                  "F; refine move (pass < 1e-7); brute (N=4e6) z (pass |z|<=4)",
                  "%.9f; %.1e %s; %.6f z %+.2f %s" % (a, abs(b - a), abs(b - a) < 1e-7,
                                                      pb, z, abs(z) <= 4))
    open("deposits/p_law_gates_claims.txt", "w", encoding="utf-8").write("\n".join(rows) + "\n")


if __name__ == "__main__":
    if sys.argv[1:] == ["gates"]:
        gates()
