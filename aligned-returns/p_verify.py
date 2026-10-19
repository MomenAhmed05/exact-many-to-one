# -*- coding: utf-8 -*-
"""p_verify.py — verification of the paired results (Section 8).

    python p_verify.py A     Track A: T2 gates (F-1..F-4), T3 maps, T4 maps
    python p_verify.py B     Track B: adversarial cells (F-6), T5 Sidak search
    python p_verify.py C     Track C: independent K = 2 route (F-5)

Each track writes deposits/p_verify_<track>_claims.txt ALWAYS. Seeds: A 10xxxx, B 20xxxx,
C 30xxxx. Tolerances were declared before the run (C-4) and are restated beside each check.
"""
import math
import sys
import time
import traceback

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import chi2, t as tdist
from numpy.polynomial.legendre import leggauss

import p_law as P
import p_killround_k5 as K5

ALPHA = 0.025
rows = []


def claim(cid, label, value):
    rows.append("CLAIM %s: %s = %s" % (cid, label, value))
    print(rows[-1], flush=True)


# ------------------------------------------------------------------ shared pieces

def c_mvt2(n, rho):
    nu = n - 1
    x, w = P.score_rule(128, 8)
    s = np.sqrt(chi2.ppf(ndtr(x), nu) / nu)
    return brentq(lambda c: float(P.phi2(c * s, c * s, rho) @ w / w.sum()) - (1 - ALPHA), 0.3, 60, xtol=1e-12)


def cstar2(n, rho, alpha=ALPHA, panels=16):
    L = P.Law2(n, rho, panels=panels)
    return brentq(lambda c: L.F(c) - (1 - alpha), 0.3, 80.0, xtol=1e-10)


def MIX():
    """Declared mixed-sign 5x5 correlation matrix (positive definite, checked)."""
    R = np.array([[1.0, 0.6, -0.3, 0.2, -0.5],
                  [0.6, 1.0, -0.2, 0.4, -0.3],
                  [-0.3, -0.2, 1.0, 0.1, 0.5],
                  [0.2, 0.4, 0.1, 1.0, 0.0],
                  [-0.5, -0.3, 0.5, 0.0, 1.0]])
    assert np.linalg.eigvalsh(R).min() > 0.05
    return R


def mvt_const(R, nu, N, rng):
    L = np.linalg.cholesky(R)
    v = []
    for _ in range(N // 500_000):
        Z = rng.standard_normal((500_000, R.shape[0])) @ L.T
        w = rng.chisquare(nu, 500_000)
        v.append(Z.max(1) / np.sqrt(w / nu))
    return float(np.quantile(np.concatenate(v), 1 - ALPHA))


def plugin_k5(R, n, NO, M, seed):
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
        # inner law at R_hat by direct data generation (valid for any nu, including nu < K)
        Lh = np.linalg.cholesky(Rh + 1e-9 * np.eye(R.shape[0]))
        Di = rng.standard_normal((M, n, R.shape[0])) @ Lh.T if n * M <= 4_000_000 else None
        if Di is None:
            Zi, si = K5.draw_stats(Rh + 1e-9 * np.eye(R.shape[0]), n, M, rng)
            mi = (Zi / si).max(1)
        else:
            mi = (math.sqrt(n) * Di.mean(1) / Di.std(1, ddof=1)).max(1)
        hits += int(tt.max() > np.quantile(mi, 1 - ALPHA))
    return hits / NO


# ------------------------------------------------------------------ Track A

def track_A():
    # F-1 known case: R = I at K = 2, deterministic
    worst = 0.0
    for n in (4, 12, 36, 120, 400):
        L = P.Law2(n, 0.0)
        for c in (0.5, 1.5, 2.3, 4.0, 6.0):
            worst = max(worst, abs(L.F(c) - tdist.cdf(c, n - 1) ** 2))
    claim("PVA-F1-K2", "max |F(rho=0) - P(t<=c)^2| over 25 (n, c) (tol 1e-12)", "%.2e pass %s" % (worst, worst <= 1e-12))
    # F-1 at K = 5, simulated
    for n in (8, 36):
        mt = K5.maxt(np.eye(5), n, 2_000_000, np.random.default_rng(100100 + n))
        q, se = K5.quantile_with_se(mt, 1 - ALPHA)
        ex = tdist.ppf((1 - ALPHA) ** (1 / 5), n - 1)
        claim("PVA-F1-K5-n%d" % n, "R=I oracle constant vs exact (|z|<=4)", "%.4f vs %.4f, z %+.2f" % (q, ex, (q - ex) / se))
    # F-2 degenerate corners rho -> +-1
    for n in (12, 36):
        for rho, lim in ((-0.999, lambda c, nu: 2 * tdist.cdf(c, nu) - 1), (0.999, lambda c, nu: tdist.cdf(c, nu))):
            L = P.Law2(n, rho, panels=32)
            d = max(abs(L.F(c) - lim(c, n - 1)) for c in (1.0, 2.0, 3.0))
            claim("PVA-F2-n%d-r%+.3f" % (n, rho), "max |F - limit| at c in {1,2,3} (tol 1e-3)", "%.2e pass %s" % (d, d <= 1e-3))
        seq = [P.Law2(n, r, panels=32).F(2.0) for r in (0.9, 0.99, 0.999)]
        claim("PVA-F2-mono-n%d" % n, "F(2) at rho 0.9, 0.99, 0.999 approaches P(t<=2) monotonically",
              "%s -> %.6f, monotone %s" % (["%.6f" % v for v in seq], tdist.cdf(2.0, n - 1),
                                           seq[0] <= seq[1] <= seq[2] <= tdist.cdf(2.0, n - 1) + 1e-12))
    # F-4 refinement over the T2 domain corners
    worst = 0.0
    for n in (4, 400):
        for rho in (-0.95, 0.0, 0.95):
            a = P.Law2(n, rho, panels=16, order=8)
            b = P.Law2(n, rho, panels=32, order=8)
            cc = P.Law2(n, rho, panels=16, order=12)
            for c in (0.5, 2.3, 6.0):
                fa = a.F(c)
                worst = max(worst, abs(b.F(c) - fa), abs(cc.F(c) - fa))
    claim("PVA-F4", "max refinement move (panels x2, order 8->12) over n in {4,400}, rho in {-.95,0,.95}, c in {.5,2.3,6} (tol 1e-12)",
          "%.2e pass %s" % (worst, worst <= 1e-12))
    # F-3 brute force at 24 points
    k = 0
    zmax = 0.0
    for n in (4, 12, 36, 120):
        for rho in (-0.9, -0.5, 0.0, 0.5, 0.9, 0.95):
            c = 2.3
            f = P.Law2(n, rho).F(c)
            pb, se = P.brute(n, rho, c, 4_000_000, 100200 + k)
            k += 1
            z = (pb - f) / se
            zmax = max(zmax, abs(z))
            claim("PVA-F3-n%d-r%+.2f" % (n, rho), "F vs brute (N=4e6)", "%.6f vs %.6f, z %+.2f" % (f, pb, z))
    claim("PVA-F3", "max |z| over 24 brute-force points (tol 4)", "%.2f pass %s" % (zmax, zmax <= 4))

    # T3 K = 2 deterministic map
    for n in (6, 8, 12, 24, 36, 60, 120):
        cells = []
        for rho in (-0.9, -0.5, 0.0, 0.5, 0.9, 0.95):
            c = c_mvt2(n, rho)
            rate = 1 - P.Law2(n, rho).F(c)
            cells.append("%+.2f:%.5f(%.3fx)" % (rho, rate, rate / ALPHA))
        claim("PVA-T3-K2-n%d" % n, "multivariate-t true rate by rho", "; ".join(cells))
    # T3 brute check at 6 points
    for j, (n, rho) in enumerate(((6, 0.5), (6, 0.95), (12, 0.9), (36, 0.5), (8, -0.5), (24, 0.9))):
        c = c_mvt2(n, rho)
        rate = 1 - P.Law2(n, rho).F(c)
        pb, se = P.brute(n, rho, c, 4_000_000, 100300 + j)
        claim("PVA-T3-BRUTE-n%d-r%+.2f" % (n, rho), "T3 rate vs brute (|z|<=4)", "%.5f vs %.5f, z %+.2f" % (rate, 1 - pb, ((1 - pb) - rate) / se))
    # T3 K = 5 simulated map
    S = dict(K5.structures(), MIX=MIX())
    for si, (name, R) in enumerate(S.items()):
        cells = []
        for n in (8, 12, 24, 36, 60):
            c = mvt_const(R, n - 1, 4_000_000, np.random.default_rng(100400 + 10 * si + n))
            mt = K5.maxt(R, n, 4_000_000, np.random.default_rng(100500 + 10 * si + n))
            r = float((mt > c).mean())
            cells.append("n%d:%.5f(%.2fx)" % (n, r, r / ALPHA))
        claim("PVA-T3-K5-%s" % name, "multivariate-t true rate by n (SE 7.8e-5)", "; ".join(cells))

    # T4 K = 2 plug-in map
    grid = np.tanh(np.linspace(np.arctanh(-0.995), np.arctanh(0.995), 81))
    for n in (6, 8, 12, 24, 36):
        sp = CubicSpline(grid, [cstar2(n, r, panels=8) for r in grid])
        cells = []
        for j, rho in enumerate((-0.5, 0.0, 0.5, 0.9)):
            rng = np.random.default_rng(100600 + 10 * n + j)
            L = np.linalg.cholesky(np.array([[1, rho], [rho, 1.0]]))
            hits = done = 0
            while done < 400_000:
                m = 100_000
                D = rng.standard_normal((m, n, 2)) @ L.T
                Dc = D - D.mean(1, keepdims=True)
                s = np.sqrt((Dc ** 2).sum(1) / (n - 1))
                rh = np.clip((Dc[:, :, 0] * Dc[:, :, 1]).sum(1) / (n - 1) / (s[:, 0] * s[:, 1]), -0.995, 0.995)
                t = math.sqrt(n) * D.mean(1) / s
                hits += int((t.max(1) > sp(rh)).sum())
                done += m
            p = hits / 400_000
            cells.append("%+.1f:%.5f(z%+.1f)" % (rho, p, (p - ALPHA) / math.sqrt(ALPHA * (1 - ALPHA) / 4e5)))
        claim("PVA-T4-K2-n%d" % n, "plug-in realised FWER by rho (N=4e5)", "; ".join(cells))
    # T4 K = 5 plug-in map
    for si, (name, R) in enumerate(S.items()):
        cells = []
        for n in (6, 8, 12, 24, 36):
            p = plugin_k5(R, n, 20_000, 20_000, 100700 + 100 * si + n)
            cells.append("n%d:%.5f(z%+.1f)" % (n, p, (p - ALPHA) / math.sqrt(ALPHA * (1 - ALPHA) / 2e4)))
        claim("PVA-T4-K5-%s" % name, "plug-in realised FWER by n (SE 1.1e-3)", "; ".join(cells))


# ------------------------------------------------------------------ Track B

def rect3_rb(R, n, cs, N, rng, nz=48):
    """Rao-Blackwellised P(|t_i| <= c, i=1..3) for each c in cs: conditional on the variance
    diagonal, integrate z1 by Gauss-Legendre and use the bivariate rectangle for (Z2, Z3) | z1."""
    nu = n - 1
    L = np.linalg.cholesky(R)
    X = rng.standard_normal((N, nu, 3)) @ L.T             # A = X'X, valid for any nu
    s = np.sqrt((X ** 2).sum(1) / nu)                      # own standard deviations (scaled)
    r12, r13, r23 = R[0, 1], R[0, 2], R[1, 2]
    v2, v3 = 1 - r12 ** 2, 1 - r13 ** 2
    cov23 = r23 - r12 * r13
    rr = cov23 / math.sqrt(v2 * v3)
    g, wg = leggauss(nz)
    out = []
    for c in cs:
        h1, h2, h3 = c * s[:, 0], c * s[:, 1], c * s[:, 2]
        acc = np.zeros(N)
        for gi, wi in zip(g, wg):
            z1 = h1 * gi                                    # maps (-1,1) to (-h1, h1)
            m2, m3 = r12 * z1, r13 * z1
            a2, b2 = (-h2 - m2) / math.sqrt(v2), (h2 - m2) / math.sqrt(v2)
            a3, b3 = (-h3 - m3) / math.sqrt(v3), (h3 - m3) / math.sqrt(v3)
            rect = (P.phi2(b2, b3, rr) - P.phi2(a2, b3, rr) - P.phi2(b2, a3, rr) + P.phi2(a2, a3, rr))
            acc += wi * h1 * np.exp(-0.5 * z1 * z1) / math.sqrt(2 * math.pi) * rect
        out.append((float(acc.mean()), float(acc.std(ddof=1) / math.sqrt(N))))
    return out


def track_B():
    # F-6 adversarial: smallest n and mixed signs are inside T3/T4 (Track A runs them);
    # here the near-singular plug-in case n <= K + 1 at K = 5 and the extreme correlations at K = 2.
    for n in (6, 7):
        p = plugin_k5(MIX(), n, 20_000, 20_000, 200100 + n)
        claim("PVB-F6-K5-MIX-n%d" % n, "plug-in realised FWER, n near K+1 (SE 1.1e-3)",
              "%.5f (z %+.1f)" % (p, (p - ALPHA) / math.sqrt(ALPHA * (1 - ALPHA) / 2e4)))
    for n in (6, 8):
        for rho in (-0.95, 0.95):
            c = c_mvt2(n, rho)
            rate = 1 - P.Law2(n, rho, panels=32).F(c)
            claim("PVB-F6-K2-n%d-r%+.2f" % (n, rho), "multivariate-t true rate at extreme rho", "%.5f (%.3fx)" % (rate, rate / ALPHA))
    # T5 Sidak search at K = 3, own denominators, two-sided
    rng = np.random.default_rng(200200)
    mats = []
    for a, b, cc in ((0.9, -0.9, -0.81), (0.9, 0.9, 0.81), (0.5, -0.5, 0.0), (0.95, 0.0, 0.0),
                     (-0.45, -0.45, -0.45), (0.7, -0.7, -0.4), (0.99, -0.5, -0.5), (0.3, 0.3, -0.9)):
        R = np.array([[1, a, b], [a, 1, cc], [b, cc, 1.0]])
        if np.linalg.eigvalsh(R).min() > 1e-3:
            mats.append(R)
    while len(mats) < 48:
        G = rng.standard_normal((3, 3))
        S = G @ G.T + 0.05 * np.eye(3)
        d = np.sqrt(np.diag(S))
        mats.append(S / np.outer(d, d))
    claim("PVB-T5-DESIGN", "correlation matrices (structured + random), n values, c values",
          "%d, (2, 3, 5, 12), (0.5, 1, 2, 3, 5)" % len(mats))
    worst = (1.0, None)
    nviol = 0
    for mi, R in enumerate(mats):
        for n in (2, 3, 5, 12):
            res = rect3_rb(R, n, (0.5, 1.0, 2.0, 3.0, 5.0), 40_000, np.random.default_rng(200300 + 10 * mi + n))
            for c, (p, se) in zip((0.5, 1.0, 2.0, 3.0, 5.0), res):
                prod = (2 * tdist.cdf(c, n - 1) - 1) ** 3
                d = p - prod
                if d < -4 * se and d < -1e-6:
                    nviol += 1
                    claim("PVB-T5-VIOL-%d" % nviol, "candidate violation (R, n, c, diff, se)",
                          "%s, %d, %.1f, %.3e, %.1e" % (np.round(R[np.triu_indices(3, 1)], 4).tolist(), n, c, d, se))
                if d / max(se, 1e-12) < worst[0]:
                    worst = (d / max(se, 1e-12), (np.round(R[np.triu_indices(3, 1)], 4).tolist(), n, c, d, se))
    claim("PVB-T5-NEARMISS", "smallest (diff / SE) over the grid, and its configuration",
          "%.2f at %s" % (worst[0], worst[1]))
    claim("PVB-T5-VERDICT", "violations beyond 4 SE and 1e-6", "%d" % nviol)


# ------------------------------------------------------------------ Track C

def track_C():
    """A second K = 2 evaluator sharing no code with p_law: Bartlett draws of the variance
    diagonal and scipy's Genz bivariate normal CDF."""
    from scipy.stats import multivariate_normal
    rng = np.random.default_rng(300100)
    zmax = 0.0
    for j in range(20):
        n = int(rng.choice([4, 6, 12, 36, 120]))
        rho = float(rng.uniform(-0.9, 0.95))
        c = float(rng.uniform(1.0, 4.0))
        nu = n - 1
        N = 200_000
        u = rng.chisquare(nu, N)
        z = rng.standard_normal(N)
        w = rng.chisquare(nu - 1, N) if nu > 1 else np.zeros(N)
        a22 = (rho * np.sqrt(u) + math.sqrt(1 - rho * rho) * z) ** 2 + (1 - rho * rho) * w
        pts = np.column_stack([c * np.sqrt(u / nu), c * np.sqrt(a22 / nu)])
        mvn = multivariate_normal([0, 0], [[1, rho], [rho, 1]])
        vals = mvn.cdf(pts)
        est, se = float(vals.mean()), float(vals.std(ddof=1) / math.sqrt(N))
        f = P.Law2(n, rho).F(c)
        z_ = (est - f) / se
        zmax = max(zmax, abs(z_))
        claim("PVC-F5-%02d" % j, "n, rho, c; route C vs T2 (|z|<=4)", "%d, %+.3f, %.3f; %.7f vs %.7f, z %+.2f" % (n, rho, c, est, f, z_))
    claim("PVC-F5", "max |z| over 20 random interior points (tol 4)", "%.2f pass %s" % (zmax, zmax <= 4))


if __name__ == "__main__":
    tr = sys.argv[1]
    t0 = time.time()
    try:
        {"A": track_A, "B": track_B, "C": track_C}[tr]()
    except Exception:
        claim("PV%s-ERROR" % tr, "the run failed", traceback.format_exc().replace("\n", " | "))
    finally:
        claim("PV%s-WALL" % tr, "seconds", "%.0f" % (time.time() - t0))
        open("deposits/p_verify_%s_claims.txt" % tr, "w", encoding="utf-8").write("\n".join(rows) + "\n")
