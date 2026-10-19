# -*- coding: utf-8 -*-
"""p_killround_k5.py — the paired kill round at K = 5 (Section 8).

DECLARATION (written before the run):

Designs: K = 5, n in {12, 36}; correlation structures
  EQ0  independence; EQ5 equicorrelated 0.5; EQ9 equicorrelated 0.9;
  BM   benchmark-dominated rho_ij = h_i h_j, h = (0.9, 0.8, 0.6, 0.4, 0.2)
       (the one-factor form a shared benchmark induces; paper 1 equation (11)).
Rules, all at alpha = 0.025 one-sided:
  ORACLE  c*(R): 0.975 quantile of max t_i under the true R, N = 4,000,000 direct draws
          (Bartlett construction of the Wishart diagonal; quantile SE by 10 batches).
  BONF    c = t_nu(1 - alpha/K).
  MVT     common-denominator multivariate t at the TRUE R (practice): 0.975 quantile of
          max Z_i / sqrt(w/nu), Z ~ N(0, R), w ~ chi2_nu, N = 4,000,000.
  PLUGIN  exact law at R_hat: per outer draw, c*(R_hat) is the 0.975 quantile of max t
          under R_hat from M = 20,000 inner Bartlett draws; outer N = 20,000 (SE 1.1e-3).
Measurements and verdicts:
  1K5  PLUGIN materially liberal (> 1.05a and > a + 3SE) in some cell -> conj 1 SURVIVES at K=5.
  3K5  MVT's true rate > 1.05a in some cell (true rate from N = 4,000,000 own-denominator
       draws at c_mvt, SE 7.8e-5) -> conj 3 SURVIVES at K=5.
  7K5  Bonferroni's power loss vs ORACLE at the effect (one candidate shifted) where the oracle
       has 80% power, N = 400,000 common random numbers -> conj 7 KILLED at K=5 if loss >= 2pp
       beyond 3SE in some cell.
Seeds 85000 + 100*design + offset.
"""
import math
import sys
import time
import traceback

import numpy as np
from scipy.optimize import brentq
from scipy.stats import t as tdist

ALPHA, K = 0.025, 5
rows = []


def claim(cid, label, value):
    rows.append("CLAIM %s: %s = %s" % (cid, label, value))
    print(rows[-1], flush=True)


def structures():
    h = np.array([0.9, 0.8, 0.6, 0.4, 0.2])
    bm = np.outer(h, h)
    np.fill_diagonal(bm, 1.0)
    eq = lambda r: np.full((K, K), r) + (1 - r) * np.eye(K)
    return {"EQ0": eq(0.0), "EQ5": eq(0.5), "EQ9": eq(0.9), "BM": bm}


def draw_stats(R, n, N, rng, shift=None):
    """Z = sqrt(n) Dbar and the sample variances s_i^2, via Z ~ N(0,R) and a Bartlett
    Wishart(nu, R) draw: A = L T T' L', s_i^2 = A_ii / nu. Returns (Z, s)."""
    nu = n - 1
    L = np.linalg.cholesky(R)
    Z = rng.standard_normal((N, K)) @ L.T
    if shift is not None:
        Z = Z + shift
    T = np.zeros((N, K, K))
    for i in range(K):
        T[:, i, i] = np.sqrt(rng.chisquare(nu - i, N))
        for j in range(i):
            T[:, i, j] = rng.standard_normal(N)
    LT = np.einsum("ij,njk->nik", L, T)
    s = np.sqrt((LT ** 2).sum(2) / nu)
    return Z, s


def maxt(R, n, N, rng, chunk=500_000):
    out = []
    done = 0
    while done < N:
        m = min(chunk, N - done)
        Z, s = draw_stats(R, n, m, rng)
        out.append((Z / s).max(1))
        done += m
    return np.concatenate(out)


def quantile_with_se(x, p, nb=10):
    qs = [np.quantile(b, p) for b in np.array_split(x, nb)]
    return float(np.quantile(x, p)), float(np.std(qs, ddof=1) / math.sqrt(nb))


def main():
    t0 = time.time()
    S = structures()
    for di, (name, R) in enumerate(S.items()):
        for n in (12, 36):
            nu = n - 1
            base = 85000 + 100 * di + n
            # ORACLE
            mt = maxt(R, n, 4_000_000, np.random.default_rng(base + 1))
            c_or, se_or = quantile_with_se(mt, 1 - ALPHA)
            c_b = tdist.isf(ALPHA / K, nu)
            # MVT at the true R, and its true rate
            rng = np.random.default_rng(base + 2)
            L = np.linalg.cholesky(R)
            vals = []
            for _ in range(8):
                Z = rng.standard_normal((500_000, K)) @ L.T
                w = rng.chisquare(nu, 500_000)
                vals.append(Z.max(1) / np.sqrt(w / nu))
            c_mvt, _ = quantile_with_se(np.concatenate(vals), 1 - ALPHA)
            rate_mvt = float((mt > c_mvt).mean())
            rate_b = float((mt > c_b).mean())
            se4 = math.sqrt(ALPHA * (1 - ALPHA) / 4_000_000)
            claim("K5-%s-n%d-CONST" % (name, n), "oracle c* (se), Bonferroni c, mvt c",
                  "%.4f (%.4f), %.4f, %.4f" % (c_or, se_or, c_b, c_mvt))
            claim("K5-%s-n%d-RATES" % (name, n), "true rate of Bonferroni / mvt (SE 7.8e-5)",
                  "%.5f / %.5f (%.2fx)" % (rate_b, rate_mvt, rate_mvt / ALPHA))
            # PLUGIN
            rng = np.random.default_rng(base + 3)
            NO, M = 20_000, 20_000
            Zo, so = draw_stats(R, n, NO, rng)
            hits = 0
            # R_hat needs the full sample covariance: rebuild from raw data for the outer draws
            rng2 = np.random.default_rng(base + 4)
            for k in range(NO):
                D = rng2.standard_normal((n, K)) @ L.T
                Sh = np.cov(D, rowvar=False)
                sd = np.sqrt(np.diag(Sh))
                Rh = Sh / np.outer(sd, sd)
                tt = math.sqrt(n) * D.mean(0) / sd
                Zi, si = draw_stats(Rh + 1e-12 * np.eye(K), n, M, rng)
                ch = np.quantile((Zi / si).max(1), 1 - ALPHA)
                hits += int(tt.max() > ch)
            p = hits / NO
            se = math.sqrt(ALPHA * (1 - ALPHA) / NO)
            claim("K5-%s-n%d-PLUGIN" % (name, n), "plug-in realised FWER (outer 20k, inner 20k; SE 1.1e-3)",
                  "%.5f (z %+.1f)" % (p, (p - ALPHA) / se))
            # power: oracle vs Bonferroni, one candidate shifted
            rng = np.random.default_rng(base + 5)
            Zp, sp = draw_stats(R, n, 400_000, rng)
            def pw(d):
                Zs = Zp.copy()
                Zs[:, 0] += d
                m = (Zs / sp).max(1)
                return float((m > c_or).mean()), float((m > c_b).mean())
            d80 = brentq(lambda d: pw(d)[0] - 0.8, 0.0, 10.0, xtol=1e-6)
            po, pb = pw(d80)
            claim("K5-%s-n%d-POWER" % (name, n), "power oracle / Bonferroni at the oracle's 80% effect, loss (pp)",
                  "%.4f / %.4f, %.2f" % (po, pb, 100 * (po - pb)))
    claim("K5-WALL", "seconds", "%.0f" % (time.time() - t0))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        claim("K5-ERROR", "the run failed", traceback.format_exc().replace("\n", " | "))
    finally:
        open("deposits/p_killround_k5_claims.txt", "w", encoding="utf-8").write("\n".join(rows) + "\n")
