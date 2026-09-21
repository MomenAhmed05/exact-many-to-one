"""W2, Phase 2A: the exact null law for ALIGNED returns with a latent factor.

THE CHANGED OBJECT.  The preceding programme assumes independent groups.  A desk
observes the benchmark and every candidate on the SAME dates, so they move
together, and the independent-groups law is then the wrong law.  The roadmap as
written formed paired differences; that turns the problem into a maximum of K
studentised statistics with different studentisers and a Wishart-coupled
denominator set, which nothing factorises.

This module takes the other route.  Keep the statistic the desk already computes
and carry the dependence in a LATENT factor:

    r_{i,t} = mu_i + beta_i f_t + eps_{i,t},     i = 0..K,  t = 1..n,

f latent Gaussian, eps independent across i.  Then conditioning on FOUR
quantities makes the K candidates exactly conditionally independent:

    fbar,   Vf = sum_t (f_t - fbar)^2,   Xbar_0,   s_0^2 .

Why it works.  Given the factor path, Xbar_i depends on it only through fbar,
and the centred residual of group i is beta_i (f - fbar) + (eps_i - eps_i-bar);
its squared norm is noncentral chi-square on n-1 df with noncentrality
beta_i^2 Vf / sigma_i^2, and by rotational symmetry of eps_i that law -- and the
joint law across i -- depends on the factor path ONLY through Vf.  eps_i-bar is
orthogonal to the centred eps_i, so Xbar_i and s_i^2 stay independent.  Nothing
couples candidate i to candidate j once (fbar, Vf) is fixed, and nothing couples
either to the benchmark once (Xbar_0, s_0^2) is fixed.

So the null law is a FOUR-dimensional outer rule over K one-dimensional inner
rules, and the outer dimension does not grow with K.  Two of the four are the
two the existing paper already integrates; the two new ones are the factor's own
mean and sum of squares.

NORMALISATION.  sigma_f = 1 (absorbed into the loadings) and sigma_0 = 1 (T is
scale invariant).  Free parameters: n, the loadings beta_0..beta_K, and the
variance ratios lam_i = sigma_i^2 / sigma_0^2.

Deposit: deposits/w2_factor.json
"""

from __future__ import annotations

import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.stats import ncx2, norm

import fcore as F
from h12_quad import chi2_over_df_nodes, normal_nodes

ABAR = F.ABAR

# outer/inner resolutions; GRID is the reference rule, FAST is for search loops
GRID = {"zf": (3, 10), "z0": (3, 10), "vf": (3, 10), "s0": 24, "si": 24}
FAST = {"zf": (2, 8), "z0": (2, 8), "vf": (2, 8), "s0": 16, "si": 16}


def gl01(m):
    """m Gauss-Legendre nodes and weights on (0, 1)."""
    x, w = leggauss(int(m))
    return 0.5 * (x + 1.0), 0.5 * w


def ncx2_nodes(df, lam, m):
    """Probability-scale nodes for a noncentral chi-square, vectorised in lam.

    Returns (nodes, weights); nodes has shape (lam.size, m).
    """
    u, w = gl01(m)
    lam = np.atleast_1d(np.asarray(lam, float))
    return ncx2.ppf(u[None, :], df, lam[:, None]), w


class AlignedFactor(object):
    """Deterministic evaluator for the aligned latent-factor null law."""

    def __init__(self, n, beta, lam, grid=None):
        self.n = int(n)
        self.beta = np.asarray(beta, float)              # length K+1
        self.lam = np.concatenate([[1.0], np.asarray(lam, float)])
        self.K = self.beta.size - 1
        assert self.lam.size == self.K + 1
        g = dict(GRID if grid is None else grid)
        self.g = g
        nu = self.n - 1.0

        self.zf, self.wzf = normal_nodes(*g["zf"])
        self.z0, self.wz0 = normal_nodes(*g["z0"])
        q, wv = chi2_over_df_nodes(nu, *g["vf"])
        self.vf, self.wvf = q * nu, wv                   # raw chi2_{n-1} nodes

        # noncentral chi-square nodes, one set per Vf node, per group
        self.s_nodes = []
        for i in range(self.K + 1):
            lamnc = self.beta[i] ** 2 * self.vf / self.lam[i]
            m = g["s0"] if i == 0 else g["si"]
            nodes, w = ncx2_nodes(nu, lamnc, m)
            self.s_nodes.append(self.lam[i] * nodes / nu)    # -> s_i^2
            if i == 0:
                self.s0_w = w
            else:
                self.s_w = w

    def noreject(self, c):
        """P(T_i <= c_i for every i) under the null, deterministic."""
        n = self.n
        cv = np.broadcast_to(np.asarray(c, float), (self.K,))
        rt = np.sqrt(n)
        total = 0.0
        wt = (self.wzf[:, None, None] * self.wz0[None, :, None])
        for j in range(self.vf.size):
            s0 = self.s_nodes[0][j]                              # (ns0,)
            acc = np.ones((self.zf.size, self.z0.size, s0.size))
            for i in range(1, self.K + 1):
                si = self.s_nodes[i][j]                          # (nsi,)
                den = np.sqrt(si[None, :] / n + s0[:, None] / n)  # (ns0,nsi)
                # Xbar_0 - beta_i fbar = ((beta_0 - beta_i) zf + z0)/sqrt(n)
                m = ((self.beta[0] - self.beta[i]) * self.zf[:, None]
                     + self.z0[None, :]) / rt                    # (nzf,nz0)
                arg = (cv[i - 1] * den[None, None, :, :]
                       + m[:, :, None, None]) / np.sqrt(self.lam[i] / n)
                acc *= (norm.cdf(arg) * self.s_w[None, None, None, :]).sum(3)
            total += self.wvf[j] * float(
                (acc * wt * self.s0_w[None, None, :]).sum())
        return total

    def fwer(self, c):
        return 1.0 - self.noreject(c)

    def critical(self, alpha=ABAR, lo=0.5, hi=12.0, tol=1e-9):
        """The common constant c with FWER = alpha."""
        for _ in range(60):
            m = 0.5 * (lo + hi)
            if self.fwer(m) > alpha:
                lo = m
            else:
                hi = m
            if hi - lo < tol:
                break
        return 0.5 * (lo + hi)


def brute(n, beta, lam, c, D, seed, batch=50000):
    """Full-panel simulation: the definition of truth for the falsifiers."""
    beta = np.asarray(beta, float)
    sig = np.sqrt(np.concatenate([[1.0], np.asarray(lam, float)]))
    K = beta.size - 1
    cv = np.broadcast_to(np.asarray(c, float), (K,))
    rng = np.random.default_rng(int(seed))
    hits = done = 0
    while done < D:
        b = min(batch, D - done)
        f = rng.normal(size=(b, n))
        eps = rng.normal(size=(b, K + 1, n)) * sig[None, :, None]
        r = beta[None, :, None] * f[:, None, :] + eps
        xb, s2 = r.mean(2), r.var(2, ddof=1)
        T = (xb[:, 1:] - xb[:, 0:1]) / np.sqrt(s2[:, 1:] / n + s2[:, 0:1] / n)
        hits += int((T > cv[None, :]).any(1).sum())
        done += b
    p = hits / D
    return p, np.sqrt(p * (1 - p) / D)


CELLS = [
    ("no loadings, independent", 36, [0.0, 0.0, 0.0], [1.0, 1.0], 2.20),
    ("equal loadings", 36, [1.0, 1.0, 1.0], [1.0, 1.0], 2.20),
    ("unequal loadings", 36, [0.8, 1.3, 0.4], [1.0, 1.0], 2.20),
    ("unequal loadings and ratios", 36, [0.8, 1.3, 0.4], [4.0, 0.25], 2.20),
    ("thirteen observations", 13, [0.8, 1.3, 0.4], [4.0, 0.25], 2.40),
    ("strong loadings", 36, [0.5, 2.5, 1.8], [1.0, 1.0], 2.20),
    ("one negative loading", 36, [1.0, -0.8, 1.5], [2.0, 0.5], 2.20),
    ("five candidates", 36, [0.9, 1.4, 0.3, 1.1, 0.6, 1.0],
     [2.25, 0.49, 1.0, 4.0, 0.64], 2.55),
]


def falsifiers(D=4000000, seed=99001):
    """Each cell: deterministic evaluator against brute force.

    The first cell is the non-vacuity anchor: at beta = 0 the law must collapse
    to the independent-groups law the existing paper computes, and at beta != 0
    it must move away from it.  A reduction that quietly ignored the factor
    would pass every other cell and fail that comparison.
    """
    rows = []
    print("%-22s %4s %10s %10s %8s %10s" %
          ("cell", "n", "brute", "exact", "z", "indep law"))
    for k, (tag, n, beta, lam, c) in enumerate(CELLS):
        A = AlignedFactor(n, beta, lam)
        pe = A.fwer(c)
        pb, eb = brute(n, beta, lam, c, D, seed + k)
        pi = AlignedFactor(n, [0.0] * len(beta), lam).fwer(c)
        z = (pe - pb) / eb
        rows.append({"cell": tag, "n": n, "beta": [float(x) for x in beta],
                     "lam": [float(x) for x in lam], "c": float(c),
                     "exact": float(pe), "brute": float(pb),
                     "brute_se": float(eb), "z": float(z),
                     "independent_law": float(pi)})
        print("%-22s %4d %10.5f %10.5f %+8.2f %10.5f"
              % (tag, n, pb, pe, z, pi))
    return rows


def main():
    t0 = time.time()
    rows = falsifiers()
    worst = max(abs(r["z"]) for r in rows)
    spread = [abs(r["exact"] - r["independent_law"]) / r["independent_law"]
              for r in rows if r["independent_law"] > 0]

    tab = []
    for tag, n, beta, lam in [
        ("equal loadings, 36 observations", 36, [1.0, 1.0, 1.0], [1.0, 1.0]),
        ("unequal loadings, 36 observations", 36, [0.8, 1.3, 0.4], [4.0, 0.25]),
        ("equal loadings, 13 observations", 13, [1.0, 1.0, 1.0], [1.0, 1.0]),
        ("no loadings, 36 observations", 36, [0.0, 0.0, 0.0], [1.0, 1.0]),
    ]:
        A = AlignedFactor(n, beta, lam)
        c = A.critical(ABAR)
        tab.append({"design": tag, "n": n, "beta": [float(x) for x in beta],
                    "lam": [float(x) for x in lam], "c_alpha": float(c),
                    "fwer_at_c": float(A.fwer(c))})
        print("critical %-18s c = %.6f" % (tag, c))

    ref = AlignedFactor(36, [0.8, 1.3, 0.4], [4.0, 0.25])
    fast = AlignedFactor(36, [0.8, 1.3, 0.4], [4.0, 0.25], grid=FAST)
    refine = abs(ref.fwer(2.20) - fast.fwer(2.20))

    out = {
        "abar": ABAR,
        "model": ("r_it = mu_i + beta_i f_t + eps_it, f latent standard "
                  "normal, eps independent across i, sigma_f = sigma_0 = 1"),
        "conditioning": ["fbar", "Vf", "Xbar_0", "s_0^2"],
        "outer_dimension": 4,
        "inner_per_candidate": 1,
        "falsifiers": rows,
        "worst_abs_z": float(worst),
        "reduction_holds": bool(worst < 4.0),
        "brute_draws": 4000000,
        "max_relative_departure_from_independent_law": float(max(spread)),
        "min_relative_departure_from_independent_law": float(min(spread)),
        "critical_table": tab,
        "refinement_gap": float(refine),
        "grid": dict((k, list(v) if isinstance(v, tuple) else v)
                     for k, v in GRID.items()),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w2_factor.json", out, "exact_in_model")
    print("\nworst |z| = %.2f over %d cells  ->  %s"
          % (worst, len(rows), "HOLDS" if out["reduction_holds"] else "FAILS"))
    print("refinement gap %.2e ; %.0f s" % (refine, out["seconds"]))


if __name__ == "__main__":
    main()
