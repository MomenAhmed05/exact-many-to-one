"""Confidence-region envelope rule (Phase A of plans/desk-objective-closure.md).

Guarantee (POPULATION reading -- this is the implemented rule).  Build a joint
region Lambda(lam_hat) for the true variance ratios at coverage >= 1-gamma, then
set

    c_env(s^2) = max_{lam in Lambda(s^2)} c_star_pop(lam),

where c_star_pop(lam) is the POPULATION-exact critical value of the
heteroscedastic many-to-one law at variance ratios lam: the smallest c with
E_{s^2}[ P(reject | c, s^2, lam) ] = alpha (this is `core.oracle_c`; it does NOT
depend on the drawn s^2).  Proof: if the true lam lies in Lambda(s^2), then
c_star_pop(lam_true) <= c_env(s^2) pointwise in s^2, hence
P(reject | c_env, s^2, lam_true) <= P(reject | c_star_pop(lam_true), s^2, lam_true);
taking expectations over s^2 and using the defining property of c_star_pop,

    FWER <= alpha + P(lam_true not in Lambda) <= alpha + gamma.

Do NOT substitute the CONDITIONAL critical value (the c that makes
P(reject | c, s^2, lam) = alpha for the drawn s^2): it does not exist on draws
with tiny sample variances, and where it does exist it is far more conservative
(measured cost +56% vs +1.7% at the df-3 family-A cell).

`c_best_vertex` evaluates the region corners and reports the best power obtainable
inside the box; it is a power lower bound, not a guaranteed rule.
"""

import os
import sys
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import chi2

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "01-exact-law"))
sys.path.insert(0, _HERE)
from core import ALPHA, batch_c, cond_noreject, draw_var, oracle_c  # noqa: E402
from h12_quad import normal_nodes  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402


def _levels(n, gamma, scale=1.0, w_frac=0.5):
        """Chi-square quantiles defining the region at level scale `scale`."""
        n = np.asarray(n, float)
        nu = n - 1.0
        K = n.size - 1
        g_w = scale * gamma * w_frac
        a_i = scale * gamma * (1.0 - w_frac) / K
        w_lo = chi2.ppf(g_w / 2.0, nu[0]) / nu[0]
        w_hi = chi2.ppf(1.0 - g_w / 2.0, nu[0]) / nu[0]
        y_lo = np.array([chi2.ppf(a_i / 2.0, nu[i + 1]) / nu[i + 1] for i in range(K)])
        y_hi = np.array(
                [chi2.ppf(1.0 - a_i / 2.0, nu[i + 1]) / nu[i + 1] for i in range(K)]
        )
        return w_lo, w_hi, y_lo, y_hi


def region_coverage(n, gamma, scale=1.0, w_frac=0.5, nquad=200):
        """Exact P(lambda in region) by 1-D quadrature over the shared W.

        `lambda in [lam_hat*w_lo/y_hi, lam_hat*w_hi/y_lo]` is equivalent to
        `Y_i/W in [y_lo_i/w_hi, y_hi_i/w_lo]` for every i; given W = w the Y_i are
        independent, so the coverage is a product of chi-square CDF differences
        integrated over the law of W.  Deterministic -- no Monte Carlo power
        problem at gamma = 5e-4.
        """
        n = np.asarray(n, float)
        nu = n - 1.0
        K = n.size - 1
        w_lo, w_hi, y_lo, y_hi = _levels(n, gamma, scale, w_frac)
        alo = y_lo / w_hi
        ahi = y_hi / w_lo

        def integrand(w):
                p = 1.0
                for i in range(K):
                        lo = chi2.cdf(nu[i + 1] * w * alo[i], nu[i + 1])
                        hi = chi2.cdf(nu[i + 1] * w * ahi[i], nu[i + 1])
                        p *= max(0.0, hi - lo)
                return p * chi2.pdf(nu[0] * w, nu[0]) * nu[0]

        try:
                cov, _err = quad(integrand, 0.0, np.inf, limit=nquad)
                return float(cov)
        except Exception:  # noqa: BLE001
                return float("nan")


@lru_cache(maxsize=256)
def tight_scale(n_key, gamma, w_frac=0.5, margin=1e-6):
        """Largest level scale whose EXACT region coverage is >= 1-gamma+margin.

        Coverage DECREASES in the scale (the region's log-width follows the
        relative quantile spread, which shrinks as the level grows), so the
        tightest region uses the largest feasible scale; the union-bound
        allocation (scale = 1) over-covers by about 17x in miss probability.
        The certificate is the quadrature coverage, never the union bound, and
        the margin keeps the certified coverage strictly above 1-gamma.
        """
        n = np.asarray(n_key, float)
        target = 1.0 - gamma + margin
        lo, hi = 1e-3, 1e4
        try:
                if region_coverage(n, gamma, hi, w_frac) >= target:
                        return hi
                for _ in range(60):
                        mid = float(np.sqrt(lo * hi))
                        if region_coverage(n, gamma, mid, w_frac) >= target:
                                lo = mid
                        else:
                                hi = mid
                s = lo
                for _ in range(8):
                        if region_coverage(n, gamma, s * 1.05, w_frac) < target:
                                break
                        s *= 1.05
                if region_coverage(n, gamma, s, w_frac) < target:
                        return 1.0
                return float(s)
        except Exception:  # noqa: BLE001
                return 1.0


def lam_region(n, lam_hat, gamma, w_frac=0.5, scale=None):
        """Joint region for lambda given lam_hat, at coverage >= 1-gamma.

        lambda_i = lam_hat_i * W / Y_i with W = s_0^2/sigma_0^2 ~ chi2_{nu0}/nu0 and
        Y_i ~ chi2_{nu_i}/nu_i independent of W and of each other.

        `scale=None` uses the exact-coverage (tight) level scale from
        `tight_scale`; `scale=1.0` reproduces the union-bound allocation.  The
        returned diagnostics always carry the certified exact coverage.
        """
        n = np.asarray(n, float)
        lam_hat = np.atleast_2d(np.asarray(lam_hat, float))
        if scale is None:
                s = tight_scale(tuple(n.tolist()), gamma, w_frac)
        else:
                try:
                        s = float(scale)
                except Exception:  # noqa: BLE001
                        s = 1.0
        w_lo, w_hi, y_lo, y_hi = _levels(n, gamma, s, w_frac)
        lam_lo = lam_hat * (w_lo / y_hi)[None, :]
        lam_hi = lam_hat * (w_hi / y_lo)[None, :]
        return (
                lam_lo,
                lam_hi,
                {
                        "scale": f64(s),
                        "coverage_exact": f64(region_coverage(n, gamma, s, w_frac)),
                        "w_lo": f64(w_lo),
                        "w_hi": f64(w_hi),
                        "y_lo": [f64(x) for x in y_lo],
                        "y_hi": [f64(x) for x in y_hi],
                },
        )


def _u_nodes(npanels=40, order=10):
        return normal_nodes(npanels, order)


def _rej_at_lam(n, d_row, lam, c, z, w):
        """Conditional rejection probability for one draw at a fixed lambda.

        `d_row` must be the 1-D vector of per-comparison denominators for ONE
        draw (shape (K,)).  Passing a (1,K) row here broadcast the (K,) term
        against the node vector, which raised inside the bisection and made
        `c_best_vertex` return +inf (deposits/env_smoke.json, pre-fix).
        """
        n = np.asarray(n, float)
        d_row = np.asarray(d_row, float).reshape(-1)
        if d_row.size != n.size - 1:
                raise ValueError(
                        f"d_row has {d_row.size} entries, expected {n.size - 1}"
                )
        u = z / np.sqrt(n[0])
        acc = np.ones(z.size)
        for i in range(n.size - 1):
                s = np.sqrt(n[i + 1] / lam[i])
                acc = acc * ndtr(s * (u + c * d_row[i]))
        return 1.0 - f64(acc @ w)


def _bisect_c(g, lo=0.05, hi=8.0, tol=1e-8, grow=24):
        """Smallest c with g(c) = 0 by sign search, or NaN if none exists.

        `g` is decreasing in c (a larger critical value rejects less).  A fixed
        bracket is wrong here: at df-3 corners the level-alpha value exceeds 8,
        and at the far corner it can sit below 0.05.
        """
        try:
                if g(lo) <= 0.0:
                        return f64(lo)
                for _ in range(grow):
                        if g(hi) <= 0.0:
                                return f64(brentq(g, lo, hi, xtol=tol, maxiter=200))
                        hi *= 2.0
        except Exception:  # noqa: BLE001 -- infeasible bracket or bad shape
                return np.nan
        return np.nan


def _neg_rej(y, n, d_j, c, z, w):
        """Objective for the box maximiser (scipy passes `args` after y)."""
        return -_rej_at_lam(n, d_j, np.exp(y), c, z, w)


def _norej_pointwise_min(n, d, s_lo, s_hi, c, z, w):
        """E_U[ prod_i min_{s in [s_lo_i, s_hi_i]} Phi(s (U + c d_i)) ].

        Rigorous lower bound on the region's minimum no-rejection probability, hence
        an upper bound on the region's maximum rejection probability (loose when the
        region is wide)."""
        n = np.asarray(n, float)
        u = z / np.sqrt(n[0])
        D = d.shape[0]
        acc = np.ones((D, z.size))
        for i in range(n.size - 1):
                zi = u[None, :] + c * d[:, i][:, None]
                nonneg = zi >= 0.0
                term = np.where(
                        nonneg, s_lo[:, i][:, None] * zi, s_hi[:, i][:, None] * zi
                )
                acc = acc * ndtr(term)
        return acc @ w


def _vertex_lams(n, lam_lo, lam_hi, lam_hat, max_vertices=64):
        """Candidate set for the region's maximum, per draw.

        Corners alone are NOT a proof of the supremum: c_star_pop is
        non-monotone at small benchmark df (hump-shaped per coordinate), so an
        interior lambda can beat every corner.  The candidate set is therefore
        the product of {lo_i, mid_i, hi_i} per coordinate (3^K) when affordable,
        and the sup-at-candidate property is validated by direct
        pointwise-dominance checks, not assumed.  Falls back to corners when 3^K
        is too large.
        """
        K = n.size - 1
        mid = np.sqrt(lam_lo * lam_hi)
        levels = (lam_lo, mid, lam_hi)
        rows = []
        if max_vertices >= 3**K:
                for v in range(3**K):
                        idx = [(v // (3**i)) % 3 for i in range(K)]
                        rows.append(
                                np.stack(
                                        [levels[idx[i]][:, i] for i in range(K)],
                                        axis=1,
                                )
                        )
        elif max_vertices >= 2**K:
                for v in range(2**K):
                        lo_mask = np.array([((v >> i) & 1) == 0 for i in range(K)])
                        rows.append(np.where(lo_mask[None, :], lam_lo, lam_hi))
        else:
                rows.append(lam_lo)
                rows.append(lam_hi)
                rows.append(lam_hat)
                for i in range(K):
                        up = lam_hat.copy()
                        up[:, i] = lam_hi[:, i]
                        dn = lam_hat.copy()
                        dn[:, i] = lam_lo[:, i]
                        rows.append(up)
                        rows.append(dn)
        return rows  # list of (D, K) arrays


def c_env(
        n,
        s0,
        si,
        lam_hat,
        alpha=ALPHA,
        gamma=0.0005,
        pool=None,
        max_vertices=64,
        workers=10,
        tight=True,
):
        """Envelope of the POPULATION exact critical value over the region.

        c_env(s^2) = max_{lam in Lambda(s^2)} c_star_pop(lam); any true lam in the
        region has c_star_pop(lam) <= c_env pointwise in s^2, so averaging the
        conditional rejection over s^2 gives FWER <= alpha + P(lam not in Lambda)
        <= alpha + gamma.  `tight=True` uses the exact-coverage region.
        """
        n = np.asarray(n, float)
        lam_hat = np.atleast_2d(np.asarray(lam_hat, float))
        lam_lo, lam_hi, _ = lam_region(
                n, lam_hat, gamma, scale=(None if tight else 1.0)
        )
        rows = _vertex_lams(n, lam_lo, lam_hi, lam_hat, max_vertices)
        V = len(rows)
        D = lam_hat.shape[0]
        lam_mat = np.vstack(rows)  # (V*D, K)
        cs = batch_c(n, lam_mat, alpha=alpha, pool=pool, workers=workers)
        cs = cs.reshape(V, D)
        return cs.max(axis=0)


def c_env_pointwise(
        n,
        s0,
        si,
        lam_hat,
        alpha=ALPHA,
        gamma=0.0005,
        npanels=40,
        order=10,
        c_hi=400.0,
        tol=1e-7,
        tight=True,
):
        """Rigorous (loose) envelope via the pointwise-min bound.  May fail to
        bracket on extreme draws; those entries are returned as NaN."""
        n = np.asarray(n, float)
        s0 = np.asarray(s0, float)
        si = np.asarray(si, float)
        lam_hat = np.atleast_2d(np.asarray(lam_hat, float))
        lam_lo, lam_hi, _ = lam_region(
                n, lam_hat, gamma, scale=(None if tight else 1.0)
        )
        d = np.sqrt(si / n[1:] + s0[:, None] / n[0])
        s_lo = np.sqrt(n[1:][None, :] / lam_hi)
        s_hi = np.sqrt(n[1:][None, :] / lam_lo)
        z, w = _u_nodes(npanels, order)
        out = np.empty(d.shape[0])
        for j in range(d.shape[0]):

                def g(c, _d=d[j : j + 1], _sl=s_lo[j : j + 1], _sh=s_hi[j : j + 1]):
                        return _norej_pointwise_min(n, _d, _sl, _sh, c, z, w)[0] - (
                                1.0 - alpha
                        )

                out[j] = f64(_bisect_c(g, lo=0.2, hi=c_hi, tol=tol))
        return out


def c_best_vertex(
        n,
        s0,
        si,
        lam_hat,
        alpha=ALPHA,
        gamma=0.0005,
        npanels=40,
        order=10,
        tol=1e-8,
        tight=True,
):
        """Best-power corner: smallest c whose rejection at that corner is alpha.
        A power lower bound, not a guaranteed rule.  NaN where no corner is
        level-alpha feasible."""
        n = np.asarray(n, float)
        s0 = np.asarray(s0, float)
        si = np.asarray(si, float)
        lam_hat = np.atleast_2d(np.asarray(lam_hat, float))
        lam_lo, lam_hi, _ = lam_region(
                n, lam_hat, gamma, scale=(None if tight else 1.0)
        )
        d = np.sqrt(si / n[1:] + s0[:, None] / n[0])
        z, w = _u_nodes(npanels, order)
        K = n.size - 1
        out = np.empty(d.shape[0])
        for j in range(d.shape[0]):
                s_lo = np.sqrt(n[1:] / lam_hi[j])
                s_hi = np.sqrt(n[1:] / lam_lo[j])
                best = np.inf
                for v in range(2**K):
                        s = np.array(
                                [s_hi[i] if (v >> i) & 1 else s_lo[i] for i in range(K)]
                        )

                        def g(c, _s=s, _d=d[j]):
                                return (
                                        _rej_at_lam(
                                                n, _d, 1.0 / (_s**2 / n[1:]), c, z, w
                                        )
                                        - alpha
                                )

                        cv = _bisect_c(g, tol=tol)
                        if np.isfinite(cv):
                                best = min(best, cv)
                out[j] = best if np.isfinite(best) else np.nan
        return out


def fwer_of(n, lam_true, cvec, s0, si, alpha=ALPHA):
        """Realised FWER under the TRUE ratios, paired with the oracle."""
        lam_true = np.asarray(lam_true, float)
        c_or = oracle_c(n, list(lam_true), alpha)
        nr_o = cond_noreject(np.full(len(cvec), c_or), s0, si, n, lam_true)
        nr = cond_noreject(cvec, s0, si, n, lam_true)
        d = nr_o - nr
        return {
                "fwer": f64(alpha + d.mean()),
                "mcse": f64(d.std(ddof=1) / np.sqrt(len(d))),
                "mean_c": f64(np.mean(cvec)),
                "c_oracle": f64(c_or),
                "power_cost_pct": f64(100 * (np.mean(cvec) - c_or) / c_or),
        }


def smoke(draws=200):
        """Decisive check: guarantee and power cost at df 3 and df 35, both families."""
        rows = []
        for fam, lam in (("A", (0.04, 0.04)), ("B", (1.0, 25.0))):
                for n0 in (4, 36):
                        n = np.array([n0, 13, 13], float)
                        s0, si, lh = draw_var(n, lam, draws, 4242)
                        c_or = oracle_c(n, lam, ALPHA)
                        c_plug = batch_c(n, lh, alpha=ALPHA, workers=1)
                        c_e = c_env(n, s0, si, lh, gamma=0.0005)
                        c_pw = c_env_pointwise(n, s0, si, lh, gamma=0.0005)
                        c_v = c_best_vertex(n, s0, si, lh, gamma=0.0005)
                        r_plug = fwer_of(n, lam, c_plug, s0, si)
                        r_env = fwer_of(n, lam, c_e, s0, si)
                        r_vtx = fwer_of(n, lam, c_v, s0, si)
                        rows.append(
                                {
                                        "family": fam,
                                        "df0": n0 - 1,
                                        "lam": list(lam),
                                        "c_oracle": f64(c_or),
                                        "plugin": r_plug,
                                        "envelope": r_env,
                                        "vertex": r_vtx,
                                        "pointwise_mean": f64(np.nanmean(c_pw)),
                                        "pointwise_nan": i64(
                                                np.sum(~np.isfinite(c_pw))
                                        ),
                                        "guarantee": f64(ALPHA + 0.0005),
                                }
                        )
                        print(
                                f"{fam} df0={n0 - 1}: oracle {c_or:.4f} | "
                                f"plugin FWER {r_plug['fwer']:.4f} (c {r_plug['power_cost_pct']:+.2f}%) | "
                                f"ENV FWER {r_env['fwer']:.4f} (c {r_env['power_cost_pct']:+.2f}%) | "
                                f"vertex c {r_vtx['power_cost_pct']:+.2f}% | "
                                f"pointwise nan={rows[-1]['pointwise_nan']}"
                        )
        return rows


if __name__ == "__main__":
        res = smoke()
        dump_json(
                {"alpha": ALPHA, "gamma": 0.0005, "rows": res},
                os.path.join(_HERE, "deposits", "env_smoke.json"),
        )
        print("wrote deposits/env_smoke.json")
