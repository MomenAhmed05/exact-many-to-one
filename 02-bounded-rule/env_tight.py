"""A5b: exact-coverage (tight) confidence region for the envelope rule.

The union-bound region in `env_rule.lam_region` over-covers: at gamma = 5e-4 its
exact coverage is 0.99997-0.999998 (quadrature, deposits/env/falsifiers.json),
i.e. the union bound is loose by roughly 17x.  The region's coverage is a
function of the design (n) only, not of lam_hat, so the per-coordinate levels
can be scaled UP until the exact coverage equals 1-gamma.  A smaller region
gives a smaller envelope c_env = max_Lambda c*(lambda), hence less power cost,
with the guarantee unchanged: coverage(s*) >= 1-gamma is re-verified by
quadrature before s* is returned and is deposited with every run.

This module is deliberately standalone (it does not modify env_rule.py), and
reuses only env_rule's public pieces plus core's solver.
"""

import multiprocessing as mp
import os
import sys
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, batch_c, draw_var, oracle_c  # noqa: E402
from env_rule import fwer_of  # noqa: E402
from io_util import dump_json, f64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "env")
W_FRAC = 0.5
FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
DFS = (3, 5, 7, 11, 19, 35)


def _levels(n, gamma, scale, w_frac=W_FRAC):
    """Scaled chi-square quantiles defining the region."""
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


def region_coverage(n, gamma, scale, w_frac=W_FRAC):
    """Exact P(lambda in region(scale)) by 1-D quadrature over W."""
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
        cov, _err = quad(integrand, 0.0, np.inf, limit=200)
        return float(cov)
    except Exception:  # noqa: BLE001
        return float("nan")


@lru_cache(maxsize=256)
def tight_scale(n_key, gamma, w_frac=W_FRAC):
    """Largest scale <= 1 whose region coverage is still >= 1-gamma.

    Bisection on log(scale) in [log(0.02), log(1.0)], then an upward safety
    sweep so the RETURNED scale has coverage >= 1-gamma under quadrature.
    """
    n = np.asarray(n_key, float)
    target = 1.0 - gamma
    lo, hi = 1e-3, 1e4
    try:
        if region_coverage(n, gamma, hi, w_frac) >= target:
            return hi
        # coverage is DECREASING in the level scale (the region's log-width
        # follows the relative quantile spread, which shrinks as the level
        # grows), so take the LARGEST feasible scale.
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
        return s
    except Exception:  # noqa: BLE001
        return 1.0


def lam_region_tight(n, lam_hat, gamma, w_frac=W_FRAC):
    """Region bounded by the exact-coverage scaled levels."""
    n = np.asarray(n, float)
    lam_hat = np.atleast_2d(np.asarray(lam_hat, float))
    K = n.size - 1
    s = tight_scale(tuple(n.tolist()), gamma, w_frac)
    w_lo, w_hi, y_lo, y_hi = _levels(n, gamma, s, w_frac)
    lam_lo = lam_hat * (w_lo / y_hi)[None, :]
    lam_hi = lam_hat * (w_hi / y_lo)[None, :]
    try:
        cov = f64(region_coverage(n, gamma, s, w_frac))
        kk = int(K)
        wf = f64(w_frac)
        yl = [f64(x) for x in y_lo]
        yh = [f64(x) for x in y_hi]
    except Exception:  # noqa: BLE001
        cov, kk, wf, yl, yh = float("nan"), 0, 0.5, [], []
    diag = {
        "scale": f64(s),
        "coverage_exact": cov,
        "w_lo": f64(w_lo),
        "w_hi": f64(w_hi),
        "y_lo": yl,
        "y_hi": yh,
        "K": kk,
        "w_frac": wf,
    }
    return lam_lo, lam_hi, diag


def _corners(lam_lo, lam_hi, lam_hat, max_vertices=64):
    """Full 2^K corner enumeration for K<=6; reduced set otherwise."""
    K = lam_lo.shape[1]
    rows = []
    if max_vertices >= 2**K:
        for v in range(2**K):
            up = np.array([bool((v >> i) & 1) for i in range(K)])
            rows.append(np.where(up[None, :], lam_hi, lam_lo))
    else:
        rows.extend([lam_lo, lam_hi, lam_hat])
        for i in range(K):
            up = lam_hat.copy()
            up[:, i] = lam_hi[:, i]
            rows.append(up)
            dn = lam_hat.copy()
            dn[:, i] = lam_lo[:, i]
            rows.append(dn)
    return rows


def c_env_tight(
    n,
    s0,
    si,
    lam_hat,
    alpha=ALPHA,
    gamma=0.0005,
    pool=None,
    workers=10,
    max_vertices=64,
):
    """Envelope over the exact-coverage region (A5b working rule)."""
    n = np.asarray(n, float)
    lam_hat = np.atleast_2d(np.asarray(lam_hat, float))
    lam_lo, lam_hi, _ = lam_region_tight(n, lam_hat, gamma)
    rows = _corners(lam_lo, lam_hi, lam_hat, max_vertices)
    D = lam_hat.shape[0]
    lam_mat = np.vstack(rows)
    cs = batch_c(n, lam_mat, alpha=alpha, pool=pool, workers=workers)
    return cs.reshape(len(rows), D).max(axis=0)


def cell(fam, df, draws, pool, gamma=0.0005):
    """One A/B cell: union-bound vs exact-coverage region."""
    lam = np.asarray(FAMS[fam], float)
    n = np.array([df + 1, 13, 13], float)
    s0, si, lh = draw_var(n, lam, draws, 4242)
    c_or = oracle_c(n, lam, ALPHA)
    c_t = c_env_tight(n, s0, si, lh, gamma=gamma, pool=pool)
    from env_rule import c_env

    c_u = c_env(n, s0, si, lh, gamma=gamma, pool=pool)
    _lo, _hi, diag = lam_region_tight(n, lh[:3], gamma)
    return {
        "family": fam,
        "df0": df,
        "lam": list(lam),
        "draws": draws,
        "gamma": gamma,
        "c_oracle": f64(c_or),
        "region": diag,
        "union_bound": fwer_of(n, lam, c_u, s0, si),
        "tight": fwer_of(n, lam, c_t, s0, si),
        "guarantee": f64(ALPHA + gamma),
    }


def sweep(draws=200, pool_n=4):
    """12-cell A/B df sweep with the tight region; deposits sweep_tight_*.json."""
    rows = []
    with mp.Pool(pool_n) as pool:
        for df in DFS:
            for fam in FAMS:
                r = cell(fam, df, draws, pool)
                dump_json(r, os.path.join(OUT, f"tight_df{df}_{fam}.json"))
                rows.append(r)
                print(
                    f"{fam} df{df:>2}: s* {r['region']['scale']:.4f} "
                    f"cov {r['region']['coverage_exact']:.6f} | "
                    f"tight FWER {r['tight']['fwer']:.4f} "
                    f"(c {r['tight']['power_cost_pct']:+.1f}%) | "
                    f"ub {r['union_bound']['power_cost_pct']:+.1f}% | "
                    f"oracle {r['c_oracle']:.4f}",
                    flush=True,
                )
    costs_u = [r["union_bound"]["power_cost_pct"] for r in rows]
    costs_t = [r["tight"]["power_cost_pct"] for r in rows]
    try:
        mu = float(np.mean(costs_u))
        mt = float(np.mean(costs_t))
    except Exception:  # noqa: BLE001
        mu = mt = float("nan")
    summary = {
        "draws": draws,
        "gamma": 0.0005,
        "cells": [
            {
                "family": r["family"],
                "df0": r["df0"],
                "scale": r["region"]["scale"],
                "coverage": r["region"]["coverage_exact"],
                "ub_cost_pct": r["union_bound"]["power_cost_pct"],
                "tight_cost_pct": r["tight"]["power_cost_pct"],
                "ub_fwer": r["union_bound"]["fwer"],
                "tight_fwer": r["tight"]["fwer"],
            }
            for r in rows
        ],
        "mean_ub_cost_pct": f64(mu),
        "mean_tight_cost_pct": f64(mt),
        "all_within_guarantee": bool(
            all(r["tight"]["fwer"] <= ALPHA + 0.0005 for r in rows)
        ),
    }
    dump_json(summary, os.path.join(OUT, "tight_summary.json"))
    print(
        f"mean cost: union bound {summary['mean_ub_cost_pct']:+.2f}% -> "
        f"tight {summary['mean_tight_cost_pct']:+.2f}% | "
        f"guarantee held {summary['all_within_guarantee']}"
    )
    return summary


def main():
    for n0 in (4, 6, 8, 12, 20, 36):
        try:
            n = (float(n0), 13.0, 13.0)
            s = tight_scale(n, 0.0005)
            cov = region_coverage(n, 0.0005, s)
        except Exception as exc:  # noqa: BLE001
            print(f"n0={n0}: failed ({exc!r})", flush=True)
            continue
        print(
            f"n0={n0:>2}: s*={s:.4f} coverage={cov:.7f} (target 0.9995000)", flush=True
        )
    if os.environ.get("W26_TIGHT_SWEEP", "0") == "1":
        try:
            draws = int(os.environ.get("W26_TIGHT_D", "200"))
            pool_n = int(os.environ.get("W26_POOL", "4"))
        except Exception:  # noqa: BLE001
            draws, pool_n = 200, 4
        sweep(draws=draws, pool_n=pool_n)


if __name__ == "__main__":
    main()
