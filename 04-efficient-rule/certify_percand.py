"""Certify the per-candidate rule over the continuum, at ANY K, one- and
two-sided.

Two things were left open and both are closed here.  Neither needed a new idea;
they needed the enclosure worked out rather than declared unavailable.

WHY THE K-DIMENSIONAL RATIO BOX COLLAPSES TO ONE DIMENSION.  For a per-candidate
threshold, conditioning on the benchmark pair (z, u_0) makes the candidates
independent and candidate i's factor depends on lambda_i alone:

    P(no reject | z, u_0) = prod_i g_i(z, u_0; lambda_i) .

The product is increasing in each factor, so a pointwise lower bound on each
factor gives a lower bound on the product.  Partition the ratio interval into
cells, take the POINTWISE MINIMUM of the factor bound across cells, and raise it
to the power K: that is a lower bound valid for every lambda in the whole cube,
including vectors whose coordinates sit in different cells.  A K-dimensional
branch-and-bound is therefore unnecessary, which is what makes this affordable.

THE FACTOR BOUNDS.  Write x = z/sqrt(n_0), p = sqrt(n_i/lambda_i),
q = sqrt(u_i + n_i u_0/(lambda_i n_0)), so that

    one-sided factor  =  Phi(x p + c q) ,
    two-sided factor  =  Phi(x p + c q) - Phi(x p - c q) .

Both p and q are decreasing in lambda_i; q is increasing in u_i and u_0.  The
one-sided factor is increasing in its argument, so its minimum over a box takes
(x p) at its minimum - which is p at its minimum if x >= 0 and at its maximum if
x < 0 - and q at its minimum.  The two-sided factor is increasing in the first
argument and decreasing in the second, so its minimum takes the first at its
minimum and the second at its maximum; both use q at its minimum, and the two
(x p) substitutions differ.  Each substitution is monotone, which
`selftest()` falsifies rather than assumes.

OMITTED MASS IS SAFE HERE, and in the opposite direction to `wb_certify.py`.
The quantity being bounded from below is a sum of non-negative terms, so
dropping the mass outside the partition can only lower it, which can only RAISE
the resulting bound on the FWER.  The bound stays conservative.

Deposit: deposits/certify_percand.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import chi2, norm

import fcore as F
import k1lab as L
import wb_certify as C

PROB_LO = 1e-9


def _cells(nu, m):
    p = C.prob_edges(m, PROB_LO, 1.0 - PROB_LO)
    q = chi2.ppf(p, nu) / nu
    q[0] = 0.0
    return q[:-1], q[1:], np.diff(p)


def _zcells(m, zmax=8.5):
    p = C.prob_edges(m, norm.cdf(-zmax), norm.cdf(zmax))
    z = norm.ppf(p)
    return z[:-1], z[1:], np.diff(p)


def factor_lower(rule, n, i, lam_lo, lam_hi, x_lo, x_hi, u0_lo, u0_hi,
                 ui_lo, ui_hi, two_sided):
    """Rigorous lower bound on candidate i's factor over the whole box.

    Shapes broadcast: x over axis 0, u_0 over axis 1, u_i over axis 2.
    """
    ni, n0 = n[i + 1], n[0]
    p_lo = np.sqrt(ni / lam_hi)
    p_hi = np.sqrt(ni / lam_lo)
    # q is decreasing in lambda, increasing in u_i and u_0
    q_lo = np.sqrt(ui_lo + ni * u0_lo / (lam_hi * n0))
    q_hi = np.sqrt(ui_hi + ni * u0_hi / (lam_lo * n0))

    # the reachable estimated ratio on this box, and the threshold's extremes
    lh_lo = lam_lo * ui_lo / np.maximum(u0_hi, 1e-300)
    lh_hi = lam_hi * ui_hi / np.maximum(u0_lo, 1e-300)
    lh_hi = np.where(np.isfinite(lh_hi), lh_hi, 1e300)
    c_min, c_max = rule.range_min_max(lh_lo.ravel(), lh_hi.ravel())
    c_min = c_min.reshape(lh_lo.shape)
    c_max = c_max.reshape(lh_lo.shape)

    xp_lo = np.minimum(x_lo * p_lo, x_lo * p_hi)
    xp_lo = np.minimum(xp_lo, np.minimum(x_hi * p_lo, x_hi * p_hi))
    xp_hi = np.maximum(x_lo * p_lo, x_lo * p_hi)
    xp_hi = np.maximum(xp_hi, np.maximum(x_hi * p_lo, x_hi * p_hi))

    a_min = xp_lo + c_min * q_lo
    if not two_sided:
        return norm.cdf(a_min)
    b_max = xp_hi - c_min * q_lo
    return np.maximum(norm.cdf(a_min) - norm.cdf(b_max), 0.0)


def fwer_upper_percand(rule, n, lam_lo, lam_hi, n_lam=64, mz=90, m0=90,
                       mi=90, two_sided=False):
    """Rigorous upper bound on the unconditional FWER, valid for EVERY ratio
    vector in [lam_lo, lam_hi]^K."""
    n = np.asarray(n, float)
    K = n.size - 1
    zl, zh, wz = _zcells(mz)
    u0l, u0h, w0 = _cells(n[0] - 1.0, m0)
    xl = zl[:, None, None] / np.sqrt(n[0])
    xh = zh[:, None, None] / np.sqrt(n[0])
    edges = np.exp(np.linspace(np.log(lam_lo), np.log(lam_hi), n_lam + 1))

    gmin = None
    for i in range(K):
        uil, uih, wi = _cells(n[i + 1] - 1.0, mi)
        acc = None
        for a, b in zip(edges[:-1], edges[1:]):
            f = factor_lower(rule, n, i, float(a), float(b), xl, xh,
                             u0l[None, :, None], u0h[None, :, None],
                             uil[None, None, :], uih[None, None, :],
                             two_sided)
            g = f @ wi                      # (nz, n0): lower bound on g_i
            acc = g if acc is None else np.minimum(acc, g)
        gmin = acc if gmin is None else gmin * acc
    if K > 1 and gmin.ndim == 2:
        pass
    lower_no_reject = float(wz @ gmin @ w0)
    return float(1.0 - lower_no_reject)


def selftest(seed=17, trials=500, two_sided=False):
    """Falsify the monotone substitutions on random boxes and interior points."""
    rng = np.random.default_rng(seed)
    bad = 0
    for _ in range(trials):
        n0 = float(rng.integers(4, 40))
        ni = float(rng.integers(4, 40))
        lam_lo = float(np.exp(rng.uniform(np.log(1e-3), np.log(1e2))))
        lam_hi = lam_lo * float(np.exp(rng.uniform(0.0, 1.5)))
        x_lo = float(rng.uniform(-3.0, 2.0))
        x_hi = x_lo + float(rng.uniform(0.0, 1.0))
        u0l = float(np.exp(rng.uniform(-1.5, 1.0)))
        u0h = u0l * float(np.exp(rng.uniform(0.0, 0.6)))
        uil = float(np.exp(rng.uniform(-1.5, 1.0)))
        uih = uil * float(np.exp(rng.uniform(0.0, 0.6)))
        c = float(rng.uniform(1.5, 5.0))

        class Flat:
            def range_min_max(self, lo, hi):
                z = np.asarray(lo, float)
                return np.full(z.shape, c), np.full(z.shape, c)

        lb = float(np.asarray(
            factor_lower(Flat(), np.array([n0, ni]), 0, lam_lo, lam_hi,
                         np.array([[[x_lo]]]), np.array([[[x_hi]]]),
                         np.array([[[u0l]]]), np.array([[[u0h]]]),
                         np.array([[[uil]]]), np.array([[[uih]]]),
                         two_sided)).ravel()[0])
        for _ in range(6):
            lam = float(np.exp(rng.uniform(np.log(lam_lo), np.log(lam_hi))))
            x = float(rng.uniform(x_lo, x_hi))
            u0 = float(np.exp(rng.uniform(np.log(u0l), np.log(u0h))))
            ui = float(np.exp(rng.uniform(np.log(uil), np.log(uih))))
            p = np.sqrt(ni / lam)
            q = np.sqrt(ui + ni * u0 / (lam * n0))
            v = norm.cdf(x * p + c * q)
            if two_sided:
                v -= norm.cdf(x * p - c * q)
            if v < lb - 1e-12:
                bad += 1
    return {"violations": int(bad), "points": int(trials * 6),
            "two_sided": bool(two_sided), "pass": bad == 0}


DESIGNS = [("K1_df11", [12.0, 13.0]), ("K2_df11", [12.0, 13.0, 13.0]),
           ("K2_df35", [36.0, 13.0, 13.0]),
           ("K3_df35", [36.0, 13.0, 13.0, 13.0]),
           ("K5_df35", [36.0, 13.0, 13.0, 13.0, 13.0, 13.0])]


def certified_level(n, two_sided, abar=F.ABAR, iters=11):
    n = np.asarray(n, float)
    n1 = np.array([n[0], n[1]])
    edges = F.log_grid(3e-3, 3e2, 59)
    mids = np.concatenate([[edges[0] / 3.0],
                           np.sqrt(edges[:-1] * edges[1:]),
                           [edges[-1] * 3.0]])
    cache = {}

    def ub(a):
        if a not in cache:
            tab = L.PluginTable(n1, a, m=41, npanel=12, order=8)
            r = C.StepRule(tab(mids), edges).build_tables()
            cache[a] = (fwer_upper_percand(r, n, 1e-3, 1e3,
                                           two_sided=two_sided), r)
        return cache[a]

    hi = abar
    if ub(hi)[0] <= abar:
        v, r = ub(hi)
        return hi, v, r, edges
    lo = 0.0008
    for _ in range(iters):
        m = 0.5 * (lo + hi)
        if ub(m)[0] <= abar:
            lo = m
        else:
            hi = m
    v, r = ub(lo)
    return lo, v, r, edges


def main():
    t0 = time.time()
    out = {"abar": F.ABAR, "lam_domain": [1e-3, 1e3],
           "selftest_one_sided": selftest(two_sided=False),
           "selftest_two_sided": selftest(two_sided=True),
           "designs": []}
    print("falsifier, one-sided:", out["selftest_one_sided"])
    print("falsifier, two-sided:", out["selftest_two_sided"])
    print("\n%-10s %-10s %9s %11s %9s %6s" %
          ("design", "sided", "alpha0", "FWER bound", "certified", "sec"))
    for name, n in DESIGNS:
        for two in (False, True):
            t1 = time.time()
            a0, bound, rule, edges = certified_level(n, two)
            rec = {"design": name, "n": list(n), "K": len(n) - 1,
                   "two_sided": bool(two), "alpha0": float(a0),
                   "fwer_upper_bound": float(bound),
                   "certified": bool(bound <= F.ABAR),
                   "c": [float(x) for x in rule.c],
                   "edges": [float(x) for x in edges],
                   "seconds": float(time.time() - t1)}
            out["designs"].append(rec)
            print("%-10s %-10s %9.5f %11.6f %9s %6.0f" %
                  (name, "two" if two else "one", a0, bound, rec["certified"],
                   rec["seconds"]))
    out["all_certified"] = all(d["certified"] for d in out["designs"])
    out["note"] = (
        "The bound is valid for every ratio VECTOR in the cube, not only for "
        "the vectors whose coordinates share a partition cell: the pointwise "
        "minimum across cells is taken before the product over candidates, "
        "which is what removes the need for a K-dimensional branch-and-bound. "
        "Mass outside the partition is dropped, which lowers a sum of "
        "non-negative terms and therefore RAISES the FWER bound: conservative, "
        "not optimistic.")
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("certify_percand.json", out,
                    "unconditional_valid_certified")
    print("\nall certified:", out["all_certified"])


if __name__ == "__main__":
    main()
