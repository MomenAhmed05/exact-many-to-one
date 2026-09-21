"""Certify the uniform FWER constraint over the CONTINUUM of ratios (K = 1).

`wb_k1.py` and `wc_frontier.py` check the constraint on a grid of lambda.  A grid
is not a proof.  This module replaces it with a rigorous two-sided bound, so the
rule can carry the theorem tier rather than the measured one:

    sup_{lambda in [L, U]}  FWER(lambda)  <=  abar .

THE BOUND.  For a threshold function c(.) that is a step function of the
estimated ratio,

    FWER(lambda) = E_{u0,u1}[ Phi( -c(lambda u1/u0) * d / sigma ) ],
    d^2 = lambda u1/n1 + u0/n0,   sigma^2 = lambda/n1 + 1/n0 .

Partition the state space into cells on the PROBABILITY scale, so each cell's
mass is an exact difference of chi-square distribution functions and the masses
sum to one with no tail term.  On a cell, and over a lambda box, the integrand is
maximised by making its argument as small as possible, which means

    c as small as possible   -> the minimum step value over the cell's
                                reachable lambda_hat range,
    d as small as possible   -> lambda_lo, u0_lo, u1_lo,
    sigma as large as possible -> lambda_hi .

Every one of those is a monotone substitution, so

    FWER(lambda) <= sum_cells  mass * Phi( -c_min * d_min / sigma_max )

holds for every lambda in the box simultaneously.  The same construction with
the inequalities reversed gives a lower bound, and the two bracket the truth, so
the enclosure can be checked for tightness rather than trusted.

The monotonicity claims are not assumed: `selftest()` falsifies them against a
direct evaluation on random cells, and the bracket is checked to contain the
grid-quadrature value at every lambda tested.

Deposit: deposits/wb_certify.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import chi2, norm

import fcore as F


PROB_LO = 1e-9


def prob_edges(m, lo=PROB_LO, hi=1.0 - PROB_LO):
    """Probability-scale cell edges, refined towards both extremes."""
    t = np.linspace(0.0, 1.0, int(m) + 1)
    # cluster towards 0 and 1 where the chi-square quantile moves fastest
    s = 0.5 * (1.0 - np.cos(np.pi * t))
    return lo + (hi - lo) * s


def state_cells(nu, m):
    p = prob_edges(m)
    q = chi2.ppf(p, nu) / nu
    q[0] = 0.0
    mass = np.diff(p)
    return q, mass


class StepRule:
    def __init__(self, c, edges):
        self.c = np.asarray(c, float)
        self.edges = np.asarray(edges, float)

    def __call__(self, lh):
        return self.c[np.searchsorted(self.edges, np.asarray(lh, float))]

    def range_min_max(self, lo, hi):
        """min and max of c over lambda_hat in [lo, hi], exactly."""
        i = np.searchsorted(self.edges, lo)
        j = np.searchsorted(self.edges, hi)
        i = np.clip(i, 0, self.c.size - 1)
        j = np.clip(j, 0, self.c.size - 1)
        lo_i = np.minimum(i, j)
        hi_i = np.maximum(i, j)
        # cumulative running extrema make the per-cell reduction O(1)
        return self._cmin[lo_i, hi_i], self._cmax[lo_i, hi_i]

    def build_tables(self):
        M = self.c.size
        self._cmin = np.empty((M, M))
        self._cmax = np.empty((M, M))
        for a in range(M):
            run_min, run_max = np.inf, -np.inf
            for b in range(a, M):
                run_min = min(run_min, self.c[b])
                run_max = max(run_max, self.c[b])
                self._cmin[a, b] = run_min
                self._cmax[a, b] = run_max
        return self


def enclose_fwer_box(rule, n, lam_lo, lam_hi, m=320):
    """Rigorous [lower, upper] bracket on FWER(lambda) for every lambda in the box."""
    n = np.asarray(n, float)
    q0, w0 = state_cells(n[0] - 1.0, m)
    q1, w1 = state_cells(n[1] - 1.0, m)
    u0_lo, u0_hi = q0[:-1], q0[1:]
    u1_lo, u1_hi = q1[:-1], q1[1:]
    W = np.outer(w0, w1)

    # lambda_hat reachable range on each cell, over the whole lambda box
    lh_lo = lam_lo * u1_lo[None, :] / np.maximum(u0_hi[:, None], 1e-300)
    lh_hi = lam_hi * u1_hi[None, :] / np.maximum(u0_lo[:, None], 1e-300)
    lh_hi = np.where(np.isfinite(lh_hi), lh_hi, 1e300)
    cmin, cmax = rule.range_min_max(lh_lo.ravel(), lh_hi.ravel())
    cmin = cmin.reshape(lh_lo.shape)
    cmax = cmax.reshape(lh_lo.shape)

    d_lo = np.sqrt(lam_lo * u1_lo[None, :] / n[1] + u0_lo[:, None] / n[0])
    d_hi = np.sqrt(lam_hi * u1_hi[None, :] / n[1] + u0_hi[:, None] / n[0])
    sig_lo = np.sqrt(lam_lo / n[1] + 1.0 / n[0])
    sig_hi = np.sqrt(lam_hi / n[1] + 1.0 / n[0])

    # The probability-scale partition spans [PROB_LO, 1 - PROB_LO] on each
    # axis, so a little mass is outside the cells.  The integrand is a
    # probability, so that mass can contribute at most its own size to the
    # truth: add it to the upper bound rather than assume it away.  With the
    # default edges it is 4e-9, which is small but is NOT nothing next to the
    # tightest certificate margin in the deposit (1.20e-8, on df35_short_cand).
    omitted = 4.0 * PROB_LO
    upper = float((W * norm.cdf(-(cmin * d_lo) / sig_hi)).sum()) + omitted
    lower = float((W * norm.cdf(-(cmax * d_hi) / sig_lo)).sum())
    return lower, upper


def certify(rule, n, lam_lo, lam_hi, abar=F.ABAR, m=320, max_boxes=4096,
            tol=None, init_split=96):
    """Branch and bound over the ratio interval until every box is decided.

    A box is accepted when its rigorous upper bound is at or below the ceiling,
    and split otherwise.  The certificate is the list of accepted boxes; if the
    budget is exhausted with boxes undecided, the routine says so rather than
    reporting a pass.
    """
    t0 = time.time()
    rule.build_tables()
    # Pre-split: the corner substitution's width grows with the box, so one
    # evaluation on the whole interval is nearly vacuous and every split after
    # it is wasted work.  Starting from a geometric partition makes each first
    # bound tight enough to decide most boxes immediately.
    e = np.exp(np.linspace(np.log(float(lam_lo)), np.log(float(lam_hi)),
                           int(init_split) + 1))
    boxes = [(float(a), float(b)) for a, b in zip(e[:-1], e[1:])]
    accepted, undecided = [], []
    worst_upper, n_eval = -1.0, 0
    while boxes and len(accepted) + len(boxes) < max_boxes:
        lo, hi = boxes.pop(0)
        low, up = enclose_fwer_box(rule, n, lo, hi, m)
        n_eval += 1
        worst_upper = max(worst_upper, up)
        if up <= abar:
            accepted.append({"lo": lo, "hi": hi, "upper": up, "lower": low})
        elif hi / lo < 1.0 + 1e-6:
            undecided.append({"lo": lo, "hi": hi, "upper": up, "lower": low})
        else:
            mid = float(np.sqrt(lo * hi))
            boxes.extend([(lo, mid), (mid, hi)])
    undecided.extend({"lo": a, "hi": b, "upper": None, "lower": None}
                     for a, b in boxes)
    ok = not undecided
    sup_upper = max((a["upper"] for a in accepted), default=float("nan"))
    return {
        "certified": bool(ok),
        "n_boxes_accepted": len(accepted),
        "n_boxes_undecided": len(undecided),
        "sup_upper_bound": float(sup_upper) if accepted else None,
        "worst_upper_seen": float(worst_upper),
        "max_bracket_width": float(max((a["upper"] - a["lower"]
                                        for a in accepted), default=0.0)),
        "undecided": undecided[:8],
        "cells_per_axis": int(m),
        "n_box_evaluations": int(n_eval),
        "seconds": float(time.time() - t0),
    }


# --------------------------------------------------------------- falsifiers
def selftest(seed=7, trials=400):
    """Falsify the monotone substitutions and the bracket.

    X2: on random cells and random lambda boxes, a direct evaluation of the
        integrand at random interior points must lie inside [lower, upper]
        corner substitution used by the bound.
    X3: the rigorous bracket must contain the quadrature value of fcore.
    """
    rng = np.random.default_rng(seed)
    bad_pt = 0
    for _ in range(trials):
        n = np.array([float(rng.integers(4, 40)), float(rng.integers(4, 40))])
        lo = float(np.exp(rng.uniform(np.log(1e-3), np.log(1e2))))
        hi = lo * float(np.exp(rng.uniform(0.0, 2.0)))
        u0l = float(np.exp(rng.uniform(-2, 1)))
        u0h = u0l * float(np.exp(rng.uniform(0, 0.7)))
        u1l = float(np.exp(rng.uniform(-2, 1)))
        u1h = u1l * float(np.exp(rng.uniform(0, 0.7)))
        c = float(rng.uniform(1.5, 5.0))
        dl = np.sqrt(lo * u1l / n[1] + u0l / n[0])
        dh = np.sqrt(hi * u1h / n[1] + u0h / n[0])
        sl = np.sqrt(lo / n[1] + 1.0 / n[0])
        sh = np.sqrt(hi / n[1] + 1.0 / n[0])
        up = norm.cdf(-(c * dl) / sh)
        low = norm.cdf(-(c * dh) / sl)
        for _ in range(6):
            lam = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            u0 = float(np.exp(rng.uniform(np.log(u0l), np.log(u0h))))
            u1 = float(np.exp(rng.uniform(np.log(u1l), np.log(u1h))))
            d = np.sqrt(lam * u1 / n[1] + u0 / n[0])
            s = np.sqrt(lam / n[1] + 1.0 / n[0])
            v = norm.cdf(-(c * d) / s)
            if not (low - 1e-12 <= v <= up + 1e-12):
                bad_pt += 1
    return {"x2_corner_substitution_violations": int(bad_pt),
            "x2_trials": int(trials * 6)}


def bracket_check(rule, n, lams, m=320):
    """X3: the rigorous bracket at a degenerate box must contain the quadrature."""
    rule.build_tables()
    out = []
    for lam in lams:
        lo, up = enclose_fwer_box(rule, n, lam, lam, m)
        q = F.fwer_exact_k1(rule, n, lam, 24, 12)
        out.append({"lam": float(lam), "lower": float(lo), "upper": float(up),
                    "quadrature": float(q),
                    "contained": bool(lo - 1e-9 <= q <= up + 1e-9),
                    "width": float(up - lo)})
    return out
