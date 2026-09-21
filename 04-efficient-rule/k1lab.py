"""The K = 1 laboratory: everything deterministic, nothing sampled.

At one candidate the whole problem is two scaled chi-squares and a normal CDF,
so the unconditional FWER of ANY threshold function is an exact 2-D quadrature
and the known-lambda oracle power has a closed form.  That makes K = 1 the right
place to settle the efficiency question before paying for general K:

  * it admits an exact uniform-FWER constraint (no Monte-Carlo slack);
  * it admits an exact power frontier;
  * it admits a genuine lower bound over a rich rule class.

WHY A LOSS IS UNAVOIDABLE (the structural statement this module quantifies).
Write rho = lambda n_0/n_1.  The conditional rejection probability of threshold c
is Phi(-c d / sigma_diff), and

    sigma_diff^2 / d^2 = (rho + 1) / (rho u_1 + u_0).

A rule attains the oracle exactly iff it can set c * d / sigma_diff to the
constant z_{1-alpha}, which needs (rho u_1 + u_0)/(rho+1).  The data supply only
rho_hat = rho u_1/u_0.  One equation, two unknowns: no threshold function of
rho_hat can equalise the conditional level across states, so some states are
over-protected and the uniform constraint binds at the least-favourable rho.
The size of that loss is what `frontier` measures.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.stats import norm

import fcore as F


# ------------------------------------------------------ population constants
def c_pop_k1(n, lam, alpha, npanel=24, order=12, lo=1e-3, hi=40.0):
    """The fixed-ratio population critical value at K = 1: P(T_1 > c) = alpha."""
    a, b = lo, hi
    for _ in range(70):
        m = 0.5 * (a + b)
        if F.fwer_exact_k1(lambda x, m=m: np.full_like(x, m), n, lam,
                           npanel, order) > alpha:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


class PluginTable:
    """c*(lam_hat) at an internal level, tabulated on a log grid and
    interpolated.  The desk computes this once per (n, level) pair."""

    def __init__(self, n, alpha, lo=1e-4, hi=1e4, m=61, npanel=20, order=10):
        self.n = np.asarray(n, float)
        self.alpha = float(alpha)
        self.x = F.log_grid(lo, hi, m)
        self.y = np.array([c_pop_k1(self.n, v, self.alpha, npanel, order)
                           for v in self.x])
        self._f = PchipInterpolator(np.log(self.x), self.y, extrapolate=True)

    def __call__(self, lh):
        lh = np.asarray(lh, float)
        return np.clip(self._f(np.log(np.clip(lh, self.x[0], self.x[-1]))),
                       1e-3, 1e3)


# ------------------------------------------------------ the uniform constraint
def sup_fwer(cfun, n, lam_grid, npanel=24, order=12):
    v = np.array([F.fwer_exact_k1(cfun, n, l, npanel, order) for l in lam_grid])
    j = int(np.argmax(v))
    return float(v[j]), float(lam_grid[j]), v


def calibrate_level(n, lam_grid, abar=F.ABAR, lo=0.002, hi=None, iters=26,
                    npanel=20, order=10, cache=None):
    """The largest internal level whose plug-in rule has sup_lambda FWER <= abar.

    This is the whole of Workstream B's first rule class: one scalar, and the
    constraint is the exact uniform one rather than a per-design calibration.
    """
    hi = float(abar if hi is None else hi)
    cache = {} if cache is None else cache

    def sup_at(a):
        if a not in cache:
            tab = PluginTable(n, a, npanel=npanel, order=order)
            cache[a] = sup_fwer(tab, n, lam_grid, npanel, order)[:2] + (tab,)
        return cache[a]

    a, b = lo, hi
    if sup_at(hi)[0] <= abar:
        s, lw, tab = sup_at(hi)
        return {"alpha0": hi, "sup_fwer": s, "lam_worst": lw, "table": tab,
                "binding": False}
    for _ in range(iters):
        m = 0.5 * (a + b)
        if sup_at(m)[0] <= abar:
            a = m
        else:
            b = m
    s, lw, tab = sup_at(a)
    return {"alpha0": a, "sup_fwer": s, "lam_worst": lw, "table": tab,
            "binding": True}


# ------------------------------------------------------------ the frontier
def frontier(n, lam_grid, alt_lams, targets=(0.5, 0.8, 0.95), abar=F.ABAR,
             npanel=24, order=12, nbins=0, cal=None):
    """Power of the calibrated plug-in rule against the known-lambda oracle.

    Reported at the effect where the ORACLE reaches each target, which is the
    only like-for-like way to price a rule in power units.
    """
    cal = cal or calibrate_level(n, lam_grid, abar, npanel=npanel, order=order)
    tab = cal["table"]
    rows = []
    for lam in alt_lams:
        for t in targets:
            dl = F.delta_for_oracle_power_k1(n, [lam], t, abar)
            p = F.power_exact_k1(tab, n, lam, dl, npanel, order)
            rows.append({"lam": float(lam), "target": float(t),
                         "delta": float(dl), "oracle_power": float(t),
                         "rule_power": float(p),
                         "power_loss_pp": float(100.0 * (t - p))})
    return {"alpha0": cal["alpha0"], "sup_fwer": cal["sup_fwer"],
            "lam_worst": cal["lam_worst"], "rows": rows,
            "max_loss_pp": max(r["power_loss_pp"] for r in rows),
            "median_loss_pp": float(np.median([r["power_loss_pp"]
                                               for r in rows]))}


# --------------------------------------------- a richer, still-valid class
class BinnedRule:
    """Piecewise-constant multiplier on the calibrated plug-in threshold.

    c(lh) = tab(lh) * exp(theta_j) for lh in bin j.  The plug-in table already
    carries the shape; theta only redistributes level across the state space,
    which is exactly the freedom the uniform constraint leaves open.
    """

    def __init__(self, tab, edges, theta=None):
        self.tab = tab
        self.edges = np.asarray(edges, float)
        self.theta = np.zeros(len(self.edges) + 1) if theta is None \
            else np.asarray(theta, float)

    def __call__(self, lh):
        lh = np.asarray(lh, float)
        j = np.searchsorted(self.edges, lh)
        return self.tab(lh) * np.exp(self.theta[j])


def optimise_binned(n, lam_grid, alt_lams, tab, edges, delta_map,
                    abar=F.ABAR, npanel=20, order=10, rounds=40, step=0.08,
                    seed=0):
    """Coordinate search for the bin multipliers that maximise mean power
    subject to the exact uniform FWER constraint.

    Deliberately a transparent local search rather than an opaque optimiser: the
    production rule has to stay inspectable, and every accepted move is one that
    both raised the objective and kept the certified-style constraint.
    """
    rng = np.random.default_rng(seed)
    rule = BinnedRule(tab, edges)

    def objective(r):
        return float(np.mean([F.power_exact_k1(r, n, lam, delta_map[lam],
                                               npanel, order)
                              for lam in alt_lams]))

    def feasible(r):
        s, _, _ = sup_fwer(r, n, lam_grid, npanel, order)
        return s <= abar, s

    ok, s0 = feasible(rule)
    if not ok:                      # shrink uniformly until feasible
        while not ok:
            rule.theta += 0.02
            ok, s0 = feasible(rule)
    best = objective(rule)
    hist = [{"round": 0, "obj": best, "sup_fwer": s0}]
    nb = len(rule.theta)
    for r in range(1, rounds + 1):
        j = int(rng.integers(nb))
        for sgn in (-1.0, +1.0):
            cand = BinnedRule(tab, edges, rule.theta.copy())
            cand.theta[j] += sgn * step
            ok, s = feasible(cand)
            if not ok:
                continue
            o = objective(cand)
            if o > best + 1e-9:
                rule, best = cand, o
                hist.append({"round": r, "bin": j, "step": sgn * step,
                             "obj": o, "sup_fwer": s})
                break
    s, lw, _ = sup_fwer(rule, n, lam_grid, npanel, order)
    return {"theta": rule.theta.tolist(), "edges": edges.tolist()
            if hasattr(edges, "tolist") else list(edges),
            "objective": best, "sup_fwer": s, "lam_worst": lw,
            "history": hist, "rule": rule}
