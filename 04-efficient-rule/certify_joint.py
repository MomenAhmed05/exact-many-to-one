"""Certify the JOINT plug-in rule over the continuum, at K = 2 and K = 3.

The joint rule reads the whole estimated ratio vector, so neither the union step
of `certify_anyk.py` nor the per-candidate factorisation applies to it directly.
It is the more powerful of the two rules, so leaving it measured had a real cost.

WHAT WAS TRIED FIRST, AND WHY IT WAS ABANDONED.  The conditional no-rejection
probability is increasing in the threshold, so replacing the joint threshold on a
product cell by its minimum over the ratios that cell can reach leaves a valid
lower bound - and with the threshold constant on the cell the product factorises
again.  That works: `git` history aside, the four-dimensional enclosure is
implemented and dominates.  It is also **loose by a factor of about two** and
converges like `m^-1` in four axes at once: 0.0612 at 32 cells per axis, 0.0558
at 80, against a true value near 0.025.  Buying a tight bound that way costs
hours.

WHAT CLOSES IT INSTEAD.  Deploy the joint rule as a STEP function on a grid of
the estimated ratios - which is what a desk runs anyway, since the table is
precomputed - and project it:

    psi_i(t)  =  min over the other coordinates of the deployed table at
                 lambda_hat_i = t .

Because the deployed rule IS the table, that minimum is EXACT rather than a
grid approximation of something else.  Rejection of candidate i requires
`T_i > c(lambda_hat) >= psi_i(lambda_hat_i)`, so

    FWER_joint  <=  sum_i  FWER_1( psi_i ) ,

and each term is a one-candidate quantity that the tight two-dimensional
enclosure of `wb_certify.py` bounds over the continuum.  The whole coupling
problem is handled by one projection.

The price is the projection gap - how far the joint threshold falls when the
other candidates' ratios are pushed to their worst - and it is measured here
rather than assumed: 3.7 to 11.2 per cent of the constant on these designs
(the deposit is the source; an earlier draft of this docstring said 1.2 to 4.7
and was stale).

Deposit: deposits/certify_joint.json
"""

from __future__ import annotations

import itertools
import time

import numpy as np

import fast
import fcore as F
import k1lab as L
import wb_certify as C
import wb_percand as W

LAM_LO, LAM_HI = 1e-3, 1e3
NBIN = 30          # bins per ratio axis for the deployed joint table
CELLS = 800


def edges_for(nbin=NBIN):
    return F.log_grid(LAM_LO * 3.0, LAM_HI / 3.0, nbin - 1)


def bin_mids(edges):
    return np.concatenate([[edges[0] / 3.0],
                           np.sqrt(edges[:-1] * edges[1:]),
                           [edges[-1] * 3.0]])


class JointStepTable:
    """The deployed joint rule: a step function of the estimated ratio vector.

    A desk precomputes this once per (design, level).  Because the deployed
    threshold is the table itself, every projection taken from it is exact.
    """

    def __init__(self, n, level, nbin=NBIN, rows=None, index=None):
        self.n = np.asarray(n, float)
        self.K = self.n.size - 1
        self.level = float(level)
        self.edges = edges_for(nbin)
        self.mids = bin_mids(self.edges)
        m = self.mids.size
        if rows is None:
            rows, index = self.tabulate(self.n, nbin)
        self.rows, self.index = rows, index
        self.grid = np.empty((m,) * self.K)
        for k, idx in enumerate(index):
            self.grid[idx] = fast.invert(rows[k], self.level)

    @staticmethod
    def tabulate(n, nbin=NBIN):
        """The level-independent half: one law evaluation per cell, in
        parallel.  Every level in a bisection then costs an inversion."""
        n = np.asarray(n, float)
        K = n.size - 1
        mids = bin_mids(edges_for(nbin))
        index = list(itertools.product(range(mids.size), repeat=K))
        lams = [[float(mids[j]) for j in idx] for idx in index]
        return fast.pmap_cdf(n, lams), index

    def __call__(self, lam_hat):
        idx = tuple(int(np.searchsorted(self.edges, float(v)))
                    for v in np.asarray(lam_hat, float).ravel())
        return float(self.grid[idx])

    def project(self, i):
        """psi_i: the exact minimum of the deployed table over the others."""
        ax = tuple(a for a in range(self.K) if a != i)
        vals = self.grid.min(axis=ax) if ax else self.grid
        return C.StepRule(np.asarray(vals, float), self.edges)

    def projection_gap(self):
        """How far the projection falls below the diagonal of the table."""
        m = self.mids.size
        diag = np.array([self.grid[(j,) * self.K] for j in range(m)])
        proj = np.array([self.project(0).c[j] for j in range(m)])
        return float(np.max(diag - proj)), float(np.max((diag - proj) / diag))


def certify_level(n, nbin=NBIN, abar=F.ABAR, iters=8, cells=CELLS):
    """Largest level at which the joint step rule certifies, via the projection."""
    n = np.asarray(n, float)
    K = n.size - 1
    cache = {}

    rows, index = JointStepTable.tabulate(n, nbin)

    def ok(a):
        if a not in cache:
            tab = JointStepTable(n, a, nbin, rows=rows, index=index)
            worst, certs = -1.0, []
            good = True
            for i in range(K):
                psi = tab.project(i)
                cert = fast.certify_parallel(
                    psi, np.array([n[0], n[i + 1]]), LAM_LO, LAM_HI,
                    abar / K, m=cells)
                certs.append(cert)
                if not cert["certified"]:
                    good = False
                else:
                    worst = max(worst, cert["sup_upper_bound"])
            cache[a] = (good, worst, certs, tab)
        return cache[a]

    hi = abar
    if ok(hi)[0]:
        return hi, ok(hi)
    lo = 0.002
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if ok(mid)[0]:
            lo = mid
        else:
            hi = mid
    return lo, ok(lo)


# (name, n, bins per ratio axis).  The table costs one law evaluation per cell,
# tabulated once and inverted per level, so the bin count can be generous at
# K = 2 and still has to fall at K = 3, where the cell count is cubic.
DESIGNS = [("K2_df35", [36.0, 13.0, 13.0], 30),
           ("K2_df11", [12.0, 13.0, 13.0], 30),
           ("K3_df35", [36.0, 13.0, 13.0, 13.0], 16),
           ("K3_df11", [12.0, 13.0, 13.0, 13.0], 16)]


def main():
    t0 = time.time()
    out = {"abar": F.ABAR, "lam_domain": [LAM_LO, LAM_HI], "nbins": NBIN,
           "designs": []}
    ak = {}
    try:
        for d in F.read_deposit("certify_anyk.json")["designs"]:
            if not d["two_sided"]:
                ak[d["design"]] = d["alpha0"]
    except Exception:  # noqa: BLE001
        pass
    meas = {}
    try:
        for d in F.read_deposit("wb_general.json")["designs"]:
            meas[d["design"]] = d["alpha0"]
    except Exception:  # noqa: BLE001
        pass

    print("%-10s %10s %12s %12s %12s %7s %6s" %
          ("design", "alpha0", "union total", "vs per-cand", "vs measured",
           "gap%", "sec"))
    for name, n, nbin in DESIGNS:
        t1 = time.time()
        nn = np.asarray(n, float)
        K = nn.size - 1
        a0, (good, worst, certs, tab) = certify_level(nn, nbin=nbin)
        total = K * worst if worst > 0 else None
        gap_abs, gap_rel = tab.projection_gap()
        # the exact law at this level, as the thing the bound must dominate
        exact = max(W.fwer_percand(
            L.PluginTable(np.array([nn[0], nn[1]]), a0, m=41, npanel=12,
                          order=8), nn, p) for p in W.lam_patterns(K))
        rec = {"design": name, "n": nn.tolist(), "K": int(K),
               "alpha0_certified_joint": float(a0),
               "per_candidate_bound": float(worst),
               "union_bound_total": float(total) if total else None,
               "certified": bool(good and total is not None
                                 and total <= F.ABAR + 1e-12),
               "projection_gap_abs": gap_abs,
               "projection_gap_pct": float(100.0 * gap_rel),
               "alpha0_certified_percand": ak.get(name),
               "alpha0_measured_joint": meas.get(name),
               "level_gain_vs_percand_x":
                   float(a0 / ak[name]) if name in ak else None,
               "level_vs_measured_pct":
                   float(100.0 * (meas[name] - a0) / meas[name])
                   if name in meas else None,
               "exact_sup_over_patterns": float(exact),
               "nbins": int(nbin),
               "boxes": [c["n_boxes_accepted"] for c in certs],
               "seconds": float(time.time() - t1)}
        out["designs"].append(rec)
        print("%-10s %10.5f %12.6f %12s %12s %7.1f %6.0f" %
              (name, a0, total or float("nan"),
               ("%.2fx" % rec["level_gain_vs_percand_x"])
               if rec["level_gain_vs_percand_x"] else "-",
               ("-%.0f%%" % rec["level_vs_measured_pct"])
               if rec["level_vs_measured_pct"] else "-",
               rec["projection_gap_pct"], rec["seconds"]))
    out["all_certified"] = all(d["certified"] for d in out["designs"])
    out["note"] = (
        "The deployed rule is the tabulated joint step function, so the "
        "projection psi_i - the minimum of the table over the other "
        "candidates' ratios - is exact rather than an approximation of some "
        "other object.  Rejection of candidate i implies T_i > psi_i, so the "
        "union of the per-candidate events bounds the FWER, and each term is a "
        "one-candidate supremum certified over the continuum by wb_certify.py.")
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("certify_joint.json", out, "unconditional_valid_certified")
    print("\nall certified:", out["all_certified"])


if __name__ == "__main__":
    main()
