"""The falsifier for the Hunt-Stein claim: let a rule cheat and see if it can.

A theorem cannot be falsified by simulation.  Its CONCLUSION can, and that is
what this does.

THE CLAIM UNDER TEST.  Section 10.8 bounds every rule measurable in
`(T, lambda_hat)` by a linear program.  The proposition extends that bound to
every test, by an invariance argument.  If the argument is being misapplied, a
rule that reads something OUTSIDE `(T, lambda_hat)` should be able to beat the
bound while still meeting the uniform size constraint.

WHAT THE EXTRA COORDINATE IS.  The sufficient statistic is
`(Xbar_0, Xbar_1, s_0^2, s_1^2)`, which is in bijection with
`(T, lambda_hat, Xbar_0, s_0^2)`.  So `Xbar_0` and `s_0^2` are exactly the two
coordinates an invariant rule throws away, one for each generator of the group.
Each is given back in turn:

    scale probe     the rule may read `s_0^2` as well as `(T, lambda_hat)`,
                    and the size constraint runs over a range of sigma_0^2;
    location probe  the rule may read `Xbar_0` as well, and the size constraint
                    runs over a range of the common mean mu.

THE CRITERION HAS TO BE THE ONE THE PROPOSITION COVERS, and the first version
of this module got that wrong in a way worth recording.  It maximised weighted
average power at ONE value of the nuisance and found gains of up to 13.9 points.
Those gains are not evidence against the Proposition.  They are evidence for
what Section 4 of `HUNT-STEIN.md` already says: at a single point of an orbit
the discarded coordinate really does carry information, and average power at a
chosen representative is not a criterion the invariance argument covers.

The Proposition is about a SUPREMUM over the orbit, so the probe takes its
objective over the orbit too:

    minimise   t
    subject to t >= beta* - power(q ; alternative a, nuisance v)
                                        for every a AND every v
               sum q mass(lambda, v, 0) <= abar    for every lambda AND v
               q nondecreasing in T,   0 <= q <= 1

Both the constraint and the objective now range over the orbit coordinate, and
the Proposition predicts the free optimum equals the flat one.  Under the null
sigma_0^2 and mu are free parameters, so a rule that leans on `s_0^2` or
`Xbar_0` must hold its size at every value of them and cannot buy worst-case
regret with a coordinate that only moves along the orbit.

THE LAW HAS TO BE THE JOINT ONE.  `Xbar_0` is NOT independent of `T`; it sits
inside `T`'s numerator, and a first draft of this module that treated them as
independent would have made the probe vacuous by construction.  Both are
therefore built from the same decomposition.  Write `A = Xbar_0 - mu` and
`B = Xbar_1 - mu - delta`, independent normals given the state.  Then

    Xbar_0 = mu + A ,        T = (B + delta - A) / d ,

so conditionally on the state and on `A` the event `T <= t` is
`B <= t d + A - delta`, which is a single normal tail.  One quadrature over `A`
therefore gives the exact joint mass in every (ratio bin, T cell, extra bin)
without a bivariate integral anywhere.

Deposit: deposits/w1_noninvariant.json
"""

from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from scipy.stats import norm

import fcore as F
from h12_quad import normal_nodes
from w1_region import ALT_LAMS, LAM_DOMAIN

ABAR = F.ABAR
ORACLE = 0.80

# The first version of this module scored weighted average power at ONE point
# of the orbit and reported this gain.  It is kept as a deposited constant so
# that the write-up can quote it without typing it, and because a falsifier
# that was itself wrong is the thing this repository records rather than
# quietly replaces.
SUPERSEDED_FIRST_VERSION_GAIN_PP = 13.868


class Plane3(object):
    """The (lambda_hat bin, T cell, extra bin) mesh and its exact masses."""

    def __init__(self, n, kind, nbins=20, nt=120, nx=8,
                 tlo=-7.0, thi=12.0, npanel=6, order=8, apanel=3, aorder=8):
        assert kind in ("scale", "location")
        self.kind = kind
        self.n = np.asarray(n, float)
        self.edges = F.log_grid(LAM_DOMAIN[0] * 3, LAM_DOMAIN[1] / 3, nbins - 1)
        self.M = self.edges.size + 1
        self.t = np.linspace(tlo, thi, int(nt) + 1)
        self.NT = int(nt) + 2

        q0, w0 = F.state_nodes(self.n[0] - 1.0, npanel, order)
        q1, w1 = F.state_nodes(self.n[1] - 1.0, npanel, order)
        U0, U1 = np.meshgrid(q0, q1, indexing="ij")
        self.u0, self.u1 = U0.ravel(), U1.ravel()
        W = np.outer(w0, w1).ravel()
        self.W = W / W.sum()

        # the quadrature over A = Xbar_0 - mu, in standard units
        self.za, self.wa = normal_nodes(apanel, aorder)

        self.NX = int(nx)
        if kind == "scale":
            self.xcut = np.log(F.log_grid(0.05, 20.0, self.NX - 1))
        else:
            self.xcut = np.linspace(-3.0, 3.0, self.NX - 1)
        self.nvar = self.M * self.NT * self.NX
        self._mono = None

    def mass(self, lam, nuis, delta=0.0):
        """Exact mass in every (ratio bin, T cell, extra bin), flattened.

        `nuis` is sigma_0^2 for the scale probe and mu for the location probe.
        """
        n0, n1 = self.n[0], self.n[1]
        v = nuis if self.kind == "scale" else 1.0
        mu = nuis if self.kind == "location" else 0.0

        s0 = v * self.u0
        s1 = v * lam * self.u1
        d = np.sqrt(s1 / n1 + s0 / n0)
        b = np.searchsorted(self.edges, lam * self.u1 / self.u0)

        sa = np.sqrt(v / n0)                 # sd of A
        sb = np.sqrt(v * lam / n1)           # sd of B
        A = sa * self.za                     # (nA,)

        if self.kind == "scale":
            ex = np.searchsorted(self.xcut, np.log(s0))          # per state
        else:
            ex = np.searchsorted(self.xcut, mu + A)              # per A node

        out = np.zeros((self.M, self.NT, self.NX))
        for ai, (a, wa) in enumerate(zip(A, self.wa)):
            # P(T <= t | state, A) = Phi((t d + A - delta) / sb)
            x = (self.t[None, :] * d[:, None] + a - delta) / sb
            cdf = norm.cdf(x)
            cells = np.empty((d.size, self.NT))
            cells[:, 0] = cdf[:, 0]
            cells[:, 1:-1] = np.diff(cdf, axis=1)
            cells[:, -1] = 1.0 - cdf[:, -1]
            cells *= (self.W * wa)[:, None]
            if self.kind == "scale":
                for j in range(self.NX):
                    m = ex == j
                    if m.any():
                        acc = np.zeros((self.M, self.NT))
                        np.add.at(acc, b[m], cells[m])
                        out[:, :, j] += acc
            else:
                acc = np.zeros((self.M, self.NT))
                np.add.at(acc, b, cells)
                out[:, :, ex[ai]] += acc
        return out.ravel()

    def monotone_rows(self):
        """q nondecreasing in the T cell, within each (ratio bin, extra bin)."""
        if self._mono is None:
            idx = np.arange(self.nvar).reshape(self.M, self.NT, self.NX)
            a = idx[:, :-1, :].ravel()
            b = idx[:, 1:, :].ravel()
            r = np.arange(a.size)
            self._mono = sp.csr_matrix(
                (np.concatenate([np.ones(a.size), -np.ones(b.size)]),
                 (np.concatenate([r, r]), np.concatenate([a, b]))),
                shape=(a.size, self.nvar))
        return self._mono

    def flat_rows(self):
        """q equal across the extra axis: the restriction to invariant rules."""
        idx = np.arange(self.nvar).reshape(self.M, self.NT, self.NX)
        a = idx[:, :, :-1].ravel()
        b = idx[:, :, 1:].ravel()
        r = np.arange(a.size)
        return sp.csr_matrix(
            (np.concatenate([np.ones(a.size), -np.ones(b.size)]),
             (np.concatenate([r, r]), np.concatenate([a, b]))),
            shape=(a.size, self.nvar))


def solve(P, alts, lam_cons, nuis_grid, flat, abar=ABAR, oracle=ORACLE):
    """Minimise the worst regret over alternatives AND over the orbit.

    Variables are (q, t) with t the epigraph variable.  A lower value is a
    better rule, so the free arm beats the flat one when its value is SMALLER.
    """
    nq = P.nvar
    rows, rhs = [], []
    for lam, delta in alts:
        for v in nuis_grid:
            m = P.mass(lam, v, delta)
            rows.append(np.concatenate([-m, [-1.0]]))     # -power - t <= -beta*
            rhs.append(-oracle)
    n_obj = len(rows)
    for l in lam_cons:
        for v in nuis_grid:
            rows.append(np.concatenate([P.mass(l, v, 0.0), [0.0]]))
            rhs.append(abar)
    A = sp.csr_matrix(np.stack(rows))
    b = np.array(rhs, float)

    Mr = P.monotone_rows()
    A = sp.vstack([A, sp.hstack([Mr, sp.csr_matrix((Mr.shape[0], 1))])],
                  format="csr")
    b = np.concatenate([b, np.zeros(Mr.shape[0])])

    kw = {}
    if flat:
        E = P.flat_rows()
        kw = {"A_eq": sp.hstack([E, sp.csr_matrix((E.shape[0], 1))]),
              "b_eq": np.zeros(E.shape[0])}
    c = np.zeros(nq + 1)
    c[-1] = 1.0
    res = linprog(c, A_ub=A, b_ub=b,
                  bounds=[(0.0, 1.0)] * nq + [(None, None)],
                  method="highs", **kw)
    q = np.clip(res.x[:nq], 0.0, 1.0)
    return {"q": q, "value": float(res.x[-1]), "status": int(res.status),
            "n_objective_rows": n_obj}


def probe(tag, n, kind, nuis_grid, ncons=21):
    t0 = time.time()
    n = np.asarray(n, float)
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], ORACLE)) for l in ALT_LAMS]
    lam_cons = F.log_grid(*LAM_DOMAIN, ncons)
    P = Plane3(n, kind)
    free = solve(P, alts, lam_cons, nuis_grid, flat=False)
    flat = solve(P, alts, lam_cons, nuis_grid, flat=True)
    # a LOWER worst regret is better, so the free arm wins when its value is
    # smaller, and the gain is the flat value minus the free one
    gain = 100.0 * (flat["value"] - free["value"])
    q = free["q"].reshape(P.M, P.NT, P.NX)
    spread = float(np.max(q.max(2) - q.min(2)))

    row = {"design": tag, "n": n.tolist(), "kind": kind,
           "extra_coordinate": "s_0^2" if kind == "scale" else "Xbar_0",
           "nuisance": "sigma_0^2" if kind == "scale" else "mu",
           "nuisance_grid": [float(x) for x in nuis_grid],
           "extra_bins": P.NX, "ratio_bins": P.M, "t_cells": P.NT,
           "free_worst_regret": free["value"],
           "flat_worst_regret": flat["value"],
           "objective_rows": free["n_objective_rows"],
           "nuisance_values": len(nuis_grid),
           "gain_pp": float(gain),
           "max_spread_across_extra_axis": spread,
           "beat_the_bound": bool(gain > 0.05),
           "lp_status_free": free["status"], "lp_status_flat": flat["status"],
           "seconds": float(time.time() - t0)}
    print("%-6s %-9s free %.5f  flat %.5f  gain %+6.3f pp  spread %.3f  "
          "(%.0f s) %s"
          % (tag, kind, free["value"], flat["value"], gain, spread,
             row["seconds"],
             "BEAT THE BOUND" if row["beat_the_bound"] else ""))
    return row


def sanity(tag, n, kind, nuis_grid, ncons=21):
    """A probe that MUST find a gain, so that finding none elsewhere means
    something.

    The same program with the size constrained at ONE nuisance value instead of
    a range.  The extra coordinate is then genuinely informative and the free
    optimum should beat the flat one.  If it does not, the mesh cannot express
    a non-invariant rule at all and the real probes prove nothing.
    """
    n = np.asarray(n, float)
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], ORACLE)) for l in ALT_LAMS]
    lam_cons = F.log_grid(*LAM_DOMAIN, ncons)
    P = Plane3(n, kind)
    one = [nuis_grid[len(nuis_grid) // 2]]
    free = solve(P, alts, lam_cons, one, flat=False)
    flat = solve(P, alts, lam_cons, one, flat=True)
    gain = 100.0 * (flat["value"] - free["value"])
    print("%-6s %-9s SANITY (size pinned at one nuisance value) gain "
          "%+6.3f pp" % (tag, kind, gain))
    return {"design": tag, "kind": kind, "gain_pp": float(gain),
            "expresses_noninvariance": bool(gain > 0.05)}


def main():
    t0 = time.time()
    designs = [("df11", [12.0, 13.0]), ("df35", [36.0, 13.0])]
    # the orbit has to be covered, not sampled: each extra bin needs a
    # nuisance value that puts mass in it, or the program finds the gap
    grids = {"scale": list(np.exp(np.linspace(-2.4, 2.4, 13))),
             "location": list(np.linspace(-2.8, 2.8, 15))}

    print("%-6s %-9s %13s %13s %12s %9s" %
          ("design", "probe", "free", "flat", "gain", "spread"))
    rows = []
    for tag, n in designs:
        for kind in ("scale", "location"):
            rows.append(probe(tag, n, kind, grids[kind]))

    print()
    sane = [sanity("df35", [36.0, 13.0], k, grids[k])
            for k in ("scale", "location")]

    beat = [r for r in rows if r["beat_the_bound"]]
    out = {
        "abar": ABAR, "oracle_power": ORACLE,
        "claim_under_test": ("the bound of Section 10.8, computed over rules "
                             "measurable in (T, lambda_hat), extends to every "
                             "test"),
        "probes": rows,
        "n_probes": len(rows),
        "n_that_beat_the_bound": len(beat),
        "worst_gain_pp": float(max(r["gain_pp"] for r in rows)),
        "superseded_first_version_gain_pp": SUPERSEDED_FIRST_VERSION_GAIN_PP,
        "superseded_first_version_error": (
            "scored weighted average power at one point of the orbit, which "
            "the Proposition does not cover; its sanity arm returned the same "
            "numbers as its real probe, which was the tell"),
        "worst_spread": float(max(r["max_spread_across_extra_axis"]
                                  for r in rows)),
        "claim_survives": bool(not beat),
        "sanity": sane,
        "sanity_shows_the_mesh_can_cheat": bool(
            all(s["expresses_noninvariance"] for s in sane)),
        "criterion": ("worst-case regret over the declared alternatives AND "
                      "over the orbit coordinate, which is the criterion the "
                      "Proposition covers; an earlier version used weighted "
                      "average power at one point of the orbit and found "
                      "gains of up to 13.9 points, which is consistent with "
                      "the Proposition rather than against it"),
        "note": ("Each probe hands the program a coordinate an invariant rule "
                 "discards, one per generator of the group, and takes both "
                 "its size constraint and its objective over the nuisance "
                 "that coordinate is sensitive to.  The flat reference is the "
                 "SAME program with equality across that axis imposed as a "
                 "constraint, so no rendering difference sits between the "
                 "two.  The sanity rows keep a single nuisance value, where "
                 "the orbit supremum degenerates and the coordinate IS "
                 "informative, so that finding no gain in the real probes is "
                 "a result rather than an inability to express one."),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_noninvariant.json", out, "frontier_measured")
    print("\n%d of %d probes beat the bound; worst gain %+.3f points"
          % (len(beat), len(rows), out["worst_gain_pp"]))
    print("the mesh can express a non-invariant rule: %s"
          % out["sanity_shows_the_mesh_can_cheat"])
    print("claim survives: %s ; %.0f s"
          % (out["claim_survives"], out["seconds"]))
    return 0 if out["claim_survives"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
