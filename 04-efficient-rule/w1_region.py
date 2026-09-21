"""W1, Question 2: what a non-threshold, possibly randomised region recovers.

Section 10 optimised over MONOTONE threshold rules, reject iff T > h(lambda_hat).
Two restrictions sit inside that: the region must be an upper set in T, and the
decision must be deterministic.  Neither is implied by validity.  This module
drops both and measures what they cost.

THE PROGRAM.  Discretise the (T, lambda_hat) plane.  Let q[b, k] in [0, 1] be
the probability of rejecting when the estimated ratio falls in bin b and T falls
in cell k.  For a true ratio lambda the state law is the same tensor quadrature
the rest of the folder uses, and conditionally on a state T is normal, so the
mass each (b, k) receives is exact:

    mass_lambda[b, k] = sum_{s : bin(s) = b} W_s
                        [ Phi(t_{k+1} d_s / sigma) - Phi(t_k d_s / sigma) ]

and under an alternative shift delta the same with t d_s replaced by
t d_s - delta.  Objective and constraints are LINEAR in q, so the best valid
region is a linear program.

THE COMPARISON HAS TO BE FAIR, and the first version of this module was not.
Rendering Section 10's threshold rule onto the T grid and dropping the cell that
straddles the threshold costs that rule about a point of power for free, which
is the same order as the effect being measured.  So monotonicity is imposed as a
CONSTRAINT instead:

    q[b, k] <= q[b, k+1]        for every bin b and cell k,

and the monotone optimum and the unrestricted optimum are then two solutions of
the same program on the same grid under the same uniform constraint.  Their
difference is the cost of monotonicity and nothing else.

Determinism is priced separately, by rounding the unrestricted solution down: a
vertex has at most (number of constraints) fractional entries, so the rounded
rule is feasible and its value is a floor on the best deterministic region.

Deposit: deposits/w1_region.json
"""

from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from scipy.stats import norm

import fcore as F

ABAR = F.ABAR
LAM_DOMAIN = (1e-3, 1e3)
ALT_LAMS = (0.04, 0.25, 1.0, 4.0, 25.0)


class Plane(object):
    """State quadrature plus the (lambda_hat bin, T cell) bookkeeping."""

    def __init__(self, n, nbins=36, nt=220, tlo=-8.0, thi=14.0,
                 npanel=9, order=8):
        self.n = np.asarray(n, float)
        self.edges = F.log_grid(LAM_DOMAIN[0] * 3, LAM_DOMAIN[1] / 3, nbins - 1)
        self.M = self.edges.size + 1
        self.t = np.linspace(tlo, thi, int(nt) + 1)
        self.NT = int(nt) + 2          # plus the two open tail cells
        q0, w0 = F.state_nodes(self.n[0] - 1.0, npanel, order)
        q1, w1 = F.state_nodes(self.n[1] - 1.0, npanel, order)
        U0, U1 = np.meshgrid(q0, q1, indexing="ij")
        self.u0, self.u1 = U0.ravel(), U1.ravel()
        W = np.outer(w0, w1).ravel()
        self.W = W / W.sum()
        self._mono = None

    def _geom(self, lam):
        lh = lam * self.u1 / self.u0
        b = np.searchsorted(self.edges, lh)
        d = np.sqrt(lam * self.u1 / self.n[1] + self.u0 / self.n[0])
        sig = np.sqrt(lam / self.n[1] + 1.0 / self.n[0])
        return b, d, sig

    def mass(self, lam, delta=0.0):
        """Exact probability mass in each (bin, T cell), shape (M * NT,)."""
        b, d, sig = self._geom(lam)
        x = (self.t[None, :] * d[:, None] - delta) / sig
        cdf = norm.cdf(x)
        cells = np.empty((d.size, self.NT))
        cells[:, 0] = cdf[:, 0]
        cells[:, 1:-1] = np.diff(cdf, axis=1)
        cells[:, -1] = 1.0 - cdf[:, -1]
        out = np.zeros((self.M, self.NT))
        np.add.at(out, b, cells * self.W[:, None])
        return out.ravel()

    def power_of(self, q, lam, delta):
        return float(q @ self.mass(lam, delta))

    def fwer_of(self, q, lam):
        return float(q @ self.mass(lam, 0.0))

    def monotone_rows(self):
        """Sparse q[b,k] - q[b,k+1] <= 0 for every bin and adjacent pair."""
        if self._mono is None:
            M, NT = self.M, self.NT
            idx = np.arange(M * NT).reshape(M, NT)
            a = idx[:, :-1].ravel()
            b = idx[:, 1:].ravel()
            r = np.arange(a.size)
            self._mono = sp.csr_matrix(
                (np.concatenate([np.ones(a.size), -np.ones(b.size)]),
                 (np.concatenate([r, r]), np.concatenate([a, b]))),
                shape=(a.size, M * NT))
        return self._mono


def solve(P, alts, lam_cons, monotone, abar=ABAR):
    obj = np.zeros(P.M * P.NT)
    for lam, delta in alts:
        obj += P.mass(lam, delta)
    obj /= len(alts)
    A = sp.csr_matrix(np.stack([P.mass(l, 0.0) for l in lam_cons]))
    b = np.full(A.shape[0], abar)
    if monotone:
        Mr = P.monotone_rows()
        A = sp.vstack([A, Mr], format="csr")
        b = np.concatenate([b, np.zeros(Mr.shape[0])])
    res = linprog(-obj, A_ub=A, b_ub=b, bounds=(0.0, 1.0), method="highs")
    q = np.clip(res.x, 0.0, 1.0)
    frac = (q > 1e-9) & (q < 1 - 1e-9)
    qdet = np.where(q > 1 - 1e-9, 1.0, 0.0)
    return {"q": q, "qdet": qdet, "n_frac": int(frac.sum()),
            "status": int(res.status), "obj": obj}


def run_design(tag, n, nt=220, nbins=36, ncons=97, nverify=193, refine=False):
    t0 = time.time()
    n = np.asarray(n, float)
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], 0.8)) for l in ALT_LAMS]
    lam_cons = F.log_grid(*LAM_DOMAIN, ncons)
    lam_ver = F.log_grid(*LAM_DOMAIN, nverify)

    P = Plane(n, nbins=nbins, nt=nt)
    rm = solve(P, alts, lam_cons, monotone=True)
    ru = solve(P, alts, lam_cons, monotone=False)

    def summarise(q):
        return (float(np.mean([P.power_of(q, l, d) for l, d in alts])),
                max(P.fwer_of(q, l) for l in lam_ver))

    pw_mono, sup_mono = summarise(rm["q"])
    pw_rand, sup_rand = summarise(ru["q"])
    pw_det, sup_det = summarise(ru["qdet"])

    if refine:
        P2 = Plane(n, nbins=nbins, nt=2 * nt)
        r2 = solve(P2, alts, lam_cons, monotone=False)
        pw_rand2 = float(np.mean([P2.power_of(r2["q"], l, d) for l, d in alts]))
    else:
        pw_rand2 = pw_rand

    oracle = 0.80
    row = {
        "design": tag, "n": n.tolist(), "nbins": nbins, "nt": nt,
        "oracle_power": oracle,
        "power_monotone": pw_mono,
        "power_randomised_region": pw_rand,
        "power_deterministic_floor": pw_det,
        "loss_pp_monotone": float(100 * (oracle - pw_mono)),
        "loss_pp_randomised": float(100 * (oracle - pw_rand)),
        "loss_pp_deterministic_floor": float(100 * (oracle - pw_det)),
        "gain_pp_from_dropping_monotonicity": float(100 * (pw_rand - pw_mono)),
        "gain_pp_deterministic": float(100 * (pw_det - pw_mono)),
        "n_fractional_cells_monotone": rm["n_frac"],
        "n_fractional_cells_unrestricted": ru["n_frac"],
        "sup_fwer_monotone": float(sup_mono),
        "sup_fwer_randomised": float(sup_rand),
        "sup_fwer_deterministic": float(sup_det),
        "overshoot_monotone": float(max(0.0, sup_mono - ABAR)),
        "overshoot_randomised": float(max(0.0, sup_rand - ABAR)),
        "t_grid_sensitivity_pp": float(100 * abs(pw_rand2 - pw_rand)),
        "lp_status_monotone": rm["status"],
        "lp_status_unrestricted": ru["status"],
        "seconds": float(time.time() - t0),
    }
    print("%-6s  monotone %6.2f pp   unrestricted %6.2f pp   det-floor "
          "%6.2f pp   gain %+5.2f pp   sup %.5f/%.5f  (%.0f s)"
          % (tag, row["loss_pp_monotone"], row["loss_pp_randomised"],
             row["loss_pp_deterministic_floor"],
             row["gain_pp_from_dropping_monotonicity"], sup_mono, sup_rand,
             row["seconds"]))
    return row


def main():
    t0 = time.time()
    designs = [("df3", [4.0, 13.0]), ("df5", [6.0, 13.0]),
               ("df7", [8.0, 13.0]), ("df11", [12.0, 13.0]),
               ("df19", [20.0, 13.0]), ("df35", [36.0, 13.0])]
    rows = [run_design(t, n, refine=(t == "df35")) for t, n in designs]
    worst_over = max(max(r["overshoot_monotone"], r["overshoot_randomised"])
                     for r in rows)
    out = {
        "abar": ABAR, "lam_domain": list(LAM_DOMAIN),
        "alt_lams": list(ALT_LAMS), "designs": rows,
        "n_constraint_grid": 97, "n_verify_grid": 193,
        "max_gain_pp_from_dropping_monotonicity":
            float(max(r["gain_pp_from_dropping_monotonicity"] for r in rows)),
        "min_gain_pp_from_dropping_monotonicity":
            float(min(r["gain_pp_from_dropping_monotonicity"] for r in rows)),
        "worst_overshoot": float(worst_over),
        "worst_overshoot_pct_of_abar": float(100 * worst_over / ABAR),
        "worst_t_grid_sensitivity_pp":
            float(max(r["t_grid_sensitivity_pp"] for r in rows)),
        "note": ("Monotonicity is imposed as a linear constraint rather than "
                 "by rendering a threshold rule onto the grid, because "
                 "rendering drops the cell that straddles the threshold and "
                 "costs that rule about a point of power for free -- the same "
                 "order as the effect being measured.  Both numbers are now "
                 "optima of the same program on the same grid under the same "
                 "uniform constraint.  The constraint is imposed at %d ratio "
                 "values and verified at %d, so a small overshoot on the "
                 "finer grid is expected and is reported rather than "
                 "shrunk away." % (97, 193)),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_region.json", out, "measured")
    print("\ncost of monotonicity and determinism: %+.2f to %+.2f points"
          % (out["min_gain_pp_from_dropping_monotonicity"],
             out["max_gain_pp_from_dropping_monotonicity"]))
    print("worst overshoot on the verification grid %.2e (%.2f%% of abar); "
          "T-grid sensitivity <= %.3f pp ; %.0f s"
          % (out["worst_overshoot"], out["worst_overshoot_pct_of_abar"],
             out["worst_t_grid_sensitivity_pp"], out["seconds"]))


if __name__ == "__main__":
    main()
