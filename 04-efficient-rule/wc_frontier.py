"""Workstreams B and C at K = 1: the exact efficiency frontier, and the rule
that attains it.

The plan separates "find a better valid rule" (B) from "show the loss is
unavoidable" (C).  At K = 1 both are the same computation, because the rule
class below is rich enough to contain every threshold function of the estimated
ratio up to bin resolution:

    c(lambda_hat) = c_j   for  lambda_hat in bin j,      j = 1..M .

Maximising power subject to the EXACT uniform constraint over that class gives
(i) the best deployable rule of this shape and (ii) an upper bound on what any
rule that sees only lambda_hat can achieve.  The gap between that optimum and
the known-lambda oracle is therefore the information cost, not a procedural
failure - which is exactly the terminal statement Workstream C asks for.

Everything is deterministic.  FWER and power are tensor quadratures over the two
scaled chi-squares, and because each state contributes through exactly one bin,
both are smooth with closed-form gradients, so the program is a small smooth NLP
rather than a search.

    FWER_m(c) = sum_s W_s Phi( -c_{b(m,s)} d_{m,s} / sigma_m )
    Power_a(c) = sum_s W_s Phi( (delta_a - c_{b(a,s)} d_{a,s}) / sigma_a )

Deposits: deposits/wc_frontier.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

import fcore as F
import k1lab as L

LAM_DOMAIN = (1e-3, 1e3)


class K1Problem:
    """Precomputed state/bin bookkeeping for one design."""

    def __init__(self, n, lam_constraint, bin_edges, npanel=20, order=10):
        self.n = np.asarray(n, float)
        self.edges = np.asarray(bin_edges, float)
        self.M = self.edges.size + 1
        q0, w0 = F.state_nodes(self.n[0] - 1.0, npanel, order)
        q1, w1 = F.state_nodes(self.n[1] - 1.0, npanel, order)
        U0, U1 = np.meshgrid(q0, q1, indexing="ij")
        self.u0, self.u1 = U0.ravel(), U1.ravel()
        W = np.outer(w0, w1).ravel()
        self.W = W / W.sum()
        self.lam_c = np.asarray(lam_constraint, float)
        self.cons = [self._slice(l) for l in self.lam_c]

    def _slice(self, lam):
        lh = lam * self.u1 / self.u0
        b = np.searchsorted(self.edges, lh)
        d = np.sqrt(lam * self.u1 / self.n[1] + self.u0 / self.n[0])
        sig = np.sqrt(lam / self.n[1] + 1.0 / self.n[0])
        return {"lam": float(lam), "bin": b, "d": d, "sig": float(sig),
                "ds": d / sig}

    # ---- FWER and its gradient ----------------------------------------
    def fwer(self, c, sl):
        x = -c[sl["bin"]] * sl["ds"]
        return float(self.W @ norm.cdf(x))

    def fwer_grad(self, c, sl):
        x = -c[sl["bin"]] * sl["ds"]
        g = -self.W * norm.pdf(x) * sl["ds"]
        out = np.zeros(self.M)
        np.add.at(out, sl["bin"], g)
        return out

    # ---- power and its gradient ---------------------------------------
    def power(self, c, sl, delta):
        x = (delta - c[sl["bin"]] * sl["d"]) / sl["sig"]
        return float(self.W @ norm.cdf(x))

    def power_grad(self, c, sl, delta):
        x = (delta - c[sl["bin"]] * sl["d"]) / sl["sig"]
        g = -self.W * norm.pdf(x) * sl["d"] / sl["sig"]
        out = np.zeros(self.M)
        np.add.at(out, sl["bin"], g)
        return out


def solve_frontier(n, alts, abar=F.ABAR, nbins=36, lam_cons=41,
                   npanel=20, order=10, c_init=None, weights=None,
                   maxiter=300):
    """Maximise weighted power subject to the exact uniform FWER constraint.

    `alts` is a list of (lam, delta) pairs declaring the alternatives the desk
    cares about.  `weights` defaults to uniform.
    """
    edges = F.log_grid(LAM_DOMAIN[0] * 3, LAM_DOMAIN[1] / 3, nbins - 1)
    lam_c = F.log_grid(*LAM_DOMAIN, lam_cons)
    P = K1Problem(n, lam_c, edges, npanel, order)
    alt_sl = [(P._slice(l), float(d)) for l, d in alts]
    w = np.ones(len(alts)) if weights is None else np.asarray(weights, float)
    w = w / w.sum()

    if c_init is None:
        c0 = np.full(P.M, 2.5)
    else:
        mid = np.concatenate([[edges[0] / 3.0],
                              np.sqrt(edges[:-1] * edges[1:]),
                              [edges[-1] * 3.0]])
        c0 = np.asarray(c_init(mid), float)

    def neg_obj(c):
        return -sum(wi * P.power(c, sl, d) for wi, (sl, d) in zip(w, alt_sl))

    def neg_jac(c):
        g = np.zeros(P.M)
        for wi, (sl, d) in zip(w, alt_sl):
            g -= wi * P.power_grad(c, sl, d)
        return g

    cons = [{"type": "ineq",
             "fun": (lambda c, s=s: abar - P.fwer(c, s)),
             "jac": (lambda c, s=s: -P.fwer_grad(c, s))}
            for s in P.cons]
    res = minimize(neg_obj, c0, jac=neg_jac, constraints=cons, method="SLSQP",
                   bounds=[(0.2, 30.0)] * P.M,
                   options={"maxiter": maxiter, "ftol": 1e-10})
    c = np.clip(res.x, 0.2, 30.0)
    sup = max(P.fwer(c, s) for s in P.cons)
    return {"c": c, "edges": edges, "problem": P, "sup_fwer": float(sup),
            "objective": float(-neg_obj(c)), "status": int(res.status),
            "message": str(res.message), "nit": int(res.nit)}


def step_rule(c, edges):
    c = np.asarray(c, float)
    edges = np.asarray(edges, float)

    def f(lh):
        return c[np.searchsorted(edges, np.asarray(lh, float))]

    return f


def main():
    designs = [("df3", [4.0, 13.0]), ("df5", [6.0, 13.0]),
               ("df7", [8.0, 13.0]), ("df11", [12.0, 13.0]),
               ("df19", [20.0, 13.0]), ("df35", [36.0, 13.0])]
    alt_lams = (0.04, 0.25, 1.0, 4.0, 25.0)
    out = {"abar": F.ABAR, "lam_domain": list(LAM_DOMAIN), "nbins": 36,
           "alt_lams": list(alt_lams), "designs": []}
    print("%-8s %8s %8s %8s %8s %8s" %
          ("design", "plugLoss", "optLoss", "gain_pp", "supFWER", "sec"))
    for name, n in designs:
        t0 = time.time()
        n = np.asarray(n, float)
        alts = [(l, F.delta_for_oracle_power_k1(n, [l], 0.8)) for l in alt_lams]

        # the calibrated plug-in of wb_k1, as the incumbent to beat
        coarse = F.log_grid(*LAM_DOMAIN, 25)
        cal = L.calibrate_level(n, coarse, F.ABAR, npanel=12, order=8, iters=10)
        tab = cal["table"]
        plug_pow = [F.power_exact_k1(tab, n, l, d, 24, 12) for l, d in alts]

        r = solve_frontier(n, alts, c_init=tab)
        rule = step_rule(r["c"], r["edges"])
        opt_pow = [F.power_exact_k1(rule, n, l, d, 24, 12) for l, d in alts]
        s_fine, lw, _ = L.sup_fwer(rule, n, F.log_grid(*LAM_DOMAIN, 97), 24, 12)

        # shrink if the fine grid finds a breach the optimiser's grid missed
        scale, tries = 1.0, 0
        while s_fine > F.ABAR and tries < 20:
            scale *= 1.004
            rule = step_rule(r["c"] * scale, r["edges"])
            s_fine, lw, _ = L.sup_fwer(rule, n, F.log_grid(*LAM_DOMAIN, 97),
                                       24, 12)
            tries += 1
        opt_pow = [F.power_exact_k1(rule, n, l, d, 24, 12) for l, d in alts]

        row = {
            "design": name, "n": n.tolist(),
            "alpha0_plugin": float(cal["alpha0"]),
            "plugin_power": [float(p) for p in plug_pow],
            "optimal_power": [float(p) for p in opt_pow],
            "plugin_loss_pp": [float(100 * (0.8 - p)) for p in plug_pow],
            "optimal_loss_pp": [float(100 * (0.8 - p)) for p in opt_pow],
            "mean_gain_pp": float(100 * (np.mean(opt_pow) - np.mean(plug_pow))),
            "sup_fwer_fine": float(s_fine), "lam_worst": float(lw),
            "shrink_scale": float(scale), "shrink_tries": int(tries),
            "valid": bool(s_fine <= F.ABAR),
            "c": [float(x) for x in (r["c"] * scale)],
            "edges": [float(x) for x in r["edges"]],
            "solver_status": r["status"], "solver_nit": r["nit"],
            "seconds": float(time.time() - t0),
        }
        out["designs"].append(row)
        print("%-8s %8.1f %8.1f %8.2f %8.5f %8.0f" %
              (name, float(np.mean(row["plugin_loss_pp"])),
               float(np.mean(row["optimal_loss_pp"])), row["mean_gain_pp"],
               s_fine, row["seconds"]))
    out["all_valid"] = all(d["valid"] for d in out["designs"])
    out["max_gain_pp"] = max(d["mean_gain_pp"] for d in out["designs"])
    out["note"] = (
        "The optimum is over piecewise-constant threshold functions of the "
        "estimated ratio with 36 bins spanning three decades.  Because every "
        "rule that sees only lambda_hat is of this form up to bin resolution, "
        "the optimal loss is an upper bound on what ANY such rule achieves, "
        "and the residual against the known-lambda oracle is the information "
        "cost of not knowing the ratio."
    )
    F.write_deposit("wc_frontier.json", out, "frontier_measured")
    print("\nall valid:", out["all_valid"], "| best mean gain %.2f pp"
          % out["max_gain_pp"])


if __name__ == "__main__":
    main()
