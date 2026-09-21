"""Step 1 of the Hunt-Stein item: the criterion the theorem actually covers.

WHY THIS MODULE EXISTS.  `w1_region.py` maximises WEIGHTED AVERAGE power subject
to the uniform size constraint.  Hunt-Stein is a statement about MINIMAX risk,
and the two are not the same criterion.  The gap matters, and it is the risk the
plan flagged.

Where the averaging argument breaks.  For an amenable group the standard device
is to average a test over the group with an invariant mean, giving an invariant
test whose power at a parameter point is the MEAN of the original test's power
over that point's ORBIT.  For a minimax criterion that is enough, because a mean
over the orbit is at least the infimum over the orbit.  For a weighted average
of powers at chosen representatives it is not: the orbit mean can sit below the
value at the representative, so nothing follows.

So the criterion is restated as worst-case REGRET against the known-ratio
oracle, which is a genuine risk with an invariant loss:

    R(phi) = sup over declared alternatives of
             [ oracle power at that alternative  -  power of phi there ] .

For an invariant test the power is constant on orbits, so R is well defined on
the invariant; for a general test the supremum runs over every point of every
orbit and therefore picks up the worst representative.  That is exactly the
quantity Hunt-Stein equates between all tests and invariant tests.

WHAT IS COMPUTED.  The same linear program as `w1_region.py`, in epigraph form:

    minimise   t
    subject to t >= beta_oracle(a) - power_a(q)   for each declared alternative
               sum q * mass_lambda <= abar        for each declared ratio
               0 <= q <= 1

which is still linear, so the optimum is still global.  Its value is the best
worst-case regret any rule measurable in (T, lambda_hat) can achieve, and by the
proposition in `HUNT-STEIN.md` it is the best any test can achieve.

Deposit: deposits/w1_minimax.json
"""

from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

import fcore as F
from w1_region import ALT_LAMS, LAM_DOMAIN, Plane

ABAR = F.ABAR
ORACLE = 0.80          # every alternative is sized so the oracle reaches this


def solve_minimax(P, alts, lam_cons, monotone, abar=ABAR, oracle=ORACLE):
    """Minimise the worst regret over the declared alternatives.

    Variables are (q, t) with t the epigraph variable, so the objective is the
    single coordinate t and every alternative contributes one row.
    """
    nq = P.M * P.NT
    masses = [P.mass(lam, delta) for lam, delta in alts]
    nulls = [P.mass(l, 0.0) for l in lam_cons]

    rows, rhs = [], []
    # t >= oracle - power_a  <=>  -power_a - t <= -oracle
    for m in masses:
        rows.append(np.concatenate([-m, [-1.0]]))
        rhs.append(-oracle)
    # size
    for m in nulls:
        rows.append(np.concatenate([m, [0.0]]))
        rhs.append(abar)
    A = sp.csr_matrix(np.stack(rows))
    b = np.array(rhs, float)

    if monotone:
        Mr = P.monotone_rows()
        Z = sp.csr_matrix((Mr.shape[0], 1))
        A = sp.vstack([A, sp.hstack([Mr, Z])], format="csr")
        b = np.concatenate([b, np.zeros(Mr.shape[0])])

    c = np.zeros(nq + 1)
    c[-1] = 1.0
    res = linprog(c, A_ub=A, b_ub=b,
                  bounds=[(0.0, 1.0)] * nq + [(None, None)], method="highs")
    q = np.clip(res.x[:nq], 0.0, 1.0)
    return {"q": q, "t": float(res.x[-1]), "status": int(res.status)}


def regret_of(P, q, alts, oracle=ORACLE):
    return [float(oracle - P.power_of(q, lam, d)) for lam, d in alts]


def run_design(tag, n, nt=220, nbins=36, ncons=97, nverify=193):
    t0 = time.time()
    n = np.asarray(n, float)
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], ORACLE)) for l in ALT_LAMS]
    lam_cons = F.log_grid(*LAM_DOMAIN, ncons)
    lam_ver = F.log_grid(*LAM_DOMAIN, nverify)

    P = Plane(n, nbins=nbins, nt=nt)
    rm = solve_minimax(P, alts, lam_cons, monotone=True)
    ru = solve_minimax(P, alts, lam_cons, monotone=False)

    # The deployed rule is NOT rendered onto the grid.  Rendering drops the
    # cell straddling the threshold and hands it a free loss, which is the
    # defect that produced a spurious 1.47-point result once already.  It is
    # evaluated by the exact quadrature instead, and the two bases differ only
    # by the T-grid resolution, which w1_region measures at 0.034 points.
    import k1lab as L
    cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                            npanel=12, order=8, iters=10)
    tab = cal["table"]
    dep_profile = [ORACLE - F.power_exact_k1(tab, n, lam, d, 24, 12)
                   for lam, d in alts]
    sup_d, _lw, _v = F.sup_fwer_k1(tab, n, lam_ver, 24, 12)

    def worst(q):
        return max(regret_of(P, q, alts))

    reg_m, reg_u, reg_d = worst(rm["q"]), worst(ru["q"]), max(dep_profile)
    sup_m = max(P.fwer_of(rm["q"], l) for l in lam_ver)
    sup_u = max(P.fwer_of(ru["q"], l) for l in lam_ver)

    row = {
        "design": tag, "n": n.tolist(), "oracle_power": ORACLE,
        "alpha0_deployed": float(cal["alpha0"]),
        "minimax_regret_monotone_pp": float(100 * reg_m),
        "minimax_regret_unrestricted_pp": float(100 * reg_u),
        "deployed_worst_regret_pp": float(100 * reg_d),
        "gap_deployed_above_floor_pp": float(100 * (reg_d - reg_m)),
        "monotonicity_costs_pp": float(100 * (reg_m - reg_u)),
        "regret_profile_monotone_pp": [100 * r
                                       for r in regret_of(P, rm["q"], alts)],
        "regret_profile_deployed_pp": [100 * r for r in dep_profile],
        "sup_fwer_monotone": float(sup_m),
        "sup_fwer_unrestricted": float(sup_u),
        "sup_fwer_deployed": float(sup_d),
        "lp_status_monotone": rm["status"],
        "lp_status_unrestricted": ru["status"],
        "seconds": float(time.time() - t0),
    }
    print("%-6s  floor %6.2f pp   unrestricted %6.2f pp   deployed %6.2f pp   "
          "gap %+5.2f pp   (%.0f s)"
          % (tag, row["minimax_regret_monotone_pp"],
             row["minimax_regret_unrestricted_pp"],
             row["deployed_worst_regret_pp"],
             row["gap_deployed_above_floor_pp"], row["seconds"]))
    return row


def main():
    t0 = time.time()
    designs = [("df3", [4.0, 13.0]), ("df5", [6.0, 13.0]),
               ("df7", [8.0, 13.0]), ("df11", [12.0, 13.0]),
               ("df19", [20.0, 13.0]), ("df35", [36.0, 13.0])]
    rows = [run_design(t, n) for t, n in designs]
    hi = [r for r in rows if int(r["design"][2:]) >= 11]
    out = {
        "abar": ABAR, "oracle_power": ORACLE,
        "lam_domain": list(LAM_DOMAIN), "alt_lams": list(ALT_LAMS),
        "criterion": ("worst-case regret against the known-ratio oracle over "
                      "the declared alternatives, which is a risk with an "
                      "invariant loss and is therefore the criterion the "
                      "Hunt-Stein reduction covers"),
        "designs": rows,
        "gap_above_floor_low_pp": float(min(r["gap_deployed_above_floor_pp"]
                                            for r in hi)),
        "gap_above_floor_high_pp": float(max(r["gap_deployed_above_floor_pp"]
                                             for r in hi)),
        "worst_monotonicity_cost_pp": float(max(
            abs(r["monotonicity_costs_pp"]) for r in rows)),
        "monotonicity_free": bool(all(abs(r["monotonicity_costs_pp"]) < 0.01
                                      for r in rows)),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_minimax.json", out, "frontier_measured")
    print("\nabove eleven degrees of freedom the deployed rule is %.2f to "
          "%.2f points above the minimax-regret floor"
          % (out["gap_above_floor_low_pp"], out["gap_above_floor_high_pp"]))
    print("monotonicity costs at most %.3f points; %.0f s"
          % (out["worst_monotonicity_cost_pp"], out["seconds"]))


if __name__ == "__main__":
    main()
