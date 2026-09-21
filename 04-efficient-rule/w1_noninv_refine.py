"""Does the falsifier's residual gain shrink when the nuisance grid refines?

THE PROBLEM THE PROBE RAN INTO.  The Proposition's size constraint runs over a
CONTINUUM of the nuisance.  A linear program can only impose it at finitely many
values, and a rule that reads the extra coordinate can always place itself in
the gaps between them.  So the free arm beats the flat one by SOMETHING at every
finite grid, and the question is not whether the gain is zero but whether it
goes to zero as the grid fills in.  If it does, the probe supports the
Proposition and the residual is discretisation.  If it plateaus, it does not.

THE SIMPLIFICATION THIS MODULE TRIED, AND WHY IT DOES NOT WORK.

A rule constant across the extra axis is a function of `(T, lambda_hat)`, so in
exact arithmetic its optimum is the two-dimensional minimax program of
`w1_minimax.py`.  Solving the flat arm there avoids the 17,000 equality rows
that made the first attempt fail numerically.

In exact arithmetic.  Not on a grid.  The free arm is an optimum on a
20 x 122 x 8 mesh and the flat arm an optimum on a 20 x 122 one, and their
discretisation errors do not cancel: the difference between them came out near
-62 points and did not move as the nuisance grid went from 3 values to 15.
Flatness in the nuisance grid is the signature, because the quantity being
measured is supposed to depend on exactly that.  So this run measures the
mismatch between two meshes and not the effect of the extra axis, its deposit
carries `VALID: false`, and the instrument is unbuilt rather than
inconclusive.

A valid version needs both arms on ONE mesh.  That is what failed numerically
in `w1_noninvariant.py`, which is where the work would have to restart, most
likely by constructing specific non-invariant rules and evaluating them with
the exact quadrature rather than optimising over a mesh at all.

The location probe is dropped.  With `n_0 = 36` the spread of `Xbar_0` about the
common mean is `1/6`, so the extra axis barely moves and the probe is weak by
construction; the scale probe carries the same argument with a coordinate that
actually varies.

Deposit: deposits/w1_noninv_refine.json
"""

from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

import fcore as F
from w1_minimax import solve_minimax
from w1_noninvariant import Plane3
from w1_region import ALT_LAMS, LAM_DOMAIN, Plane

ABAR = F.ABAR
ORACLE = 0.80


def free_arm(P, alts, lam_cons, nuis, abar=ABAR, oracle=ORACLE):
    """Worst regret over alternatives and the orbit, rule free in the extra."""
    nq = P.nvar
    rows, rhs = [], []
    for lam, delta in alts:
        for v in nuis:
            rows.append(np.concatenate([-P.mass(lam, v, delta), [-1.0]]))
            rhs.append(-oracle)
    for l in lam_cons:
        for v in nuis:
            rows.append(np.concatenate([P.mass(l, v, 0.0), [0.0]]))
            rhs.append(abar)
    A = sp.csr_matrix(np.stack(rows))
    b = np.array(rhs, float)
    Mr = P.monotone_rows()
    A = sp.vstack([A, sp.hstack([Mr, sp.csr_matrix((Mr.shape[0], 1))])],
                  format="csr")
    b = np.concatenate([b, np.zeros(Mr.shape[0])])
    c = np.zeros(nq + 1)
    c[-1] = 1.0
    res = linprog(c, A_ub=A, b_ub=b,
                  bounds=[(0.0, 1.0)] * nq + [(None, None)], method="highs")
    if res.x is None:
        return None, int(res.status)
    return float(res.x[-1]), int(res.status)


def main():
    t0 = time.time()
    n = np.array([36.0, 13.0])
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], ORACLE)) for l in ALT_LAMS]
    lam_cons = F.log_grid(*LAM_DOMAIN, 15)

    # the flat arm: a rule constant across the extra axis is a function of
    # (T, lambda_hat), so its optimum is the two-dimensional program
    P2 = Plane(n, nbins=20, nt=120)
    rm = solve_minimax(P2, alts, lam_cons, monotone=True)
    flat = rm["t"]
    print("flat arm (two-dimensional minimax): %.6f" % flat)

    rows = []
    print("\n%8s %12s %12s %10s %8s" %
          ("nuisance", "free", "flat", "gain pp", "sec"))
    for m in (3, 5, 9, 15):
        nuis = list(np.exp(np.linspace(-2.4, 2.4, m)))
        P3 = Plane3(n, "scale", nbins=20, nt=120, nx=8)
        t1 = time.time()
        val, st = free_arm(P3, alts, lam_cons, nuis)
        sec = time.time() - t1
        if val is None:
            print("%8d   LP FAILED (status %d)" % (m, st))
            rows.append({"nuisance_values": m, "failed": True,
                         "lp_status": st})
            continue
        gain = 100.0 * (flat - val)
        rows.append({"nuisance_values": m, "failed": False,
                     "free_worst_regret": val, "flat_worst_regret": flat,
                     "gain_pp": float(gain), "lp_status": st,
                     "seconds": float(sec)})
        print("%8d %12.6f %12.6f %10.3f %8.0f" % (m, val, flat, gain, sec))

    good = [r for r in rows if not r["failed"]]
    shrinking = bool(len(good) >= 2
                     and all(good[i]["gain_pp"] <= good[i - 1]["gain_pp"] + 1e-9
                             for i in range(1, len(good))))
    # The two arms are not comparable and the run is recorded as invalid.
    # `free` is the optimum on a 20 x 122 x 8 mesh; `flat` is the optimum on a
    # 20 x 122 mesh solved by a different module.  Their discretisation errors
    # do not cancel, and the difference between them, about -62 points and flat
    # in the nuisance grid, is that mismatch rather than any effect of the
    # extra axis.  A valid version needs both arms on ONE mesh, which is what
    # failed numerically in w1_noninvariant.py, so the instrument is unbuilt
    # rather than merely inconclusive.
    valid = bool(good and abs(good[0]["gain_pp"]) < 5.0)
    out = {
        "VALID": valid,
        "invalid_because": (None if valid else
                            "the free and flat arms are optima on different "
                            "meshes, so their difference is a discretisation "
                            "mismatch and not a measurement of the extra "
                            "axis; the gain is flat in the nuisance grid, "
                            "which is the signature of that"),
        "abar": ABAR, "oracle_power": ORACLE, "n": n.tolist(),
        "probe": "scale",
        "flat_worst_regret": float(flat),
        "rows": rows,
        "n_solved": len(good),
        "gain_at_coarsest_pp": float(good[0]["gain_pp"]) if good else None,
        "gain_at_finest_pp": float(good[-1]["gain_pp"]) if good else None,
        "monotonically_shrinking_IGNORE_IF_INVALID": shrinking,
        "note": ("The size constraint of the Proposition runs over a "
                 "continuum; a program can impose it only at finitely many "
                 "nuisance values, and a rule reading the extra coordinate "
                 "can sit in the gaps.  So a positive gain at any finite grid "
                 "is expected, and the question is whether it goes to zero as "
                 "the grid fills in.  The flat arm is solved on the "
                 "two-dimensional mesh, since a rule constant across the "
                 "extra axis is a function of (T, lambda_hat); that is exact "
                 "and it removes the equality rows whose size made the first "
                 "attempt fail."),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_noninv_refine.json", out, "frontier_measured")
    if good:
        print("\ngain %.3f -> %.3f points as the grid goes %d -> %d; "
              "monotonically shrinking: %s"
              % (good[0]["gain_pp"], good[-1]["gain_pp"],
                 good[0]["nuisance_values"], good[-1]["nuisance_values"],
                 shrinking))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
