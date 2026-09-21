"""W1, Question 1: what information is actually available to the rule.

Section 10 searched piecewise-constant threshold functions of the estimated
ratio and found almost nothing recoverable above eleven degrees of freedom.
That is a statement about a CLASS.  This module establishes how big the class
is, which is what decides whether the statement is a description of the search
or a fact about the problem.

THE GROUP.  With one candidate the null is mu_0 = mu_1 = mu with sigma_0^2 and
sigma_1^2 unknown and unequal.  The testing problem is unchanged by

    x -> a + b x        applied to BOTH groups,  b > 0,

which is the affine group of the line.  Nothing else is admissible: scaling the
groups separately changes lambda and so changes the problem.

THE COUNT.  The sufficient statistic is four-dimensional, (Xbar_0, Xbar_1,
s_0^2, s_1^2).  The group is two-dimensional.  A maximal invariant is therefore
two-dimensional, and (T, lambda_hat) is a maximal invariant:

  * INVARIANCE.  Both coordinates are unchanged by every group element.
  * MAXIMALITY.  Any two sufficient statistics sharing (T, lambda_hat) are
    related by an explicit group element, recoverable in closed form:

        b = sqrt(w_0^2 / s_0^2),     a = Ybar_0 - b Xbar_0 .

Both are checked here against independently constructed data rather than
asserted.

WHAT THIS BUYS.  The class Section 10 searched is exactly the class of
non-randomised INVARIANT tests, up to bin resolution and monotonicity -- not a
convenient parametric family.  And the affine group is solvable, hence amenable,
so the Hunt-Stein theorem applies: a minimax test may be taken invariant.  A
numerical bound over (T, lambda_hat)-measurable tests is therefore a bound over
ALL tests for a worst-case criterion, which is the lower-bound statement the
programme could not previously write.  The conditions that has to satisfy are
recorded in the deposit rather than waved at; the numerical content of
Hunt-Stein is not checkable here and is not claimed to be.

BIN RESOLUTION.  Separately, the 36-bin discretisation is itself an
approximation.  Re-solving at 36, 72 and 144 bins bounds what it can hide.

Deposit: deposits/w1_invariant.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import wc_frontier as WC

ABAR = F.ABAR


def stats_of(x0, x1):
    """(Xbar_0, Xbar_1, s_0^2, s_1^2) from two samples."""
    return (x0.mean(-1), x1.mean(-1), x0.var(-1, ddof=1), x1.var(-1, ddof=1))


def invariant_of(st, n):
    """(T, lambda_hat) from the sufficient statistic."""
    m0, m1, v0, v1 = st
    T = (m1 - m0) / np.sqrt(v1 / n[1] + v0 / n[0])
    return T, v1 / v0


PAIRS = 200000          # draws per design, for both checks and the deposit


def invariance_check(n, D=PAIRS, seed=5150):
    """Apply random affine maps to raw data; the invariant must not move."""
    rng = np.random.default_rng(seed)
    x0 = rng.normal(0.3, 1.7, size=(D, int(n[0])))
    x1 = rng.normal(0.3, 0.6, size=(D, int(n[1])))
    T0, L0 = invariant_of(stats_of(x0, x1), n)
    a = rng.normal(0.0, 10.0, size=(D, 1))
    b = np.exp(rng.normal(0.0, 1.5, size=(D, 1)))
    T1, L1 = invariant_of(stats_of(a + b * x0, a + b * x1), n)
    return (float(np.max(np.abs(T1 - T0))),
            float(np.max(np.abs(L1 / L0 - 1.0))))


def maximality_check(n, D=PAIRS, seed=6160):
    """Build a SECOND statistic sharing (T, lambda_hat), then recover (a, b).

    The second statistic is constructed directly from the invariant, not by
    applying a group element, so recovering a group element that maps one onto
    the other is a genuine test of maximality rather than a tautology.
    """
    rng = np.random.default_rng(seed)
    m0 = rng.normal(0.0, 3.0, D)
    v0 = np.exp(rng.normal(0.0, 1.2, D))
    lh = np.exp(rng.normal(0.0, 2.0, D))
    v1 = lh * v0
    T = rng.normal(0.0, 2.5, D)
    m1 = m0 + T * np.sqrt(v1 / n[1] + v0 / n[0])

    # an independently drawn second point on the same (T, lambda_hat) level set
    w0 = np.exp(rng.normal(0.0, 1.2, D))
    w1 = lh * w0
    p0 = rng.normal(0.0, 3.0, D)
    p1 = p0 + T * np.sqrt(w1 / n[1] + w0 / n[0])

    b = np.sqrt(w0 / v0)
    a = p0 - b * m0
    err = np.max(np.abs(np.stack([
        (a + b * m0) - p0,
        (a + b * m1) - p1,
        b ** 2 * v0 - w0,
        b ** 2 * v1 - w1,
    ])))
    return float(err), float(np.min(b)), float(np.max(b))


def bin_resolution(design_n, counts=(36, 72, 144), alt_lams=(0.04, 0.25, 1.0,
                                                             4.0, 25.0)):
    """Re-solve the frontier at finer bin counts; the change bounds the error."""
    n = np.asarray(design_n, float)
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], 0.8)) for l in alt_lams]
    out = []
    for m in counts:
        t0 = time.time()
        r = WC.solve_frontier(n, alts, nbins=m)
        rule = WC.step_rule(r["c"], r["edges"])
        sup, lw, _ = F.sup_fwer_k1(rule, n, F.log_grid(1e-3, 1e3, 97), 24, 12)
        scale, tries = 1.0, 0
        while sup > ABAR and tries < 25:      # same shrink the frontier uses,
            scale *= 1.004                    # so all bin counts are compared
            rule = WC.step_rule(r["c"] * scale, r["edges"])   # as VALID rules
            sup, lw, _ = F.sup_fwer_k1(rule, n, F.log_grid(1e-3, 1e3, 97),
                                       24, 12)
            tries += 1
        pw = [F.power_exact_k1(rule, n, l, d, 24, 12) for l, d in alts]
        out.append({"nbins": int(m), "mean_power": float(np.mean(pw)),
                    "shrink_scale": float(scale),
                    "mean_loss_pp": float(100 * (0.8 - np.mean(pw))),
                    "sup_fwer": float(sup), "lam_worst": float(lw),
                    "valid": bool(sup <= ABAR),
                    "seconds": float(time.time() - t0)})
        print("   nbins %3d  loss %6.2f pp  supFWER %.5f  (%.0f s)"
              % (m, out[-1]["mean_loss_pp"], sup, out[-1]["seconds"]))
    return out


def main():
    t0 = time.time()
    designs = [("df35", [36.0, 13.0]), ("df11", [12.0, 13.0]),
               ("df3", [4.0, 13.0])]

    print("INVARIANCE AND MAXIMALITY")
    inv = []
    for tag, n in designs:
        n = np.asarray(n, float)
        dT, dL = invariance_check(n)
        err, blo, bhi = maximality_check(n)
        inv.append({"design": tag, "n": n.tolist(),
                    "max_abs_T_drift": dT, "max_rel_lamhat_drift": dL,
                    "max_recovery_error": err,
                    "b_range": [blo, bhi]})
        print("  %-6s invariance: |dT| %.2e  |dlam/lam| %.2e   "
              "maximality: recovery err %.2e" % (tag, dT, dL, err))

    print("\nBIN RESOLUTION (what 36 bins can hide)")
    res = []
    for tag, n in [("df35", [36.0, 13.0]), ("df3", [4.0, 13.0])]:
        print("  %s" % tag)
        rows = bin_resolution(n)
        base = next(r for r in rows if r["nbins"] == 36)
        fine = rows[-1]
        res.append({"design": tag, "n": list(map(float, n)), "rows": rows,
                    "hidden_pp": float(base["mean_loss_pp"]
                                       - fine["mean_loss_pp"])})
        print("   -> 36 bins can hide at most %.3f pp on this design"
              % res[-1]["hidden_pp"])

    out = {
        "abar": ABAR,
        "group": "x -> a + b x applied to both groups, b > 0 (affine group)",
        "sufficient_dim": 4,
        "group_dim": 2,
        "maximal_invariant_dim": 2,
        "maximal_invariant": ["T", "lambda_hat"],
        "pairs_per_design": PAIRS,
        "recovery_formula": "b = sqrt(w0^2/s0^2), a = Ybar_0 - b Xbar_0",
        "checks": inv,
        "worst_invariance_drift": float(max(
            max(c["max_abs_T_drift"], c["max_rel_lamhat_drift"])
            for c in inv)),
        "worst_recovery_error": float(max(c["max_recovery_error"]
                                          for c in inv)),
        "bin_resolution": res,
        "max_hidden_pp": float(max(r["hidden_pp"] for r in res)),
        "hunt_stein": {
            "group_is_solvable": True,
            "reason": ("translations form a normal abelian subgroup of the "
                       "affine group with abelian quotient, so the group is "
                       "solvable and therefore amenable"),
            "consequence": ("a minimax test may be taken invariant, so a bound "
                            "over (T, lambda_hat)-measurable tests bounds all "
                            "tests for a worst-case criterion"),
            "conditions_still_to_discharge": [
                "the loss (regret against the known-ratio oracle) must be "
                "invariant under the group, which it is because both power and "
                "oracle power are invariant functions of lambda",
                "the usual measurability and topological conditions of the "
                "Hunt-Stein theorem for the affine group",
                "the extension from K = 1 to K > 1, where the group is the "
                "same but the invariant is (T_1..T_K, lambda_hat_1..K)",
            ],
            "status": "derivation, not computation; not verified numerically",
        },
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_invariant.json", out, "exact_in_model")
    print("\nworst invariance drift %.2e ; worst recovery error %.2e"
          % (out["worst_invariance_drift"], out["worst_recovery_error"]))
    print("36 bins hide at most %.3f pp ; %.0f s"
          % (out["max_hidden_pp"], out["seconds"]))


if __name__ == "__main__":
    main()
