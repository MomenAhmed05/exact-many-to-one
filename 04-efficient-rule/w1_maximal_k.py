"""Step 4 of the Hunt-Stein item: the maximal invariant at any candidate count.

`w1_invariant.py` establishes the one-candidate case.  The proposition needs the
general one, and the claim is that nothing changes shape.

THE COUNT.  With K candidates the sufficient statistic is
(Xbar_0..Xbar_K, s_0^2..s_K^2), which is 2K + 2 coordinates.  On the boundary of
the null every mean is equal, so the parameter is one mean and K + 1 variances,
which is K + 2.  The group is the affine group applied to every sample at once,
which is two-dimensional.  A maximal invariant is therefore (2K + 2) - 2 = 2K
dimensional, and

    (T_1, .., T_K, lambda_hat_1, .., lambda_hat_K)

has exactly that many coordinates.

THE RECOVERY, which is what makes it MAXIMAL rather than merely invariant.
Given two sufficient statistics sharing the whole invariant, the group element
that maps one onto the other is determined by the BENCHMARK alone,

    b = sqrt(w_0^2 / s_0^2),      a = Ybar_0 - b Xbar_0 ,

and every candidate coordinate then follows without being used to fit anything:

    w_i^2 = (s_i^2 / s_0^2) w_0^2 = b^2 s_i^2                    from lambda_hat
    Ybar_i - Ybar_0 = b (Xbar_i - Xbar_0)                        from T_i
        because the denominator of T_i scales by exactly b,
    so  Ybar_i = a + b Xbar_i .

That is K + 1 pairs of coordinates recovered from a two-parameter element fitted
on one pair, so the check has 2K free coordinates of slack and is not a
tautology.  It is run here at K = 1, 2, 3 and 5.

Deposit: deposits/w1_maximal_k.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F

PAIRS = 200000


def invariant_of(m, v, n):
    """(T_1..T_K, lambda_hat_1..lambda_hat_K) from the sufficient statistic.

    `m` and `v` are (D, K+1); `n` is (K+1,).
    """
    d = np.sqrt(v[:, 1:] / n[None, 1:] + v[:, 0:1] / n[0])
    T = (m[:, 1:] - m[:, 0:1]) / d
    return T, v[:, 1:] / v[:, 0:1]


def invariance_check(n, D=PAIRS, seed=5150):
    """Random affine maps applied to raw data; the invariant must not move."""
    n = np.asarray(n, int)
    rng = np.random.default_rng(int(seed))
    K = n.size - 1
    xs = [rng.normal(0.4, 0.7 + 0.5 * i, size=(D, n[i])) for i in range(K + 1)]
    m = np.stack([x.mean(1) for x in xs], 1)
    v = np.stack([x.var(1, ddof=1) for x in xs], 1)
    T0, L0 = invariant_of(m, v, n.astype(float))

    a = rng.normal(0.0, 10.0, size=(D, 1))
    b = np.exp(rng.normal(0.0, 1.5, size=(D, 1)))
    ys = [a + b * x for x in xs]
    m2 = np.stack([y.mean(1) for y in ys], 1)
    v2 = np.stack([y.var(1, ddof=1) for y in ys], 1)
    T1, L1 = invariant_of(m2, v2, n.astype(float))
    return (float(np.max(np.abs(T1 - T0))),
            float(np.max(np.abs(L1 / L0 - 1.0))))


def maximality_check(n, D=PAIRS, seed=6160):
    """Two statistics built to share the invariant; recover the group element.

    Neither point is produced from the other, so the element has to be found
    rather than remembered, and it is fitted on the benchmark pair only.
    """
    n = np.asarray(n, float)
    K = int(n.size - 1)
    rng = np.random.default_rng(int(seed))

    # a point, drawn freely
    m = np.empty((D, K + 1))
    v = np.empty((D, K + 1))
    m[:, 0] = rng.normal(0.0, 3.0, D)
    v[:, 0] = np.exp(rng.normal(0.0, 1.2, D))
    lh = np.exp(rng.normal(0.0, 2.0, size=(D, K)))
    T = rng.normal(0.0, 2.5, size=(D, K))
    v[:, 1:] = lh * v[:, 0:1]
    d = np.sqrt(v[:, 1:] / n[None, 1:] + v[:, 0:1] / n[0])
    m[:, 1:] = m[:, 0:1] + T * d

    # a second point on the same level set, built from the invariant directly
    p = np.empty((D, K + 1))
    w = np.empty((D, K + 1))
    p[:, 0] = rng.normal(0.0, 3.0, D)
    w[:, 0] = np.exp(rng.normal(0.0, 1.2, D))
    w[:, 1:] = lh * w[:, 0:1]
    d2 = np.sqrt(w[:, 1:] / n[None, 1:] + w[:, 0:1] / n[0])
    p[:, 1:] = p[:, 0:1] + T * d2

    # the element, from the benchmark pair alone
    b = np.sqrt(w[:, 0] / v[:, 0])
    a = p[:, 0] - b * m[:, 0]

    err_m = np.max(np.abs((a[:, None] + b[:, None] * m) - p))
    err_v = np.max(np.abs((b[:, None] ** 2) * v - w))
    return float(max(err_m, err_v)), float(b.min()), float(b.max())


def dimension_count(K):
    p = K + 1
    return {"K": int(K), "groups": int(p),
            "sufficient_dim": int(2 * p),
            "null_parameter_dim": int(p + 1),
            "group_dim": 2,
            "maximal_invariant_dim": int(2 * p - 2),
            "invariant_coordinates": int(2 * K),
            "counts_agree": bool(2 * p - 2 == 2 * K)}


def main():
    t0 = time.time()
    designs = [("K1", [36.0, 13.0]),
               ("K2", [36.0, 13.0, 13.0]),
               ("K3", [36.0, 13.0, 20.0, 9.0]),
               ("K5", [36.0, 13.0, 20.0, 9.0, 30.0, 11.0])]
    rows = []
    print("%-5s %4s %14s %14s %16s" %
          ("case", "K", "|dT|", "|dlam/lam|", "recovery error"))
    for tag, n in designs:
        n = np.asarray(n, float)
        K = int(n.size - 1)
        dT, dL = invariance_check(n.astype(int))
        err, blo, bhi = maximality_check(n)
        c = dimension_count(K)
        rows.append(dict(c, design=tag, n=n.tolist(),
                         max_abs_T_drift=dT, max_rel_lamhat_drift=dL,
                         max_recovery_error=err,
                         b_range=[blo, bhi]))
        print("%-5s %4d %14.2e %14.2e %16.2e" % (tag, K, dT, dL, err))

    out = {
        "abar": F.ABAR,
        "group": "x -> a + b x, b > 0, applied to every sample at once",
        "pairs_per_design": PAIRS,
        "cases": rows,
        "worst_invariance_drift": float(max(
            max(r["max_abs_T_drift"], r["max_rel_lamhat_drift"])
            for r in rows)),
        "worst_recovery_error": float(max(r["max_recovery_error"]
                                          for r in rows)),
        "all_counts_agree": bool(all(r["counts_agree"] for r in rows)),
        "max_K": int(max(r["K"] for r in rows)),
        "recovery_formula": ("b = sqrt(w_0^2/s_0^2), a = Ybar_0 - b Xbar_0, "
                             "fitted on the benchmark pair alone"),
        "note": ("The element is fitted on one pair of coordinates and then "
                 "checked against all 2K + 2, so at K = 5 there are ten "
                 "coordinates of slack.  A formula that merely re-expressed "
                 "the invariant would not survive that."),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_maximal_k.json", out, "exact_in_model")
    print("\nworst drift %.2e, worst recovery error %.2e, counts agree at "
          "every K up to %d"
          % (out["worst_invariance_drift"], out["worst_recovery_error"],
             out["max_K"]))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
