"""Are the latent loadings knowable at all, and from what?

This settles a boundary that was stated two different ways and can only be one.

The rule of Section 10 sees the data through `T` and `lambda_hat` and nothing
else, so the loadings are not a function of what the RULE sees.  That is true at
every candidate count and it is what the paper says.

Whether they are knowable from the DATA is a different question with a different
answer, and the difference is the three-indicator rule of factor analysis:

    Sigma = beta beta' + diag(sigma^2)

    K = 1, two series   3 distinct covariance entries against 4 parameters,
                        so the model is UNDER-identified and a declared box is
                        the only route;
    K >= 2, three or    beta_i^2 = Sigma_ij Sigma_ik / Sigma_jk for any two
    more series         other indices j, k, so the loadings are identified up
                        to a common sign, exactly, from the covariance alone.

The consequence for the programme: a confidence region for the loadings, of the
kind Sections 9.2 and 9.3 build for the ratio, is constructible at two or more
candidates and is not constructible at one.

One structural fact makes that region cheap to use.  For a normal sample the
sample mean and the sample covariance are independent, so a region built from
the covariance carries no information about the numerator of the statistic, and
the error allocation between the region's miss probability and the conditional
level composes the way Section 9.3's does.  That independence is measured here
rather than assumed.

Deposit: deposits/w2_identify.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F


def population_sigma(beta, lam):
    beta = np.asarray(beta, float)
    return np.outer(beta, beta) + np.diag(np.asarray(lam, float))


def three_indicator(S):
    """|beta_i| from the covariance, by the three-indicator identity."""
    p = S.shape[0]
    out = []
    for i in range(p):
        j, k = [q for q in range(p) if q != i][:2]
        out.append(np.sqrt(abs(S[i, j] * S[i, k] / S[j, k])))
    return np.array(out)


def count_identification(K):
    """Free parameters against distinct covariance entries."""
    p = K + 1
    params = 2 * p                      # p loadings and p residual variances
    entries = p * (p + 1) // 2
    return {"K": int(K), "series": int(p), "parameters": int(params),
            "covariance_entries": int(entries),
            "identified": bool(entries >= params and p >= 3)}


def sampling_check(beta, lam, n, D, seed):
    """How well the loadings come back from a SAMPLE covariance, and whether
    the sample mean and the sample covariance are independent in practice."""
    beta = np.asarray(beta, float)
    sig = np.sqrt(np.asarray(lam, float))
    p = beta.size
    rng = np.random.default_rng(int(seed))
    errs = np.empty(D)
    corr = np.empty(D)
    xb = np.empty((D, p))
    s01 = np.empty(D)
    for d in range(D):
        f = rng.normal(size=n)
        r = np.outer(f, beta) + rng.normal(size=(n, p)) * sig[None, :]
        S = np.cov(r, rowvar=False)
        est = three_indicator(S)
        errs[d] = np.max(np.abs(est - np.abs(beta)))
        xb[d] = r.mean(0)
        s01[d] = S[0, 1]
    # independence of Xbar and S: the correlation between the first mean and an
    # off-diagonal of the covariance, which is exactly zero in the model
    corr = float(np.corrcoef(xb[:, 0], s01)[0, 1])
    se = 1.0 / np.sqrt(D - 3.0)
    return {"n": int(n), "draws": int(D),
            "median_abs_error": float(np.median(errs)),
            "p90_abs_error": float(np.quantile(errs, 0.9)),
            "corr_xbar_with_cov": corr,
            "corr_z": float(corr / se),
            "independent_at_3se": bool(abs(corr / se) < 3.0)}


def main():
    t0 = time.time()
    print("IDENTIFICATION BY COUNTING")
    counts = [count_identification(k) for k in (1, 2, 3, 5)]
    for c in counts:
        print("  K=%d  %d series  %d parameters vs %d entries  ->  %s"
              % (c["K"], c["series"], c["parameters"],
                 c["covariance_entries"],
                 "identified" if c["identified"] else "under-identified"))

    print("\nEXACT RECOVERY FROM THE POPULATION COVARIANCE")
    exact = []
    for beta, lam, tag in (([1.0, -0.8, 1.5], [1.0, 2.0, 0.5], "K=2"),
                           ([0.9, 1.4, 0.3, 1.1], [1.0, 1.5, 0.7, 2.0], "K=3"),
                           ([0.5, 2.5, 1.8, 0.4, 1.2, 0.9],
                            [1.0, 1.0, 2.0, 0.5, 1.5, 0.8], "K=5")):
        S = population_sigma(beta, lam)
        got = three_indicator(S)
        err = float(np.max(np.abs(got - np.abs(beta))))
        exact.append({"case": tag, "beta": [float(b) for b in beta],
                      "max_abs_error": err})
        print("  %-4s max |error| %.2e" % (tag, err))

    print("\nFROM A SAMPLE, AND THE INDEPENDENCE THE ALLOCATION NEEDS")
    samp = []
    for n in (36, 120, 360):
        r = sampling_check([1.0, -0.8, 1.5], [1.0, 2.0, 0.5], n, 4000,
                           seed=1000 + n)
        samp.append(r)
        print("  n=%-4d median |error| %.3f   corr(Xbar, S) %+.4f (%.1f se)  %s"
              % (n, r["median_abs_error"], r["corr_xbar_with_cov"],
                 r["corr_z"], "independent" if r["independent_at_3se"] else
                 "NOT INDEPENDENT"))

    out = {
        "counts": counts,
        "identified_from_K": 2,
        "under_identified_at_K1": True,
        "exact_recovery": exact,
        "worst_population_error": float(max(e["max_abs_error"]
                                            for e in exact)),
        "sampling": samp,
        "worst_corr_z": float(max(abs(r["corr_z"]) for r in samp)),
        "mean_and_covariance_independent": bool(
            all(r["independent_at_3se"] for r in samp)),
        "what_the_rule_sees": ("the rule of Section 10 reads the data only "
                               "through T and lambda_hat, so the loadings are "
                               "not a function of what it sees at any "
                               "candidate count; that is a separate statement "
                               "from whether they are knowable from the data"),
        "consequence": ("a confidence region for the loadings is constructible "
                        "at two or more candidates and is not constructible at "
                        "one, where a declared box is the only route"),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w2_identify.json", out, "exact_in_model")
    print("\nidentified from K = %d upward; worst population error %.2e; "
          "mean and covariance independent: %s"
          % (out["identified_from_K"], out["worst_population_error"],
             out["mean_and_covariance_independent"]))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
