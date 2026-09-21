"""M6, part 3: factor exposure, by the same route as serial dependence.

Serial dependence turned out not to be a boundary but a DECLARED-NUISANCE
pre-processing step: whiten with a declared coefficient and the iid certificate
applies unchanged.  Factor exposure has exactly the same shape, and this module
establishes it and then measures what happens when the declaration is replaced
by an estimate.

THE CONSTRUCTION.  For `r_{i,t} = alpha_i + beta_i f_t + eps_{i,t}` with a
Gaussian factor and Gaussian residuals, the factor-residual series

    e_{i,t} = r_{i,t} - beta_i f_t

is iid Gaussian whenever `beta_i` is DECLARED rather than estimated.  The
estimand changes - a promotion now means a higher factor alpha, not a higher raw
mean - and with that changed estimand the problem is again the iid problem this
repository solves, so every rule and every certificate transfers.

WHERE IT STOPS.  Estimating the loadings from the same returns being tested is a
new nuisance parameter, and it is not free: the residual variance estimate loses
a degree of freedom, and the residuals are no longer independent of the
estimated loading.  That is measured here rather than assumed away, together
with the version that keeps the guarantee - loadings from a pre-period that is
not reused for the performance test.

Deposit: deposits/factor.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L

ABAR = F.ABAR


def simulate(n0, n1, beta0, beta1, lam, R, rng, sig_f=1.0, alpha_shift=0.0,
             pre=0):
    """Returns for benchmark and candidate under a one-factor model.

    `pre` extra periods are generated first and returned separately, so a
    loading estimated on them is independent of the test window.
    """
    T = max(n0, n1) + pre
    f = rng.normal(0.0, sig_f, (R, T))
    e0 = rng.normal(0.0, 1.0, (R, T))
    e1 = rng.normal(0.0, np.sqrt(lam), (R, T))
    r0 = beta0 * f + e0
    r1 = alpha_shift + beta1 * f + e1
    return f, r0, r1


def ols_beta(y, x):
    """Slope of y on x, row-wise, no intercept removal beyond centring."""
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    return (xc * yc).sum(axis=1) / (xc * xc).sum(axis=1)


def welch(b, c, ddof_extra=0):
    n0, n1 = b.shape[1], c.shape[1]
    s0 = b.var(axis=1, ddof=1 + ddof_extra)
    s1 = c.var(axis=1, ddof=1 + ddof_extra)
    d = np.sqrt(s1 / n1 + s0 / n0)
    t = (c.mean(axis=1) - b.mean(axis=1)) / d
    return t, s1 / s0


def fwer_of(t, lh, tab):
    return float((t > tab(lh)).mean())


def main():
    t0 = time.time()
    n0 = n1 = 36
    n = np.array([float(n0), float(n1)])
    R = 40000
    lam = 1.0
    out = {"abar": ABAR, "n": n.tolist(), "draws": R, "lam": lam, "rows": []}

    cal = L.calibrate_level(n, F.log_grid(1e-3, 1e3, 25), ABAR, npanel=12,
                            order=8, iters=11)
    tab = cal["table"]
    out["alpha0_iid"] = float(cal["alpha0"])
    # residualising with an estimated loading costs a degree of freedom, so the
    # rule for that arm must be the one built for n-1 observations
    n_res = np.array([float(n0 - 1), float(n1 - 1)])
    cal_res = L.calibrate_level(n_res, F.log_grid(1e-3, 1e3, 25), ABAR,
                               npanel=12, order=8, iters=11)
    tab_res = cal_res["table"]
    out["alpha0_iid_one_df_less"] = float(cal_res["alpha0"])
    print("alpha0 (n=36) %.5f   alpha0 (n=35, for the residual arm) %.5f"
          % (cal["alpha0"], cal_res["alpha0"]))

    print("\n%-34s %9s %9s %7s" % ("arm", "FWER", "se", "z"))
    for label, b0, b1 in (("no exposure (beta = 0, 0)", 0.0, 0.0),
                          ("same sign (1.0, 0.8)", 1.0, 0.8),
                          ("opposite signs (1.0, -0.8)", 1.0, -0.8),
                          ("benchmark unexposed (0, 1.0)", 0.0, 1.0)):
        rng = np.random.default_rng(31 + int(10 * (b0 + 3 * b1)))
        f, r0, r1 = simulate(n0, n1, b0, b1, lam, R, rng)
        f_t, r0_t, r1_t = f[:, :n0], r0[:, :n0], r1[:, :n1]

        # (a) raw returns, no adjustment: the estimand is the raw mean
        t, lh = welch(r0_t, r1_t)
        raw = fwer_of(t, lh, tab)

        # (b) declared loadings: residualise with the TRUE beta
        e0 = r0_t - b0 * f_t
        e1 = r1_t - b1 * f[:, :n1]
        t, lh = welch(e0, e1)
        declared = fwer_of(t, lh, tab)

        # (c) estimated loadings, in-sample: a plug-in, not a theorem
        bh0 = ols_beta(r0_t, f_t)
        bh1 = ols_beta(r1_t, f[:, :n1])
        e0h = r0_t - bh0[:, None] * f_t
        e1h = r1_t - bh1[:, None] * f[:, :n1]
        t, lh = welch(e0h, e1h, ddof_extra=1)
        est_correct_df = fwer_of(t, lh, tab_res)
        t, lh = welch(e0h, e1h, ddof_extra=0)
        est_wrong_df = fwer_of(t, lh, tab)

        for nm, v in (("raw returns", raw), ("declared loadings", declared),
                      ("estimated loadings, df corrected", est_correct_df),
                      ("estimated loadings, df NOT corrected",
                       est_wrong_df)):
            se = float(np.sqrt(v * (1 - v) / R))
            z = (v - ABAR) / se
            out["rows"].append({"exposure": label, "arm": nm, "fwer": v,
                                "mcse": se, "z": float(z),
                                "within_ceiling": bool(v <= ABAR + 2 * se)})
            print("%-34s %9.5f %9.5f %+7.2f   [%s]"
                  % (label[:20] + " | " + nm[:12], v, se, z, nm))

    print("\npre-period loadings: estimated on data NOT reused for the test")
    rng = np.random.default_rng(777)
    f, r0, r1 = simulate(n0, n1, 1.0, -0.8, lam, R, rng, pre=60)
    fp, f_t = f[:, :60], f[:, 60:60 + n0]
    bh0 = ols_beta(r0[:, :60], fp)
    bh1 = ols_beta(r1[:, :60], fp)
    e0 = r0[:, 60:60 + n0] - bh0[:, None] * f_t
    e1 = r1[:, 60:60 + n1] - bh1[:, None] * f[:, 60:60 + n1]
    t, lh = welch(e0, e1)
    v = fwer_of(t, lh, tab)
    se = float(np.sqrt(v * (1 - v) / R))
    out["pre_period_loadings"] = {"fwer": v, "mcse": se,
                                  "z": float((v - ABAR) / se),
                                  "pre_periods": 60,
                                  "within_ceiling": bool(v <= ABAR + 2 * se)}
    print("    FWER %.5f (se %.5f, z %+.2f)" % (v, se, (v - ABAR) / se))

    out["status"] = (
        "With DECLARED loadings the factor-residual series is iid Gaussian and "
        "the certified rule applies unchanged to the factor-alpha estimand; "
        "that is the same structure as the whitening result for serial "
        "dependence.  With loadings estimated IN SAMPLE the guarantee is a "
        "plug-in: the degree of freedom lost to the regression must be given "
        "back to the rule, and the arm that does not do so is measured here to "
        "show the size of that mistake.  Loadings from a pre-period that is "
        "not reused restore the exact statement.")
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("factor.json", out, "unconditional_candidate_measured")
    print("\nwrote factor.json")


if __name__ == "__main__":
    main()
