"""M7, part 1: an independent cross-check of every decisive quantity.

A blind replication by another party is the plan's M7 gate and cannot be run
by the author, who has seen the code.  What can be run, and is worth
more than nothing, is a re-derivation of each decisive number by a route that
shares no code path with the one that produced it.  Where the two routes agree,
a defect would have to be present in both to survive.

  X-A  the deposited threshold values against `HeteroMaxT`, the paper's own
       exact-law implementation.  My tables are built by bisection on MY
       quadrature; this checks them against an implementation written for a
       different purpose by a different route.

  X-B  the unconditional FWER by Monte Carlo over the variance estimates with
       the EXACT conditional rejection inside, against the tensor quadrature.
       The two share no integration code: one samples the state law, the other
       integrates it on a probability-scale Gauss-Legendre grid.

  X-C  the power at the declared alternatives, the same way.

  X-D  the certified upper bound against a direct Monte Carlo of the same
       quantity: the bound must sit ABOVE the sampled value at every ratio, and
       a bound that sat below one would be a refutation, not a discrepancy.

Deposit: deposits/crosscheck.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import norm

import fcore as F
import wb_certify as C
from h12_hetero_exact import HeteroMaxT

RNG_SEED = 424242


def mc_fwer_k1(rule, n, lam, R=200000, seed=0):
    """Unconditional FWER by sampling the state law, exact inner.

    Shares no integration code with `fcore.fwer_exact_k1`: the states are drawn
    rather than quadratured, and the conditional rejection is the closed form.
    """
    n = np.asarray(n, float)
    rng = np.random.default_rng(seed)
    u0 = rng.chisquare(n[0] - 1.0, R) / (n[0] - 1.0)
    u1 = rng.chisquare(n[1] - 1.0, R) / (n[1] - 1.0)
    lh = lam * u1 / u0
    c = rule(lh)
    d = np.sqrt(lam * u1 / n[1] + u0 / n[0])
    sig = np.sqrt(lam / n[1] + 1.0 / n[0])
    p = norm.cdf(-(c * d) / sig)
    return float(p.mean()), float(p.std(ddof=1) / np.sqrt(R))


def mc_power_k1(rule, n, lam, delta, R=200000, seed=0):
    n = np.asarray(n, float)
    rng = np.random.default_rng(seed)
    u0 = rng.chisquare(n[0] - 1.0, R) / (n[0] - 1.0)
    u1 = rng.chisquare(n[1] - 1.0, R) / (n[1] - 1.0)
    lh = lam * u1 / u0
    c = rule(lh)
    d = np.sqrt(lam * u1 / n[1] + u0 / n[0])
    sig = np.sqrt(lam / n[1] + 1.0 / n[0])
    p = norm.cdf((delta - c * d) / sig)
    return float(p.mean()), float(p.std(ddof=1) / np.sqrt(R))


def main():
    t0 = time.time()
    dep = F.read_deposit("wb_certified_rule.json")
    out = {"abar": dep["abar"], "seed": RNG_SEED, "xa": [], "xb": [], "xc": [],
           "xd": []}

    print("X-A  deposited thresholds vs HeteroMaxT (independent implementation)")
    worst_a = 0.0
    for d in dep["designs"][:6]:
        n = np.asarray(d["n"], float)
        edges = np.asarray(d["edges"], float)
        mids = np.concatenate([[edges[0] / 3.0],
                               np.sqrt(edges[:-1] * edges[1:]),
                               [edges[-1] * 3.0]])
        # check a spread of bins rather than all sixty
        idx = np.linspace(5, len(mids) - 6, 7).astype(int)
        diffs = []
        for j in idx:
            lam_hat = float(mids[j])
            ref = float(HeteroMaxT(n, np.sqrt([1.0, lam_hat]),
                                   grid=F.GRID).critical_value(d["alpha0"]))
            diffs.append(abs(ref - d["c"][j]))
        w = float(max(diffs))
        worst_a = max(worst_a, w)
        out["xa"].append({"design": d["design"], "bins_checked": len(idx),
                          "worst_abs_diff": w})
        print("    %-16s worst |mine - HeteroMaxT| = %.2e" % (d["design"], w))
    out["xa_worst"] = worst_a

    print("\nX-B  FWER: Monte Carlo (exact inner) vs tensor quadrature")
    worst_z = 0.0
    for d in dep["designs"][:6]:
        n = np.asarray(d["n"], float)
        rule = C.StepRule(d["c"], d["edges"])
        for lam in (0.04, 1.0, 25.0):
            q = F.fwer_exact_k1(rule, n, lam, 24, 12)
            m, se = mc_fwer_k1(rule, n, lam, seed=RNG_SEED + int(lam * 7))
            z = (m - q) / se if se > 0 else 0.0
            worst_z = max(worst_z, abs(z))
            out["xb"].append({"design": d["design"], "lam": float(lam),
                              "quadrature": float(q), "monte_carlo": float(m),
                              "mcse": float(se), "z": float(z)})
        print("    %-16s worst |z| so far %.2f" % (d["design"], worst_z))
    out["xb_worst_abs_z"] = float(worst_z)

    print("\nX-C  power: Monte Carlo vs tensor quadrature")
    worst_zc = 0.0
    for d in dep["designs"][:4]:
        n = np.asarray(d["n"], float)
        rule = C.StepRule(d["c"], d["edges"])
        for lam in (0.04, 1.0, 25.0):
            dl = F.delta_for_oracle_power_k1(n, [lam], 0.8, dep["abar"])
            q = F.power_exact_k1(rule, n, lam, dl, 24, 12)
            m, se = mc_power_k1(rule, n, lam, dl,
                                seed=RNG_SEED + 99 + int(lam * 7))
            z = (m - q) / se if se > 0 else 0.0
            worst_zc = max(worst_zc, abs(z))
            out["xc"].append({"design": d["design"], "lam": float(lam),
                              "quadrature": float(q), "monte_carlo": float(m),
                              "mcse": float(se), "z": float(z)})
    out["xc_worst_abs_z"] = float(worst_zc)
    print("    worst |z| %.2f" % worst_zc)

    print("\nX-D  the certified upper bound must dominate a direct sample")
    bad = []
    for d in dep["designs"][:6]:
        n = np.asarray(d["n"], float)
        rule = C.StepRule(d["c"], d["edges"])
        ub = d["sup_upper_bound"]
        for lam in (0.01, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0):
            m, se = mc_fwer_k1(rule, n, lam, R=100000,
                               seed=RNG_SEED + 7 + int(lam * 11))
            # a breach is only a refutation if it is beyond sampling error
            if m - 3.0 * se > ub:
                bad.append({"design": d["design"], "lam": float(lam),
                            "mc": float(m), "bound": float(ub)})
            out["xd"].append({"design": d["design"], "lam": float(lam),
                              "monte_carlo": float(m), "mcse": float(se),
                              "bound": float(ub),
                              "dominated": bool(m - 3.0 * se <= ub)})
    out["xd_violations"] = bad
    out["xd_all_dominated"] = not bad
    print("    %d ratio points checked, %d violations beyond 3 se"
          % (len(out["xd"]), len(bad)))

    out["verdict"] = {
        "xa_pass": bool(worst_a < 5e-3),
        "xb_pass": bool(worst_z < 4.0),
        "xc_pass": bool(worst_zc < 4.0),
        "xd_pass": bool(not bad),
    }
    out["all_pass"] = all(out["verdict"].values())
    out["scope"] = (
        "This is a cross-check, NOT a blind replication: it was written by the "
        "same author in the same session as the code it checks, and it reuses "
        "the deposited thresholds rather than re-deriving them from a "
        "specification.  What it establishes is that the quadrature, the "
        "tabulation and the enclosure are not all wrong in the same direction, "
        "which is the failure mode a single implementation cannot detect.  The "
        "blind version is the kit in replicate/, and it has not been executed.")
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("crosscheck.json", out, "frontier_measured")
    print("\nverdict:", out["verdict"], "| all pass:", out["all_pass"])


if __name__ == "__main__":
    main()
