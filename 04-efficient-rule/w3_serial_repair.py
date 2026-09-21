"""W3: why tier B and tier C fail for serial dependence, and the repair.

`w3_tiers.py` measures the three validity tiers.  Tier A holds exactly, as
Section 11 says it does.  Tier B and tier C do not, and the size of the failure
is not a rounding error: at a lag-one coefficient of 0.6 and a common mean of 2,
estimating the coefficient on the evaluation period takes the realised
familywise error rate to about eleven times its ceiling.  Even a 360-observation
calibration period leaves it above.

THE MECHANISM.  The Cochrane-Orcutt transform maps a series with mean mu to one
with mean mu(1 - phi).  With phi DECLARED, both groups are multiplied by the
same constant, so under the null the difference of the transformed means is
still zero whatever mu is, and that is why tier A is exact.

With phi ESTIMATED SEPARATELY IN EACH GROUP the two multipliers differ, and the
numerator acquires

    mu (phi_hat_benchmark - phi_hat_candidate),

which is zero in expectation but not zero in any sample, and whose spread grows
linearly in the common mean.  That is a bias in the numerator that the
denominator does not know about, so the statistic inflates.  It is the same
shape as the Prais-Winsten defect the preceding programme found, arriving by a
different route: there the multiplier varied with the group's LENGTH, here it
varies with the group's ESTIMATE.

THE REPAIR.  Pool.  Estimate one coefficient from both groups and apply it to
both.  The multiplier is then identical by construction and mu cancels exactly,
whatever the estimate happens to be.  The residual error -- imperfect whitening,
which biases the variance estimate -- remains, and is what the mu = 0 rows
isolate.

This module measures separate against pooled at both tiers, and separates the
two error sources by running each at mu = 0 and mu = 2.

Deposit: deposits/w3_serial_repair.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L

ABAR = F.ABAR
LAM_DOMAIN = (1e-3, 1e3)
DRAWS = 400000


def ar1_paths(D, n, phi, sig, rng):
    e = rng.normal(size=(D, n))
    x = np.empty((D, n))
    x[:, 0] = e[:, 0]
    for t in range(1, n):
        x[:, t] = phi * x[:, t - 1] + np.sqrt(1.0 - phi ** 2) * e[:, t]
    return sig * x


def co_whiten(x, phi):
    phi = np.asarray(phi, float)
    if phi.ndim == 1:
        phi = phi[:, None]
    return (x[:, 1:] - phi * x[:, :-1]) / np.sqrt(1.0 - phi ** 2)


def acf1_parts(x):
    """Numerator and denominator of the lag-one autocorrelation, per row.

    Returned separately so they can be POOLED across groups before dividing,
    which is what makes one shared estimate rather than an average of two.
    """
    z = x - x.mean(1, keepdims=True)
    return (z[:, 1:] * z[:, :-1]).sum(1), (z * z).sum(1)


def welch(b, c):
    n0, n1 = b.shape[1], c.shape[1]
    s0, s1 = b.var(1, ddof=1), c.var(1, ddof=1)
    return (c.mean(1) - b.mean(1)) / np.sqrt(s1 / n1 + s0 / n0), s1 / s0


def cell(tab, n, phi, mu, tier, pooled, cal_len=0, lam=1.0, D=DRAWS, seed=1,
         batch=80000):
    rng = np.random.default_rng(int(seed))
    n0, n1 = int(n[0]), int(n[1])
    hits = done = 0
    while done < D:
        m = min(batch, D - done)
        b = ar1_paths(m, n0 + 1, phi, 1.0, rng)
        c = ar1_paths(m, n1 + 1, phi, np.sqrt(lam), rng)
        if tier == "A":
            pb = pc = phi
        else:
            if tier == "B":
                sb = ar1_paths(m, cal_len, phi, 1.0, rng)
                sc = ar1_paths(m, cal_len, phi, np.sqrt(lam), rng)
            else:
                sb, sc = b, c
            nb, db = acf1_parts(sb)
            nc, dc = acf1_parts(sc)
            if pooled:
                # ONE coefficient from both groups: the scales differ, so the
                # parts are standardised before pooling
                p = np.clip((nb / db + nc / dc) * 0.5, -0.95, 0.95)
                pb = pc = p
            else:
                pb = np.clip(nb / db, -0.95, 0.95)
                pc = np.clip(nc / dc, -0.95, 0.95)
        t, lh = welch(co_whiten(b + mu, pb), co_whiten(c + mu, pc))
        hits += int((t > tab(lh)).sum())
        done += m
    p = hits / D
    return p, float(np.sqrt(p * (1 - p) / D))


def main():
    t0 = time.time()
    out = {"abar": ABAR, "draws": DRAWS, "rows": []}
    designs = [("n36_36", [36.0, 36.0]), ("n36_13", [36.0, 13.0])]
    print("%-8s %5s %5s %-8s %-9s %9s %8s" %
          ("design", "phi", "mu", "tier", "estimate", "FWER", "z"))
    for dtag, n in designs:
        n = np.asarray(n, float)
        cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                                npanel=12, order=8, iters=10)
        tab = cal["table"]
        k = 0
        for phi in (0.2, 0.4, 0.6):
            for mu in (0.0, 2.0, 5.0):
                for tier, cl in (("B", 120), ("C", 0)):
                    for pooled in (False, True):
                        k += 1
                        p, se = cell(tab, n, phi, mu, tier, pooled, cl,
                                     seed=k + 17)
                        z = (p - ABAR) / se
                        lab = ("%s(%d)" % (tier, cl)) if tier == "B" else tier
                        out["rows"].append({
                            "design": dtag, "phi": float(phi), "mu": float(mu),
                            "tier": tier, "cal_len": int(cl),
                            "pooled": bool(pooled), "fwer": p, "mcse": se,
                            "z": float(z), "breach": bool(z > 3.0)})
                        print("%-8s %5.1f %5.1f %-8s %-9s %9.5f %+8.1f %s"
                              % (dtag, phi, mu, lab,
                                 "pooled" if pooled else "separate", p, z,
                                 "BREACH" if z > 3.0 else ""))

    def sub(**kw):
        r = out["rows"]
        for k, v in kw.items():
            r = [x for x in r if x[k] == v]
        return r

    sep = sub(pooled=False)
    poo = sub(pooled=True)
    sep_mu0 = [r for r in sep if r["mu"] == 0.0]
    poo_mu0 = [r for r in poo if r["mu"] == 0.0]
    sep_mu5 = [r for r in sep if r["mu"] == 5.0]
    poo_mu5 = [r for r in poo if r["mu"] == 5.0]

    # the signature of the repair: pooled must be FLAT in the common mean
    flat = {}
    for tier in ("B", "C"):
        vals = [max(r["fwer"] for r in poo if r["tier"] == tier
                    and r["mu"] == m) for m in (0.0, 2.0, 5.0)]
        flat[tier] = float(max(vals) - min(vals))

    out.update({
        "separate_worst_fwer": float(max(r["fwer"] for r in sep)),
        "pooled_worst_fwer": float(max(r["fwer"] for r in poo)),
        "separate_worst_ratio": float(max(r["fwer"] for r in sep) / ABAR),
        "pooled_worst_ratio": float(max(r["fwer"] for r in poo) / ABAR),
        "separate_worst_at_mu0": float(max(r["fwer"] for r in sep_mu0)),
        "pooled_worst_at_mu0": float(max(r["fwer"] for r in poo_mu0)),
        "separate_worst_at_mu5": float(max(r["fwer"] for r in sep_mu5)),
        "pooled_worst_at_mu5": float(max(r["fwer"] for r in poo_mu5)),
        "pooled_spread_over_mu_tierB": flat["B"],
        "pooled_spread_over_mu_tierC": flat["C"],
        "separate_breaches": int(sum(1 for r in sep if r["breach"])),
        "pooled_breaches": int(sum(1 for r in poo if r["breach"])),
        "n_per_arm": len(sep),
        "repair_removes_the_mean_effect": bool(
            max(flat.values()) < 0.004
            and max(r["fwer"] for r in poo)
            < max(r["fwer"] for r in sep) * 0.5),
        "residual_error_source": (
            "imperfect whitening at an estimated coefficient, which biases the "
            "variance estimate; it is what the mu = 0 rows isolate and pooling "
            "does not remove it"),
        "seconds": float(time.time() - t0),
    })
    F.write_deposit("w3_serial_repair.json", out, "measured")
    print("\nseparate: worst %.5f (%.1fx), %d of %d breach"
          % (out["separate_worst_fwer"], out["separate_worst_ratio"],
             out["separate_breaches"], out["n_per_arm"]))
    print("pooled:   worst %.5f (%.1fx), %d of %d breach"
          % (out["pooled_worst_fwer"], out["pooled_worst_ratio"],
             out["pooled_breaches"], out["n_per_arm"]))
    print("pooled spread over the common mean: B %.5f  C %.5f"
          % (flat["B"], flat["C"]))
    print("repair removes the mean effect: %s ; %.0f s"
          % (out["repair_removes_the_mean_effect"], out["seconds"]))


if __name__ == "__main__":
    main()
