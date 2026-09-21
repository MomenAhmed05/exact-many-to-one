"""W3: the three validity tiers for a declared-nuisance transformation.

Section 11 treats serial correlation and factor exposure as nuisances a desk
DECLARES and removes before the statistic is formed.  That is honest, and it is
also the easiest of three cases.  This module separates them and measures each,
because "declared" is a word doing a great deal of work and the difference
between the tiers is the difference between a certificate and a hope.

  TIER A  parameters known.  The transform is exact, the transformed groups are
          iid, and every certificate transfers at the shortened count.  This is
          the tier Section 11 establishes.

  TIER B  parameters estimated on a DISJOINT calibration period, then frozen.
          The transform is now random but independent of the evaluation data,
          so the conditional interpretation survives; what it costs is a
          function of the calibration length, measured here.

  TIER C  parameters estimated on the evaluation period itself.  The transform
          is a function of the data being tested.  This is the realistic case
          and it is the one with no guarantee.

The AR(1) transform is the repaired one: drop the first observation rather than
rescale it, which is what makes the whitened mean vector constant and is the
defect the preceding programme's falsifier located.

Deposit: deposits/w3_tiers.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L

ABAR = F.ABAR
LAM_DOMAIN = (1e-3, 1e3)
DRAWS = 300000


def ar1_paths(D, n, phi, sig, rng, burn=0):
    """Stationary AR(1), length n + burn, unconditional variance sig^2."""
    m = n + burn
    e = rng.normal(size=(D, m))
    x = np.empty((D, m))
    x[:, 0] = e[:, 0]
    for t in range(1, m):
        x[:, t] = phi * x[:, t - 1] + np.sqrt(1.0 - phi ** 2) * e[:, t]
    return sig * x


def co_whiten(x, phi):
    """Cochrane-Orcutt: y_t = x_t - phi x_{t-1}, t >= 2; y_1 discarded.

    `phi` may be a scalar or one value per row.
    """
    phi = np.asarray(phi, float)
    if phi.ndim == 1:
        phi = phi[:, None]
    y = x[:, 1:] - phi * x[:, :-1]
    return y / np.sqrt(1.0 - phi ** 2)


def phi_hat(x):
    """Lag-one autocorrelation, per row, clipped to a stationary range."""
    z = x - x.mean(1, keepdims=True)
    num = (z[:, 1:] * z[:, :-1]).sum(1)
    den = (z * z).sum(1)
    return np.clip(num / np.maximum(den, 1e-12), -0.95, 0.95)


def welch(b, c):
    n0, n1 = b.shape[1], c.shape[1]
    s0, s1 = b.var(1, ddof=1), c.var(1, ddof=1)
    t = (c.mean(1) - b.mean(1)) / np.sqrt(s1 / n1 + s0 / n0)
    return t, s1 / s0


def serial_tier(tab, n, phi, lam, tier, cal_len=0, D=DRAWS, seed=0,
                mu=0.0, batch=75000):
    """Realised FWER under AR(1) at one validity tier."""
    rng = np.random.default_rng(int(seed))
    n0, n1 = int(n[0]), int(n[1])
    hits = done = 0
    while done < D:
        m = min(batch, D - done)
        # one extra observation so the dropped one leaves n0 / n1 behind
        b = ar1_paths(m, n0 + 1, phi, 1.0, rng)
        c = ar1_paths(m, n1 + 1, phi, np.sqrt(lam), rng)
        if tier == "A":
            pb, pc = phi, phi
        elif tier == "B":
            cb = ar1_paths(m, cal_len, phi, 1.0, rng)
            cc = ar1_paths(m, cal_len, phi, np.sqrt(lam), rng)
            pb, pc = phi_hat(cb), phi_hat(cc)
        else:                                     # tier C
            pb, pc = phi_hat(b), phi_hat(c)
        bw = co_whiten(b + mu, pb)
        cw = co_whiten(c + mu, pc)
        t, lh = welch(bw, cw)
        hits += int((t > tab(lh)).sum())
        done += m
    p = hits / D
    return p, float(np.sqrt(p * (1 - p) / D))


def factor_tier(tab, n, beta, lam, tier, cal_len=0, D=DRAWS, seed=0,
                batch=75000):
    """Realised FWER under an OBSERVED one-factor model at one tier.

    The factor is observed here (Section 11.2's setting); the latent case is
    w2_factor.py.  Tier A residualises with the true loading, tier B with a
    loading fitted on a disjoint pre-period, tier C with one fitted on the
    evaluation returns themselves.
    """
    rng = np.random.default_rng(int(seed))
    n0, n1 = int(n[0]), int(n[1])
    nn = max(n0, n1)
    b0, b1 = beta

    def fit(y, f):
        fz = f - f.mean(1, keepdims=True)
        yz = y - y.mean(1, keepdims=True)
        return (yz * fz).sum(1) / np.maximum((fz * fz).sum(1), 1e-12)

    hits = done = 0
    while done < D:
        m = min(batch, D - done)
        f = rng.normal(size=(m, nn))
        r0 = b0 * f[:, :n0] + rng.normal(size=(m, n0))
        r1 = b1 * f[:, :n1] + np.sqrt(lam) * rng.normal(size=(m, n1))
        if tier == "A":
            e0 = r0 - b0 * f[:, :n0]
            e1 = r1 - b1 * f[:, :n1]
        elif tier == "B":
            fp = rng.normal(size=(m, cal_len))
            p0 = b0 * fp + rng.normal(size=(m, cal_len))
            p1 = b1 * fp + np.sqrt(lam) * rng.normal(size=(m, cal_len))
            h0, h1 = fit(p0, fp), fit(p1, fp)
            e0 = r0 - h0[:, None] * f[:, :n0]
            e1 = r1 - h1[:, None] * f[:, :n1]
        else:
            h0, h1 = fit(r0, f[:, :n0]), fit(r1, f[:, :n1])
            e0 = r0 - h0[:, None] * f[:, :n0]
            e1 = r1 - h1[:, None] * f[:, :n1]
        t, lh = welch(e0, e1)
        hits += int((t > tab(lh)).sum())
        done += m
    p = hits / D
    return p, float(np.sqrt(p * (1 - p) / D))


def main():
    t0 = time.time()
    out = {"abar": ABAR, "draws": DRAWS, "serial": [], "factor": []}
    for dtag, n in [("n36_36", [36.0, 36.0]), ("n36_13", [36.0, 13.0])]:
        n = np.asarray(n, float)
        cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                                npanel=12, order=8, iters=10)
        tab = cal["table"]
        print("\n%s  alpha0 = %.5f" % (dtag, cal["alpha0"]))

        print(" SERIAL (AR(1), Cochrane-Orcutt)")
        for phi in (0.2, 0.4, 0.6):
            for mu in (0.0, 2.0):
                for tier, cl in [("A", 0), ("B", 36), ("B", 120), ("B", 360),
                                 ("C", 0)]:
                    p, se = serial_tier(tab, n, phi, 1.0, tier, cl, mu=mu,
                                        seed=abs(int(100 * phi + 10 * mu + cl)) + 1)
                    z = (p - ABAR) / se
                    out["serial"].append({
                        "design": dtag, "phi": float(phi), "mu": float(mu),
                        "tier": tier, "cal_len": int(cl), "fwer": p,
                        "mcse": se, "z": float(z), "breach": bool(z > 3.0)})
                    lab = tier if tier != "B" else "B(%d)" % cl
                    print("   phi %.1f mu %.1f tier %-6s FWER %.5f  z %+6.1f %s"
                          % (phi, mu, lab, p, z, "BREACH" if z > 3.0 else ""))

        print(" FACTOR (observed one-factor)")
        for beta in ([1.0, 1.0], [0.5, 1.5], [1.0, -0.5]):
            for tier, cl in [("A", 0), ("B", 36), ("B", 120), ("B", 360),
                             ("C", 0)]:
                p, se = factor_tier(tab, n, beta, 1.0, tier, cl,
                                    seed=abs(int(1000 * beta[1] + cl)) + 1)
                z = (p - ABAR) / se
                out["factor"].append({
                    "design": dtag, "beta": [float(x) for x in beta],
                    "tier": tier, "cal_len": int(cl), "fwer": p, "mcse": se,
                    "z": float(z), "breach": bool(z > 3.0)})
                lab = tier if tier != "B" else "B(%d)" % cl
                print("   beta %-11s tier %-6s FWER %.5f  z %+6.1f %s"
                      % (str(beta), lab, p, z, "BREACH" if z > 3.0 else ""))

    for key in ("serial", "factor"):
        rows = out[key]
        for tier in ("A", "B", "C"):
            sub = [r for r in rows if r["tier"] == tier]
            out["%s_tier%s_worst_fwer" % (key, tier)] = \
                float(max(r["fwer"] for r in sub))
            out["%s_tier%s_breaches" % (key, tier)] = \
                int(sum(1 for r in sub if r["breach"]))
        for cl in (36, 120, 360):
            sub = [r for r in rows if r["tier"] == "B" and r["cal_len"] == cl]
            out["%s_tierB_cal%d_worst_fwer" % (key, cl)] = \
                float(max(r["fwer"] for r in sub))
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("w3_tiers.json", out, "measured")
    print("\nserial   A %.5f  B(360) %.5f  C %.5f"
          % (out["serial_tierA_worst_fwer"],
             out["serial_tierB_cal360_worst_fwer"],
             out["serial_tierC_worst_fwer"]))
    print("factor   A %.5f  B(360) %.5f  C %.5f"
          % (out["factor_tierA_worst_fwer"],
             out["factor_tierB_cal360_worst_fwer"],
             out["factor_tierC_worst_fwer"]))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
