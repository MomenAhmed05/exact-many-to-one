"""W3 stress tests: deliberately break every assumption and report where.

The certificate is a normal-theory object.  It assumes independent Gaussian
groups with variance estimates of the usual chi-square form.  A desk deploying
it is entitled to know which violations it survives, which it does not, and by
how much -- measured, not asserted, and with the direction of the breach stated
because a conservative failure and a liberal one are not the same news.

Every cell runs the SAME deployed rule: the calibrated table at the internal
level that makes the iid rule uniformly valid over the declared ratio domain.
Only the data-generating process changes.

The battery, each declared before it was run:

  serial          AR(1) from strongly negative to strongly positive, untreated,
                  and treated with the repaired (first-observation-dropped)
                  whitening at a declared coefficient;
  garch           GARCH(1,1) innovations, stationary, unit unconditional
                  variance, so only the CONDITIONAL heteroskedasticity differs;
  heavy tails     Student t on 3 and 5 degrees of freedom, standardised;
  skew            standardised centred chi-square on 3 df;
  regime          variance doubles at the midpoint of the sample;
  short variance  the variance estimated on a trailing subwindow while the mean
                  uses the whole sample, which is what a rolling risk model does;
  factor          a second, undeclared factor left in after the declared one is
                  residualised, and residual correlation left after removal.

Deposit: deposits/w3_stress.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L

ABAR = F.ABAR
LAM_DOMAIN = (1e-3, 1e3)
DRAWS = 400000
BREACH_SE = 3.0          # declared: a breach is 3 standard errors above abar


# ----------------------------------------------------------------- generators
def g_normal(D, n, sig, rng):
    return sig * rng.normal(size=(D, n))


def g_t(df):
    def f(D, n, sig, rng):
        x = rng.standard_t(df, size=(D, n))
        return sig * x / np.sqrt(df / (df - 2.0))
    return f


def g_skew(df=3):
    def f(D, n, sig, rng):
        x = rng.chisquare(df, size=(D, n))
        return sig * (x - df) / np.sqrt(2.0 * df)
    return f


def g_garch(alpha=0.10, beta=0.85):
    """GARCH(1,1) with unit unconditional variance."""
    omega = 1.0 - alpha - beta

    def f(D, n, sig, rng):
        z = rng.normal(size=(D, n))
        h = np.ones(D)
        x = np.empty((D, n))
        for t in range(n):
            x[:, t] = np.sqrt(h) * z[:, t]
            h = omega + alpha * x[:, t] ** 2 + beta * h
        return sig * x
    return f


def g_regime(ratio=2.0):
    """Variance steps up at the midpoint; unconditional variance held at sig^2."""
    def f(D, n, sig, rng):
        h = n // 2
        v = np.concatenate([np.ones(h), np.full(n - h, ratio)])
        v = v / v.mean()
        return sig * rng.normal(size=(D, n)) * np.sqrt(v)[None, :]
    return f


def g_ar1(phi, repair=False):
    """Stationary AR(1); `repair` applies the Cochrane-Orcutt transform."""
    def f(D, n, sig, rng):
        m = n + 1 if repair else n
        e = rng.normal(size=(D, m))
        x = np.empty((D, m))
        x[:, 0] = e[:, 0]
        for t in range(1, m):
            x[:, t] = phi * x[:, t - 1] + np.sqrt(1.0 - phi ** 2) * e[:, t]
        x = sig * x
        if not repair:
            return x
        y = x[:, 1:] - phi * x[:, :-1]          # drop y_1, keep n observations
        return y / np.sqrt(1.0 - phi ** 2)      # hold the innovation scale
    return f


# ------------------------------------------------------------------- the rule
def welch(b, c, short=0):
    n0, n1 = b.shape[1], c.shape[1]
    if short:
        s0 = b[:, -short:].var(1, ddof=1)
        s1 = c[:, -short:].var(1, ddof=1)
    else:
        s0 = b.var(1, ddof=1)
        s1 = c.var(1, ddof=1)
    t = (c.mean(1) - b.mean(1)) / np.sqrt(s1 / n1 + s0 / n0)
    return t, s1 / s0


def measure(tab, n, lam, gen, D=DRAWS, seed=0, short=0, batch=100000):
    rng = np.random.default_rng(seed)
    hits = done = 0
    while done < D:
        m = min(batch, D - done)
        b = gen(m, int(n[0]), 1.0, rng)
        c = gen(m, int(n[1]), np.sqrt(lam), rng)
        t, lh = welch(b, c, short)
        hits += int((t > tab(lh)).sum())
        done += m
    p = hits / D
    return p, float(np.sqrt(p * (1 - p) / D))


# --------------------------------------------------------------- factor cases
def factor_cells(tab, n, lam, D=DRAWS, seed=7000):
    """Undeclared second factor, and residual correlation after removal."""
    rng = np.random.default_rng(seed)
    out = []
    for tag, b2, rho in [("second undeclared factor", 0.7, 0.0),
                         ("residual correlation 0.3", 0.0, 0.3),
                         ("undeclared factor and residual correlation",
                          0.7, 0.3)]:
        hits = done = 0
        while done < D:
            m = min(100000, D - done)
            n0, n1 = int(n[0]), int(n[1])
            nn = max(n0, n1)
            g2 = rng.normal(size=(m, nn))                 # the undeclared factor
            w = rng.normal(size=(m, nn))                  # the shared residual
            e0 = rng.normal(size=(m, nn))
            e1 = rng.normal(size=(m, nn))
            r0 = b2 * g2 + np.sqrt(rho) * w + np.sqrt(1 - rho) * e0
            r1 = np.sqrt(lam) * (b2 * g2 + np.sqrt(rho) * w
                                 + np.sqrt(1 - rho) * e1)
            t, lh = welch(r0[:, :n0], r1[:, :n1])
            hits += int((t > tab(lh)).sum())
            done += m
        p = hits / D
        se = float(np.sqrt(p * (1 - p) / D))
        out.append({"cell": tag, "fwer": p, "mcse": se,
                    "z": float((p - ABAR) / se),
                    "breach": bool((p - ABAR) / se > BREACH_SE)})
    return out


def main():
    t0 = time.time()
    designs = [("n36_36", [36.0, 36.0]), ("n36_13", [36.0, 13.0])]
    lams = (1.0, 4.0)
    battery = [
        ("normal (control)", g_normal, 0),
        ("t5", g_t(5), 0),
        ("t3", g_t(3), 0),
        ("skew chi2_3", g_skew(3), 0),
        ("garch(.10,.85)", g_garch(), 0),
        ("regime x2", g_regime(2.0), 0),
        ("ar1 phi=-0.4 raw", g_ar1(-0.4), 0),
        ("ar1 phi=+0.2 raw", g_ar1(0.2), 0),
        ("ar1 phi=+0.4 raw", g_ar1(0.4), 0),
        ("ar1 phi=+0.4 repaired", g_ar1(0.4, repair=True), 0),
        ("ar1 phi=+0.6 repaired", g_ar1(0.6, repair=True), 0),
        ("short variance window", g_normal, 12),
    ]

    out = {"abar": ABAR, "draws": DRAWS, "breach_se": BREACH_SE,
           "designs": []}
    for dtag, n in designs:
        n = np.asarray(n, float)
        cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                                npanel=12, order=8, iters=10)
        tab = cal["table"]
        print("\n%s  alpha0 = %.5f" % (dtag, cal["alpha0"]))
        rows = []
        for k, (tag, gen, short) in enumerate(battery):
            for lam in lams:
                p, se = measure(tab, n, lam, gen, seed=1000 * k + int(lam),
                                short=short)
                z = (p - ABAR) / se
                rows.append({"cell": tag, "lam": float(lam), "fwer": p,
                             "mcse": se, "z": float(z),
                             "breach": bool(z > BREACH_SE)})
                print("  %-24s lam %-4.1f  FWER %.5f (%.5f)  z %+6.1f  %s"
                      % (tag, lam, p, se, z, "BREACH" if z > BREACH_SE else ""))
        frows = factor_cells(tab, n, 1.0)
        for r in frows:
            print("  %-24s lam 1.0   FWER %.5f (%.5f)  z %+6.1f  %s"
                  % (r["cell"], r["fwer"], r["mcse"], r["z"],
                     "BREACH" if r["breach"] else ""))
        out["designs"].append({
            "design": dtag, "n": n.tolist(), "alpha0": float(cal["alpha0"]),
            "rows": rows, "factor_rows": frows})

    allr = [r for d in out["designs"] for r in d["rows"] + d["factor_rows"]]
    br = [r for r in allr if r["breach"]]
    out["n_cells"] = len(allr)
    out["n_breaches"] = len(br)
    out["worst_fwer"] = float(max(r["fwer"] for r in allr))
    out["worst_cell"] = max(allr, key=lambda r: r["fwer"])["cell"]
    out["worst_ratio_to_abar"] = float(out["worst_fwer"] / ABAR)
    out["breaching_cells"] = sorted(set(r["cell"] for r in br))
    out["surviving_cells"] = sorted(set(r["cell"] for r in allr
                                        if not r["breach"]))
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("w3_stress.json", out, "diagnostic")
    print("\n%d of %d cells breach; worst %.5f (%.2fx abar) at '%s'"
          % (len(br), len(allr), out["worst_fwer"],
             out["worst_ratio_to_abar"], out["worst_cell"]))
    print("breaching:", ", ".join(out["breaching_cells"]) or "none")
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
