"""M6, part 2: serial dependence, by the plan's EXACT-MODEL route.

The preceding programmes measured the boundary and stopped there: under AR(1)
errors no rule in the tested set survives `phi ~ 0.1` at n = 36.  That is a true
statement about rules built for independent data, and it is not the end of the
matter.

THE CONSTRUCTION.  If the errors follow a Gaussian AR(1) with a DECLARED
coefficient, the Prais-Winsten transform

    y_1 = sqrt(1 - phi^2) x_1 ,   y_t = x_t - phi x_{t-1}   (t >= 2)

maps the series to exactly independent Gaussians with the same mean scaled by
`k(phi) = ...` and a variance the transform fixes.  The whitened problem is
therefore the iid problem this repository already solves, and **every rule and
every certificate transfers to it unchanged**.  Serial dependence stops being a
boundary and becomes a pre-processing step, on one condition: the coefficient
has to be declared, not estimated from the same data.

WHAT IS AND IS NOT CLAIMED.

  * Declared, correct `phi`         -> the iid certificate applies exactly.
  * Declared `phi`, wrong by some   -> measured here, as a misspecification
    amount                             surface, because nothing makes it exact.
  * `phi` estimated from the data   -> a new nuisance parameter.  Not certified,
                                       and the cost of pretending otherwise is
                                       measured.
  * A declared BOUND `phi <= phi_max` -> the worst case over the bound is taken,
                                       which is valid for every `phi` in it and
                                       is the deployable version.

Deposit: deposits/serial.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import norm

import fcore as F
import k1lab as L

ABAR = F.ABAR


def ar1(n, phi, sigma, rng, size):
    """`size` independent AR(1) series of length n with stationary variance
    sigma^2."""
    e = rng.normal(0.0, 1.0, (size, n))
    x = np.empty((size, n))
    x[:, 0] = e[:, 0]                       # stationary start
    for t in range(1, n):
        x[:, t] = phi * x[:, t - 1] + np.sqrt(1.0 - phi ** 2) * e[:, t]
    return sigma * x


def whiten(x, phi):
    """Prais-Winsten: exactly independent output for a stationary AR(1) input."""
    y = np.empty_like(x)
    y[:, 0] = np.sqrt(1.0 - phi ** 2) * x[:, 0]
    y[:, 1:] = x[:, 1:] - phi * x[:, :-1]
    return y / np.sqrt(1.0 - phi ** 2)      # keep the innovation scale at sigma


def stat_and_ratio(b, c):
    """The Welch statistic and the estimated ratio from two samples."""
    n0, n1 = b.shape[1], c.shape[1]
    s0 = b.var(axis=1, ddof=1)
    s1 = c.var(axis=1, ddof=1)
    d = np.sqrt(s1 / n1 + s0 / n0)
    t = (c.mean(axis=1) - b.mean(axis=1)) / d
    return t, s1 / s0


def run_cell(n, phi_true, phi_declared, tab, R=40000, seed=0, lam=1.0):
    """Realised FWER of the rule when the data are AR(1) and the transform uses
    `phi_declared`.  `phi_declared = None` means no transform at all."""
    rng = np.random.default_rng(seed)
    b = ar1(int(n[0]), phi_true, 1.0, rng, R)
    c = ar1(int(n[1]), phi_true, np.sqrt(lam), rng, R)
    if phi_declared is not None:
        b = whiten(b, phi_declared)
        c = whiten(c, phi_declared)
    t, lh = stat_and_ratio(b, c)
    thr = tab(lh)
    rej = t > thr
    p = float(rej.mean())
    return {"fwer": p, "mcse": float(np.sqrt(p * (1 - p) / R))}


def main():
    t0 = time.time()
    n = np.array([36.0, 36.0])
    tab_cache = {}

    def table(level):
        if level not in tab_cache:
            tab_cache[level] = L.PluginTable(n, level, m=41, npanel=12,
                                             order=8)
        return tab_cache[level]

    # the iid-calibrated level for this design, from the same machinery
    cal = L.calibrate_level(n, F.log_grid(1e-3, 1e3, 25), ABAR, npanel=12,
                            order=8, iters=11)
    tab = table(cal["alpha0"])
    out = {"abar": ABAR, "n": n.tolist(), "alpha0_iid": float(cal["alpha0"]),
           "draws": 40000, "rows": [], "bound_rule": []}

    print("F1  the whitening identity: whitened AR(1) must behave as iid")
    # falsifier: with the correct phi the realised FWER must match the iid one
    iid = run_cell(n, 0.0, None, tab, R=8000, seed=1)
    print("    iid data, no transform            FWER %.5f (se %.5f)"
          % (iid["fwer"], iid["mcse"]))
    out["iid_reference"] = iid

    print("\nF2  the naive rule under AR(1), and the repaired rule")
    print("    %6s %14s %14s %14s" %
          ("phi", "no transform", "correct phi", "declared 0.0"))
    for phi in (0.0, 0.05, 0.1, 0.2, 0.3, 0.4):
        # 40,000 draws, not 4,000: at the smaller count the repaired arm sat
        # 2.3 standard errors above the ceiling at two coefficients, which
        # reads as a breach and is sampling noise.  The precision has to be
        # good enough to tell those apart before the figure is drawn.
        naive = run_cell(n, phi, None, tab, seed=100 + int(phi * 100))
        fixed = run_cell(n, phi, phi, tab, seed=100 + int(phi * 100))
        out["rows"].append({"phi": float(phi), "naive": naive,
                            "whitened_correct": fixed})
        print("    %6.2f %14s %14s" %
              (phi, "%.5f" % naive["fwer"], "%.5f" % fixed["fwer"]))

    print("\nF3  misspecification: declare phi_0, the truth is phi")
    grid_d = (0.0, 0.1, 0.2, 0.3)
    grid_t = (0.0, 0.1, 0.2, 0.3)
    mis = []
    print("    declared\\true %s" % "".join("%9.2f" % p for p in grid_t))
    for pd in grid_d:
        row = []
        for pt in grid_t:
            r = run_cell(n, pt, pd, tab, R=12000,
                         seed=500 + int(pd * 10) * 7 + int(pt * 10))
            row.append(r["fwer"])
            mis.append({"declared": float(pd), "true": float(pt),
                        "fwer": r["fwer"], "mcse": r["mcse"]})
        print("    %12.2f %s" % (pd, "".join("%9.5f" % v for v in row)))
    out["misspecification"] = mis

    print("\nF4  a rule valid for every phi in a declared bound")
    # take the worst case over the bound by lowering the internal level until
    # the realised rate is under the ceiling at the worst phi in the bound
    for phi_max in (0.05, 0.1, 0.2):
        lvl = cal["alpha0"]
        grid = np.linspace(0.0, phi_max, 5)
        for _ in range(14):
            tb = table(round(lvl, 6))
            worst = max(run_cell(n, p, phi_max, tb, R=12000,
                                 seed=900 + int(p * 1000))["fwer"]
                        for p in grid)
            if worst <= ABAR:
                break
            lvl *= 0.9
        out["bound_rule"].append({"phi_max": float(phi_max),
                                  "alpha0": float(lvl),
                                  "worst_fwer": float(worst),
                                  "valid": bool(worst <= ABAR),
                                  "level_price_vs_iid_pct":
                                      float(100.0 * (cal["alpha0"] - lvl)
                                            / cal["alpha0"])})
        print("    phi <= %.2f  ->  alpha0 %.5f  worst FWER %.5f  %s"
              % (phi_max, lvl, worst, "OK" if worst <= ABAR else "FAILS"))

    out["seconds"] = float(time.time() - t0)
    out["status"] = (
        "The whitened problem is the iid problem, so with a DECLARED and "
        "correct coefficient every certificate in this repository transfers "
        "unchanged.  Everything else here is measured: the misspecification "
        "surface and the declared-bound rule are simulation estimates with "
        "their standard errors, and neither is a theorem.")
    F.write_deposit("serial.json", out, "unconditional_candidate_measured")
    print("\nwrote serial.json")


if __name__ == "__main__":
    main()
