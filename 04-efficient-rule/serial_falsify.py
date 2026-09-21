"""Does the whitened problem really reduce to the independent one?

Section 11.1 claims the Prais-Winsten transform returns the problem Sections 3
and 4 solve, so every certificate transfers.  Table 25 supports that, but it was
run at ONE cell: mean zero and equal group sizes.  The transform's own mean
vector is

    mu * k(phi),    k = (sqrt(1 - phi^2), 1 - phi, ..., 1 - phi),

which is NOT constant in t.  The deployed statistic takes the plain mean of the
whitened series, so it estimates `mu * mean(k)`, and `mean(k)` depends on the
group's LENGTH.  Two predictions follow, neither of them visible in the cell
Table 25 ran:

  * with n0 == n1 the two groups share mean(k), the difference of means is
    unbiased for zero under the null whatever mu is, and the rule should hold
    the ceiling exactly;
  * with n0 != n1 the numerator carries a bias proportional to mu, so the rule
    should drift off the ceiling as mu grows.

Each design is calibrated at its own (n0, n1) here, so a deviation cannot be an
artefact of running one design's table on another's data.  That confound is
what made the first version of this script uninterpretable.

Deposit: deposits/serial_falsify.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L
import serial as S

ABAR = F.ABAR
DRAWS = 40000


def kvec(n, phi):
    """The transform's mean vector, up to mu."""
    k = np.full(int(n), 1.0 - phi)
    k[0] = np.sqrt(1.0 - phi ** 2)
    return k


def whiten_co(x, phi):
    """Cochrane-Orcutt: drop the first observation, keep the rest.

    y_t = x_t - phi x_{t-1} for t >= 2.  The mean vector is then
    mu (1 - phi) at every t, so the transformed sample is a location model on
    n - 1 observations and the group's length does not enter it.
    """
    return (x[:, 1:] - phi * x[:, :-1]) / np.sqrt(1.0 - phi ** 2)


def run(n0, n1, phi, mu, tab, R=DRAWS, seed=0, co_tab=None):
    """Realised FWER of the deployed rule on whitened AR(1) with mean mu."""
    rng = np.random.default_rng(seed)
    b = S.ar1(int(n0), phi, 1.0, rng, R) + mu
    c = S.ar1(int(n1), phi, 1.0, rng, R) + mu
    t, lh = S.stat_and_ratio(S.whiten(b, phi), S.whiten(c, phi))
    p = float((t > tab(lh)).mean())
    m0 = float(kvec(n0, phi).mean())
    m1 = float(kvec(n1, phi).mean())
    out = {"n": [float(n0), float(n1)], "phi": float(phi), "mu": float(mu),
           "fwer": p, "mcse": float(np.sqrt(p * (1 - p) / R)),
           "mean_k": [m0, m1], "numerator_bias_per_mu": m1 - m0}
    if co_tab is not None:
        tc, lhc = S.stat_and_ratio(whiten_co(b, phi), whiten_co(c, phi))
        pc = float((tc > co_tab(lhc)).mean())
        out["fwer_co"] = pc
        out["mcse_co"] = float(np.sqrt(pc * (1 - pc) / R))
    return out


def main():
    t0 = time.time()
    out = {"abar": ABAR, "draws": DRAWS, "rows": [],
           "note": ("each design is calibrated at its own (n0, n1); the "
                    "statistic is the deployed one, the plain mean and "
                    "variance of the whitened series")}
    print("%-12s %5s %6s %10s %10s %10s %s"
          % ("n", "phi", "mu", "plain", "repaired", "se",
             "bias in the numerator per mu"))
    for n0, n1 in ((36, 36), (36, 13)):
        n = np.array([float(n0), float(n1)])
        cal = L.calibrate_level(n, F.log_grid(1e-3, 1e3, 25), ABAR,
                                npanel=12, order=8, iters=11)
        tab = L.PluginTable(n, cal["alpha0"], m=41, npanel=12, order=8)
        nco = np.array([float(n0 - 1), float(n1 - 1)])
        calc = L.calibrate_level(nco, F.log_grid(1e-3, 1e3, 25), ABAR,
                                 npanel=12, order=8, iters=11)
        co_tab = L.PluginTable(nco, calc["alpha0"], m=41, npanel=12, order=8)
        for phi in (0.0, 0.2, 0.4):
            for mu in (0.0, 0.5, 2.0):
                r = run(n0, n1, phi, mu, tab, seed=11 + int(phi * 10) + n1,
                        co_tab=co_tab)
                out["rows"].append(r)
                print("%-12s %5.2f %6.2f %10.5f %10.5f %10.5f   %+.4f"
                      % ("[%d, %d]" % (n0, n1), phi, mu, r["fwer"],
                         r["fwer_co"], r["mcse"],
                         r["numerator_bias_per_mu"]))

    bal = [r for r in out["rows"] if r["n"][0] == r["n"][1]]
    unb = [r for r in out["rows"] if r["n"][0] != r["n"][1]]
    out["balanced_worst_abs_dev"] = float(
        max(abs(r["fwer"] - ABAR) for r in bal))
    out["balanced_worst_z"] = float(
        max(abs(r["fwer"] - ABAR) / r["mcse"] for r in bal))
    out["unbalanced_worst_abs_dev"] = float(
        max(abs(r["fwer"] - ABAR) for r in unb))
    out["unbalanced_worst_z"] = float(
        max(abs(r["fwer"] - ABAR) / r["mcse"] for r in unb))
    out["unbalanced_drift_with_mu"] = float(
        max(r["fwer"] for r in unb if r["mu"] > 0)
        - min(r["fwer"] for r in unb if r["mu"] == 0))
    out["repair"] = "Cochrane-Orcutt: drop the first observation, so the " \
                    "transformed mean vector is constant at 1 - phi and the " \
                    "group's length does not enter it"
    out["unbalanced_repaired_worst_abs_dev"] = float(
        max(abs(r["fwer_co"] - ABAR) for r in unb))
    out["unbalanced_repaired_worst_z"] = float(
        max(abs(r["fwer_co"] - ABAR) / r["mcse_co"] for r in unb))
    out["repair_works"] = bool(out["unbalanced_repaired_worst_z"] < 3.0)
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("serial_falsify.json", out, "frontier_measured")
    print("\nbalanced   worst |dev| %.5f  (%.1f se)"
          % (out["balanced_worst_abs_dev"], out["balanced_worst_z"]))
    print("unbalanced worst |dev| %.5f  (%.1f se)"
          % (out["unbalanced_worst_abs_dev"], out["unbalanced_worst_z"]))
    print("unbalanced, repaired  %.5f  (%.1f se)   repair works: %s"
          % (out["unbalanced_repaired_worst_abs_dev"],
             out["unbalanced_repaired_worst_z"], out["repair_works"]))


if __name__ == "__main__":
    main()
