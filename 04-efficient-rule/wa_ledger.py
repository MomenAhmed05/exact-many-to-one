"""Workstream A: the efficiency ledger, with every object labelled.

The plan's first rule is that no comparison may be called a price unless the two
objects condition on the same data and target the same rejection event.  At
K = 1 every object below is exact in closed form, so the ledger carries no Monte
Carlo error at all:

  oracle_conditional      c = z_{1-abar} * sigma(lambda_true) / d .  Needs the
                          truth; the ceiling on what any rule can do.
  plugin_measured         the exact fixed-ratio constant at lambda_hat, at the
                          full ceiling.  Not valid: its sup over lambda exceeds
                          abar (measured here).
  region_sup_conditional  the 03-certified-supremum construction, specialised to K = 1: build an
                          exact 1 - gamma interval for lambda and take the
                          supremum of the conditional oracle constant over it.
                          Because c_cond is increasing in lambda, that supremum
                          is attained at the interval's upper endpoint exactly,
                          so the region rule can be evaluated here with no
                          branch-and-bound and no approximation.
  uncond_certified        this programme's rule: the plug-in shape at the
                          largest internal level whose unconditional FWER is
                          certified at or below abar over the whole ratio
                          domain.

Reporting all four on one draw set, at one set of alternatives, is the ledger.

Deposit: deposits/wa_ledger.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import f as fdist
from scipy.stats import norm

import fcore as F

GAMMA = 5e-4                      # the region rule's miss budget, as in 03-certified-supremum


def region_interval(lam_hat, n, gamma=GAMMA):
    """Exact 1 - gamma interval for lambda from lambda_hat.

    lambda_hat / lambda = u_1/u_0 ~ F(nu_1, nu_0), so inverting the two tails
    gives an interval whose coverage is exactly 1 - gamma.  No quadrature.
    """
    n = np.asarray(n, float)
    nu1, nu0 = n[1] - 1.0, n[0] - 1.0
    lo_q = fdist.ppf(gamma / 2.0, nu1, nu0)
    hi_q = fdist.ppf(1.0 - gamma / 2.0, nu1, nu0)
    return lam_hat / hi_q, lam_hat / lo_q


def c_region(u, n, gamma=GAMMA, abar=F.ABAR):
    """The region rule's conditional threshold at K = 1, exactly.

    c_cond(lambda) = z_{1-alpha0} sigma(lambda)/d is increasing in lambda, so the
    supremum over the admitted interval sits at its upper endpoint.
    """
    n = np.asarray(n, float)
    u = np.asarray(u, float)
    alpha0 = abar - gamma
    lam_hat = u[1] / u[0]                     # sigma cancels: this is S1^2/S0^2
    _, hi = region_interval(lam_hat, n, gamma)
    d = np.sqrt(u[1] / n[1] + u[0] / n[0])    # placeholder; caller rescales
    return alpha0, hi, d


def ledger_design(name, n, lam_true_list, certified_rule, plugin_tab,
                  gamma=GAMMA, abar=F.ABAR, npanel=24, order=12):
    n = np.asarray(n, float)
    z = norm.ppf(1.0 - abar)
    z0 = norm.ppf(1.0 - (abar - gamma))
    q0, w0 = F.state_nodes(n[0] - 1.0, npanel, order)
    q1, w1 = F.state_nodes(n[1] - 1.0, npanel, order)
    U0, U1 = np.meshgrid(q0, q1, indexing="ij")
    W = np.outer(w0, w1)
    W = W / W.sum()
    nu1, nu0 = n[1] - 1.0, n[0] - 1.0
    hi_q = fdist.ppf(1.0 - gamma / 2.0, nu1, nu0)
    lo_q = fdist.ppf(gamma / 2.0, nu1, nu0)

    rows = []
    for lam in lam_true_list:
        d = np.sqrt(lam * U1 / n[1] + U0 / n[0])
        sig = np.sqrt(lam / n[1] + 1.0 / n[0])
        lh = lam * U1 / U0

        c_or = z * sig / d                                   # oracle conditional
        lam_hi = lh / lo_q                                   # region upper end
        c_reg = z0 * np.sqrt(lam_hi / n[1] + 1.0 / n[0]) / d
        c_plug = plugin_tab(lh.ravel()).reshape(lh.shape)
        c_new = certified_rule(lh.ravel()).reshape(lh.shape)

        def fw(c):
            return float((W * norm.cdf(-(c * d) / sig)).sum())

        def pw(c, delta):
            return float((W * norm.cdf((delta - c * d) / sig)).sum())

        d80 = F.delta_for_oracle_power_k1(n, [lam], 0.8, abar)
        d50 = F.delta_for_oracle_power_k1(n, [lam], 0.5, abar)
        ent = {"lam": float(lam), "delta80": float(d80)}
        for lab, c in (("oracle_conditional", c_or),
                       ("plugin_measured", c_plug),
                       ("region_sup_conditional", c_reg),
                       ("uncond_certified", c_new)):
            ent[lab] = {
                "mean_c": float((W * c).sum()),
                "median_c": float(np.median(c)),
                "max_c": float(c.max()),
                "ratio_to_oracle_mean": float((W * (c / c_or)).sum()),
                "fwer": fw(c),
                "power_at_oracle80": pw(c, d80),
                "power_at_oracle50": pw(c, d50),
                "loss_pp_at_80": float(100.0 * (0.8 - pw(c, d80))),
            }
        rows.append(ent)
    return {"design": name, "n": n.tolist(), "gamma": gamma, "abar": abar,
            "rows": rows}


def main():
    import k1lab as L
    import wb_certified_rule as R

    designs = [("df3", [4.0, 13.0]), ("df5", [6.0, 13.0]),
               ("df7", [8.0, 13.0]), ("df11", [12.0, 13.0]),
               ("df19", [20.0, 13.0]), ("df35", [36.0, 13.0])]
    lams = (0.04, 1.0, 25.0)
    cert = None
    try:
        cert = F.read_deposit("wb_certified_rule.json")
    except Exception:  # noqa: BLE001
        pass
    out = {"abar": F.ABAR, "gamma": GAMMA, "lams": list(lams), "designs": []}
    print("%-7s %6s | %-22s %8s %8s %8s %8s" %
          ("design", "lam", "object", "meanC", "xOracle", "FWER", "loss80"))
    for name, n in designs:
        nn = np.asarray(n, float)
        edges = R.step_edges()
        if cert:
            rec = [d for d in cert["designs"] if d["design"] == name]
        else:
            rec = []
        if rec:
            rule = __import__("wb_certify").StepRule(rec[0]["c"],
                                                     rec[0]["edges"])
            a0 = rec[0]["alpha0"]
        else:
            a0 = L.calibrate_level(nn, F.log_grid(1e-3, 1e3, 25), F.ABAR,
                                   npanel=12, order=8, iters=9)["alpha0"]
            rule = R.rule_from_level(nn, a0, edges)
        tab = L.PluginTable(nn, F.ABAR, m=41, npanel=12, order=8)
        r = ledger_design(name, nn, lams, rule, tab)
        r["alpha0_certified"] = float(a0)
        out["designs"].append(r)
        for row in r["rows"]:
            for lab in ("oracle_conditional", "plugin_measured",
                        "region_sup_conditional", "uncond_certified"):
                e = row[lab]
                print("%-7s %6g | %-22s %8.3f %8.2f %8.5f %8.1f" %
                      (name, row["lam"], lab, e["mean_c"],
                       e["ratio_to_oracle_mean"], e["fwer"],
                       e["loss_pp_at_80"]))
    F.write_deposit("wa_ledger.json", out, "frontier_measured")
    print("\nwrote wa_ledger.json")


if __name__ == "__main__":
    main()
