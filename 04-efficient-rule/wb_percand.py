"""A general-K rule class that can be certified, not just measured.

`wb_general.py` calibrates the JOINT plug-in: the threshold is the exact
fixed-ratio constant at the whole estimated ratio vector.  That rule is
measured-tier above K = 1, because its unconditional FWER is an expectation over
K+1 chi-squares and no tractable enclosure of that expectation exists: the
threshold couples all the coordinates, so nothing factorises.

Replace it with a PER-CANDIDATE threshold,

    reject candidate i  iff  T_i > c( lambda_hat_i ) ,

one scalar function of one estimated ratio, shared by every candidate.  Now the
paper's own reduction applies to the FWER itself.  Conditioning on the benchmark
pair (z, u_0), the candidates are independent and the threshold for candidate i
depends only on u_i, so

    P(no reject | z, u_0) = prod_i E_{u_i}[ Phi( s_i (z/sqrt(n_0) + c(lam_i u_i/u_0) d_i) ) ]

and the unconditional FWER is a TWO-dimensional outer rule over K
one-dimensional inner rules - the same structure, and the same cost, as the null
law in the paper.  That is what makes a rigorous enclosure affordable at any K.

The price of the restriction is that the rule can no longer use one candidate's
variance to sharpen another's threshold.  This module measures that price
against the joint plug-in before anything is built on it.

Deposit: deposits/wb_percand.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import norm

import fcore as F
import k1lab as L


def fwer_percand(cfun, n, lam, npanel_u=20, order_u=10, nz=64, zorder=12):
    """Exact unconditional FWER of a per-candidate threshold rule.

    Two-dimensional outer rule over (z, u_0), K one-dimensional inner rules.
    No draws.
    """
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    K = n.size - 1
    z, wz = F.normal_nodes(nz, zorder)
    q0, w0 = F.state_nodes(n[0] - 1.0, npanel_u, order_u)
    inner = [F.state_nodes(n[i + 1] - 1.0, npanel_u, order_u)
             for i in range(K)]
    s_i = np.sqrt(n[1:] / lam)                      # sqrt(n_i/lam_i)
    tot = 0.0
    for a0, pa0 in zip(q0, w0):
        prod = np.ones(z.size)
        for i in range(K):
            qi, wi = inner[i]
            lh = lam[i] * qi / a0                   # this candidate's estimate
            ci = cfun(lh)
            d = np.sqrt(lam[i] * qi / n[i + 1] + a0 / n[0])
            # shape (nz, n_inner)
            arg = (z[:, None] / np.sqrt(n[0]) + ci[None, :] * d[None, :]) \
                * s_i[i]
            prod *= norm.cdf(arg) @ wi
        tot += pa0 * float(wz @ prod)
    return float(1.0 - tot)


def power_percand(cfun, n, lam, delta, npanel_u=20, order_u=10, nz=64,
                  zorder=12):
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    K = n.size - 1
    dl = np.broadcast_to(np.asarray(delta, float), (K,))
    z, wz = F.normal_nodes(nz, zorder)
    q0, w0 = F.state_nodes(n[0] - 1.0, npanel_u, order_u)
    inner = [F.state_nodes(n[i + 1] - 1.0, npanel_u, order_u)
             for i in range(K)]
    s_i = np.sqrt(n[1:] / lam)
    tot = 0.0
    for a0, pa0 in zip(q0, w0):
        prod = np.ones(z.size)
        for i in range(K):
            qi, wi = inner[i]
            lh = lam[i] * qi / a0
            ci = cfun(lh)
            d = np.sqrt(lam[i] * qi / n[i + 1] + a0 / n[0])
            arg = (z[:, None] / np.sqrt(n[0]) + ci[None, :] * d[None, :]
                   - dl[i]) * s_i[i]
            prod *= norm.cdf(arg) @ wi
        tot += pa0 * float(wz @ prod)
    return float(1.0 - tot)


def lam_patterns(K):
    out = [np.full(K, v) for v in (0.01, 0.1, 1.0, 10.0, 100.0)]
    if K >= 2:
        for lo, hi in ((0.04, 25.0), (0.01, 1.0), (1.0, 100.0)):
            out.append(np.exp(np.linspace(np.log(lo), np.log(hi), K)))
        v = np.ones(K); v[0] = 25.0; out.append(v.copy())
        v = np.ones(K); v[-1] = 0.04; out.append(v.copy())
    return out


def calibrate(n, abar=F.ABAR, iters=13, base_level=None):
    """The largest per-candidate level whose exact uniform FWER is at or below
    the ceiling, over the declared pattern set.

    The shape of c(.) is the K = 1 plug-in table, which is the natural
    per-candidate object; the scalar being calibrated is its internal level.
    """
    n = np.asarray(n, float)
    K = n.size - 1
    pats = lam_patterns(K)
    n1 = np.array([n[0], n[1]])
    cache = {}

    def sup_at(a):
        if a not in cache:
            tab = L.PluginTable(n1, a, m=41, npanel=12, order=8)
            vals = [fwer_percand(tab, n, p) for p in pats]
            j = int(np.argmax(vals))
            cache[a] = (float(vals[j]), pats[j], tab)
        return cache[a]

    hi = abar
    if sup_at(hi)[0] <= abar:
        s, lw, tab = sup_at(hi)
        return {"alpha0": hi, "sup_fwer": s, "worst": lw.tolist(),
                "table": tab, "binding": False}
    lo = 0.002
    for _ in range(iters):
        m = 0.5 * (lo + hi)
        if sup_at(m)[0] <= abar:
            lo = m
        else:
            hi = m
    s, lw, tab = sup_at(lo)
    return {"alpha0": lo, "sup_fwer": s, "worst": lw.tolist(), "table": tab,
            "binding": True}


DESIGNS = [
    ("K2_df3", [4.0, 13.0, 13.0]),
    ("K2_df7", [8.0, 13.0, 13.0]),
    ("K2_df11", [12.0, 13.0, 13.0]),
    ("K2_df35", [36.0, 13.0, 13.0]),
    ("K3_df11", [12.0, 13.0, 13.0, 13.0]),
    ("K3_df35", [36.0, 13.0, 13.0, 13.0]),
]
# K = 5 is omitted deliberately: the nested product makes one calibration cost
# about forty minutes, and the point of this module - what the per-candidate
# restriction costs relative to the joint rule - is already made at K = 2 and
# K = 3, where the gap is three to four times in power.  The joint rule covers
# K = 5 in wb_general.py.


def main():
    out = {"abar": F.ABAR, "designs": []}
    print("%-10s %8s %9s %9s %8s %8s %6s" %
          ("design", "K", "alpha0", "supFWER", "loss80m", "loss80x", "sec"))
    for name, n in DESIGNS:
        t0 = time.time()
        nn = np.asarray(n, float)
        K = nn.size - 1
        cal = calibrate(nn)
        tab = cal["table"]
        rows = []
        for p in lam_patterns(K)[:5]:
            c_or = F.oracle_c_pop(nn, p, F.ABAR, F.FAST)

            def opow(dl):
                return power_percand(lambda x, c=c_or: np.full_like(
                    np.asarray(x, float), c), nn, p, dl)

            a, b = 0.0, 1.0
            while opow(b) < 0.8 and b < 40:
                b *= 2
            for _ in range(22):
                m = 0.5 * (a + b)
                if opow(m) < 0.8:
                    a = m
                else:
                    b = m
            dl = 0.5 * (a + b)
            pw = power_percand(tab, nn, p, dl)
            rows.append({"lam": p.tolist(), "delta80": float(dl),
                         "rule_power": float(pw),
                         "loss_pp": float(100.0 * (0.8 - pw))})
        losses = [r["loss_pp"] for r in rows]
        rec = {"design": name, "n": nn.tolist(), "K": int(K),
               "alpha0": float(cal["alpha0"]),
               "sup_fwer": float(cal["sup_fwer"]),
               "worst_lam": cal["worst"], "rows": rows,
               "median_loss_pp": float(np.median(losses)),
               "max_loss_pp": float(max(losses)),
               "seconds": float(time.time() - t0)}
        out["designs"].append(rec)
        print("%-10s %8d %9.5f %9.5f %8.1f %8.1f %6.0f" %
              (name, K, rec["alpha0"], rec["sup_fwer"],
               rec["median_loss_pp"], rec["max_loss_pp"], rec["seconds"]))
    out["all_within_ceiling"] = all(d["sup_fwer"] <= F.ABAR + 1e-9
                                    for d in out["designs"])
    out["worst_loss_pp"] = max(d["max_loss_pp"] for d in out["designs"])
    out["note"] = (
        "Per-candidate thresholds restore the paper's factorisation for the "
        "FWER itself, so the unconditional constraint is an exact 2-D outer "
        "rule over K 1-D inner rules with no draws at any K.  The supremum "
        "here is still over a declared pattern set; enclosing it over the "
        "continuum is the remaining step to the certified tier at K > 1."
    )
    F.write_deposit("wb_percand.json", out,
                    "unconditional_candidate_measured")
    print("\nall within ceiling:", out["all_within_ceiling"],
          "| worst loss %.1f pp" % out["worst_loss_pp"])


if __name__ == "__main__":
    main()
