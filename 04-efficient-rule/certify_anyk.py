"""Certify the per-candidate rule over the continuum at ANY K, one- and
two-sided.  This closes both items that were previously reported as open.

THE ARGUMENT, in one line.  For a per-candidate threshold the rejection events
are one per candidate, so

    1 - prod_i g_i  <=  sum_i (1 - g_i)          pointwise in (z, u_0),

and taking expectations,

    FWER(lambda_1..lambda_K)  <=  sum_i  FWER_1(lambda_i) ,

where `FWER_1` is the ONE-CANDIDATE unconditional error rate for that
candidate's design.  The right-hand side separates completely, so

    sup over the cube  <=  K * sup_lambda FWER_1(lambda) ,

and the inner supremum is exactly the quantity `wb_certify.py` already encloses
tightly over the continuum.  Certifying `sup_lambda FWER_1 <= abar / K` therefore
certifies the K-candidate rule, at any K, with no new enclosure.

Two-sided follows immediately from the exact factor of two established in
`twosided.py`: certify `sup_lambda FWER_1 <= abar / (2K)`.

WHY THIS ROUTE AND NOT A K-DIMENSIONAL ENCLOSURE.  A direct three-dimensional
enclosure over (z, u_0, u_i) was written first and works - it is retained in
`certify_percand.py` and its substitutions are falsified - but it is loose: at
K = 2 it returns 0.0655 where the true supremum is 0.0246, a factor of 2.7, and
tightening it costs cells in three axes at once.  The union step costs far less.
Measured against the exact law at the ratios that actually bind:

    K = 2   slack 5.4e-05 .. 2.6e-03
    K = 3   slack 4.3e-04 .. 6.0e-03
    K = 5   slack 1.4e-03 .. 1.4e-02   (the largest at an all-equal ratio, which
                                        is far below nominal and does not bind)

so the price of the union step is about half a per cent of the level at K = 2 and
about ten per cent at K = 5 - against a factor of 2.7 for the direct enclosure.
The looser construction is kept and reported rather than deleted, because a
reader is entitled to see which route was chosen and why.

Deposit: deposits/certify_anyk.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L
import wb_certify as C
import wb_percand as W

LAM_LO, LAM_HI = 1e-3, 1e3
NBIN = 60
CELLS = 800


def step_edges(nbin=NBIN):
    return F.log_grid(LAM_LO * 3.0, LAM_HI / 3.0, nbin - 1)


def rule_at(n1, level, edges):
    tab = L.PluginTable(np.asarray(n1, float), level, m=41, npanel=12, order=8)
    mids = np.concatenate([[edges[0] / 3.0],
                           np.sqrt(edges[:-1] * edges[1:]),
                           [edges[-1] * 3.0]])
    return C.StepRule(tab(mids), edges)


def certify_for(n1, target, iters=10, cells=CELLS):
    """Largest internal level whose ONE-CANDIDATE FWER certifies at `target`."""
    edges = step_edges()
    cache = {}

    def ok(a):
        if a not in cache:
            r = rule_at(n1, a, edges)
            cache[a] = (C.certify(r, np.asarray(n1, float), LAM_LO, LAM_HI,
                                  target, m=cells, max_boxes=640), r)
        return cache[a]

    hi = target
    if ok(hi)[0]["certified"]:
        cert, rule = ok(hi)
        return hi, cert, rule, edges
    lo = 0.0004
    for _ in range(iters):
        m = 0.5 * (lo + hi)
        if ok(m)[0]["certified"]:
            lo = m
        else:
            hi = m
    cert, rule = ok(lo)
    return lo, cert, rule, edges


DESIGNS = [("K2_df35", [36.0, 13.0, 13.0]),
           ("K3_df35", [36.0, 13.0, 13.0, 13.0]),
           ("K5_df35", [36.0, 13.0, 13.0, 13.0, 13.0, 13.0]),
           ("K2_df11", [12.0, 13.0, 13.0]),
           ("K3_df11", [12.0, 13.0, 13.0, 13.0])]


def main():
    t0 = time.time()
    out = {"abar": F.ABAR, "lam_domain": [LAM_LO, LAM_HI], "designs": []}
    print("%-10s %-6s %3s %10s %12s %10s %10s %6s" %
          ("design", "sided", "K", "alpha0", "1-cand bound", "K x bound",
           "certified", "sec"))
    for name, n in DESIGNS:
        n = np.asarray(n, float)
        K = n.size - 1
        n1 = np.array([n[0], n[1]])
        for two in (False, True):
            t1 = time.time()
            div = (2 * K) if two else K
            target = F.ABAR / div
            a0, cert, rule, edges = certify_for(n1, target)
            one = cert["sup_upper_bound"]
            tot = (2.0 if two else 1.0) * K * one if one is not None else None
            # exact check against the real law at the declared patterns
            tab = L.PluginTable(n1, a0, m=41, npanel=12, order=8)
            pats = W.lam_patterns(K)
            if two:
                import twosided as T
                exact = max(T.fwer_2s_percand(tab, n, p) for p in pats)
            else:
                exact = max(W.fwer_percand(tab, n, p) for p in pats)
            rec = {"design": name, "n": n.tolist(), "K": int(K),
                   "two_sided": bool(two),
                   "target_per_candidate": float(target),
                   "alpha0": float(a0),
                   "one_candidate_bound": float(one) if one else None,
                   "union_bound_total": float(tot) if tot else None,
                   "certified": bool(tot is not None and tot <= F.ABAR + 1e-12),
                   "exact_sup_over_patterns": float(exact),
                   "union_slack": float(tot - exact) if tot else None,
                   "boxes": cert["n_boxes_accepted"],
                   "c": [float(x) for x in rule.c],
                   "edges": [float(x) for x in edges],
                   "seconds": float(time.time() - t1)}
            out["designs"].append(rec)
            print("%-10s %-6s %3d %10.5f %12.6f %10.6f %10s %6.0f" %
                  (name, "two" if two else "one", K, a0, one, tot,
                   rec["certified"], rec["seconds"]))
    out["all_certified"] = all(d["certified"] for d in out["designs"])
    out["worst_union_slack"] = max(d["union_slack"] for d in out["designs"]
                                   if d["union_slack"] is not None)
    out["theorem"] = (
        "For a per-candidate threshold rule, sup over the ratio cube of the "
        "unconditional FWER is at most K times the supremum of the "
        "one-candidate unconditional FWER (2K times, two-sided).  The "
        "one-candidate supremum is certified over the continuum by the "
        "enclosure of wb_certify.py.  Therefore the rule is certified at any K "
        "and for either side, with no K-dimensional enclosure.")
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("certify_anyk.json", out, "unconditional_valid_certified")
    print("\nall certified:", out["all_certified"],
          "| worst union slack %.2e" % out["worst_union_slack"])


if __name__ == "__main__":
    main()
