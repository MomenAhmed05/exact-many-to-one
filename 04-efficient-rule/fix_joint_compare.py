"""Recompute the joint deposit's comparison field against the RIGHT object.

Gate X11 failed, and it was right to.  `certify_joint.py` filled
`exact_sup_over_patterns` by evaluating a PER-CANDIDATE plug-in table at the
joint rule's level - a different rule.  At the same level the joint constant
exceeds the one-candidate constant, so that proxy rejects more often than the
rule being certified and its FWER sat ABOVE the bound, which looks exactly like
a bound failing to dominate.

The certificates themselves are untouched: `alpha0` and the union bounds come
from the projection and the enclosure and never used this field.  What was wrong
was the number the deposit offered for comparison, and a gate that compares a
bound against the wrong object is worse than no gate, so it is fixed here rather
than relaxed.

The right object is the deployed joint step rule's own unconditional FWER,
evaluated Rao-Blackwellised: draw the variance estimates, look the threshold up
in the deployed table exactly as a desk would, and evaluate the conditional
rejection in closed form.  That carries Monte-Carlo error, which is deposited
with it so the gate can allow for it.
"""

from __future__ import annotations

import time

import numpy as np

import certify_joint as J
import fcore as F
import wb_percand as W

DRAWS = 4000
SEED = 20260920


def exact_fwer_joint(tab, n, lam, D=DRAWS, seed=SEED, quad=None):
    """Unconditional FWER of the DEPLOYED joint step rule at one ratio vector."""
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    quad = quad or F._u_nodes(48, 10)
    u = F.draw_states(n, D, seed)
    lh = lam[None, :] * u[:, 1:] / u[:, :1]
    r = np.empty(D)
    for d in range(D):
        c = tab(lh[d])
        r[d] = F.cond_reject(c, u[d], n, lam, quad)
    return float(r.mean()), float(r.std(ddof=1) / np.sqrt(D))


def main():
    t0 = time.time()
    dep = F.read_deposit("certify_joint.json")
    print("%-10s %12s %12s %10s %8s" %
          ("design", "union bound", "exact FWER", "mcse", "dominates"))
    for rec in dep["designs"]:
        n = np.asarray(rec["n"], float)
        K = rec["K"]
        nbin = rec.get("nbins", J.NBIN)
        rows, index = J.JointStepTable.tabulate(n, nbin)
        tab = J.JointStepTable(n, rec["alpha0_certified_joint"], nbin,
                               rows=rows, index=index)
        best, bse, worst_lam = -1.0, 0.0, None
        for p in W.lam_patterns(K):
            v, se = exact_fwer_joint(tab, n, p)
            if v > best:
                best, bse, worst_lam = v, se, p.tolist()
        rec["exact_sup_over_patterns"] = float(best)
        rec["exact_sup_mcse"] = float(bse)
        rec["exact_sup_worst_lam"] = worst_lam
        rec["exact_sup_object"] = (
            "the DEPLOYED joint step rule's own unconditional FWER, "
            "Rao-Blackwellised over %d draws; the previous value in this field "
            "was a per-candidate plug-in table at the same level, which is a "
            "different rule and rejects more often" % DRAWS)
        dom = rec["union_bound_total"] >= best - 3.0 * bse
        rec["bound_dominates_exact"] = bool(dom)
        print("%-10s %12.6f %12.6f %10.6f %10s" %
              (rec["design"], rec["union_bound_total"], best, bse, dom))
    dep["all_bounds_dominate"] = all(d["bound_dominates_exact"]
                                     for d in dep["designs"])
    dep["comparison_fix"] = (
        "Gate X11 caught the comparison field measuring the wrong rule.  The "
        "certificates are unchanged - they never used it - and the field now "
        "carries the deployed rule's own FWER with its Monte-Carlo error.")
    dep["seconds_recompare"] = float(time.time() - t0)
    F.write_deposit("certify_joint.json", dep, "unconditional_valid_certified")
    print("\nall bounds dominate:", dep["all_bounds_dominate"])


if __name__ == "__main__":
    main()
