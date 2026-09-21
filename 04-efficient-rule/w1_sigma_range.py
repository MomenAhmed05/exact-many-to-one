"""Has the supremum over sigma_0^2 actually converged?

THE DEFECT THIS EXISTS TO MEASURE.  The Proposition's size constraint runs over
every sigma_0^2 in (0, infinity).  The probes impose it on a finite sweep, and
`w1_hypotheses_directed.py` swept exp(-6) to exp(6).  Two rules that differ only
in their guess v0 are the SAME rule at relocated sigma_0^2, so they must score
identically once the sweep is wide enough to contain both.  On the e^6 sweep
they did not: deproj@0.3 and deproj@1.0 came out 0.39 points apart at
gamma = 0.25, which is the same order as the gap the Proposition concerns.

So the discrepancy between two guesses is a DIRECT read on whether the sweep is
wide enough, and it costs nothing extra to measure: it is the difference between
two numbers the probe already computes.  This module widens the sweep and
watches that difference.  If it shrinks toward zero the supremum has converged
and the probe's numbers can be quoted; if it plateaus the sweep is not the
binding issue and the instrument needs rethinking rather than widening.

WHICH WAY THE BIAS RUNS.  Truncating a supremum UNDERSTATES each rule's worst
FWER, so the calibrating shift comes out too small and the rule is left too
liberal, which lowers its measured regret.  A truncated sweep therefore
FLATTERS the non-invariant rules.  They lost anyway on the narrow sweep, so the
qualitative finding is conservative; what a narrow sweep forbids is quoting the
regrets precisely.

Deposit: deposits/w1_sigma_range.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L
from w1_hypotheses_directed import (ABAR, ALT_LAMS, LAM_DOMAIN, ORACLE,
                                    rule_fwer, rule_power, tilt)

RANGES = (6.0, 10.0, 14.0, 18.0)
GAMMAS = (0.25, 1.0)
GUESSES = ("deproj@1.0", "deproj@0.3")


def main():
    t0 = time.time()
    n = np.array([36.0, 13.0])
    cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                            npanel=12, order=8, iters=10)
    tab = cal["table"]
    lam_grid = F.log_grid(*LAM_DOMAIN, 41)
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], ORACLE)) for l in ALT_LAMS]

    rows = []
    print("%6s %6s %8s %12s %12s %10s"
          % ("range", "gamma", "n_v", "deproj@1.0", "deproj@0.3", "gap pp"))
    for R in RANGES:
        m = int(4 * R + 1)                     # hold the node density fixed
        v_grid = np.exp(np.linspace(-R, R, m))

        def worst_fwer(cf):
            return max(rule_fwer(cf, n, l, v)
                       for l in lam_grid for v in v_grid)

        def worst_regret(cf):
            return max(ORACLE - rule_power(cf, n, l, v, d)
                       for l, d in alts for v in v_grid)

        for g in GAMMAS:
            got = {}
            for kind in GUESSES:
                r = tilt(tab, g, kind)
                lo, hi = -2.0, 40.0
                for _ in range(26):
                    mid = 0.5 * (lo + hi)
                    if worst_fwer(lambda lh, s0, t=mid, rr=r:
                                  rr(lh, s0) + t) > ABAR:
                        lo = mid
                    else:
                        hi = mid
                cr = (lambda lh, s0, t=hi, rr=r: rr(lh, s0) + t)
                got[kind] = (100.0 * worst_regret(cr), float(hi),
                             float(worst_fwer(cr)))
            gap = abs(got[GUESSES[0]][0] - got[GUESSES[1]][0])
            rows.append({
                "range_log_sigma2": float(R), "n_v": m, "gamma": float(g),
                "regret_guess_1_pp": got["deproj@1.0"][0],
                "regret_guess_0p3_pp": got["deproj@0.3"][0],
                "shift_guess_1": got["deproj@1.0"][1],
                "shift_guess_0p3": got["deproj@0.3"][1],
                "fwer_guess_1": got["deproj@1.0"][2],
                "fwer_guess_0p3": got["deproj@0.3"][2],
                "gap_pp": float(gap)})
            print("%6.0f %6.2f %8d %12.3f %12.3f %10.3f"
                  % (R, g, m, got["deproj@1.0"][0], got["deproj@0.3"][0], gap))

    # the gap must shrink as the sweep widens, at every gamma separately
    conv = {}
    for g in GAMMAS:
        sub = [r for r in rows if r["gamma"] == g]
        gaps = [r["gap_pp"] for r in sub]
        conv["gamma_%.2f" % g] = {
            "gaps_pp": gaps,
            "first": gaps[0], "last": gaps[-1],
            "shrank": bool(gaps[-1] < gaps[0]),
            "monotone": bool(all(gaps[i] <= gaps[i - 1] + 1e-9
                                 for i in range(1, len(gaps)))),
            "converged": bool(gaps[-1] < 0.02)}
    allconv = all(c["converged"] for c in conv.values())
    out = {
        "abar": ABAR, "oracle_power": ORACLE, "n": n.tolist(),
        "ranges": list(RANGES), "gammas": list(GAMMAS),
        "rows": rows, "convergence": conv,
        "sup_has_converged": allconv,
        "widest_range": float(RANGES[-1]),
        "worst_gap_at_widest_pp": float(max(
            r["gap_pp"] for r in rows
            if r["range_log_sigma2"] == RANGES[-1])),
        "worst_gap_at_narrowest_pp": float(max(
            r["gap_pp"] for r in rows
            if r["range_log_sigma2"] == RANGES[0])),
        "bias_direction": ("a truncated sweep understates the worst FWER, so "
                           "it undercalibrates and FLATTERS the non-invariant "
                           "rules; they lost anyway, so the qualitative "
                           "finding is conservative"),
        "note": ("Two guesses v0 describe the same rule at relocated "
                 "sigma_0^2, so the gap between their scores is zero when the "
                 "sweep is wide enough to contain both.  That makes the gap a "
                 "direct read on convergence of the supremum, computed from "
                 "numbers the probe already produces."),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_sigma_range.json", out, "frontier_measured")
    for g in GAMMAS:
        c = conv["gamma_%.2f" % g]
        print("gamma %.2f: gap %.3f -> %.3f points as the sweep goes e^%.0f "
              "to e^%.0f; monotone %s, converged %s"
              % (g, c["first"], c["last"], RANGES[0], RANGES[-1],
                 c["monotone"], c["converged"]))
    print("the supremum has converged: %s" % allconv)
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
