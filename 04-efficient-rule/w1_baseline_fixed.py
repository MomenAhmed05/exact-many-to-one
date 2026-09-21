"""The comparison baseline does not move when the sweep widens.

WHY THIS NEEDS ITS OWN DEPOSIT.  `w1_sigma_range.py` shows every non-invariant
rule getting WORSE as the sweep over sigma_0^2 widens.  That trend only means
something if the thing they are being compared against holds still.  If both
arms drifted with the sweep, a rising loss would say nothing.

The invariant rule's threshold is a function of lambda_hat alone, and the law of
lambda_hat does not involve sigma_0^2, so its score should be EXACTLY
independent of the sweep rather than approximately so.  That is a prediction
with no free parameter, and it is worth measuring rather than asserting: it is
the Hunt-Stein mechanism stated as arithmetic, since widening the sweep
punishes every rule that reads the scale and leaves the deployed rule untouched.

Deposit: deposits/w1_baseline_fixed.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L
from w1_hypotheses_directed import (ABAR, ALT_LAMS, LAM_DOMAIN, ORACLE,
                                    rule_fwer, rule_power)

RANGES = ((6.0, 25), (10.0, 41), (14.0, 57), (18.0, 73))


def main():
    t0 = time.time()
    n = np.array([36.0, 13.0])
    cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                            npanel=12, order=8, iters=10)
    tab = cal["table"]
    lam = F.log_grid(*LAM_DOMAIN, 41)
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], ORACLE)) for l in ALT_LAMS]

    def flat(lh, s0):
        return tab(lh)

    rows = []
    print("%6s %6s %16s %16s" % ("range", "n_v", "worst FWER", "worst regret"))
    for R, m in RANGES:
        v = np.exp(np.linspace(-R, R, m))
        f = max(rule_fwer(flat, n, l, x) for l in lam for x in v)
        r = max(ORACLE - rule_power(flat, n, l, x, d) for l, d in alts
                for x in v)
        rows.append({"range_log_sigma2": float(R), "n_v": int(m),
                     "worst_fwer": float(f), "worst_regret_pp": float(100 * r)})
        print("%6.0f %6d %16.10f %16.6f" % (R, m, f, 100 * r))

    fw = [r["worst_fwer"] for r in rows]
    rg = [r["worst_regret_pp"] for r in rows]
    spread_f = float(max(fw) - min(fw))
    spread_r = float(max(rg) - min(rg))
    out = {
        "abar": ABAR, "oracle_power": ORACLE, "n": n.tolist(),
        "rows": rows,
        "worst_fwer": fw[0],
        "worst_regret_pp": rg[0],
        "fwer_spread_across_sweeps": spread_f,
        "regret_spread_across_sweeps_pp": spread_r,
        "exactly_fixed": bool(spread_f == 0.0 and spread_r == 0.0),
        "narrowest_range": float(RANGES[0][0]),
        "widest_range": float(RANGES[-1][0]),
        "orders_of_magnitude_spanned": float(
            2 * RANGES[-1][0] / np.log(10.0)),
        "note": ("The deployed rule reads lambda_hat only, and the law of "
                 "lambda_hat does not involve sigma_0^2, so widening the sweep "
                 "cannot move it.  Measured rather than asserted, because it "
                 "is what makes the rising losses in w1_sigma_range a "
                 "comparison and not a common drift."),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_baseline_fixed.json", out, "exact_in_model")
    print("spread across sweeps: FWER %.3e, regret %.3e points"
          % (spread_f, spread_r))
    print("exactly fixed across %.0f orders of magnitude in sigma_0^2: %s"
          % (out["orders_of_magnitude_spanned"], out["exactly_fixed"]))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
