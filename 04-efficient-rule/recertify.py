"""Re-certify the deposited K = 1 rules after the tail-mass correction.

`wb_certify.enclose_fwer_box` previously summed over a probability-scale
partition of [1e-9, 1 - 1e-9] per axis and reported that sum as an upper bound,
which silently dropped 4e-9 of mass.  The integrand is a probability, so the
dropped mass can contribute at most its own size; it is now added to the upper
bound.  The tightest margin in the deposit was 1.6e-8, so the correction cannot
be waved away as beneath notice, and every certificate has to be re-checked
against the corrected bound.

This re-runs only the final certification of each deposited rule - not the
bisection that found its level - and rewrites the deposit's certificate fields.
A rule that no longer certifies is reported as such and its level is lowered.

Deposit: rewrites deposits/wb_certified_rule.json in place.
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import wb_certify as C
import wb_certified_rule as R


def main():
    dep = F.read_deposit("wb_certified_rule.json")
    abar = dep["abar"]
    print("%-16s %12s %12s %10s %8s" %
          ("design", "old supUpper", "new supUpper", "certified", "sec"))
    changed = 0
    for d in dep["designs"]:
        t0 = time.time()
        n = np.asarray(d["n"], float)
        rule = C.StepRule(d["c"], d["edges"])
        cert = C.certify(rule, n, R.LAM_LO, R.LAM_HI, abar,
                         m=dep["cells_per_axis"], max_boxes=640)
        old = d["sup_upper_bound"]
        if not cert["certified"]:
            # lower the level until the corrected bound certifies again
            a0 = d["alpha0"]
            for _ in range(8):
                a0 *= 0.995
                rule = R.rule_from_level(n, a0, np.asarray(d["edges"], float))
                cert = C.certify(rule, n, R.LAM_LO, R.LAM_HI, abar,
                                 m=dep["cells_per_axis"], max_boxes=640)
                if cert["certified"]:
                    break
            d["alpha0"] = float(a0)
            d["c"] = [float(x) for x in rule.c]
            d["rows"] = R.price(rule, n)
            l80 = [r["loss_pp"] for r in d["rows"] if r["target"] == 0.8]
            d["median_loss_pp_at_80"] = float(np.median(l80))
            d["max_loss_pp_at_80"] = float(max(l80))
            changed += 1
        d["certificate"] = {k: v for k, v in cert.items() if k != "undecided"}
        d["certified"] = bool(cert["certified"])
        d["sup_upper_bound"] = cert["sup_upper_bound"]
        print("%-16s %12.9f %12.9f %10s %8.0f" %
              (d["design"], old, cert["sup_upper_bound"] or float("nan"),
               cert["certified"], time.time() - t0))
    dep["all_certified"] = all(x["certified"] for x in dep["designs"])
    dep["worst_loss_pp_at_80"] = max(x["max_loss_pp_at_80"]
                                     for x in dep["designs"])
    dep["tail_mass_correction"] = {
        "omitted_probability_mass": 4.0 * C.PROB_LO,
        "note": ("the probability-scale partition leaves this much mass "
                 "outside its cells; it is added to every upper bound because "
                 "the integrand is a probability and can be 1 there"),
        "designs_relevelled": int(changed),
    }
    F.write_deposit("wb_certified_rule.json", dep,
                    "unconditional_valid_certified")
    print("\nall certified after the correction:", dep["all_certified"],
          "| designs re-levelled:", changed)


if __name__ == "__main__":
    main()
