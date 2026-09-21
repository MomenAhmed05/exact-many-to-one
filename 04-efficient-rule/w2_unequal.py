"""W2, Phase 2B: unequal histories, and what the incumbent's prehistory buys.

A promotion gate has an asymmetry the clinical-trial model does not: the
incumbent has a long track record and the candidate a short one.  Aligning them
means discarding the part of the benchmark's history the candidate does not
cover, and the obvious objection is that this throws away information.

The question is whether it throws away TEST-RELEVANT information, which is a
different question and has a different answer.  Three predeclared designs:

  PAIRED-WINDOW   use only the overlap, for both groups.  Exact, and the
                  aligned-factor law of w2_factor.py applies directly.
  TWO-WINDOW      use the overlap for both, plus the benchmark's prehistory to
                  sharpen its VARIANCE estimate only.  Requires that the
                  benchmark's variance is constant across the two windows; it
                  does NOT require the mean to be, because the mean is only ever
                  taken on the overlap.
  FULL-HISTORY    use the benchmark's whole record for both its mean and its
                  variance.  This is what a desk does by default, and it is the
                  one with a stationarity assumption on the MEAN.

The measurement that matters: the prehistory enters the statistic only through
s_0^2, whose sampling error is second order in the numerator's noise once the
overlap is long.  So the power the prehistory buys should saturate, and the
point at which it saturates is the operational answer to "is it worth assuming
stationarity to use it".

A stationarity break is then injected deliberately -- the benchmark's mean
shifts before the overlap begins -- to price what FULL-HISTORY costs when the
assumption it needs is false.

Deposit: deposits/w2_unequal.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L

ABAR = F.ABAR
LAM_DOMAIN = (1e-3, 1e3)
DRAWS = 300000


def run(design, n_over, n_pre, lam, beta, shift, tab_of, D=DRAWS, seed=0,
        delta=0.0, batch=75000):
    """Realised FWER (delta = 0) or power at one design.

    `beta` is the common latent-factor loading contrast; the benchmark carries
    loading 0 and the candidate `beta`, since only the contrast is identified.
    `shift` is a mean break applied to the benchmark's PREHISTORY only.
    """
    rng = np.random.default_rng(int(seed))
    hits = done = 0
    while done < D:
        m = min(batch, D - done)
        f = rng.normal(size=(m, n_over))
        b_over = rng.normal(size=(m, n_over))
        b_pre = rng.normal(size=(m, n_pre)) + shift if n_pre else None
        c_over = (np.sqrt(lam) * rng.normal(size=(m, n_over))
                  + beta * f + delta)

        if design == "paired":
            bm, bv, nb = b_over.mean(1), b_over.var(1, ddof=1), n_over
        elif design == "two_window":
            allb = np.concatenate([b_pre, b_over], 1)
            # mean from the overlap, variance from everything
            bm, bv, nb = b_over.mean(1), allb.var(1, ddof=1), n_over
        else:                                       # full_history
            allb = np.concatenate([b_pre, b_over], 1)
            bm, bv, nb = allb.mean(1), allb.var(1, ddof=1), n_pre + n_over

        cm, cv = c_over.mean(1), c_over.var(1, ddof=1)
        t = (cm - bm) / np.sqrt(cv / n_over + bv / nb)
        tab = tab_of(nb)
        hits += int((t > tab(cv / bv)).sum())
        done += m
    p = hits / D
    return p, float(np.sqrt(p * (1 - p) / D))


def main():
    t0 = time.time()
    n_over = 24
    pres = (0, 12, 36, 120, 360)
    cache = {}

    def tab_of(nb):
        key = int(nb)
        if key not in cache:
            n = np.array([float(key), float(n_over)])
            cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                                    npanel=12, order=8, iters=10)
            cache[key] = cal["table"]
        return cache[key]

    # the alternative: sized so the paired design lands near 50% power, which
    # is where differences between designs are most visible
    delta = 0.62

    rows = []
    print("%-13s %5s %6s %9s %9s %8s" %
          ("design", "pre", "shift", "FWER", "power", "z(FWER)"))
    for design in ("paired", "two_window", "full_history"):
        for n_pre in pres:
            if design == "paired" and n_pre != 0:
                continue
            if design != "paired" and n_pre == 0:
                continue
            for shift in (0.0, 0.5, -0.5):
                if design == "paired" and shift != 0.0:
                    continue
                p, se = run(design, n_over, n_pre, 1.0, 0.0, shift, tab_of,
                            seed=abs(int(7 * n_pre + 100 * shift)) + 1)
                pw, _ = run(design, n_over, n_pre, 1.0, 0.0, shift, tab_of,
                            seed=abs(int(7 * n_pre + 100 * shift)) + 6, delta=delta)
                z = (p - ABAR) / se
                rows.append({"design": design, "n_pre": int(n_pre),
                             "shift": float(shift), "fwer": p, "mcse": se,
                             "z": float(z), "power": float(pw),
                             "breach": bool(z > 3.0)})
                print("%-13s %5d %6.1f %9.5f %9.4f %+8.1f %s"
                      % (design, n_pre, shift, p, pw, z,
                         "BREACH" if z > 3.0 else ""))

    base = next(r for r in rows if r["design"] == "paired")
    clean = [r for r in rows if r["shift"] == 0.0]
    broke = [r for r in rows if r["shift"] != 0.0]
    fh = [r for r in clean if r["design"] == "full_history"]
    tw = [r for r in clean if r["design"] == "two_window"]

    # where does the prehistory stop buying power?
    gains = sorted([(r["n_pre"], 100 * (r["power"] - base["power"]))
                    for r in fh])
    sat = "it does not saturate by %d" % max(pres)
    for i in range(1, len(gains)):
        if abs(gains[i][1] - gains[i - 1][1]) < 0.5:
            sat = str(gains[i - 1][0])
            break

    out = {
        "abar": ABAR, "n_overlap": n_over, "prehistory_lengths": list(pres),
        "delta": delta, "rows": rows,
        "paired_fwer": base["fwer"], "paired_power": base["power"],
        "two_window_max_power_gain_pp":
            float(max(100 * (r["power"] - base["power"]) for r in tw)),
        "full_history_max_power_gain_pp":
            float(max(100 * (r["power"] - base["power"]) for r in fh)),
        "full_history_power_gain_by_pre":
            dict((str(k), float(v)) for k, v in gains),
        "saturates_at_prehistory": str(sat),
        "two_window_worst_fwer": float(max(r["fwer"] for r in tw)),
        "full_history_worst_fwer": float(max(r["fwer"] for r in fh)),
        "broken_stationarity_worst_fwer": float(max(r["fwer"] for r in broke)),
        "broken_stationarity_worst_design":
            max(broke, key=lambda r: r["fwer"])["design"],
        "broken_ratio_to_abar":
            float(max(r["fwer"] for r in broke) / ABAR),
        "n_breaches": int(sum(1 for r in rows if r["breach"])),
        "note": ("The two-window design uses the prehistory for the benchmark's "
                 "VARIANCE only and so needs constant variance but not a "
                 "constant mean; full-history uses it for the mean too and "
                 "needs both.  The shift rows break exactly that assumption "
                 "and price it in BOTH directions: a prehistory ABOVE the "
                 "later mean makes the test powerless, one BELOW it makes the "
                 "test liberal, and only the second is dangerous."),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w2_unequal.json", out, "measured")
    print("\nprehistory buys at most %.2f pp (two-window) / %.2f pp "
          "(full history); saturates at %s"
          % (out["two_window_max_power_gain_pp"],
             out["full_history_max_power_gain_pp"],
             out["saturates_at_prehistory"]))
    print("with the mean break: worst FWER %.5f (%.2fx abar) on '%s'"
          % (out["broken_stationarity_worst_fwer"],
             out["broken_ratio_to_abar"],
             out["broken_stationarity_worst_design"]))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
