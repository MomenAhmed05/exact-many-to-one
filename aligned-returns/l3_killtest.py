"""L-3 kill test: does the EQUAL-WINDOW clean reduction survive UNEQUAL windows?

Built from the instrument of paper 1's `04-efficient-rule/w2_unequal.py`
(read-only; imported conceptually, never edited). That instrument already carries a
shared factor on the candidate series. This driver changes one thing: instead of
lengthening a benchmark PREHISTORY, it lets the two series be estimated on two
windows of DIFFERENT lengths that share `k` dates, which is L-3's configuration.

The falsifier, stated before running (after a review of the overlapping-samples
literature): the equal-window clean reduction is
the claim that one critical value, calibrated on the equal-window design, controls
size on ANY design. So:

  1. calibrate a critical value on the EQUAL, FULLY OVERLAPPING design;
  2. transfer it unchanged to unequal designs;
  3. measure realised size (rejection frequency under the null).

A departure of more than 3 Monte Carlo standard errors from the nominal level is
the predeclared failure criterion. Controls: c=0 (no factor coupling) and k=0
(disjoint windows). A zero MEAN signed statistic is deliberately NOT the criterion:
the null is symmetric, so wrong tails pass a mean check unnoticed.

Seeds are recorded. All output is CLAIM lines; no number is typed by hand.
"""

from __future__ import annotations

import json
import math
import pathlib

import numpy as np

ALPHA = 0.05
DRAWS = 200_000
BATCH = 50_000
LAM = 1.0
N_CAL = 24          # equal-window calibration design
# unequal designs: (n0, n1, k) with k = shared dates
DESIGNS = [
    ("equal-full-overlap", N_CAL, N_CAL, N_CAL),
    ("equal-half-overlap", N_CAL, N_CAL, N_CAL // 2),
    ("unequal-n1-longer",  12,    24,    12),
    ("unequal-n0-longer",  24,    12,    12),
    ("unequal-disjoint",   24,    24,     0),
    ("unequal-short",       4,    24,     4),
]
Z = 3.0              # predeclared failure threshold, in Monte Carlo standard errors


def _simulate(m, n0, n1, k, b0, b1, lam, seed):
    """Draw dates, share the factor on the overlap, return the statistic."""
    rng = np.random.default_rng(int(seed))
    # W0 = dates [0, n0); W1 = last k of W0 plus n1-k fresh dates
    span = n0 + max(0, n1 - k)
    f = rng.normal(size=(m, span))
    e0 = rng.normal(size=(m, span))
    e1 = rng.normal(size=(m, span))
    r0 = b0 * f + e0
    r1 = b1 * f + math.sqrt(lam) * e1
    w1 = slice(n0 - k, span)
    w0 = slice(0, n0)
    a, b = r0[:, w0], r1[:, w1]
    m0, s0 = a.mean(1), a.var(1, ddof=1)
    m1, s1 = b.mean(1), b.var(1, ddof=1)
    return (m1 - m0) / np.sqrt(s1 / n1 + s0 / n0)


def tail_level(n0, n1, k, b0, b1, crit, seed, draws=DRAWS):
    """Rejection frequency under the null, with its Monte Carlo standard error."""
    hits = done = 0
    while done < draws:
        m = min(BATCH, draws - done)
        t = _simulate(m, n0, n1, k, b0, b1, LAM, seed + done)
        hits += int((t > crit).sum())
        done += m
    p = hits / draws
    return p, math.sqrt(p * (1 - p) / draws)


def calibrate(n, b0, b1, seed, alpha=ALPHA):
    """The 1-alpha quantile of the statistic under the null, on the equal design."""
    t = _simulate(DRAWS, n, n, n, b0, b1, LAM, seed)
    return float(np.quantile(t, 1.0 - alpha))


def main():
    out = {"alpha": ALPHA, "draws": DRAWS, "lam": LAM, "n_cal": N_CAL,
           "z_threshold": Z, "rows": []}

    # --- calibration on the equal, fully overlapping design, at the coupled loadings
    B0, B1 = -1.5, 1.5                      # c = b0*b1 < 0, the sheet of interest
    crit_coupled = calibrate(N_CAL, B0, B1, seed=10_001)
    crit_uncoupled = calibrate(N_CAL, 0.0, 0.0, seed=10_002)

    for name, n0, n1, k in DESIGNS:
        for cname, b0, b1, crit in (
            ("coupled(c<0)", B0, B1, crit_coupled),
            ("uncoupled(c=0)", 0.0, 0.0, crit_uncoupled),
        ):
            p, se = tail_level(n0, n1, k, b0, b1, crit,
                               seed=abs(hash((name, cname))) % 100_000)
            z = (p - ALPHA) / se if se > 0 else float("nan")
            out["rows"].append({
                "design": name, "arm": cname, "n0": n0, "n1": n1, "k": k,
                "size": p, "mcse": se, "z_vs_alpha": z,
                "breach": bool(abs(z) > Z),
            })

    coupled = [r for r in out["rows"] if r["arm"].startswith("coupled")]
    base = next(r for r in coupled if r["design"] == "equal-full-overlap")

    # the reduction's claim: transferring the equal-window critical value changes nothing.
    # so compare each design's size against the EQUAL-WINDOW design's realised size.
    for r in out["rows"]:
        r["z_vs_equal_design"] = ((r["size"] - base["size"]) / r["mcse"]
                                  if r["mcse"] > 0 else float("nan"))
        r["survives_reduction"] = bool(abs(r["z_vs_equal_design"]) <= Z)

    out["critical_value_coupled"] = crit_coupled
    out["critical_value_uncoupled"] = crit_uncoupled
    out["equal_design_realised_size"] = base["size"]
    out["n_breaches_vs_equal"] = int(sum(1 for r in out["rows"]
                                          if not r["survives_reduction"]))
    out["n_breaches_vs_alpha"] = int(sum(1 for r in out["rows"] if r["breach"]))
    out["worst_coupled_z_vs_equal"] = float(max(
        abs(r["z_vs_equal_design"]) for r in coupled))
    out["verdict_reduction_survives"] = bool(out["n_breaches_vs_equal"] == 0)

    p = pathlib.Path("deposits/l3_killtest.json")
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    return out


if __name__ == "__main__":
    o = main()
    print(f"alpha={ALPHA} draws={DRAWS}  "
          f"equal-design realised size={o['equal_design_realised_size']:.5f}")
    hdr = (f"{'design':22} {'arm':14} {'n0':>3} {'n1':>3} {'k':>3}  "
           f"{'size':>8} {'mcse':>8} {'z(vs a)':>7} {'z(vs eq)':>8}  survives")
    print(hdr)
    for r in o["rows"]:
        surv = "yes" if r["survives_reduction"] else "NO"
        print(f"{r['design']:22} {r['arm']:14} {r['n0']:3d} {r['n1']:3d} "
              f"{r['k']:3d}  {r['size']:8.5f} {r['mcse']:8.5f} "
              f"{r['z_vs_alpha']:+7.1f} {r['z_vs_equal_design']:+8.1f}  {surv}")
    print(f"\nreduction survives: {o['verdict_reduction_survives']}   "
          f"(breaches vs equal design: {o['n_breaches_vs_equal']})")
