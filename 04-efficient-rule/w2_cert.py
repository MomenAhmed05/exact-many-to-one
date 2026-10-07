"""W2, Phase 2A continued: can the aligned-factor law carry a CERTIFICATE?

`w2_factor.py` gives the exact null law when the loadings and the variance
ratios are known.  That is the oracle tier.  A desk does not know them, and the
whole point of the preceding programme was that an exact law used as a plug-in
is not a guarantee.  So the question here is whether the construction that
worked for the variance ratio -- declare a domain, reduce the internal level,
take the supremum over the domain -- survives a higher-dimensional nuisance.

WHERE THE LIBERALITY COMES FROM, which the first version of this module could
not see.  Write the loadings as beta_0 for the benchmark and beta_1 for the
candidate.  Then

    Var(Xbar_1 - Xbar_0)  contains  (beta_1 - beta_0)^2 / n
    E[s_1^2/n + s_0^2/n]  contains  (beta_1^2 + beta_0^2) / n .

The numerator carries the SQUARED CONTRAST and the denominator carries the SUM
OF SQUARES.  Those are equal only when one loading is zero.  With loadings of
the same sign the contrast is smaller and the test is conservative; with
loadings of OPPOSITE SIGN the contrast is larger and the test is liberal.

The first version of this module set beta_0 = 0 on the argument that only the
contrast is identified.  That is false, and it placed the entire scan on the
exact boundary between the two regimes, which is why it found nothing.  A
long-short candidate against a directional benchmark is precisely the
opposite-sign case, so this is a configuration a desk is likely to be in.

THE NUISANCE IS THEREFORE THREE-DIMENSIONAL per candidate: (lam, beta_0,
beta_1).  Only lam is estimable from what the rule sees.  The loadings can only
be declared, and the guarantee is uniform over the declared box.

Deposit: deposits/w2_cert.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
import k1lab as L
from w2_factor import GRID, AlignedFactor

# Resolution ladder for the box scan.  MID was used for the first version of
# this scan and is kept for the record: at the small-ratio corners of the box it
# is not converged (at lam = 0.04 it returns negative rates when both loadings
# are 1.5, and overstates the rate at the argsup corner by 7.7e-04).  FINE is the
# resolution of record; FINER is used only to check that FINE has converged.
MID = {"zf": (2, 10), "z0": (2, 10), "vf": (2, 10), "s0": 20, "si": 20}
FINE = {"zf": (4, 14), "z0": (4, 14), "vf": (4, 14), "s0": 40, "si": 40}
FINER = {"zf": (5, 16), "z0": (5, 16), "vf": (5, 16), "s0": 48, "si": 48}
# declared before the run: FINER at the FINE constant must sit within this of
# ABAR at the argsup and the argmin point (about 1e-4 in the constant), and the
# full-panel brute force at the argsup must sit within BRUTE_Z standard errors
LADDER_TOL = 4e-6
BRUTE_D, BRUTE_Z = 4000000, 3.0


def _critical_fine(p):
    """The known-parameter constant at one box point, at FINE resolution."""
    a = AlignedFactor(36, [p[1], p[2]], [p[0]], grid=FINE)
    return a.critical(ABAR, lo=0.5, hi=4.0, tol=1e-7)

ABAR = F.ABAR
LAM_DOMAIN = (1e-3, 1e3)


def aligned_fwer_of_table(tab, n, lam, b0, b1, D=300000, seed=1,
                          batch=100000):
    """Realised FWER of a lambda_hat-indexed table on aligned factor data."""
    rng = np.random.default_rng(int(seed))
    n0, n1 = int(n[0]), int(n[1])
    nn = max(n0, n1)
    hits = done = 0
    while done < D:
        m = min(batch, D - done)
        f = rng.normal(size=(m, nn))
        x0 = rng.normal(size=(m, n0)) + b0 * f[:, :n0]
        x1 = np.sqrt(lam) * rng.normal(size=(m, n1)) + b1 * f[:, :n1]
        s0 = x0.var(1, ddof=1)
        s1 = x1.var(1, ddof=1)
        t = (x1.mean(1) - x0.mean(1)) / np.sqrt(s1 / n1 + s0 / n0)
        hits += int((t > tab(s1 / s0)).sum())
        done += m
    p = hits / D
    return p, float(np.sqrt(p * (1 - p) / D))


def main():
    t0 = time.time()
    n = np.array([36.0, 36.0])

    cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                            npanel=12, order=8, iters=10)
    tab = cal["table"]
    print("iid-certified internal level alpha0 = %.5f" % cal["alpha0"])

    print("\n1. THE PRICE OF IGNORING THE FACTOR (iid rule on aligned data)")
    rows = []
    k = 0
    for lam in (0.25, 1.0, 4.0):
        for b0, b1 in ((0.0, 1.0), (1.0, 1.0), (1.0, 0.2),
                       (1.0, -0.5), (1.0, -1.0), (1.5, -1.5)):
            k += 1
            p, se = aligned_fwer_of_table(tab, n, lam, b0, b1, seed=k + 3)
            z = (p - ABAR) / se
            rows.append({"lam": float(lam), "b0": float(b0), "b1": float(b1),
                         "same_sign": bool(b0 * b1 >= 0),
                         "fwer": p, "mcse": se, "z": float(z),
                         "liberal": bool(z > 3.0)})
            print("   lam %-5.2f b0 %-5.1f b1 %-5.1f  FWER %.5f (%.5f)  "
                  "z %+8.1f  %s"
                  % (lam, b0, b1, p, se, z, "LIBERAL" if z > 3.0 else ""))

    lib = [r for r in rows if r["liberal"]]
    worst = max(rows, key=lambda r: r["fwer"])
    opp = [r for r in rows if not r["same_sign"]]
    same = [r for r in rows if r["same_sign"]]

    print("\n2. EXACT ALIGNED LAW OVER THE DECLARED BOX (grid supremum)")
    bgrid = np.linspace(-1.5, 1.5, 5)
    lgrid = F.log_grid(0.04, 25.0, 7)
    pts = [(float(l), float(a), float(b))
           for l in lgrid for a in bgrid for b in bgrid]
    print("   box: %d points in (lam, beta_0, beta_1)" % len(pts))

    # The rate is increasing in the constant at every point, so the smallest
    # constant holding over the whole grid is the largest of the per-point
    # known-parameter constants; each is root-found at FINE, in parallel.
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=11) as ex:
        oc = np.array(list(ex.map(_critical_fine, pts)))
    assert oc.min() > 0.5 + 1e-6 and oc.max() < 4.0 - 1e-6, \
        "a per-point constant sits on the root-finder's bracket"
    i_max, i_min = int(oc.argmax()), int(oc.argmin())
    arg, arg_min = pts[i_max], pts[i_min]
    c_cert = float(oc[i_max])
    sup_at_cert = AlignedFactor(int(n[0]), [arg[1], arg[2]], [arg[0]],
                                grid=FINE).fwer(c_cert)
    print("   constant holding over the grid c = %.5f  (rate %.6f at lam=%.3f "
          "b0=%.2f b1=%.2f)" % (c_cert, sup_at_cert, arg[0], arg[1], arg[2]))

    print("\n   convergence: FINER at the FINE constant, argsup and argmin")
    ladder = []
    for tag, p, c in (("argsup", arg, c_cert), ("argmin", arg_min,
                                                  float(oc[i_min]))):
        r = AlignedFactor(int(n[0]), [p[1], p[2]], [p[0]], grid=FINER).fwer(c)
        ladder.append({"point": tag, "lam": p[0], "b0": p[1], "b1": p[2],
                       "c_fine": c, "fwer_finer_at_c_fine": float(r),
                       "gap": float(abs(r - ABAR)),
                       "pass": bool(abs(r - ABAR) <= LADDER_TOL)})
        print("   %-6s FINER rate %.8f  gap %.2e  %s"
              % (tag, r, abs(r - ABAR),
                 "PASS" if abs(r - ABAR) <= LADDER_TOL else "FAIL"))
    from w2_factor import brute
    bp, bse = brute(int(n[0]), [arg[1], arg[2]], [arg[0]], c_cert, BRUTE_D,
                    seed=99701)
    bp, bse = float(bp), float(bse)
    bz = float((bp - ABAR) / bse)
    print("   brute force at the argsup: %.6f (se %.6f, z %+.2f) %s"
          % (bp, bse, bz, "PASS" if abs(bz) <= BRUTE_Z else "FAIL"))

    print("\n3. THE PRICE AGAINST A KNOWN-PARAMETER ORACLE")
    print("   oracle constants %.4f to %.4f" % (oc.min(), oc.max()))

    out = {
        "abar": ABAR, "n": n.tolist(),
        "iid_alpha0": float(cal["alpha0"]),
        "mechanism": ("the numerator carries (beta_1 - beta_0)^2 and the "
                      "denominator carries beta_1^2 + beta_0^2; equal only "
                      "when one loading is zero, conservative when the "
                      "loadings share a sign, liberal when they do not"),
        "price_rows": rows,
        "n_liberal_cells": len(lib),
        "n_cells": len(rows),
        "n_opposite_sign_cells": len(opp),
        "n_opposite_sign_liberal": int(sum(1 for r in opp if r["liberal"])),
        "n_same_sign_liberal": int(sum(1 for r in same if r["liberal"])),
        "worst_fwer_iid_rule_on_aligned_data": float(worst["fwer"]),
        "worst_cell": {"lam": worst["lam"], "b0": worst["b0"],
                       "b1": worst["b1"]},
        "worst_ratio_to_abar": float(worst["fwer"] / ABAR),
        "worst_same_sign_fwer": float(max(r["fwer"] for r in same)),
        "grid_gap_mid_vs_reference": float(abs(
            AlignedFactor(int(n[0]), [0.0, 1.0], [1.0], grid=MID).fwer(2.2)
            - AlignedFactor(int(n[0]), [0.0, 1.0], [1.0],
                            grid=GRID).fwer(2.2))),
        "box": {"lam": [float(lgrid[0]), float(lgrid[-1])],
                "beta": [float(bgrid[0]), float(bgrid[-1])],
                "points": len(pts)},
        "certified_constant_over_box": float(c_cert),
        "sup_at_certified_constant": float(sup_at_cert),
        "argsup": {"lam": arg[0], "b0": arg[1], "b1": arg[2]},
        "argmin": {"lam": arg_min[0], "b0": arg_min[1], "b1": arg_min[2]},
        "resolution": {"of_record": FINE, "convergence_check": FINER,
                       "superseded": MID,
                       "ladder_tol": LADDER_TOL, "ladder": ladder,
                       "ladder_pass": bool(all(x["pass"] for x in ladder))},
        "brute_at_argsup": {"draws": BRUTE_D, "seed": 99701, "rate": float(bp),
                            "se": bse, "z": bz,
                            "pass": bool(abs(bz) <= BRUTE_Z)},
        "oracle_constant_min": float(oc.min()),
        "oracle_constant_max": float(oc.max()),
        "certificate_price_vs_worst_oracle_x": float(c_cert / oc.max()),
        "certificate_price_vs_best_oracle_x": float(c_cert / oc.min()),
        "note": ("The supremum is a GRID maximum over a %d-point box, not a "
                 "continuum certificate.  The branch-and-bound construction "
                 "of Section 9.11 applies in principle because the "
                 "branch-fixed integrand is entire in the same way, but the "
                 "enclosure is three-dimensional per candidate rather than "
                 "one and has not been built." % len(pts)),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w2_cert.json", out, "measured")
    print("\n%d of %d cells liberal (%d of %d opposite-sign, %d of %d "
          "same-sign); worst %.5f (%.2fx)"
          % (len(lib), len(rows), out["n_opposite_sign_liberal"], len(opp),
             out["n_same_sign_liberal"], len(same),
             out["worst_fwer_iid_rule_on_aligned_data"],
             out["worst_ratio_to_abar"]))
    print("certificate costs %.2fx the worst oracle constant, %.2fx the best"
          % (out["certificate_price_vs_worst_oracle_x"],
             out["certificate_price_vs_best_oracle_x"]))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
