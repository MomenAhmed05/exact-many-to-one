# -*- coding: utf-8 -*-
"""l2_challenge_v5.py — L-2, a challenge to the certified supremum at a finer grid.

Challenges 2B's certified supremum [L, U] from outside its construction: a rigorous
bracket of c*(theta) at every point of the 1053-point grid P1, nested over the
175-point construction grid P0 of paper 1's w2_cert.py, using the 2B rate engine
(formal01_f8a_battery.cstar_bracket). AlignedFactor is not used anywhere.

Verdict (declared before any run):
  (D) dominance  max_{P1} u(theta) <= U
  (A) anchor     max_{P1} u(theta) - max_{P0} l(theta) <= tol = 1e-4
  STANDS iff (D) and (A).
Falsifiers: F-1 t_70 containment at (1,0,0); F-2 corner bracket inside [L, U] and equal
to 2B's F-3 interval to 1e-12; F-4 point Monte Carlo only if the argmax is not the corner.

    python l2_challenge_v5.py [--time-only]
Writes deposits/l2_challenge_v5_claims.txt ALWAYS (P-5), including on failure.
"""
import hashlib
import math
import pathlib
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
OUT = ROOT / "deposits" / "l2_challenge_v5_claims.txt"
TOL = F(1, 10000)
F2_REF = (2.717339803515846, 2.7173398035159595)   # 2B's recorded F-3 gridmax interval
rows = []


def claim(cid, label, value):
    rows.append("CLAIM %s: %s = %s" % (cid, label, value))
    print(rows[-1], flush=True)


def grids():
    """P0 exactly as paper 1's w2_cert.py builds it; P1 = P0 plus midpoints."""
    lam0 = list(np.exp(np.linspace(np.log(0.04), np.log(25.0), 7)))
    b0 = list(np.linspace(-1.5, 1.5, 5))
    lam1 = []
    for i, l in enumerate(lam0):
        lam1.append(l)
        if i + 1 < len(lam0):
            lam1.append(math.sqrt(l * lam0[i + 1]))     # geometric midpoint
    b1 = []
    for i, b in enumerate(b0):
        b1.append(b)
        if i + 1 < len(b0):
            b1.append(0.5 * (b + b0[i + 1]))
    P0 = [(l, a, b) for l in lam0 for a in b0 for b in b0]
    P1 = [(l, a, b) for l in lam1 for a in b1 for b in b1]
    return P0, P1


def bracket(pt):
    import formal01_f8a_battery as B                     # imported in the worker
    th = tuple(F(x) for x in pt)                          # exact value of each float
    lo, hi = B.cstar_bracket(th)
    return pt, F(lo), F(hi)


def main():
    t0 = time.time()
    import formal01_f8a_battery as B
    U, L = F(B.UC), F(B.LC)
    claim("L2V5-CONSTANTS", "U, L, tol (exact)", "%s, %s, %s" % (U, L, TOL))
    P0, P1 = grids()
    s0, s1 = set(P0), set(P1)
    claim("L2V5-GRIDS", "|P0|, |P1|, P0 subset of P1 exactly",
          "%d, %d, %s" % (len(s0), len(s1), s0 <= s1))
    assert len(s0) == 175 and len(s1) == 1053 and s0 <= s1

    # F-1: t_70 containment at (1,0,0)
    import mpmath as mp
    mp.mp.dps = 50
    t70 = mp.findroot(lambda x: mp.betainc(35, mp.mpf(1) / 2, 0, 70 / (70 + x * x),
                                           regularized=True) / 2 - mp.mpf(1) / 40, 1.99)
    _, lo, hi = bracket((1.0, 0.0, 0.0))
    ok1 = lo <= F(str(t70)) <= hi
    claim("L2V5-F1", "t_70 upper-1/40 quantile (mpmath) inside bracket at (1,0,0)",
          "q=%s bracket=[%.15f, %.15f] contained=%s" % (mp.nstr(t70, 20), lo, hi, ok1))

    # F-2: the corner
    corner = (0.04000000000000001, -1.5, 1.5)
    assert corner in s0
    _, clo, chi = bracket(corner)
    ok2 = (L <= clo and chi <= U and abs(float(clo) - F2_REF[0]) <= 1e-12
           and abs(float(chi) - F2_REF[1]) <= 1e-12)
    claim("L2V5-F2", "corner bracket, inside [L,U], equal to 2B F-3 interval to 1e-12",
          "[%.16f, %.16f] pass=%s" % (clo, chi, ok2))

    if "--time-only" in sys.argv:
        t = time.time()
        for pt in P1[:6]:
            bracket(pt)
        claim("L2V5-TIMING", "seconds per bracket, serial, 6 points",
              "%.2f" % ((time.time() - t) / 6))
        return 0

    with ProcessPoolExecutor(max_workers=11) as ex:
        res = list(ex.map(bracket, P1, chunksize=4))
    by = {pt: (lo, hi) for pt, lo, hi in res}
    umax_pt = max(P1, key=lambda p: by[p][1])
    umax = by[umax_pt][1]
    lmax0 = max(by[p][0] for p in P0)
    wide = max(by[p][1] - by[p][0] for p in P1)
    D = umax <= U
    A = umax - lmax0 <= TOL
    claim("L2V5-BRACKETS", "points bracketed, widest bracket", "%d, %.3e" % (len(by), wide))
    claim("L2V5-ARGMAX", "P1 point with the largest upper bracket", str(umax_pt))
    claim("L2V5-D", "max_P1 u <= U", "max_P1 u = %.16f, U = %s, holds = %s" % (umax, U, D))
    claim("L2V5-A", "max_P1 u - max_P0 l <= tol",
          "difference = %.3e, tol = %s, holds = %s" % (umax - lmax0, TOL, A))
    claim("L2V5-MIN", "smallest lower bracket on P1 and its point",
          "%.6f at %s" % (min(by[p][0] for p in P1), min(P1, key=lambda p: by[p][0])))
    claim("L2V5-F4", "independent point check needed (argmax off the corner)",
          str(umax_pt != corner))
    verdict = "STANDS" if (D and A and ok1 and ok2) else "CONTAMINATED"
    claim("L2V5-VERDICT", "L-2 verdict", verdict)
    claim("L2V5-WALL", "seconds", "%.0f" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    try:
        rc = main()
    except Exception:
        claim("L2V5-ERROR", "the run failed", traceback.format_exc().replace("\n", " | "))
        rc = 1
    finally:
        rows.append("CLAIM L2V5-SHA: sha256 of this script = "
                    + hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest())
        OUT.write_text("\n".join(rows) + "\n", encoding="utf-8")
    sys.exit(rc)
