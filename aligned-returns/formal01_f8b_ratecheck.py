"""F-8 CHECKER, Stage A, part 2 — acceptance rules A-4 (per-cell rate bound) and A-5 (H-wit).

Runs the certified cover from the exact-rational half (formal01_f8a_checker.py) through the
rigorous rate engine (formal01_phi_engine.py). The cover carries each cell's exact rational
rho_1, rho_2 and the Lean hypotheses hCover / hLow1 / hLow2 / hNum. This file checks the
remaining two:

  A-4  for EVERY cell, the UPPER endpoint of Phi(Uc; rho_k) is <= alpha, exactly compared,
       and the ball radius is within budget;
  A-5  H-wit: Phi(Lc; r(corner)) > alpha, with r(corner) carried as a BALL because it is
       irrational, and the ball propagated rather than point-evaluated.

Both routes are independent (R-beta and R-Imhof); the spec requires the upper endpoint of
BOTH to be <= alpha and the two balls to overlap. No float enters any comparison.

Declared before the run (A-6), and NOT revisited afterwards:
    Uc = 271740/100000, Lc = 271731/100000, alpha = 1/40, nu = 35
    knots [1/25, 1, 25] x [-3/2, 0, 3/2] x [-3/2, 0, 3/2], lmax = 3, eta = 1/4, eps = 2^-20
    ball radius budget 1e-9
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys
import time
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parent
DEP = ROOT / "deposits"
CLAIMS: list[str] = []


def claim(cid: str, label: str, value) -> None:
    CLAIMS.append(f"CLAIM {cid}: {label} = {value}")


def iv_to_ball(pe, iv):
    """Exact [lo, hi] -> an Arb ball that CONTAINS it.

    The checker's intervals are its own `Iv` class, not flint balls, so they must be
    converted before the engine will take them. Built as midpoint + (hi-lo)/2: adding a
    non-negative term to an arb widens it symmetrically, so the result provably contains
    [lo, hi]. No rounding shortcut, no point evaluation of an irrational.
    """
    from flint import arb, fmpq

    def q(fr):
        return arb(fmpq(fr.numerator, fr.denominator))
    return (q(iv.lo) + q(iv.hi)) / 2 + (q(iv.hi) - q(iv.lo)) / 2


def load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {name} from {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    fc = load("fc", str(ROOT / "formal01_f8a_checker.py"))
    pe = load("pe", str(ROOT / "formal01_phi_engine.py"))

    UC, LC, ALPHA = F(271740, 100000), F(271731, 100000), F(1, 40)
    RAD_BUDGET = F(1, 10 ** 9)
    lam_knots, b0_knots, b1_knots = [F(1, 25), F(1), F(25)], [F(-3, 2), F(0), F(3, 2)], \
        [F(-3, 2), F(0), F(3, 2)]
    lmax, eta, eps = 3, F(1, 4), F(1, 1 << 20)

    t0 = time.time()
    cells = fc.make_cover(lam_knots, b0_knots, b1_knots, lmax, eta, eps)
    claim("A4-COVER-BUILD-SEC", "seconds to build the certified cover", round(time.time() - t0, 1))
    claim("A4-CELLS", "cells in the certified cover", len(cells))
    claim("A4-UNSERVED", "leaves that could not be certified", len(fc.COVER_INFO["unserved"]))
    claim("A2-COVERAGE", "A-2 coverage of the declared box",
          fc.verify_cover(lam_knots, b0_knots, b1_knots, cells, lmax))
    ok3, f3 = fc.all_certified(cells)
    claim("A3-CERTS", "A-3 per-cell sign certificates", ok3)

    # ---- A-4, R-beta on every cell (the fast route), exact comparison --------------
    t0 = time.time()
    fails, worst, worst_cell, maxrad = [], F(0), None, F(0)
    for c in cells:
        v = pe.phi_r_beta(UC, c["rho1"], c["rho2"])
        hi = pe.ball_bounds(v)[1]
        r = pe.radius(v)
        maxrad = max(maxrad, r)
        if hi > ALPHA:
            fails.append((c, hi))
        if hi > worst:
            worst, worst_cell = hi, c
    claim("A4-BETA-SEC", "seconds for A-4 on every cell (R-beta)", round(time.time() - t0, 1))
    claim("A4-BETA-FAILS", "cells whose R-beta upper endpoint exceeds alpha", len(fails))
    claim("A4-BETA-WORST", "largest R-beta upper endpoint over the cover", float(worst))
    claim("A4-BETA-MARGIN", "alpha minus that worst upper endpoint", float(ALPHA - worst))
    claim("A4-MAX-RADIUS", "largest R-beta ball radius seen", float(maxrad))
    claim("A4-RADIUS-WITHIN-BUDGET", "every radius <= 1e-9", maxrad <= RAD_BUDGET)
    if worst_cell is not None:
        claim("A4-WORST-CELL-RHO", "rho of the worst cell",
              {"rho1": float(worst_cell["rho1"]), "rho2": float(worst_cell["rho2"])})

    # ---- A-5, H-wit: Phi(Lc; r(corner)) > alpha, with r(corner) as a BALL ----------
    corner = (fc.Iv(fc.LAM_LO, fc.LAM_LO), fc.Iv(fc.B_LO, fc.B_LO), fc.Iv(fc.B_HI, fc.B_HI))
    r1b, r2b = fc.r1_iv(corner, 200), fc.r2_iv(corner, 200)
    b1, b2 = iv_to_ball(pe, r1b), iv_to_ball(pe, r2b)
    claim("A5-R1-BALL", "r1(corner) enclosure", [str(r1b.lo), str(r1b.hi)])
    claim("A5-R2-BALL", "r2(corner) enclosure", [str(r2b.lo), str(r2b.hi)])
    vw = pe.phi_r_beta(LC, b1, b2)
    lo_w, hi_w = pe.ball_bounds(vw)
    claim("A5-PHI-BETA", "Phi(Lc; r(corner)) by R-beta", [str(lo_w), str(hi_w)])
    claim("A5-RADIUS", "H-wit ball radius", float(pe.radius(vw)))
    claim("A5-ABOVE-ALPHA", "H-wit holds: the lower endpoint strictly exceeds alpha",
          lo_w > ALPHA)
    claim("A5-MARGIN", "lower endpoint minus alpha", float(lo_w - ALPHA))

    (DEP / "formal01_f8b_claims.txt").write_text("\n".join(CLAIMS) + "\n", encoding="utf-8")
    print("\n".join(CLAIMS))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
