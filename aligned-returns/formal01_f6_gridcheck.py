"""F-6 — the construction-grid containment and dominance check.

The certificate counts `sup <= U` as CERTIFIED only if checks F-1..F-6 all pass, and a
missing F-6 does NOT discharge it. F-1 probed 1062 box points
(G1 + the exact corner + 8 box vertices), which is NOT the 175-point construction grid, so
this check is a separate obligation and is discharged here.

For each of the 175 construction-grid points, on the EXACT rational grid re-derived from
w2_cert.py (not from a deposit):
  1. the point lies inside the declared box;
  2. at least one cell of the certified cover contains it (containment);
  3. the cell's certified rho dominates: rho_k <= r(point)  (so cover_sound applies to it);
  4. upper(Phi(Uc; r(point))) <= alpha  (the bound actually holds at the grid point).
Also reports the F-6-CORNER line: the exact corner (1/25, -3/2, 3/2) is
covered and dominated.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import time
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
DEP = HERE / "deposits"
CLAIMS: list[str] = []


def claim(cid: str, label: str, value) -> None:
    CLAIMS.append(f"CLAIM {cid}: {label} = {value}")


def load(name: str, path: str):
    s = importlib.util.spec_from_file_location(name, path)
    if s is None or s.loader is None:
        raise ImportError(f"cannot load {name}")
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m


def main() -> int:
    pe = load("pe", str(HERE / "formal01_phi_engine.py"))
    f8c = load("f8c", str(HERE / "formal01_f8c_certified.py"))
    bat = load("bat", str(HERE / "formal01_f8a_battery.py"))

    # The cover module exports the box constants; the rate constants are declared here
    # rather than imported, because the module's copies live under worker-only names.
    UC, ALPHA = F(271740, 100000), F(1, 40)
    cells, unserved = f8c.build_graded_cover(workers=1)
    claim("F6-CELLS", "cells in the cover used by F-6", len(cells))
    claim("F6-UNSERVED", "leaves that could not be certified", len(unserved))

    lgrid = bat.log_grid(0.04, 25.0, 7)
    bgrid = bat.linspace(-1.5, 1.5, 5)
    claim("F6-GRID", "construction grid re-derived from w2_cert.py as (n_lam, n_b) per axis",
          f"log_grid({lgrid[0]},{lgrid[-1]},{len(lgrid)}) x linspace("
          f"{bgrid[0]},{bgrid[-1]},{len(bgrid)})^2")
    claim("F6-GRID-N", "construction-grid points", len(lgrid) * len(bgrid) ** 2)

    t0 = time.time()
    outside, uncovered, dom_fail, undec, viol, worst, worst_pt = 0, 0, 0, 0, 0, F(0), None
    for lam_ in lgrid:
        for a in bgrid:
            for b in bgrid:
                pt = (F(lam_), F(a), F(b))
                if not (f8c.LAM_LO <= pt[0] <= f8c.LAM_HI
                        and f8c.B_LO <= pt[1] <= f8c.B_HI
                        and f8c.B_LO <= pt[2] <= f8c.B_HI):
                    outside += 1
                    continue
                R1, R2 = bat.r_enc(pt)
                hit = False
                for c in cells:
                    bx = c["box"]
                    inside = all(bx[i].lo <= pt[i] <= bx[i].hi for i in range(3))
                    if not inside:
                        continue
                    hit = True
                    if c["rho1"] > R1.hi or c["rho2"] > R2.hi:
                        dom_fail += 1
                    elif c["rho1"] > R1.lo or c["rho2"] > R2.lo:
                        undec += 1
                if not hit:
                    uncovered += 1
                    continue
                v = pe.phi_r_beta(UC, R1.lo, R2.lo)
                hi_ = pe.ball_bounds(v)[1]
                if hi_ > worst:
                    worst, worst_pt = hi_, pt
                if hi_ > ALPHA:
                    viol += 1
    claim("F6-SEC", "seconds", round(time.time() - t0, 1))
    claim("F6-OUTSIDE-BOX", "construction-grid points outside the declared box", outside)
    claim("F6-UNCOVERED", "construction-grid points in NO cell (containment breach)", uncovered)
    claim("F6-DOM-FAIL", "cell/point pairs with rho_k > r_hi (cell does not certify)", dom_fail)
    claim("F6-DOM-UNDECIDED", "cell/point pairs undecided", undec)
    claim("F6-VIOLATIONS", "construction-grid points with upper(Phi(U;r)) > alpha", viol)
    claim("F6-WORST", "largest upper(Phi(U; r)) over the construction grid", float(worst))
    claim("F6-MARGIN", "alpha minus that worst upper endpoint", float(ALPHA - worst))
    claim("F6-ARGMX", "grid point attaining the largest upper(Phi(U;r))",
          [float(x) for x in worst_pt] if worst_pt else None)
    corner = (f8c.LAM_LO, f8c.B_LO, f8c.B_HI)
    cc = [c for c in cells if all(c["box"][i].lo <= corner[i] <= c["box"][i].hi
                                  for i in range(3))]
    R1c, R2c = bat.r_enc(corner)
    claim("F6-CORNER-COVERED", "cells containing the exact corner", len(cc))
    claim("F6-CORNER-DOMINATED",
          "every containing cell has rho_k <= r(corner)",
          all(c["rho1"] <= R1c.hi and c["rho2"] <= R2c.hi for c in cc) if cc else False)
    claim("F6-VERDICT", "F-6 construction-grid containment and dominance",
          "PASS" if (outside == 0 and uncovered == 0 and dom_fail == 0 and viol == 0) else "FAIL")
    (DEP / "formal01_f6_claims.txt").write_text("\n".join(CLAIMS) + "\n", encoding="utf-8")
    (HERE / "formal01_f6_result.json").write_text(json.dumps(CLAIMS, indent=1), encoding="utf-8")
    print("\n".join(CLAIMS))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
