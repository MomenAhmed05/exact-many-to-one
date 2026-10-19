"""F-8 CHECKER, Stage A, part 1 — the exact-rational half.

Implements acceptance rules A-2 (coverage) and A-3 (per-cell sign certificates), as
declared, and the falsifiers that exercise them with NO numerical integration:
FK-1 (dropped leaf), FK-2 (moved breakpoint), FK-3 (inflated rho_1), FK-4 (swapped rhos),
FK-10 (interior minimum that vertex values miss).

Why interval extension and not Bernstein: the Bernstein range-enclosure lemma would need
its own justification, so it is not used here. Natural interval extension is
sound by the DEFINITION of interval arithmetic -- each operation's interval extension
contains the true result for every input drawn from the inputs' intervals -- so no imported
lemma is relied upon. It is conservative relative to Bernstein and can fail where Bernstein
succeeds; that is recorded as a finding, never used as a licence to fall back on the lemma.

Every model formula below is transcribed from the APPROVED Lean statements in
lean/W4lean/Unit2BFormalR1.lean (aF, g00, g11, g01, trG, detG, charQ, r1, r2, Box), so a
divergence between the checker and the theorem is a build failure here rather than a silent
pass at the gate.
"""

from __future__ import annotations

import pathlib
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parent
DEP = ROOT / "deposits"
CLAIMS: list[str] = []


def claim(cid: str, label: str, value) -> None:
    CLAIMS.append(f"CLAIM {cid}: {label} = {value}")


# ------------------------------------------------------------- exact intervals ----

class Iv:
    """A closed exact-rational interval [lo, hi]. Invariant: lo <= hi."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        if hi is None:
            hi = lo
        if lo > hi:
            lo, hi = hi, lo
        self.lo, self.hi = lo, hi

    def __add__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        return Iv(self.lo - o.hi, self.hi - o.lo)

    def __mul__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        c = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv(min(c), max(c))

    def sq(self):
        if self.lo >= 0:
            return Iv(self.lo * self.lo, self.hi * self.hi)
        if self.hi <= 0:
            return Iv(self.hi * self.hi, self.lo * self.lo)
        return Iv(F(0), max(self.lo * self.lo, self.hi * self.hi))

    def scale(self, k):
        """Multiply by the exact rational k (k >= 0)."""
        return Iv(self.lo * k, self.hi * k)

    def __truediv__(self, o):
        """Exact interval division. Refuses a divisor straddling zero rather than
        silently returning an unbounded enclosure."""
        o = o if isinstance(o, Iv) else Iv(o)
        if o.contains_zero():
            raise ZeroDivisionError(f"interval divisor {o} contains 0; result unbounded")
        c = (self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi)
        return Iv(min(c), max(c))

    def __repr__(self) -> str:
        return f"[{self.lo}, {self.hi}]"

    def le_zero(self) -> bool:
        return self.hi <= 0

    def ge_zero(self) -> bool:
        return self.lo >= 0

    def contains_zero(self) -> bool:
        return self.lo <= 0 <= self.hi


def isqrt_interval(x: Iv, prec: int = 60) -> Iv:
    """Rigorous enclosure of sqrt([x.lo, x.hi]): [floor(sqrt(x.lo)), ceil(sqrt(x.hi))] on
    the 2^-prec grid, by exact integer square roots.

    The previous version bisected with the single test `mid*mid <= x.lo`, so BOTH endpoints
    converged on sqrt(x.lo) and the upper endpoint never saw x.hi: isqrt_interval([1,4])
    returned [1,1]. On a point interval that is harmless, but on a cell it made r1_iv.hi
    fall below the true r1 (0.5050090 against 0.5050254 at the lambda-min/b0-min corner of
    width 2^-11) and r2_iv.lo rise above the true r2. Each endpoint is now rounded outward
    from ITS OWN radicand endpoint.
    """
    from math import isqrt
    if x.hi < 0:
        raise ValueError("negative radicand")
    xlo = max(x.lo, F(0))
    xhi = F(x.hi)
    xlo = F(xlo)
    lo = isqrt((xlo.numerator << (2 * prec)) // xlo.denominator)   # floor(sqrt(xlo) 2^p)
    m = -((-xhi.numerator << (2 * prec)) // xhi.denominator)       # ceil(xhi 4^p)
    hi = isqrt(m)
    if hi * hi < m:
        hi += 1                                                     # ceil(sqrt(m)) >= sqrt(xhi 4^p)
    return Iv(F(lo, 1 << prec), F(hi, 1 << prec))


# ------------------------------------------- model, exactly as the approved Lean ----
# theta = (lam, b0, b1);  Box = {1/25 <= lam <= 25, -3/2 <= b0 <= 3/2, -3/2 <= b1 <= 3/2}
LAM_LO, LAM_HI = F(1, 25), F(25)
B_LO, B_HI = F(-3, 2), F(3, 2)
BOX = (Iv(LAM_LO, LAM_HI), Iv(B_LO, B_HI), Iv(B_LO, B_HI))


def aF_pt(lam, b0, b1):
    return (b1 - b0) ** 2 + 1 + lam


def trG_pt(lam, b0, b1):
    return (b0 * b0 + 1) + (b1 * b1 + lam)


def detG_pt(lam, b0, b1):
    return (b0 * b0 + 1) * (b1 * b1 + lam) - (b0 * b1) ** 2


def charQ_pt(m, lam, b0, b1):
    return m * m - trG_pt(lam, b0, b1) * m + detG_pt(lam, b0, b1)


def _iv(fn, box):
    """Interval extension of a point function on a box, by the given closed form."""
    raise NotImplementedError


def aF_iv(b: tuple[Iv, Iv, Iv]) -> Iv:
    d = b[2] - b[1]                      # b1 - b0
    return d.sq() + Iv(1) + b[0]


def trG_iv(b):
    return b[1].sq() + Iv(1) + b[2].sq() + b[0]


def detG_iv(b):
    return (b[1].sq() + Iv(1)) * (b[2].sq() + b[0]) - (b[1] * b[2]).sq()


def charQ_iv(m: Iv, b) -> Iv:
    return m.sq() - trG_iv(b) * m + detG_iv(b)


def r1_iv(b, prec=60):
    """Enclosure of r1 = mu1/aF, mu1 = (trG + sqrt(disc))/2."""
    tr, det = trG_iv(b), detG_iv(b)
    disc = tr.sq() - det.scale(4)
    sq = isqrt_interval(disc, prec)
    mu1 = (tr + sq).scale(F(1, 2))
    return mu1 / aF_iv(b)


def r2_iv(b, prec=60):
    """Enclosure of r2 = mu2/aF.

    Uses mu2 = 2*detG/(tr + sqrt(tr^2 - 4*detG)) rather than (tr - sqrt(...))/2. Both
    are algebraically equal, but the first does not cancel: on a whole cell the second
    subtracts two nearly equal quantities whose intervals are independent, which widens
    the enclosure by ~2x and destroys the hLow2 certificate at every refinement level.
    """
    tr, det = trG_iv(b), detG_iv(b)
    disc = tr.sq() - det.scale(4)
    sq = isqrt_interval(disc, prec)
    denom = tr + sq
    mu2 = det.scale(2) / denom
    return mu2 / aF_iv(b)


def hlow2_holds(rho2, b, prec=200):
    """The two hLow2 conditions over the whole cell, by exact interval extension."""
    if rho2 <= 0:
        return False
    lin = (m_iv(rho2, b)) - trG_iv(b).scale(F(1, 2))
    if not lin.le_zero():
        return False
    q = charQ_iv(m_iv(rho2, b), b)
    return q.ge_zero()


def m_iv(rho, b):
    """rho * aF(theta) as an interval over the cell."""
    return Iv(rho) * aF_iv(b)


def largest_certified_rho2(b, prec=200, steps=200):
    """Largest certified rho2 for hLow2 (a CONJUNCTION, so its admissible set is an
    interval [0, rho2_inf]).  Same candidate-then-descent shape as rho1: try r2_iv's
    upper endpoint, its lower endpoint, and rho_lin, then a short geometric descent."""
    r2b = r2_iv(b, prec)
    lin = (trG_iv(b).scale(F(1, 2))) / aF_iv(b)
    best = None
    for cand in (r2b.hi, r2b.lo, lin.lo):
        if cand > 0 and hlow2_holds(cand, b, prec):
            best = cand if best is None else max(best, cand)
    if best is not None:
        return best
    rho = r2b.hi
    for _ in range(steps):
        rho = rho / 2
        if hlow2_holds(rho, b, prec):
            return rho
    return None


# ------------------------------------------------- A-3: per-cell sign certificates -


def cert_cell(rho1: F, rho2: F, b: tuple[Iv, Iv, Iv]) -> dict:
    """Check hLow1 and hLow2 over the whole cell by exact rational interval extension.

    hLow1 : rho1*aF <= trG/2   OR   charQ(rho1*aF) <= 0
    hLow2 : rho2*aF <= trG/2   AND  charQ(rho2*aF) >= 0
    """
    if rho1 <= 0 or rho2 <= 0:
        return {"ok": False, "why": "rho not strictly positive (Real.sqrt is 0 on negatives)"}
    lin1 = (Iv(rho1) * aF_iv(b)) - trG_iv(b).scale(F(1, 2))
    q1 = charQ_iv(Iv(rho1) * aF_iv(b), b)
    lin2 = (Iv(rho2) * aF_iv(b)) - trG_iv(b).scale(F(1, 2))
    q2 = charQ_iv(Iv(rho2) * aF_iv(b), b)
    low1_ok = lin1.le_zero() or q1.le_zero()
    low1_branch = ("linear" if lin1.le_zero() else
                   "quadratic" if q1.le_zero() else "none")
    low2_ok = lin2.le_zero() and q2.ge_zero()
    return {"ok": low1_ok and low2_ok, "low1": low1_ok, "low1_branch": low1_branch,
            "low2": low2_ok, "lin1": lin1, "q1": q1, "lin2": lin2, "q2": q2}


# -------------------------------------------------- A-2: coverage of the box --------


def cells_of_block(lam_knots, b0_knots, b1_knots, cells, ib):
    return [c for c in cells if c["block"] == ib]


def eight_children(box: tuple[Iv, Iv, Iv]):
    b = box
    l0, h0, l1, h1, l2, h2 = b[0].lo, b[0].hi, b[1].lo, b[1].hi, b[2].lo, b[2].hi
    out = []
    for a, b in ((l0, (l0 + h0) / 2), ((l0 + h0) / 2, h0)):
        for c, d in ((l1, (l1 + h1) / 2), ((l1 + h1) / 2, h1)):
            for e, f in ((l2, (l2 + h2) / 2), ((l2 + h2) / 2, h2)):
                out.append((Iv(a, b), Iv(c, d), Iv(e, f)))
    return out


def verify_block(leaves, block, lmax, hard_cap=64):
    """Top-down reconstruction: a node is covered iff it is a leaf, or all 8 children are.

    Also rejects a leaf that is a strict ancestor of another leaf (no leaf is an ancestor
    of another), which is the FK-1 failure class.

    The recursion is NOT stopped at the producer's declared `lmax`: the producer may
    refine deeper than that to certify a cell whose enclosure is too wide (an algebraic,
    Phi-free rule), and a verifier that refuses to descend past the declared depth would
    reject a complete cover. `hard_cap` exists only to bound a pathological descent; the
    cell set is finite, so a genuine dyadic refinement terminates well inside it.
    """
    keys = {c["key"] for c in leaves}

    def key(box):
        return (str(box[0].lo), str(box[0].hi), str(box[1].lo), str(box[1].hi),
                str(box[2].lo), str(box[2].hi))

    # no leaf may strictly contain another
    boxes = [(c["key"], c["box"]) for c in leaves]
    for (k1, b1), (k2, b2) in ((a, b) for i, a in enumerate(boxes) for b in boxes[i + 1:]):
        if k1 == k2:
            continue
        for ax, bx in zip(b1, b2, strict=True):
            if not (ax.lo <= bx.lo and bx.hi <= ax.hi and (ax.lo, ax.hi) != (bx.lo, bx.hi)):
                break
        else:
            return False, "a leaf strictly contains another leaf"

    def rec(box, depth):
        k = key(box)
        if k in keys:
            return True
        if depth >= hard_cap:
            return False
        return all(rec(ch, depth + 1) for ch in eight_children(box))

    if not rec(block, 0):
        return False, "block not fully covered by the submitted leaves"
    return True, "ok"


def verify_cover(lam_knots, b0_knots, b1_knots, cells, lmax):
    """A-2 (i) blocks partition the box; (ii) each block's leaves form a complete dyadic
    partition; (iii) regenerating from the parameter vector reproduces the list exactly."""
    for name, knots in (("lam", lam_knots), ("b0", b0_knots), ("b1", b1_knots)):
        if len(knots) < 2 or any(knots[i] >= knots[i + 1] for i in range(len(knots) - 1)):
            return False, f"{name} breakpoints not strictly increasing"
    if (lam_knots[0], lam_knots[-1]) != (LAM_LO, LAM_HI):
        return False, "lam breakpoints do not span the declared box"
    if (b0_knots[0], b0_knots[-1]) != (B_LO, B_HI):
        return False, "b0 breakpoints do not span the declared box"
    if (b1_knots[0], b1_knots[-1]) != (B_LO, B_HI):
        return False, "b1 breakpoints do not span the declared box"
    n = 0
    for i in range(len(lam_knots) - 1):
        for j in range(len(b0_knots) - 1):
            for k in range(len(b1_knots) - 1):
                block = (Iv(lam_knots[i], lam_knots[i + 1]),
                         Iv(b0_knots[j], b0_knots[j + 1]),
                         Iv(b1_knots[k], b1_knots[k + 1]))
                ok, why = verify_block(cells_of_block(lam_knots, b0_knots, b1_knots,
                                                      cells, (i, j, k)), block, lmax)
                if not ok:
                    return False, f"block {(i, j, k)}: {why}"
                n += 1
    return True, f"{n} blocks, all fully covered"


def _serve(box, depth, cap, prec):
    """Serve a box: try to certify rho2 at this depth; if no positive certified rho2
    exists, split 8 ways and recurse. This is an ALGEBRAIC refinement rule -- it is driven
    by the exact rational sign conditions only, never by the rate, so the cover still does
    not depend on any numerical acceptance."""
    rho2 = largest_certified_rho2(box, prec)
    if rho2 is not None and rho2 > 0:
        return rho2, depth
    if depth >= cap:
        return None, depth
    for ch in eight_children(box):
        r2, d2 = _serve(ch, depth + 1, cap, prec)
        if r2 is None:
            return None, depth
    return None, depth


def make_cell(lam_knots, b0_knots, b1_knots, blk, box, rho1, rho2):
    return {"block": blk, "box": box, "rho1": rho1, "rho2": rho2,
            "key": (str(box[0].lo), str(box[0].hi), str(box[1].lo), str(box[1].hi),
                    str(box[2].lo), str(box[2].hi))}


# ------------------------------------------------------- cover construction --------
# Static scheme: a declared rational tensor partition, dyadic 8-way
# refinement inside each block, split iff depth < lmax and width_w(N) > eta*(d(N)+eps),
# where d(N) is the weighted inf-norm distance from N to the critical set. Phi NEVER
# influences a split. Here the weight is the sup-norm in all three coordinates.


def box_width(b):
    return max(hh.hi - hh.lo for hh in b)


def dist_to_critical(b):
    """Weighted inf-norm distance from the cell to the critical set
    C = {(1/25,-3/2,3/2), (1/25,3/2,-3/2)}."""
    def coord_d(c, corner):
        return max(abs(c.lo - corner), abs(c.hi - corner))
    return max(min(coord_d(b[0], LAM_LO),
                    coord_d(b[1], B_LO), coord_d(b[2], B_HI)),
               min(coord_d(b[0], LAM_LO),
                   coord_d(b[1], B_HI), coord_d(b[2], B_LO)))


def build_block_leaves(block, lmax, eta, eps, w=None):
    """Return the leaves of one block under the declared static rule."""
    w = w or (F(1), F(1), F(1))
    leaves = []
    stack = [(block, 0)]
    while stack:
        node, depth = stack.pop()
        if depth >= lmax:
            leaves.append(node)
            continue
        width_w = max(w[i] * (node[i].hi - node[i].lo) for i in range(3))
        d = dist_to_critical(node)
        if width_w <= eta * (d + eps):
            leaves.append(node)
            continue
        for ch in eight_children(node):
            stack.append((ch, depth + 1))
    return leaves


def make_cover(lam_knots, b0_knots, b1_knots, lmax, eta, eps, prec=200, extra_levels=8):
    """Build the static cover, then REFINE any leaf that cannot be certified.

    The refine-never-drop rule is the direct lesson of the barred producer, whose
    `cap_split` discarded half its children and therefore never covered the box. A leaf
    that cannot be served is replaced by its eight children; it is NEVER dropped. If a
    branch is still unserved at the depth cap it is reported in `unserved` and the run is
    OPEN, and A-2 will reject the cover.
    """
    leaves = []
    for i in range(len(lam_knots) - 1):
        for j in range(len(b0_knots) - 1):
            for k in range(len(b1_knots) - 1):
                block = (Iv(lam_knots[i], lam_knots[i + 1]),
                         Iv(b0_knots[j], b0_knots[j + 1]),
                         Iv(b1_knots[k], b1_knots[k + 1]))
                blk = (i, j, k)
                for bx in build_block_leaves(block, lmax, eta, eps):
                    leaves.append((blk, bx, 0))
    cells, unserved = [], []
    pending = leaves
    for _ in range(extra_levels):
        nxt = []
        for blk, box, depth in pending:
            rho1 = largest_certified_rho1(box, prec)
            rho2 = largest_certified_rho2(box, prec)
            if rho1 is not None and rho2 is not None and rho1 > 0 and rho2 > 0:
                cells.append(make_cell(lam_knots, b0_knots, b1_knots, blk, box,
                                       rho1, rho2))
            else:
                nxt.append((blk, box, depth + 1))
        if not nxt:
            # every box in `pending` was served in this round; clear it so the final
            # pass does NOT serve them a second time. Leaving it set produced 624
            # duplicate keys, which made FK-1 (drop one leaf) vacuous: removing a
            # duplicate changes nothing, so the coverage check passed on a cover that
            # had lost a leaf. That is the exact class of defect FK-1 exists to catch,
            # reproduced inside the checker meant to catch it.
            pending = []
            break
        pending = [(blk, ch, depth + 1) for blk, box, depth in nxt
                   for ch in eight_children(box)]
    # Final pass: serve whatever is still pending. A leaf that still cannot be served is
    # NEVER dropped: it is reported and A-2 rejects the cover, making the run OPEN.
    for blk, box, _depth in pending:
        rho1 = largest_certified_rho1(box, prec)
        rho2 = largest_certified_rho2(box, prec)
        if rho1 is not None and rho2 is not None and rho1 > 0 and rho2 > 0:
            cells.append(make_cell(lam_knots, b0_knots, b1_knots, blk, box,
                                   rho1, rho2))
        else:
            unserved.append(box)
    COVER_INFO["unserved"] = unserved
    COVER_INFO["n_leaves_initial"] = len(leaves)
    COVER_INFO["n_cells"] = len(cells)
    return cells


COVER_INFO: dict = {}


def all_certified(cells, prec=200):
    """A-3 for every cell. Returns (ok, first_failure)."""
    for c in cells:
        res = cert_cell(c["rho1"], c["rho2"], c["box"])
        if not res["ok"]:
            return False, (c, res)
        if c["rho1"] <= 0 or c["rho2"] <= 0:
            return False, (c, {"why": "rho not positive"})
    return True, None


# ============================================================== Stage A harness ====

def declare_parameters():
    """The parameter vector, DECLARED BEFORE any run (A-6).
    Phi never influences a split."""
    lam_knots = [F(1, 25), F(1), F(25)]
    b0_knots = [F(-3, 2), F(0), F(3, 2)]
    b1_knots = [F(-3, 2), F(0), F(3, 2)]
    lmax, eta, eps = 3, F(1, 2), F(1, 1 << 20)
    return lam_knots, b0_knots, b1_knots, lmax, eta, eps


def critical_cell(cells):
    for c in cells:
        b = c["box"]
        if (b[0].lo <= LAM_LO <= b[0].hi and b[1].lo <= B_LO <= b[1].hi
                and b[2].lo <= B_HI <= b[2].hi):
            return c
    return None


def fq(rho, lam, b0, b1):
    """f(theta) = charQ(rho*aF(theta), theta), the quantity hLow2 must be non-negative on."""
    return charQ_pt(rho * aF_pt(lam, b0, b1), lam, b0, b1)


def fk10_search(box, rho_grid):
    """Find rho with f >= 0 at all 8 VERTICES but < 0 at an interior point. A-3 must
    REJECT such a cell: a vertex-only check would wrongly accept it."""
    verts = [tuple(Iv(box[i].hi if (n >> i) & 1 else box[i].lo) for i in range(3))
             for n in range(8)]
    interior = [tuple(Iv(box[i].lo + (box[i].hi - box[i].lo) * F(num, den))
                      for i in range(3))
                for num in (1, 3, 5) for den in (4, 6, 8)]
    for rho in rho_grid:
        if rho <= 0:
            continue
        vmin = min(fq(rho, v[0].lo, v[1].lo, v[2].lo) for v in verts)
        if vmin < 0:
            continue
        imin = min(fq(rho, v[0].lo, v[1].lo, v[2].lo) for v in interior)
        if imin < 0:
            return rho, vmin, imin
    return None




def interior_point_violates(box, rho, axis=0):
    """Look for an interior point where charQ(rho*aF) < 0 while the centre is >= 0."""
    ctr = [Iv((x.lo + x.hi) / 2) for x in box]
    for num in (1, 3, 7):
        for den in (4, 8):
            off = (box[axis].hi - box[axis].lo) * F(num, den)
            pt = list(ctr)
            pt[axis] = Iv((box[axis].lo + box[axis].hi) / 2 + off)
            for sgn in (1, -1):
                q = ctr[axis] if sgn > 0 else pt[axis]
                p2 = list(ctr)
                p2[axis] = Iv(q.lo - off if sgn > 0 else q.lo + off)
                val = charQ_pt(rho * aF_pt(*[x.lo for x in ctr]), ctr[0].lo, ctr[1].lo,
                               ctr[2].lo)
                if val >= 0:
                    # check a real interior point for negativity
                    for t in (F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(7, 8)):
                        ip = list(ctr)
                        ip[axis] = Iv(ctr[axis].lo + (box[axis].hi - box[axis].lo) * t)
                        v = charQ_pt(rho * aF_pt(ip[0].lo, ip[1].lo, ip[2].lo),
                                     ip[0].lo, ip[1].lo, ip[2].lo)
                        if v < 0:
                            return True
    return False


def main() -> int:
    lam_knots, b0_knots, b1_knots, lmax, eta, eps = declare_parameters()
    cells = make_cover(lam_knots, b0_knots, b1_knots, lmax, eta, eps)
    claim("A3-CELLS", "cells in the declared static cover", len(cells))

    ok2, why2 = verify_cover(lam_knots, b0_knots, b1_knots, cells, lmax)
    claim("A2-COVERAGE", "A-2 coverage verified", ok2)
    claim("A2-DETAIL", "A-2 outcome", why2)

    ok3, fail3 = all_certified(cells)
    claim("A3-CERTS", "A-3 per-cell sign certificates verified for every cell", ok3)
    if fail3 is not None:
        c, res = fail3
        claim("A3-FIRST-FAIL", "first cell failing A-3",
              {"rho1": str(c["rho1"]), "rho2": str(c["rho2"]),
               "why": res.get("why") or f"low1={res.get('low1')}/{res.get('low1_branch')},"
                                         f" low2={res.get('low2')}"})
    else:
        claim("A3-FIRST-FAIL", "first cell failing A-3", "none")

    crit = critical_cell(cells)
    claim("A3-CRITICAL-CELL-FOUND", "the cell containing the critical point exists",
          crit is not None)

    # ---- FK-1: drop one leaf -> A-2 must fail
    if cells:
        dropped = cells[:-1]
        ok, why = verify_cover(lam_knots, b0_knots, b1_knots, dropped, lmax)
        claim("FK1-DETECTED", "FK-1 dropped leaf is rejected by A-2", not ok)
        claim("FK1-WHY", "FK-1 reason", why)

    # ---- FK-2: move one breakpoint by 2^-60 -> A-2 must fail
    bk = list(b0_knots)
    bk[1] = bk[1] + F(1, 1 << 60)
    ok, why = verify_cover(bk, b0_knots, b1_knots, cells, lmax)
    claim("FK2-DETECTED", "FK-2 moved breakpoint is rejected by A-2", not ok)
    claim("FK2-WHY", "FK-2 reason", why)

    # ---- FK-3: inflate rho1 at the critical cell -> A-3 must fail
    if crit is not None:
        bad = [dict(c) for c in cells]
        idx = cells.index(crit)
        centre = [Iv((x.lo + x.hi) / 2) for x in crit["box"]]
        bad[idx] = dict(crit, rho1=r1_iv(centre).hi + F(1, 10 ** 6))
        ok, _ = all_certified(bad)
        claim("FK3-DETECTED", "FK-3 inflated rho_1 is rejected by A-3", not ok)

    # ---- FK-4: swap rho1 and rho2 -> A-3 must fail somewhere
    swapped = [dict(c, rho1=c["rho2"], rho2=c["rho1"]) for c in cells]
    ok, fail = all_certified(swapped)
    claim("FK4-DETECTED", "FK-4 swapped rhos are rejected by A-3", not ok)
    claim("FK4-CELLS-AFFECTED", "cells where swapping changes the verdict",
          "some" if not ok else "none")

    # ---- FK-10: vertex-positive but interior-negative -> A-3 must NOT accept
    fk10 = None
    for c in cells:
        b = c["box"]
        base = r2_iv(b, 400).hi
        grid = [base + F(k, 1 << 20) * base for k in range(0, 4000, 1)]
        fk10 = fk10_search(b, grid)
        if fk10:
            rho, vmin, imin = fk10
            res = cert_cell(c["rho1"], rho, b)
            claim("FK10-DETECTED", "FK-10 vertex-positive/interior-negative is REJECTED by A-3",
                  not res["ok"])
            claim("FK10-EVIDENCE", "FK-10 vertex min and interior min of charQ(rho*aF)",
                  {"vertex_min": str(vmin), "interior_min": str(imin)})
            break
    if not fk10:
        claim("FK10-DETECTED", "FK-10 vertex-positive/interior-negative is REJECTED by A-3",
              "no such cell found in the searched cover")

    (DEP / "formal01_f8a_claims.txt").write_text("\n".join(CLAIMS) + "\n", encoding="utf-8")
    print("\n".join(CLAIMS))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# ---------------------------------------------------- larger certified rho1 ---------
# hLow1 admits rho1 up to r1 via the QUADRATIC branch (charQ(rho1*aF) <= 0, i.e.
# rho1*aF in [mu2, mu1]). The linear branch only reaches (mu1+mu2)/2, which is far
# below r1 and makes Phi(Uc; rho) exceed alpha -- caught by acceptance rule A-4.

def hlow1_holds(rho1, b, prec=200):
    """The hLow1 disjunction over the whole cell, by exact interval extension."""
    if rho1 <= 0:
        return False
    lin = m_iv(rho1, b) - trG_iv(b).scale(F(1, 2))
    if lin.le_zero():
        return True
    return charQ_iv(m_iv(rho1, b), b).le_zero()


def largest_certified_rho1(b, prec=200, steps=200):
    """Largest certified rho1 for hLow1, with the disjoint admissible set handled
    ANALYTICALLY rather than by ladder search.

    hLow1 is a DISJUNCTION, so the admissible set is
        [0, rho_lin]  UNION  [rho2_sup, rho1_inf]
    where rho_lin = inf_cell (trG/2)/aF.  A ladder or a bisection is only valid under a
    monotonicity that does not hold here, and a 800-step ladder costs ~7.4e6 interval
    evaluations over the cover.  Instead the two candidate maxima are evaluated directly:
    rho_lin (always admissible via the first disjunct) and r1_iv(b).lo (a certified lower
    bound for min_cell r1, admissible via the second when the certificate holds), with a
    short geometric descent as a fallback.  O(few) evaluations per cell.
    """
    r1b = r1_iv(b, prec)
    lin = (trG_iv(b).scale(F(1, 2))) / aF_iv(b)
    cands = [r1b.hi, r1b.lo, lin.lo]
    best = None
    for cand in cands:
        if cand > 0 and hlow1_holds(cand, b, prec):
            best = cand if best is None else max(best, cand)
    if best is not None:
        return best
    rho = r1b.hi
    for _ in range(steps):
        rho = rho / 2
        if hlow1_holds(rho, b, prec):
            return rho
    return None
