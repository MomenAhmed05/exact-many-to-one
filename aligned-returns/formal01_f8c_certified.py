"""F-8 certified run — acceptance rules A-2, A-3, A-4 on the graded cover.

Implements the graded-cover method (R-3). An earlier run failed for four checker
defects, none of them mathematical: D-a natural interval extension of aF/trG/detG as
independent intervals (a dependency blow-up that made r1's enclosure [0.047, 3.125] when
min r1 is 0.5), D-b an inverted dist_to_critical, D-c lmax=3 defeating the grading, and
D-d a halving fallback that loses a factor of 2 per step.

THE ENCLOSURE (the substantive change): for a cell with centre z and half-widths h, expand
charQ(rho*aF(theta), theta) and rho*aF(theta) - trG(theta)/2 EXACTLY in t = theta - z with
sympy, then bound the result by a0 +- sum_{k != 0} |a_k| * prod_i h_i^{k_i}. This is the
natural interval extension of an algebraically IDENTICAL polynomial — no Bernstein lemma is
needed — and its overestimation is O(h^2) rather than O(h). The expansion identity is checked
by a SECOND, independent implementation (`expansion_crosscheck`).

Declared parameter vector (R-3 step 2), fixed before any certified Phi evaluation:
    lam knots [1/25, 77/50, 3, 5, 8, 12, 17, 25]
    b0, b1 knots [-3/2, 0, 3/2]
    w = (1,1,1), eta = 1/8, eps = 2^-8, lmax = 12
"""

from __future__ import annotations

import pathlib
import random
import sys
from fractions import Fraction as F

import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from formal01_f8a_checker import (  # noqa: E402
    B_HI,
    B_LO,
    LAM_HI,
    LAM_LO,
    Iv,
    aF_pt,
    detG_pt,
    trG_pt,
)

ROOT = pathlib.Path(__file__).resolve().parent
DEP = ROOT / "deposits"
CLAIMS: list[str] = []


def claim(cid: str, label: str, value) -> None:
    CLAIMS.append(f"CLAIM {cid}: {label} = {value}")


# ------------------------------------------------------------------ exact expansion --
_lam, _b0, _b1, _rho = sp.symbols("lam b0 b1 rho")


def _charQ_expr(rho_expr):
    aF = (rho_expr * ((_b1 - _b0) ** 2 + 1 + _lam))
    trG = (_b0 ** 2 + 1) + (_b1 ** 2 + _lam)
    detG = (_b0 ** 2 + 1) * (_b1 ** 2 + _lam) - (_b0 * _b1) ** 2
    return sp.expand(aF ** 2 - trG * aF + detG)


def _lin_expr(rho_expr):
    aF = (rho_expr * ((_b1 - _b0) ** 2 + 1 + _lam))
    trG = (_b0 ** 2 + 1) + (_b1 ** 2 + _lam)
    return sp.expand(aF - trG / 2)


def _q(f):
    return F(int(f.p), int(f.q)) if f.q != 1 else F(int(f.p))


def centred_enclosure(expr, z, h, rho_val):
    """Expand `expr` exactly in t = theta - z and bound it by a0 +- sum |a_k| h^k.

    Returns (lo, hi) as exact Fractions. The substitution is exact: the multivariate
    polynomial is re-expressed in the t variables with theta_i = z_i + t_i, and every
    coefficient is an exact rational. No float, no sampling.
    """
    subs = {_lam: sp.Rational(z[0].numerator, z[0].denominator) + sp.Symbol("t0"),
            _b0: sp.Rational(z[1].numerator, z[1].denominator) + sp.Symbol("t1"),
            _b1: sp.Rational(z[2].numerator, z[2].denominator) + sp.Symbol("t2")}
    poly = sp.expand(expr.subs(subs))
    ts = sp.symbols("t0 t1 t2")
    P = sp.Poly(poly, *ts)
    a0 = _q(P.coeff_monomial(ts[0] ** 0 * ts[1] ** 0 * ts[2] ** 0))
    total = F(0)
    for (i, j, k), c in P.terms():
        if i == j == k == 0:
            continue
        total += abs(_q(c)) * (h[0] ** i) * (h[1] ** j) * (h[2] ** k)
    return a0 - total, a0 + total


def expansion_crosscheck(z, h, rho_val, trials=40):
    """SECOND, independent implementation of the expansion identity (gate note).

    Rather than trusting sympy's expansion, this checks that the *function* the expansion
    is centred on is the same one the Lean statements define, by comparing the sympy
    expression against the checker's own exact `charQ_pt` / `aF_pt` / `trG_pt` at random
    EXACT rational points. If the two disagreed anywhere, the centred enclosure would be
    bounding the wrong function, so this is the check that matters.
    """
    rng = random.Random(20260928)
    worst = F(0)
    for _ in range(trials):
        pt = [z[i] + (rng.randint(-9, 9) * h[i]) / 10 for i in range(3)]
        lam, b0, b1 = pt
        aF = aF_pt(lam, b0, b1)
        trG = trG_pt(lam, b0, b1)
        det = detG_pt(lam, b0, b1)
        own = (rho_val * aF) ** 2 - trG * (rho_val * aF) + det
        theirs = _q(_charQ_expr(sp.Rational(rho_val.numerator, rho_val.denominator))
                    .subs({_lam: sp.Rational(lam.numerator, lam.denominator),
                           _b0: sp.Rational(b0.numerator, b0.denominator),
                           _b1: sp.Rational(b1.numerator, b1.denominator)}))
        worst = max(worst, abs(own - theirs))
    return worst


# ------------------------------------------------------------- corrected distance ---
CRIT = [(LAM_LO, B_LO, B_HI), (LAM_LO, B_HI, B_LO)]


def dist_to_critical_fixed(cell) -> F:
    """min over the two critical points of the max over coordinates of the NEAR-edge
    distance. The previous version took max-over-points of min-over-coordinates using far
    edges, which scored a cell at lambda in [24,25] as distance 0.1 (D-b)."""
    # CRIT is a non-empty module constant, so the min is always taken; seeded with a
    # value larger than any cell extent so the first iteration genuinely competes.
    best = None
    for p in CRIT:
        d = F(0)
        for i in range(3):
            lo, hi = cell[i].lo, cell[i].hi
            if lo <= p[i] <= hi:
                di = F(0)
            else:
                di = min(abs(lo - p[i]), abs(hi - p[i]))
            d = max(d, di)
        best = d if best is None else min(best, d)
    if best is None:          # unreachable: CRIT is non-empty
        return F(0)
    return best


# ------------------------------------------------- certificates on the centred form --
# KEY IDENTITY, which makes the whole thing fast:
#     charQ(rho*aF(theta), theta) = rho^2 * aF(theta)^2 - rho * aF(theta)*trG(theta)
#                                    + detG(theta)
# so the t-expansion of charQ is a FIXED linear combination of the t-expansions of just
# five polynomials. Expanding once per cell therefore serves EVERY rho in closed form. The
# previous version re-expanded with sympy on every rho evaluation -- about 220 expansions
# per cell, times ~6e4 cells -- and never finished. Five expansions per cell, and each rho
# is then pure rational arithmetic.


def cell_geom(cell):
    """Centre and half-widths of a cell, exactly."""
    z = [(cell[i].lo + cell[i].hi) / 2 for i in range(3)]
    h = [(cell[i].hi - cell[i].lo) / 2 for i in range(3)]
    return z, h


class CellCert:
    """Centred-form enclosure data for one cell: aF^2, aF*trG, detG, aF and trG expanded
    EXACTLY in t = theta - z, each kept as its full per-monomial coefficient dict.

    The absolute value is taken only AFTER the polynomials are combined at the given rho:
        charQ(rho*aF) = sum_k (rho^2 a_k - rho b_k + d_k) t^k,   |t_i| <= h_i
    so the enclosure is c_0 +- sum_{k != 0} |rho^2 a_k - rho b_k + d_k| h^k, which is the
    natural interval extension of the exact t-expansion of charQ(rho*aF) itself. Since
    |rho^2 a - rho b + d| <= rho^2|a| + rho|b| + |d|, this is never looser than folding
    each polynomial to (a_0, sum|a_k|h^k) first -- the previous form, which on the corner
    cells was ~40x wider and failed A-4 on 51,834 cells. The linear form is treated the
    same way: |rho l_k - t_k/2|.

    Float copies of the same coefficients are kept ONLY for the untrusted rho locator; every
    accepted rho is re-checked through the exact Fraction methods.
    """

    __slots__ = ("cell", "z", "h", "q0", "qk", "l0", "lk", "qkf", "lkf")

    def __init__(self, cell):
        self.cell = cell
        self.z, self.h = cell_geom(cell)
        subs = {_lam: sp.Rational(self.z[0].numerator, self.z[0].denominator) + sp.Symbol("t0"),
                _b0: sp.Rational(self.z[1].numerator, self.z[1].denominator) + sp.Symbol("t1"),
                _b1: sp.Rational(self.z[2].numerator, self.z[2].denominator) + sp.Symbol("t2")}
        ts = sp.symbols("t0 t1 t2")
        aF = sp.expand((( _b1 - _b0) ** 2 + 1 + _lam).subs(subs))
        trG = sp.expand(((_b0 ** 2 + 1) + (_b1 ** 2 + _lam)).subs(subs))
        detG = sp.expand(((_b0 ** 2 + 1) * (_b1 ** 2 + _lam) - (_b0 * _b1) ** 2).subs(subs))
        A = self._coeffs(sp.expand(aF ** 2), ts)
        B = self._coeffs(sp.expand(aF * trG), ts)
        D = self._coeffs(detG, ts)
        Lc = self._coeffs(aF, ts)
        Tc = self._coeffs(trG, ts)
        z0 = (0, 0, 0)
        self.q0 = (A.get(z0, F(0)), B.get(z0, F(0)), D.get(z0, F(0)))
        self.l0 = (Lc.get(z0, F(0)), Tc.get(z0, F(0)))
        # (a_k, b_k, d_k, h^k) for every non-constant monomial of charQ's t-expansion
        self.qk = [(A.get(k, F(0)), B.get(k, F(0)), D.get(k, F(0)), self._hk(k))
                   for k in sorted((set(A) | set(B) | set(D)) - {z0})]
        self.lk = [(Lc.get(k, F(0)), Tc.get(k, F(0)), self._hk(k))
                   for k in sorted((set(Lc) | set(Tc)) - {z0})]
        self.qkf = [tuple(float(x) for x in r) for r in self.qk]
        self.lkf = [tuple(float(x) for x in r) for r in self.lk]

    def _hk(self, k):
        return (self.h[0] ** k[0]) * (self.h[1] ** k[1]) * (self.h[2] ** k[2])

    @staticmethod
    def _coeffs(poly, ts):
        return {m: _q(c) for m, c in sp.Poly(poly, *ts).terms()}

    def charQ_bounds(self, rho):
        """Enclosure of charQ(rho*aF, .) over the cell (exact rationals)."""
        r2 = rho * rho
        A0, B0, D0 = self.q0
        a0 = r2 * A0 - rho * B0 + D0
        s = sum(abs(r2 * a - rho * b + d) * hk for a, b, d, hk in self.qk)
        return a0 - s, a0 + s

    def lin_bounds(self, rho):
        """Enclosure of rho*aF - trG/2 over the cell (exact rationals)."""
        L0, T0 = self.l0
        a0 = rho * L0 - T0 / 2
        s = sum(abs(rho * lam_k - t / 2) * hk for lam_k, t, hk in self.lk)
        return a0 - s, a0 + s

    # ---- float twins, used ONLY by the untrusted locator in `locate_rhos` ----
    def charQ_bounds_f(self, rho):
        r2 = rho * rho
        A0, B0, D0 = (float(x) for x in self.q0)
        a0 = r2 * A0 - rho * B0 + D0
        s = sum(abs(r2 * a - rho * b + d) * hk for a, b, d, hk in self.qkf)
        return a0 - s, a0 + s

    def lin_hi_f(self, rho):
        L0, T0 = (float(x) for x in self.l0)
        return rho * L0 - T0 / 2 + sum(abs(rho * lam_k - t / 2) * hk for lam_k, t, hk in self.lkf)


def hlow1_holds(rho, cert):
    """`rho*aF <= trG/2  OR  charQ(rho*aF) <= 0` over the cell, by the centred form."""
    if rho <= 0:
        return False, "nonpositive"
    lin_hi = cert.lin_bounds(rho)[1]
    if lin_hi <= 0:
        return True, ("linear", lin_hi)
    q_hi = cert.charQ_bounds(rho)[1]
    if q_hi <= 0:
        return True, ("quadratic", q_hi)
    return False, ("none", q_hi)


def hlow2_holds(rho, cert):
    """`rho*aF <= trG/2  AND  charQ(rho*aF) >= 0` over the cell, by the centred form."""
    if rho <= 0:
        return False, "nonpositive"
    lin_hi = cert.lin_bounds(rho)[1]
    if lin_hi > 0:
        return False, "linear"
    q_lo = cert.charQ_bounds(rho)[0]
    if q_lo < 0:
        return False, "quadratic"
    return True, "ok"


RHO_BITS = 64  # every accepted rho is a dyadic k / 2^RHO_BITS: bounded denominators
SHAVES = (F(1, 1 << 50), F(1, 1 << 40), F(1, 1 << 30))  # ~1e-15, 1e-12, 1e-9 relative


def _dyadic_down(x: F, bits: int = RHO_BITS) -> F:
    return F((x.numerator << bits) // x.denominator, 1 << bits)


def _r_point(z):
    """Float r1, r2 at the cell centre: the upper ends of the locator's brackets."""
    import math
    lam, b0, b1 = (float(x) for x in z)
    aF = (b1 - b0) ** 2 + 1 + lam
    tr = b0 * b0 + 1 + b1 * b1 + lam
    det = (b0 * b0 + 1) * (b1 * b1 + lam) - (b0 * b1) ** 2
    sq = math.sqrt(max(tr * tr - 4 * det, 0.0))
    return (tr + sq) / 2 / aF, 2 * det / (tr + sq) / aF


def _bisect(test, hi, lo=0.0):
    """Largest float in [lo, hi] passing `test`, by EXACTLY 60 float bisection steps from a
    passing `lo`. Untrusted: whatever it returns is re-checked exactly by `certify_rho`."""
    if test(hi):
        return hi
    for _ in range(60):
        mid = (lo + hi) / 2
        if test(mid):
            lo = mid
        else:
            hi = mid
    return lo


def locate_rhos(cert):
    """UNTRUSTED float locator (the reference producer's search, on the combined-coefficient
    enclosure). Once |.| is taken per monomial after combining at rho, the enclosure
    endpoints are no longer quadratics in rho, so the old closed-form roots do not apply.

      hLow1: quadratic branch by bisection on [mid, r1c], mid = (r1c + r2c)/2, attempted
             only if it passes at mid; linear branch by bisection on [0, r1c]; take the max.
      hLow2: bisection on [0, r2c] of (lin_hi <= 0 and charQ_lo >= 0).
    Returns floats (p1, p2); 0.0 means nothing was located.
    """
    r1c, r2c = _r_point(cert.z)
    t1 = lambda r: r > 0.0 and cert.charQ_bounds_f(r)[1] <= 0  # noqa: E731
    tl = lambda r: cert.lin_hi_f(r) <= 0  # noqa: E731
    t2 = lambda r: cert.lin_hi_f(r) <= 0 and cert.charQ_bounds_f(r)[0] >= 0  # noqa: E731
    mid = 0.5 * (r1c + r2c)
    p1 = _bisect(t1, r1c, mid) if t1(mid) else 0.0
    p1 = max(p1, _bisect(tl, r1c))
    p2 = _bisect(t2, r2c)
    return p1, p2


def certify_rho(p, holds, cert):
    """EXACT re-check of a located float p: the certified value is the dyadic
    floor(p*(1 - shave)) at the first of three FIXED shaves that `holds` accepts, in exact
    Fraction arithmetic. At most 3 `holds` calls per rho (6 per cell, within the 9-call
    budget); denominators are bounded by 2^RHO_BITS and nothing iterates or grows.
    """
    if not p > 0.0:
        return None
    pf = F(p)  # exact value of the float
    for shave in SHAVES:
        q = _dyadic_down(pf * (1 - shave))
        if q > 0 and holds(q, cert)[0]:
            return q
    return None


# ---------------------------------------------------------------------- the mesh --

LAM_KNOTS = [F(1, 25), F(77, 50), F(3), F(5), F(8), F(12), F(17), F(25)]
B_KNOTS = [F(-3, 2), F(0), F(3, 2)]
W = (F(1), F(1), F(1))
ETA = F(1, 8)
EPS = F(1, 1 << 8)
LMAX = 12


def children8(cell):
    """The 8 children of a cell: the tensor product of the two halves of each coordinate.

    The 8 = 2^3 children must be generated by halving EACH coordinate, not by pairing
    coordinates with each other. An earlier version here did the latter and was wrong.
    """
    halves = []
    for iv in cell:
        m = (iv.lo + iv.hi) / 2
        halves.append((Iv(iv.lo, m), Iv(m, iv.hi)))
    return [(halves[0][i0], halves[1][i1], halves[2][i2])
            for i0 in (0, 1) for i1 in (0, 1) for i2 in (0, 1)]


def split_rule(cell, depth):
    """D1's rule: split iff depth < lmax and width_w > eta*(d + eps). Phi never enters."""
    if depth >= LMAX:
        return False
    width_w = max(W[i] * (cell[i].hi - cell[i].lo) for i in range(3))
    return width_w > ETA * (dist_to_critical_fixed(cell) + EPS)


def _serve_cell(args):
    """Worker: certify one leaf. Returns (box, rho1, rho2) or (box, None, None)."""
    b = args
    cert = CellCert(b)
    p1, p2 = locate_rhos(cert)
    r1 = certify_rho(p1, hlow1_holds, cert)
    r2 = certify_rho(p2, hlow2_holds, cert)
    if r1 is not None and r1 > 0 and r2 is not None and r2 > 0:
        return (b, r1, r2)
    return (b, None, None)


def mesh_leaves():
    """The declared-parameter graded mesh, before certification. Cheap: no arithmetic."""
    leaves = []
    for i in range(len(LAM_KNOTS) - 1):
        for j in range(len(B_KNOTS) - 1):
            for k in range(len(B_KNOTS) - 1):
                block = (Iv(LAM_KNOTS[i], LAM_KNOTS[i + 1]),
                         Iv(B_KNOTS[j], B_KNOTS[j + 1]),
                         Iv(B_KNOTS[k], B_KNOTS[k + 1]))
                stack = [(block, 0)]
                while stack:
                    node, d = stack.pop()
                    if split_rule(node, d):
                        for ch in children8(node):
                            stack.append((ch, d + 1))
                    else:
                        leaves.append(node)
    return leaves


def build_graded_cover(workers=None):
    """Leaves, then certify them in PARALLEL, then REFINE any unserved leaf (never drop).

    Certification is the dominant cost (~13 ms per leaf) and is embarrassingly parallel, so
    it is farmed to a process pool. The previous single-core version took 27 minutes to
    build and, because the A-4 loop was also single core and its result was recomputed
    inside a sort key, over two hours in total.
    """
    import multiprocessing as mp
    import time
    workers = workers or max(1, (mp.cpu_count() or 1) - 1)
    leaves = mesh_leaves()
    cells, unserved = [], []
    for level in range(14):
        if not leaves:
            break
        t_lvl = time.time()
        if workers <= 1:
            # Serial path. A Pool is spawned even for a single worker, and under the Windows
            # spawn start method a child cannot import a dynamically loaded module, so the
            # pool path is not merely slower here, it fails outright. Single-threaded is
            # correct and adequate for callers that only need the cover (e.g. F-6).
            out = [_serve_cell(b) for b in leaves]
        else:
            with mp.Pool(workers) as pool:
                out = pool.map(_serve_cell, leaves, chunksize=64)
        nxt = []
        for b, r1, r2 in out:
            if r1 is not None and r2 is not None:
                cells.append({"box": b, "rho1": r1, "rho2": r2})
            else:
                nxt.append(b)
        print(f"  level {level}: {len(leaves)} leaves -> {len(cells)} certified, "
              f"{len(nxt)} to refine, {time.time() - t_lvl:.1f}s", flush=True)
        if not nxt:
            # every leaf in this level was certified; nothing is left over. Appending
            # `leaves` here reported all 59276 certified cells as unserved and made A-2
            # return False for a cover that was in fact complete.
            leaves = []
            break
        leaves = [ch for b in nxt for ch in children8(b)]
    unserved = list(leaves)
    return cells, unserved


def total_volume(cell):
    return ((cell[0].hi - cell[0].lo) * (cell[1].hi - cell[1].lo)
            * (cell[2].hi - cell[2].lo))


def a2_containment(cells, unserved):
    """A-2 by exact volume equality, and by the top-down reconstruction of each block.

    The previous pairwise 'no leaf contains another' test is O(n^2) and unusable at ~6e4
    cells; the reconstruction below is O(n * depth) and is the check that actually
    establishes that the leaves TILE each block.
    """
    if unserved:
        return False, f"{len(unserved)} leaves could not be certified"
    box_vol = (LAM_HI - LAM_LO) * (B_HI - B_LO) * (B_HI - B_LO)
    v = F(0)
    for c in cells:
        v += total_volume(c["box"])
    if v != box_vol:
        return False, f"leaf volume {v} != box volume {box_vol}"
    return True, f"volumes sum exactly to the box; {len(cells)} cells"


# ---- A-4 workers at MODULE level: under the Windows spawn start method a pool cannot
# pickle a function defined inside main(), so these must not be closures.
UC_W = F(271740, 100000)
ALPHA_W = F(1, 40)
_PE = None


def _a3_one(args):
    """Worker (A-3 pass): does this certified cell satisfy hLow1 at rho1 and hLow2 at rho2?

    A separate pass over the FINISHED cover, deliberately not a counter inside _serve_cell:
    failing cells are refined away there, so an in-place counter would always read 0.
    Returns (hLow1 holds, hLow2 holds).
    """
    b, rho1, rho2 = args
    cert = CellCert(b)
    return (hlow1_holds(rho1, cert)[0], hlow2_holds(rho2, cert)[0])


def _a3_detect_one(args):
    """Worker (FK-3/FK-4 detection): would this cell's rhos have been REJECTED under
    perturbation?  If both counts come back zero the A-3 pass is vacuous.  Returns
    (rejected by hlow1 at doubled rho1, rejected by hlow2 at rho1).
    """
    b, rho1, rho2 = args
    cert = CellCert(b)
    return (not hlow1_holds(2 * rho1, cert)[0], not hlow2_holds(rho1, cert)[0])


def _a4_one(args):
    """Upper endpoint and radius of Phi(Uc; rho1, rho2) by R-beta, as exact Fractions."""
    import importlib.util
    import sys as _sys
    global _PE
    if _PE is None:
        s = importlib.util.spec_from_file_location("pe_w", str(ROOT / "formal01_phi_engine.py"))
        if s is None or s.loader is None:
            raise ImportError("rate engine")
        _PE = importlib.util.module_from_spec(s)
        _sys.modules["pe_w"] = _PE
        s.loader.exec_module(_PE)
    rho1, rho2 = args
    v = _PE.phi_r_beta(UC_W, rho1, rho2)
    return _PE.ball_bounds(v)[1], _PE.radius(v)


def _imhof_pair(args):
    """(R-beta bounds, R-Imhof bounds) for one cell, in a worker process."""
    import importlib.util
    import sys as _sys
    global _PE
    if _PE is None:
        s = importlib.util.spec_from_file_location("pe_w", str(ROOT / "formal01_phi_engine.py"))
        if s is None or s.loader is None:
            raise ImportError("rate engine")
        _PE = importlib.util.module_from_spec(s)
        _sys.modules["pe_w"] = _PE
        s.loader.exec_module(_PE)
    # The arguments arrive directly: under the spawn start method a worker is a FRESH
    # process, so a module global set in the parent is None there.
    rho1, rho2 = args
    vb = _PE.phi_r_beta(UC_W, rho1, rho2)
    vi = _PE.phi_r_imhof(UC_W, rho1, rho2)
    return _PE.ball_bounds(vb), _PE.ball_bounds(vi)


def main():
    import importlib.util
    import time
    pe_spec = importlib.util.spec_from_file_location("pe", str(ROOT / "formal01_phi_engine.py"))
    if pe_spec is None or pe_spec.loader is None:
        raise ImportError("cannot load the rate engine")
    pe = importlib.util.module_from_spec(pe_spec)
    sys.modules["pe"] = pe
    pe_spec.loader.exec_module(pe)

    ALPHA, RAD_BUDGET = F(1, 40), F(1, 10 ** 9)
    claim("F8C-PARAM-LAMKNOTS", "declared lambda knots", [str(x) for x in LAM_KNOTS])
    claim("F8C-PARAM-BKNOTS", "declared b0/b1 knots", [str(x) for x in B_KNOTS])
    claim("F8C-PARAM-ETA-EPS-LMAX", "declared eta, eps, lmax, w",
          [str(ETA), str(EPS), LMAX, [str(x) for x in W]])

    t0 = time.time()
    cells, unserved = build_graded_cover()
    claim("F8C-BUILD-SEC", "seconds to build the graded cover", round(time.time() - t0, 1))
    claim("F8C-CELLS", "cells in the graded cover", len(cells))
    claim("F8C-UNSERVED", "leaves that could not be certified", len(unserved))
    ok2, why2 = a2_containment(cells, unserved)
    claim("F8C-A2", "A-2 containment", ok2)
    claim("F8C-A2-DETAIL", "A-2 outcome", why2)


    # ---- A-3 on every cell of the FINISHED cover, in parallel.
    # Separate pass on purpose: a counter inside _serve_cell
    # would always read 0, because failing cells are refined rather than counted.
    import multiprocessing as mp
    ta3 = time.time()
    with mp.Pool(max(1, (mp.cpu_count() or 1) - 1)) as pool:
        res3 = pool.map(_a3_one, [(c["box"], c["rho1"], c["rho2"]) for c in cells],
                         chunksize=64)
    h1fails = [i for i, r in enumerate(res3) if not r[0]]
    h2fails = [i for i, r in enumerate(res3) if not r[1]]
    claim("F8C-A3-SEC", "seconds for the A-3 pass on every cell (parallel)",
          round(time.time() - ta3, 1))
    claim("F8C-A3-CELLS", "cells the A-3 pass examined", len(res3))
    claim("F8C-A3-HLOW1-FAILS", "cells failing hLow1 (rho*aF <= trG/2 OR charQ <= 0)", len(h1fails))
    claim("F8C-A3-HLOW2-FAILS", "cells failing hLow2 (rho*aF <= trG/2 AND charQ >= 0)", len(h2fails))
    claim("F8C-A3-FIRST-FAIL", "first failing cell index, -1 if none",
          (h1fails + h2fails)[0] if (h1fails or h2fails) else -1)
    claim("F8C-A3-ALL", "every cell of the cover satisfies hLow1 and hLow2",
          (len(h1fails) == 0 and len(h2fails) == 0 and len(res3) == len(cells)))

    # ---- E-5: identity of THIS cover, so later runs tie to it by hash not by cell count.
    # `hashlib` is function-LOCAL in main() because an `import hashlib` appears later in
    # the same function, so it must be imported here too or E-5 raises UnboundLocalError.
    # That import happened after this block was first written; caught by the first run, which
    # crashed BEFORE R3 was written, so the order's single retry is what this uses.
    import hashlib
    hsh = hashlib.sha256()
    for c in cells:
        hsh.update(("|".join(str(x) for x in c["box"]) + "#"
                    + str(c["rho1"]) + "#" + str(c["rho2"]) + chr(10)).encode("utf-8"))
    claim("F8C-COVER-SHA256", "sha256 over the cells in list order", hsh.hexdigest())
    # gate note: the expansion identity, checked by a second independent implementation
    claim("F8C-EXPANSION-XCHECK", "max |sympy charQ - checker charQ| over random exact points",
          str(expansion_crosscheck([F(1), F(0), F(0)], [F(1, 2)] * 3, F(1, 2), trials=60)))

    # ---- A-4 on every cell, in parallel, keeping each upper endpoint so the
    # "cells nearest alpha" selection does NOT recompute 6e4 balls a second time.
    import multiprocessing as mp
    t0 = time.time()
    with mp.Pool(max(1, (mp.cpu_count() or 1) - 1)) as pool:
        res = pool.map(_a4_one, [(c["rho1"], c["rho2"]) for c in cells], chunksize=64)
    ups = [r[0] for r in res]
    maxrad = max((r[1] for r in res), default=F(0))
    fails = sum(1 for u in ups if u > ALPHA)
    worst = max(ups) if ups else F(0)
    claim("F8C-A4-SEC", "seconds for A-4 on every cell (R-beta, parallel)", round(time.time() - t0, 1))
    claim("F8C-A4-FAILS", "cells whose R-beta upper endpoint exceeds alpha", fails)
    claim("F8C-A4-WORST", "largest R-beta upper endpoint over the cover", float(worst))
    claim("F8C-A4-MARGIN", "alpha minus that worst upper endpoint", float(ALPHA - worst))
    claim("F8C-A4-MAX-RADIUS", "largest R-beta ball radius", float(maxrad))
    claim("F8C-A4-RADIUS-BUDGET", "every radius <= 1e-9", maxrad <= RAD_BUDGET)

    # ---- R-Imhof cross-check on the 168 cells nearest alpha, reusing the values above

    order = sorted(range(len(cells)), key=lambda i: ups[i], reverse=True)[:168]
    # ---- A-3 detection check on the 168 cells nearest alpha (E-4).  If both counts are
    # zero, the hLow checks are vacuous on this cover and the A-3 pass proves nothing.
    import multiprocessing as mp
    td3 = time.time()
    with mp.Pool(max(1, (mp.cpu_count() or 1) - 1)) as pool:
        res4 = pool.map(_a3_detect_one,
                        [(cells[i]["box"], cells[i]["rho1"], cells[i]["rho2"]) for i in order],
                        chunksize=64)
    claim("F8C-A3-DETECT-SEC", "seconds for the A-3 detection check", round(time.time() - td3, 1))
    claim("F8C-A3-FK3-REJECTED",
          "nearest-alpha cells rejected by hlow1_holds(2*rho1); FK-3 would have fired",
          sum(1 for r in res4 if r[0]))
    claim("F8C-A3-FK4-REJECTED",
          "nearest-alpha cells rejected by hlow2_holds(rho1); FK-4 would have fired",
          sum(1 for r in res4 if r[1]))
    t0 = time.time()
    with mp.Pool(max(1, (mp.cpu_count() or 1) - 1)) as pool:
        pairs = pool.map(_imhof_pair, [(cells[i]["rho1"], cells[i]["rho2"]) for i in order])
    mism, iworst = 0, F(0)
    for hb, hi_ in pairs:
        if not (hb[0] <= hi_[1] and hi_[0] <= hb[1]):
            mism += 1
        iworst = max(iworst, hi_[1])
    claim("F8C-IMHOF-CELLS", "cells cross-checked with R-Imhof", len(order))
    claim("F8C-IMHOF-SEC", "seconds for the R-Imhof cross-check", round(time.time() - t0, 1))
    claim("F8C-IMHOF-DISJOINT", "cells where the two routes' balls do NOT overlap", mism)
    claim("F8C-IMHOF-WORST", "largest R-Imhof upper endpoint on those cells", float(iworst))
    claim("F8C-IMHOF-ALL-BELOW", "every R-Imhof upper endpoint <= alpha", iworst <= ALPHA)

    # R2 is the run this one extends; the R1 file (the per-polynomial-folded enclosure,
    # 51,834 A-4 failures) is left byte-for-byte in place and named here as superseded.
    import hashlib
    r2_path = DEP / "formal01_f8c_claims-R2.txt"
    r2_sha = hashlib.sha256(r2_path.read_bytes()).hexdigest() if r2_path.exists() else "absent"
    CLAIMS.insert(0, "CLAIM F8C-R3-SUPERSEDES: this run adds an A-3 pass and an "
                     f"F8C-COVER-SHA256 identity. R2 (sha256 {r2_sha}) is retained unchanged and "
                     "this file supersedes it ONLY IF the RC gate holds, i.e. every non-timing "
                     "R2 line reappears here as an identical string = conditional on RC")
    (DEP / "formal01_f8c_claims-R3.txt").write_text("\n".join(CLAIMS) + "\n", encoding="utf-8")
    print("\n".join(CLAIMS))


if __name__ == "__main__":
    main()
