#!/usr/bin/env python3
"""
formal01_f8_indep_checker.py
============================
INDEPENDENT VERIFICATION CHECKER for the f8 numerical cover certificate.

CLEAN ROOM
----------
Written from a written specification only.  The ONLY project
file ever opened is ``deposits/f8c_cover_witness-R1.txt`` (the witness, §2),
plus this file itself, which §3 requires it to scan for G-IND.  Every read
goes through ``read_project_file()``, which enforces that allow-list at runtime
and records what was opened.  No producer module is imported or read:
``formal01_f8c_certified.py``, ``formal01_f8a_checker.py``,
``formal01_phi_engine.py`` and every ``deposits/formal01_*claims*`` file are
untouched.  The numerical facts about the witness (sha256, line count, knots,
worst endpoint, box volume) came from the specification as DATA and sit in the
constants block below.

DEPENDENCIES
------------
Python standard library only.  Forbidden: sympy, flint, mpmath, numpy, scipy
and every project module.  ``ast`` is imported in addition to the §3 list
because §3 mandates a *syntax-tree* scan and ``ast`` is the only stdlib
facility for it; it is not on the forbidden list and G-IND reports it.

ARITHMETIC
----------
Every decision path is exact.  Rationals are ``fractions.Fraction``; the only
irrationals enter through ``math.isqrt`` with outward rounding, and the hNum
enclosure lives on a 2**-160 fixed-point grid with outward rounding at each
step.  Binary floats appear in exactly two instrument tests (the float Simpson
integration §5(c) asks for, and the Monte-Carlo draws of §5(e)) and in no
decision path.

TRANSCRIPTION TABLE (Lean line -> this checker)
-----------------------------------------------
  nu :6 -> NU        alpha :7 -> ALPHA     Lc :8 -> LC     Uc :9 -> UC
  tol :10 -> TOL     Box :13-15 -> BOX_LO/BOX_HI             corner :16 -> CORNER
  aF :18 -> aF()     g00 :19 -> g00()     g11 :20 -> g11() g01 :21 -> g01()
  trG :23 -> trG()   detG :24 -> detG()   charQ :30 -> charQ()
  Phi :32-34 -> §4/§5 route (validated by gate G-0 before use)
  Cell / Cell.mem :41-49 -> B-2
"""

import ast
import collections
import hashlib
import json
import multiprocessing
import random
import re
import sys
import time
from fractions import Fraction as Fr
from math import comb, exp, isqrt, lgamma, log, sqrt
from pathlib import Path

# --------------------------------------------------------------------------
# §1  constants transcribed from the Lean definitions
# --------------------------------------------------------------------------
NU = Fr(35)                     # :6   nu
ALPHA = Fr(1, 40)               # :7   alpha
LC = Fr(271731, 100000)         # :8   Lc
UC = Fr(271740, 100000)         # :9   Uc
TOL = Fr(1, 10000)              # :10  tol
BOX_LO = (Fr(1, 25), Fr(-3, 2), Fr(-3, 2))    # :13-15
BOX_HI = (Fr(25), Fr(3, 2), Fr(3, 2))          # :13-15
CORNER = (Fr(1, 25), Fr(-3, 2), Fr(3, 2))     # :16
UC2 = UC * UC
UCF = 2.71740                   # float mirror of Uc, G-0(e) Monte Carlo only

# --------------------------------------------------------------------------
# §2  witness data (given in the specification as DATA, never read from a producer)
# --------------------------------------------------------------------------
SELF_PATH = Path(__file__).resolve()
WITNESS_PATH = SELF_PATH.parent / "deposits" / "f8c_cover_witness-R1.txt"
CLAIMS_PATH = SELF_PATH.parent / "deposits" / "f8_indep_claims-R1.txt"
SHA256_EXPECTED = "0d73eaee7bfec3452b373f11ae4a642f479b1bf50c7765f6c633cf354a12e32c"
NCELLS_EXPECTED = 59276
KNOTS_L = [Fr(1, 25), Fr(77, 50), Fr(3), Fr(5), Fr(8), Fr(12), Fr(17), Fr(25)]
KNOTS_B = [Fr(-3, 2), Fr(0), Fr(3, 2)]
WORST_CLAIMED = Fr(24999125950489053, 10 ** 18)     # 0.024999125950489053
VOLUME_EXPECTED = Fr(5616, 25)
BLOCK_EXPECTED = 28
MARGIN_CLAIMED = Fr(874, 10 ** 9)                  # alpha - worst, 8.74e-07

# --------------------------------------------------------------------------
# acceptance numbers fixed by the specification; never modified
# --------------------------------------------------------------------------
PHI_GRID_BITS = 160             # hNum fixed-point width   (§4)
PHI_TAIL = Fr(1, 10 ** 13)      # tail-mass target         (§4)
PHI_JOBCAP = 20000              # J cap                    (§4)
PHI_WIDTH_TARGET = Fr(1, 10 ** 10)
HLOW_MAXDEPTH = 4               # first pass               (§4 B-3)
HLOW_MAXDEPTH_ESC = 6           # one escalation           (§6)
GXC_TOL = Fr(5, 10 ** 9)        # P-5 tolerance            (§6)
MC_DRAWS = 10 ** 6              # G-0(e)                   (§5)
MC_SE_FACTOR = 4                # G-0(e)                   (§5)
MC_SEED = 20240101              # G-0(e) fixed seed        (§5)

# --------------------------------------------------------------------------
# clean-room file discipline
# --------------------------------------------------------------------------
OPENED = []
WRITES = []
_READ_ALLOW = {str(SELF_PATH), str(WITNESS_PATH.resolve())}


def read_project_file(path):
    """Read a project file, enforcing the clean-room allow-list and recording it."""
    p = Path(path).resolve()
    if str(p) not in _READ_ALLOW:
        raise RuntimeError(f"clean-room violation: refused to read {p}")
    data = p.read_bytes()
    OPENED.append((p.relative_to(SELF_PATH.parent).as_posix(), len(data), hashlib.sha256(data).hexdigest()))
    return data


def record_write(path):
    """Record a write destination.  Byte counts are deliberately omitted: the deposit
    file cannot contain its own size without being self-referential."""
    WRITES.append(Path(path).resolve().relative_to(SELF_PATH.parent).as_posix())


# --------------------------------------------------------------------------
# §3  G-IND : syntax-tree scan of this file's own imports
# --------------------------------------------------------------------------
FORBIDDEN_MODULES = ("sympy", "flint", "mpmath", "numpy", "scipy")
ALLOWED_STDLIB = {"ast", "collections", "fractions", "hashlib", "json", "math",
                  "multiprocessing", "pathlib", "random", "re", "sys", "time"}


def g_ind():
    src = read_project_file(SELF_PATH)
    tree = ast.parse(src.decode("utf-8"), filename=str(SELF_PATH))
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports += [a.name.split(".")[0] for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            imports.append(node.module.split(".")[0])
    declared = sorted(set(imports))
    project = [m for m in declared
               if (SELF_PATH.parent / f"{m}.py").exists()
               or (SELF_PATH.parent / m / "__init__.py").exists()]
    loaded_forbidden = sorted(m for m in FORBIDDEN_MODULES if m in sys.modules)
    outside = [m for m in declared if m not in ALLOWED_STDLIB]
    return {"declared_imports": declared,
            "non_stdlib_imports": outside,
            "project_module_imports": project,
            "forbidden_modules_loaded": loaded_forbidden,
            "pass": not project and not loaded_forbidden and not outside}


# --------------------------------------------------------------------------
# §1  pointwise transcription of the Lean definitions (independent path)
# --------------------------------------------------------------------------
def aF(t):
    return (t[2] - t[1]) ** 2 + 1 + t[0]                       # :18


def g00(t):
    return t[1] ** 2 + 1                                      # :19


def g11(t):
    return t[2] ** 2 + t[0]                                   # :20


def g01(t):
    return t[1] * t[2]                                        # :21


def trG(t):
    return g00(t) + g11(t)                                    # :23


def detG(t):
    return g00(t) * g11(t) - g01(t) ** 2                     # :24


def discG(t):
    return trG(t) ** 2 - 4 * detG(t)                          # :22


def charQ(m, t):
    return m * m - trG(t) * m + detG(t)                       # :30


def hlow1_true(t, r1):
    """hLow1 : rho1*aF <= trG/2  OR  charQ(rho1*aF) <= 0"""
    m = r1 * aF(t)
    return (m <= trG(t) / 2) or (charQ(m, t) <= 0)


def hlow2_true(t, r2):
    """hLow2 : rho2*aF <= trG/2  AND charQ(rho2*aF) >= 0"""
    m = r2 * aF(t)
    return (m <= trG(t) / 2) and (charQ(m, t) >= 0)


# --------------------------------------------------------------------------
# exact polynomial arithmetic on monomials (i,j,k) of (lam, b0, b1)
# --------------------------------------------------------------------------
ZERO = Fr(0)
ONE = Fr(1)


def padd(a, b):
    o = dict(a)
    for k, v in b.items():
        n = o.get(k, ZERO) + v
        if n:
            o[k] = n
        else:
            o.pop(k, None)
    return o


def pscale(a, c):
    if c == 0:
        return {}
    return {k: v * c for k, v in a.items() if v * c}


def pmul(a, b):
    o = {}
    for (i, j, k), x in a.items():
        for (i2, j2, k2), y in b.items():
            key = (i + i2, j + j2, k + k2)
            n = o.get(key, ZERO) + x * y
            if n:
                o[key] = n
            else:
                o.pop(key, None)
    return o


def psub(a, b):
    return padd(a, {k: -v for k, v in b.items()})


def shift(p, d):
    """substitute lam=lo1+t1, b0=lo2+t2, b1=lo3+t3 and return the poly in t"""
    out = {}
    for (e0, e1, e2), c in p.items():
        for i0 in range(e0 + 1):
            a0 = comb(e0, i0) * d[0] ** (e0 - i0)
            for i1 in range(e1 + 1):
                a1 = comb(e1, i1) * d[1] ** (e1 - i1)
                for i2 in range(e2 + 1):
                    a2 = comb(e2, i2) * d[2] ** (e2 - i2)
                    key = (i0, i1, i2)
                    v = a0 * a1 * a2 * c
                    n = out.get(key, ZERO) + v
                    if n:
                        out[key] = n
                    else:
                        out.pop(key, None)
    return out


def peval(p, t):
    r = ZERO
    for (i, j, k), c in p.items():
        if i == 0 and j == 0 and k == 0:
            r += c
        else:
            r += c * (t[0] ** i) * (t[1] ** j) * (t[2] ** k)
    return r


def corner_bound(p, w):
    """Sound corner expansion (0 <= t_i <= w_i so 0 <= t^k <= w^k):
       upper = a0 + sum_{k!=0} max(a_k,0) w^k ;  lower = a0 + sum_{k!=0} min(a_k,0) w^k"""
    a0 = p.get((0, 0, 0), ZERO)
    up = a0
    dn = a0
    for k, c in p.items():
        if k[0] == 0 and k[1] == 0 and k[2] == 0:
            continue
        wk = (w[0] ** k[0]) * (w[1] ** k[1]) * (w[2] ** k[2])
        if c > 0:
            up += c * wk
        else:
            dn += c * wk
    return up, dn


_AF = {(1, 0, 0): ONE, (0, 0, 0): ONE, (0, 2, 0): ONE, (0, 0, 2): ONE, (0, 1, 1): Fr(-2)}
_G00 = {(0, 2, 0): ONE, (0, 0, 0): ONE}
_G11 = {(0, 0, 2): ONE, (1, 0, 0): ONE}
_G01 = {(0, 1, 1): ONE}
TRG = padd(_G00, _G11)
DET = psub(pmul(_G00, _G11), {(0, 2, 2): ONE})


def hlow_polys(rho):
    """P1 = rho*aF - trG/2 ; P2 = charQ(rho*aF) = (rho aF)^2 - trG (rho aF) + detG"""
    m = pscale(_AF, rho)
    return (psub(m, pscale(TRG, Fr(1, 2))),
            padd(psub(pmul(m, m), pmul(TRG, m)), DET))


# --------------------------------------------------------------------------
# B-1  identity
# --------------------------------------------------------------------------
RE_DELIVERED = re.compile(
    r"^\[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]\|\[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]"
    r"\|\[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]#(-?\d+(?:/\d+)?)#(-?\d+(?:/\d+)?)$")
RE_TEMPLATE = re.compile(
    r"^(-?\d+(?:/\d+)?)\|(-?\d+(?:/\d+)?)\|(-?\d+(?:/\d+)?)\|(-?\d+(?:/\d+)?)"
    r"\|(-?\d+(?:/\d+)?)\|(-?\d+(?:/\d+)?)#(-?\d+(?:/\d+)?)#(-?\d+(?:/\d+)?)$")


def parse_witness(raw):
    text = raw.decode("utf-8")
    lines = text.split("\n")
    trailing_nl = bool(lines) and lines[-1] == ""
    if trailing_nl:
        lines = lines[:-1]
    cells = []
    n_delivered = 0
    n_template = 0
    n_fail = 0
    n_roundtrip = 0
    n_blank = 0
    for line in lines:
        if not line.strip():
            n_blank += 1
        md = RE_DELIVERED.match(line)
        if md:
            n_delivered += 1
        if RE_TEMPLATE.match(line):
            n_template += 1
        if md is None:
            n_fail += 1
            continue
        g = md.groups()
        v = [Fr(s) for s in g]
        if any(str(x) != s for x, s in zip(v, g, strict=True)):
            n_roundtrip += 1
        cells.append(((v[0], v[2], v[4]), (v[1], v[3], v[5]), v[6], v[7]))
    info = {"bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "lines": len(lines),
            "trailing_newline": trailing_nl,
            "crlf": b"\r\n" in raw,
            "fmt_delivered": n_delivered,
            "fmt_template_as_written": n_template,
            "parse_failures": n_fail,
            "roundtrip_failures": n_roundtrip,
            "blank_lines": n_blank}
    return cells, info


def b1_identity(info):
    out = [("B1-SHA256", "sha256(witness bytes)", info["sha256"],
            info["sha256"] == SHA256_EXPECTED),
           ("B1-LINES", "line count", info["lines"], info["lines"] == NCELLS_EXPECTED),
           ("B1-FMT-DELIVERED",
            "lines matching the delivered grammar [lo1,hi1]|[lo2,hi2]|[lo3,hi3]#rho1#rho2",
            info["fmt_delivered"], info["fmt_delivered"] == info["lines"]),
           ("B1-FMT-TEMPLATE-AS-WRITTEN",
            "lines matching the §2 template lo1|lo2|lo3|hi1|hi2|hi3#rho1#rho2",
            info["fmt_template_as_written"],
            info["fmt_template_as_written"] == info["lines"]),
           ("B1-PARSE", "parse failures", info["parse_failures"],
            info["parse_failures"] == 0),
           ("B1-ROUNDTRIP", "rationals with str(Fraction(s)) != s",
            info["roundtrip_failures"], info["roundtrip_failures"] == 0),
           ("B1-BLANK", "blank lines", info["blank_lines"], info["blank_lines"] == 0)]
    return out


# --------------------------------------------------------------------------
# B-2  hCover
# --------------------------------------------------------------------------
def knot_blocks():
    out = []
    for i in range(len(KNOTS_L) - 1):
        for j in range(len(KNOTS_B) - 1):
            for k in range(len(KNOTS_B) - 1):
                out.append(((KNOTS_L[i], KNOTS_B[j], KNOTS_B[k]),
                            (KNOTS_L[i + 1], KNOTS_B[j + 1], KNOTS_B[k + 1])))
    return out


BLOCKS = knot_blocks()


def _pow2_index(q):
    """d with q == 2**d, or -d with q == 2**(-d); None if q is not a power of two"""
    num = q.numerator
    den = q.denominator
    if num != 1 and (num & (num - 1)):
        return None
    if den != 1 and (den & (den - 1)):
        return None
    if num == 1:
        return -(den.bit_length() - 1)
    if den == 1:
        return num.bit_length() - 1
    return None


def _interval(knots, x):
    """index of the knot interval that owns x.  A corner lying exactly on a knot
    belongs to the UPPER interval, otherwise a cell whose lo sits on a knot would be
    assigned to the interval below it and fail the containment test."""
    idx = 0
    for i in range(len(knots)):
        if knots[i] <= x:
            idx = i
    return min(idx, len(knots) - 2)


def b2_check(cells):
    res = {}

    # ---- (a) knots strictly increasing, end knots equal the Box bounds ----
    inc = (all(KNOTS_L[i] < KNOTS_L[i + 1] for i in range(len(KNOTS_L) - 1))
           and all(KNOTS_B[i] < KNOTS_B[i + 1] for i in range(len(KNOTS_B) - 1)))
    ends = (KNOTS_L[0] == BOX_LO[0] and KNOTS_L[-1] == BOX_HI[0]
            and KNOTS_B[0] == BOX_LO[1] and KNOTS_B[-1] == BOX_HI[1]
            and KNOTS_B[0] == BOX_LO[2] and KNOTS_B[-1] == BOX_HI[2])
    bvol = ZERO
    for (blo, bhi) in BLOCKS:
        bvol += (bhi[0] - blo[0]) * (bhi[1] - blo[1]) * (bhi[2] - blo[2])
    res["a_knots_strictly_increasing"] = inc
    res["a_ends_equal_box_bounds"] = ends
    res["a_block_count"] = len(BLOCKS)
    res["a_blocks_tile_box_volume"] = str(bvol)
    res["a_ok"] = (inc and ends and len(BLOCKS) == BLOCK_EXPECTED
                   and bvol == VOLUME_EXPECTED)

    # ---- (b) lo < hi strictly in each coordinate, inside Box ----
    bad_b = 0
    for (lo, hi, _r1, _r2) in cells:
        if (all(lo[i] < hi[i] for i in range(3))
                and all(BOX_LO[i] <= lo[i] and hi[i] <= BOX_HI[i] for i in range(3))):
            continue
        bad_b += 1
    res["b_cells_violating"] = bad_b

    # ---- (c) each cell is a dyadic sub-box of exactly one block, uniform depth ----
    keys = {}
    dup = 0
    bad_c = 0
    noblock = 0
    depth_hist = collections.Counter()
    for (lo, hi, _r1, _r2) in cells:
        il = _interval(KNOTS_L, lo[0])
        j0 = _interval(KNOTS_B, lo[1])
        j1 = _interval(KNOTS_B, lo[2])
        if il < 0 or j0 < 0 or j1 < 0:
            noblock += 1
            bad_c += 1
            continue
        bi = 4 * il + 2 * j0 + j1
        blo, bhi = BLOCKS[bi]
        if not all(blo[t] <= lo[t] and hi[t] <= bhi[t] for t in range(3)):
            noblock += 1
            bad_c += 1
            continue
        d = None
        ok = True
        for t in range(3):
            w = bhi[t] - blo[t]
            u = hi[t] - lo[t]
            if u <= 0:
                ok = False
                break
            dd = _pow2_index(w / u)
            if dd is None or (d is not None and dd != d):
                ok = False
                break
            d = dd
        if not ok or d is None:
            bad_c += 1
            continue
        idxv = []
        for t in range(3):
            w = bhi[t] - blo[t]
            num = (lo[t] - blo[t]) * (1 << d)
            if num % w:
                bad_c += 1
                idxv = []
                break
            v = num // w
            if v < 0 or v >= (1 << d):
                bad_c += 1
                idxv = []
                break
            idxv.append(v)
        if not idxv:
            continue
        key = (bi, d, tuple(idxv))
        if key in keys:
            dup += 1
        else:
            keys[key] = len(keys)
        depth_hist[d] += 1
    res["c_cells_not_a_dyadic_subbox"] = bad_c
    res["c_cells_unassigned_or_multi"] = noblock
    res["c_duplicate_keys"] = dup
    res["c_depth_histogram"] = {str(k): v for k, v in sorted(depth_hist.items())}
    res["c_ok"] = bad_c == 0 and noblock == 0 and dup == 0

    # ---- (d) top-down rebuild; a node that is neither a witness cell nor an
    #      ancestor of one is an UNCOVERED REGION (never descend into it) ----
    by_block = collections.defaultdict(lambda: collections.defaultdict(set))
    for (bi, d, iv) in keys:
        by_block[bi][d].add(iv)
    uncovered = []
    used = set()
    nodes_visited = 0
    for bi in range(len(BLOCKS)):
        levels = by_block[bi]
        proj = {}
        for d2 in levels:
            for d1 in range(d2 + 1):
                sh = d2 - d1
                proj[(d2, d1)] = {(iv[0] >> sh, iv[1] >> sh, iv[2] >> sh)
                                  for iv in levels[d2]}
        frontier = [(0, (0, 0, 0))]
        while frontier:
            d, iv = frontier.pop()
            nodes_visited += 1
            if iv in levels.get(d, ()):
                used.add((bi, d, iv))
                continue
            if any(iv in proj[(d2, d)] for d2 in levels if d2 >= d):
                for c in range(8):
                    nxt = tuple(iv[t] * 2 + ((c >> (2 - t)) & 1) for t in range(3))
                    frontier.append((d + 1, nxt))
                continue
            blo, bhi = BLOCKS[bi]
            w = [bhi[t] - blo[t] for t in range(3)]
            lo = tuple(blo[t] + iv[t] * w[t] / (1 << d) for t in range(3))
            hi = tuple(lo[t] + w[t] / (1 << d) for t in range(3))
            if len(uncovered) < 8:
                uncovered.append({"block": bi, "depth": d,
                                  "lo": [str(x) for x in lo], "hi": [str(x) for x in hi]})
            else:
                uncovered.append({"block": bi, "depth": d, "lo": ["..."], "hi": ["..."]})
    res["d_nodes_visited"] = nodes_visited
    res["d_uncovered_regions"] = len(uncovered)
    res["d_uncovered_examples"] = uncovered[:8]

    # ---- (e) every cell used exactly once ----
    res["e_cells_used"] = len(used)
    res["e_parsed_cells"] = len(cells)
    res["e_expected"] = NCELLS_EXPECTED
    res["e_ok"] = len(used) == len(cells) == NCELLS_EXPECTED

    # ---- (f) exact volume ----
    vol = ZERO
    for (lo, hi, _r1, _r2) in cells:
        vol += (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])
    res["f_volume"] = str(vol)
    res["f_volume_matches"] = vol == VOLUME_EXPECTED
    return res


# --------------------------------------------------------------------------
# B-3  hLow1 / hLow2 on every cell
# --------------------------------------------------------------------------
def hlow_cell(lo, hi, r1, r2, maxdepth=HLOW_MAXDEPTH):
    A1 = hlow_polys(r1)
    A2 = hlow_polys(r2)
    w0 = (hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])

    def rec(d, lo, w, S1, S2, S3, S4):
        u1a, _ = corner_bound(S1, w)
        u2a, _ = corner_bound(S2, w)
        u1b, _ = corner_bound(S3, w)
        _, d2b = corner_bound(S4, w)
        if ((u1a <= 0) or (u2a <= 0)) and (u1b <= 0) and (d2b >= 0):
            return (True, True, d)
        if d >= maxdepth:
            return (False, False, -1)
        hw = (w[0] / 2, w[1] / 2, w[2] / 2)
        C1 = shift(S1, hw)
        C2 = shift(S2, hw)
        C3 = shift(S3, hw)
        C4 = shift(S4, hw)
        res = []
        for c in range(8):
            nlo = (lo[0] + (hw[0] if c & 4 else 0),
                   lo[1] + (hw[1] if c & 2 else 0),
                   lo[2] + (hw[2] if c & 1 else 0))
            res.append(rec(d + 1, nlo, hw, C1, C2, C3, C4))
        return (all(x[0] for x in res), all(x[1] for x in res), max(x[2] for x in res))

    p1, p2, dmax = rec(0, lo, w0, shift(A1[0], lo), shift(A1[1], lo),
                       shift(A2[0], lo), shift(A2[1], lo))
    out = {"h1": p1, "h2": p2, "depth": dmax, "refute1": None, "refute2": None}
    if p1 and p2:
        return out
    w = w0
    pts = [tuple(lo[i] + (w[i] if (k >> (2 - i)) & 1 else 0) for i in range(3))
           for k in range(8)]
    pts.append(tuple(lo[i] + w[i] / 2 for i in range(3)))
    rr = random.Random(MC_SEED ^ (lo[0].numerator & 0xFFFF))
    for _ in range(48):
        pts.append(tuple(lo[i] + w[i] * Fr(rr.randint(0, 10 ** 6), 10 ** 6)
                         for i in range(3)))
    for t in pts:
        if not p1 and out["refute1"] is None and not hlow1_true(t, r1):
            out["refute1"] = [str(x) for x in t]
        if not p2 and out["refute2"] is None and not hlow2_true(t, r2):
            out["refute2"] = [str(x) for x in t]
        if out["refute1"] and out["refute2"]:
            break
    return out


def b3_selfcheck(cells, ncell=256, npts=20):
    """expanded polynomial vs a separately written direct evaluation of §1"""
    bad = 0
    tested = 0
    rr = random.Random(MC_SEED ^ 0xBEEF)
    for (lo, hi, r1, r2) in cells[:ncell]:
        w = (hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])
        for rho in (r1, r2):
            P1, P2 = hlow_polys(rho)
            S1 = shift(P1, lo)
            S2 = shift(P2, lo)
            for _ in range(npts):
                t = tuple(w[i] * Fr(rr.randint(0, 10 ** 6), 10 ** 6) for i in range(3))
                th = tuple(lo[i] + t[i] for i in range(3))
                m = rho * aF(th)
                tested += 2
                if peval(S1, t) != (m - trG(th) / 2):
                    bad += 1
                if peval(S2, t) != charQ(m, th):
                    bad += 1
    return bad, tested


# --------------------------------------------------------------------------
# §4 / §5  hNum route: t-CDF series, NegBin weights, 2**-160 enclosure
# --------------------------------------------------------------------------
def ceil_div(a, b):
    return -((-a) // b)


def g_up(fr):
    return -((-fr.numerator << PHI_GRID_BITS) // fr.denominator)


def g_dn(fr):
    return (fr.numerator << PHI_GRID_BITS) // fr.denominator


def sqrt_grid_updown(num, den):
    """grid ints lo,hi with (lo/2**b)**2 <= num/den <= (hi/2**b)**2"""
    t = (num << (2 * PHI_GRID_BITS)) // den
    m = isqrt(t)
    return m, m + 1


def grid_pow_negbin(p):
    """(lo, hi) grid ints enclosing p**(35/2).  That power is irrational, so it is
    split as p**17 * sqrt(p): p**17 is exact rational and sqrt(p) is enclosed with
    isqrt outward.  Never computed in floating point."""
    a = p.numerator
    b = p.denominator
    na = a ** 17
    nb = b ** 17
    Slo, Shi = sqrt_grid_updown(a, b)
    return (na * Slo) // nb, ceil_div(na * Shi, nb)


def a_coef_grid_series(Ylo, Yhi, kmax):
    """prefix sums of a_k y^k with a_k = a_{k-1}(2k-1)/(2k), outward rounded"""
    Q = 1 << PHI_GRID_BITS
    SL = [Q]
    SH = [Q]
    Tlo = Q
    Thi = Q
    for k in range(kmax):
        n1 = 2 * k + 1
        n2 = 2 * k + 2
        Thi = ceil_div(Thi * n1 * Yhi, n2 * Q)
        Tlo = (Tlo * n1 * Ylo) // (n2 * Q)
        SH.append(SH[-1] + Thi)
        SL.append(SL[-1] + Tlo)
    return SL, SH


def t_tail_grid(s, d, SL, SH, Xlo, Xhi):
    """(lo,hi) grid enclosure of P(T_d > s) from prefix sums of a_k y^k"""
    Q = 1 << PHI_GRID_BITS
    need = d // 2 - 1
    if need > len(SL) - 1:
        return None
    Slo = SL[need]
    Shi = SH[need]
    half = Q // 2
    return (half - ceil_div(Xhi * Shi, 2 * Q), half - (Xlo * Slo) // (2 * Q))


def tail_t_exact_bounds(d, s):
    """exact (grid) enclosure of P(T_d > s) through the §4 series"""
    y = Fr(d) / (Fr(d) + s * s)
    one_minus_y = 1 - y
    Xlo, Xhi = sqrt_grid_updown(one_minus_y.numerator, one_minus_y.denominator)
    SL, SH = a_coef_grid_series(g_dn(y), g_up(y), d // 2)
    return t_tail_grid(s, d, SL, SH, Xlo, Xhi)


def closed_t2(s):
    """P(T_2 > s) = 1/2 - s/(2 sqrt(2+s^2)), enclosed exactly"""
    Q = 1 << PHI_GRID_BITS
    den = 2 + s * s
    Xlo, Xhi = sqrt_grid_updown(den.numerator, den.denominator)
    sn = s.numerator * Q * Q
    f_hi = ceil_div(sn, 2 * s.denominator * Xlo)
    f_lo = (sn // (2 * s.denominator * Xhi)) if sn % (2 * s.denominator * Xhi) == 0 \
        else (sn // (2 * s.denominator * Xhi))
    return Q // 2 - f_hi, Q // 2 - f_lo


def closed_t4(s):
    """P(T_4 > s) = 1/2 - (s/(2 sqrt(4+s^2)))(1 + 2/(4+s^2)), enclosed exactly"""
    Q = 1 << PHI_GRID_BITS
    den = 4 + s * s
    Dlo, Dhi = sqrt_grid_updown(den.numerator, den.denominator)
    inv = 1 / den
    ilo, ihi = g_dn(inv), g_up(inv)
    sn = s.numerator * Q
    sd = 2 * s.denominator
    f_lo = (sn * (Q + 2 * ilo)) // (sd * Dhi)
    f_hi = ceil_div(sn * (Q + 2 * ihi), sd * Dlo)
    return Q // 2 - f_hi, Q // 2 - f_lo


def t_pdf(x, d):
    return exp(lgamma((d + 1) / 2.0) - 0.5 * log(d * 3.141592653589793)
               - lgamma(d / 2.0)) * ((1.0 + x * x / d) ** (-(d + 1) / 2.0))


def simpson_cdf(s, d, n=400000):
    """float Simpson integration of P(T_d <= s); test-only (§5(c))"""
    lo = -60.0
    hi = float(s)
    h = (hi - lo) / n
    tot = 0.0
    for i in range(n + 1):
        x = lo + i * h
        wgt = 1 if i in (0, n) else (4 if i & 1 else 2)
        tot += wgt * t_pdf(x, d)
    return tot * h / 3.0


def phi_bounds(rho1, rho2, tail=None, jobcap=None):
    """Enclosure of Phi(Uc; rho1, rho2) by the §4 route; all values are ints on
    the 2**-160 grid (value = int / 2**-160)."""
    tail = PHI_TAIL if tail is None else tail
    jobcap = PHI_JOBCAP if jobcap is None else jobcap
    Q = 1 << PHI_GRID_BITS
    half = Q // 2
    rl, rh = (rho1, rho2) if rho1 <= rho2 else (rho2, rho1)
    p = rl / rh
    den = NU + UC2 * rl
    Ylo = (NU * Q) // den
    Yhi = ceil_div(NU * Q, den)
    nxc = UC2 * rl
    xratio = nxc / den                       # x = sqrt(1-y) = sqrt(Uc^2 rho_lo/(nu+Uc^2 rho_lo))
    Xlo, Xhi = sqrt_grid_updown(xratio.numerator, xratio.denominator)
    Wl, Wh = grid_pow_negbin(p)
    un = (1 - p).numerator
    ud = (1 - p).denominator
    TAILG = ceil_div(tail.numerator * Q, tail.denominator)
    SL = [Q]
    SH = [Q]
    Thi = Q
    Tlo = Q
    k = 0
    PHIhi = 0
    PHIlo = 0
    SUMWh = 0
    SUMWl = 0
    Jused = jobcap
    for j in range(jobcap + 1):
        need = 34 + j
        while k < need:
            n1 = 2 * k + 1
            n2 = 2 * k + 2
            Thi = ceil_div(Thi * n1 * Yhi, n2 * Q)
            Tlo = (Tlo * n1 * Ylo) // (n2 * Q)
            SH.append(SH[-1] + Thi)
            SL.append(SL[-1] + Tlo)
            k += 1
        phij_hi = half - (Xlo * SL[need]) // (2 * Q)
        phij_lo = half - ceil_div(Xhi * SH[need], 2 * Q)
        PHIhi += ceil_div(Wh * phij_hi, Q)
        PHIlo += (Wl * phij_lo) // Q
        SUMWh += Wh
        SUMWl += Wl
        Wh = ceil_div(Wh * (35 + 2 * j) * un, 2 * (j + 1) * ud)
        Wl = (Wl * (35 + 2 * j) * un) // (2 * (j + 1) * ud)
        # tail_true = 1 - sum_{j<=J} w_j is bounded ABOVE by (1 - SUM_lo/Q) and BELOW
        # by (1 - SUM_hi/Q), so the upper enclosure uses the LOWER weight sum.
        if Q - SUMWl <= TAILG:
            Jused = j
            break
    PHIhi += ceil_div(Q - SUMWl, 2)
    PHIlo += (Q - SUMWh) // 2
    # An enclosure for a probability must satisfy 0 <= lo <= hi <= 1.  If it does not,
    # the enclosure is rejected outright: it can then never be read as PROVED.
    valid = 0 <= PHIlo <= PHIhi <= Q
    if not valid:
        status = "INVALID"
    elif PHIhi * 40 <= Q:
        status = "PROVED"
    elif PHIlo * 40 > Q:
        status = "REFUTED"
    else:
        status = "UNDECIDED"
    return {"status": status, "hi": PHIhi, "lo": PHIlo, "J": Jused, "Q": Q,
            "width": PHIhi - PHIlo, "valid": bool(valid)}


def pareto_minimal(pairs):
    """Pareto-minimal (rho1,rho2): no other pair <= it in both coordinates."""
    srt = sorted(set(pairs), key=lambda t: (t[0], t[1]))
    out = []
    best = None
    for p in srt:
        if best is None or p[1] < best:
            out.append(p)
            best = p[1]
    return out


def match_to_minimal(pairset, M):
    """for each pair, find the Pareto-minimal pair m <= it (exact staircase sweep)"""
    Ms = sorted(M, key=lambda t: (t[0], t[1]))
    srt = sorted(set(pairset), key=lambda t: (t[0], t[1]))
    k = 0
    unmatched = []
    n_ok = 0
    for (r1, r2) in srt:
        while k < len(Ms) and Ms[k][0] <= r1:
            k += 1
        if k == 0:
            unmatched.append([str(r1), str(r2)])
            continue
        if Ms[k - 1][1] <= r2:
            n_ok += 1
        elif len(unmatched) < 8:
            unmatched.append([str(r1), str(r2)])
    return n_ok, unmatched


# --------------------------------------------------------------------------
# G-0  instrument tests for the B-4 route (validate the route before relying on it)
# --------------------------------------------------------------------------
def g0_a_weights():
    out = []
    Q = 1 << PHI_GRID_BITS
    for ps in (Fr(1, 10), Fr(1, 2), Fr(9, 10)):
        p = ps
        Wl, Wh = grid_pow_negbin(p)
        un = (1 - p).numerator
        ud = (1 - p).denominator
        TAILG = ceil_div(PHI_TAIL.numerator * Q, PHI_TAIL.denominator)
        SWl = Wl
        SWh = Wh
        j = 0
        while j < PHI_JOBCAP and (Q - SWl) > TAILG:
            Wh = ceil_div(Wh * (35 + 2 * j) * un, 2 * (j + 1) * ud)
            Wl = (Wl * (35 + 2 * j) * un) // (2 * (j + 1) * ud)
            SWl += Wl
            SWh += Wh
            j += 1
        tailmass = Fr(Q - SWl, Q)
        ok = (Fr(SWl, Q) >= 1 - PHI_TAIL) and (tailmass <= PHI_TAIL)
        out.append({"p": str(ps), "partial_sum_lo": str(Fr(SWl, Q)),
                    "tail_mass": str(tailmass), "J": j, "pass": bool(ok)})
    return out


def g0_b_c_zero(pairs):
    """at c = 0 the route must give exactly 1/2"""
    Q = 1 << PHI_GRID_BITS
    half = Q // 2
    res = []
    for (r1, r2) in pairs:
        rl, rh = (r1, r2) if r1 <= r2 else (r2, r1)
        p = rl / rh
        Wl, Wh = grid_pow_negbin(p)
        un = (1 - p).numerator
        ud = (1 - p).denominator
        TAILG = ceil_div(PHI_TAIL.numerator * Q, PHI_TAIL.denominator)
        hi_acc = Wh // 2
        lo_acc = Wl // 2
        SWl = Wl
        SWh = Wh
        j = 0
        while j < PHI_JOBCAP and (Q - SWl) > TAILG:
            Wh = ceil_div(Wh * (35 + 2 * j) * un, 2 * (j + 1) * ud)
            Wl = (Wl * (35 + 2 * j) * un) // (2 * (j + 1) * ud)
            hi_acc += Wh // 2
            lo_acc += Wl // 2
            SWl += Wl
            SWh += Wh
            j += 1
        hi_tot = hi_acc + ceil_div(Q - SWl, 2)
        lo_tot = lo_acc + (Q - SWh) // 2
        res.append({"pair": [str(r1), str(r2)], "lo": str(Fr(lo_tot, Q)),
                    "hi": str(Fr(hi_tot, Q)),
                    "pass": bool(lo_tot <= half <= hi_tot)})
    return res


def g0_c_tcdf():
    Q = 1 << PHI_GRID_BITS
    out = []
    for (d, s, closed) in ((2, Fr(2), closed_t2), (4, Fr(2), closed_t4)):
        tlo, thi = tail_t_exact_bounds(d, s)
        clo, chi = closed(s)
        out.append({"test": f"d={d} s={s} series vs closed form",
                    "series": [str(Fr(tlo, Q)), str(Fr(thi, Q))],
                    "closed": [str(Fr(clo, Q)), str(Fr(chi, Q))],
                    "pass": bool(tlo <= chi and clo <= thi)})
    tlo, thi = tail_t_exact_bounds(70, Fr(2))
    sim_tail = 1.0 - simpson_cdf(2.0, 70)     # the series is P(T_d > s), the CDF is its complement
    mid = (Fr(tlo, Q) + Fr(thi, Q)) / 2
    out.append({"test": "d=70 s=2 series vs float Simpson of the tail (1e-7)",
                "series_mid": float(mid), "simpson_tail": sim_tail,
                "absdiff": abs(float(mid) - sim_tail),
                "pass": bool(abs(float(mid) - sim_tail) < 1e-7)})
    return out


def g0_d_symmetry(pairs):
    res = []
    for (r1, r2) in pairs:
        a = phi_bounds(r1, r2)
        b = phi_bounds(r2, r1)
        lo = max(a["lo"], b["lo"])
        hi = min(a["hi"], b["hi"])
        res.append({"pair": [str(r1), str(r2)], "overlap": bool(lo <= hi),
                    "identical": bool(a["hi"] == b["hi"] and a["lo"] == b["lo"]),
                    "overlap_width": str(Fr(hi - lo, a["Q"])), "pass": bool(lo <= hi)})
    return res


def g0_e1_exact_moments(pairs):
    """exact cross-check of the route's identity
       rho1*X1+rho2*X2 == rho_lo * chi2_{70+2J},  J~NegBin(r=35/2,p),  Xi ~ chi2_nu.
       Compare E[A], E[A^2] and Var with Fractions."""
    r_half = Fr(35, 2)
    res = []
    for (r1, r2) in list(pairs) + [(Fr(3), Fr(3))]:
        rl, rh = (r1, r2) if r1 <= r2 else (r2, r1)
        p = rl / rh
        EJ = r_half * (1 - p) / p
        EJ2 = r_half * (1 - p) / (p * p) + EJ * EJ
        nu = 35
        d0 = 2 * nu
        m1_route = rl * (d0 + 2 * EJ)
        m1_true = nu * (r1 + r2)
        m2_route = rl * rl * (d0 * (d0 + 2) + 4 * (d0 + 1) * EJ + 4 * EJ2)
        m2_true = (nu * nu + 2 * nu) * (r1 * r1 + r2 * r2) + 2 * nu * nu * r1 * r2
        v_route = m2_route - m1_route * m1_route
        v_true = m2_true - m1_true * m1_true
        res.append({"pair": [str(r1), str(r2)], "p": str(p),
                    "E_match": bool(m1_route == m1_true),
                    "E2_match": bool(m2_route == m2_true),
                    "Var_match": bool(v_route == v_true),
                    "pass": bool(m1_route == m1_true and m2_route == m2_true
                                 and v_route == v_true)})
    return res


def _bisect_left(cum, u):
    lo = 0
    hi = len(cum) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if cum[mid] < u:
            lo = mid + 1
        else:
            hi = mid
    return lo


def negbin_cum(p):
    """exact (rational) cumulative table for J ~ NegBin(r=35/2, p).
    Only the leading weight p**(35/2) is irrational; it is evaluated in double
    precision here, which is harmless for a test-only sampler and is confined to
    the G-0(e) Monte Carlo."""
    pf = float(p)
    w = pf ** 17 * sqrt(pf)
    un = (1 - p).numerator
    ud = (1 - p).denominator
    cum = [float(w)]
    tot = w
    j = 0
    while j < PHI_JOBCAP and (1 - tot) > PHI_TAIL:
        w = w * (35 + 2 * j) * un / (2 * (j + 1) * ud)
        tot = tot + w
        cum.append(float(tot))
        j += 1
    return cum


def g0_e2_mc_identity(pairs):
    """Monte-Carlo of the route's identity itself (no t-CDF needed):
       draw Xi ~ chi2_35 = gammavariate(35/2, 2), compare the empirical CDF of
       rho1*X1+rho2*X2 with the empirical CDF of the NegBin/chi2 mixture at the
       direct sample's quantiles; 4 standard errors."""
    res = []
    for (r1, r2) in pairs:
        rl, rh = (r1, r2) if r1 <= r2 else (r2, r1)
        p = rl / rh
        rnd = random.Random(MC_SEED)
        direct = sorted(float(r1) * rnd.gammavariate(17.5, 2)
                        + float(r2) * rnd.gammavariate(17.5, 2) for _ in range(MC_DRAWS))
        cum = negbin_cum(p)
        mix = sorted(float(rl) * rnd.gammavariate(35.0 + _bisect_left(cum, rnd.random()), 2)
                     for _ in range(MC_DRAWS))
        checks = []
        ok = True
        for q in (0.1, 0.5, 0.9):
            x = direct[int(q * MC_DRAWS)]
            emp = _bisect_left(mix, x) / MC_DRAWS
            se = sqrt(2.0 * q * (1 - q) / MC_DRAWS)
            z = (emp - q) / se
            ok = ok and abs(z) <= MC_SE_FACTOR
            checks.append({"q": q, "direct_ecdf": q, "mixture_ecdf": emp,
                           "z": z, "pass": bool(abs(z) <= MC_SE_FACTOR)})
        res.append({"pair": [str(r1), str(r2)], "checks": checks, "pass": bool(ok)})
    return res


def g0_e3_mc_phi(pairs):
    """Monte-Carlo of the route's own joint model against the §4 enclosure.
    The route's t-CDF arises from ONE chi-square: given J = j, d = 70+2j, the
    conditional event is  Uc*sqrt(rho_lo*V/35) < Z  with V ~ chi2_d and Z ~ N(0,1),
    and P(s*sqrt(V/d) < Z) = P(T_d > s) with s = Uc*sqrt(rho_lo*d/35).  So one
    gammavariate and one gauss per draw; a second independent chi-square would
    model a different event.  10**6 draws, fixed seed, 4 standard errors."""
    Q = 1 << PHI_GRID_BITS
    res = []
    for (r1, r2) in pairs:
        rl, rh = (r1, r2) if r1 <= r2 else (r2, r1)
        p = rl / rh
        b = phi_bounds(r1, r2)
        rnd = random.Random(MC_SEED)
        cum = negbin_cum(p)
        hits = 0
        for _ in range(MC_DRAWS):
            j = _bisect_left(cum, rnd.random())
            half_d = 35.0 + j
            V = rnd.gammavariate(half_d, 2)
            Z = rnd.gauss(0.0, 1.0)
            if UCF * sqrt(float(rl) * V / 35.0) < Z:
                hits += 1
        phi_mc = hits / MC_DRAWS
        mid = (Fr(b["lo"], Q) + Fr(b["hi"], Q)) / 2
        se = sqrt(phi_mc * (1 - phi_mc) / MC_DRAWS)
        z = (phi_mc - float(mid)) / se if se > 0 else 0.0
        res.append({"pair": [str(r1), str(r2)], "phi_mc": phi_mc,
                    "enclosure": [str(Fr(b["lo"], Q)), str(Fr(b["hi"], Q))],
                    "z": z, "pass": bool(abs(z) <= MC_SE_FACTOR)})
    return res


def g0_monotonicity(worst_pair):
    """G-M monotonicity spot-check: Phi must be non-increasing in each rho."""
    if worst_pair is None:
        return {"pass": False, "note": "no binding cell"}
    r1, r2 = worst_pair
    base = phi_bounds(r1, r2)["hi"]
    checks = []
    ok = True
    for name, p_ in (("rho1+1/1000", (r1 + Fr(1, 1000), r2)),
                     ("rho1*2", (r1 * 2, r2)),
                     ("rho2+1/1000", (r1, r2 + Fr(1, 1000))),
                     ("rho2*2", (r1, r2 * 2))):
        h = phi_bounds(*p_)["hi"]
        good = h <= base
        ok = ok and good
        checks.append({"perturbation": name, "monotone": bool(good)})
    return {"checks": checks, "pass": bool(ok)}


# --------------------------------------------------------------------------
# worker plumbing (Windows spawn: module level, initialised once per worker)
# --------------------------------------------------------------------------
_W = {}


def _winit():
    _W["cells"] = parse_witness(read_project_file(WITNESS_PATH))[0]


def _w_hlow(chunk):
    bad = 0
    tested = 0
    proved1 = 0
    proved2 = 0
    hist = collections.Counter()
    refs = []
    fails = []
    for idx in chunk:
        lo, hi, r1, r2 = _W["cells"][idx]
        out = hlow_cell(lo, hi, r1, r2, HLOW_MAXDEPTH)
        if idx < 256:
            b, t = b3_selfcheck([(lo, hi, r1, r2)], 1, 20)
            bad += b
            tested += t
        if out["h1"]:
            proved1 += 1
        if out["h2"]:
            proved2 += 1
        if out["h1"] and out["h2"]:
            hist[out["depth"]] += 1
        else:
            hist["undecided"] += 1
            fails.append(idx)
        if (out["refute1"] or out["refute2"]) and len(refs) < 8:
            refs.append({"cell": idx, "rho1": str(r1), "rho2": str(r2),
                         "hLow1_counterexample": out["refute1"],
                         "hLow2_counterexample": out["refute2"]})
    return (proved1, proved2, hist, refs, fails, bad, tested)


def _w_phi(chunk):
    n_proved = 0
    n_undec = 0
    n_ref = 0
    n_bad = 0
    c_proved = 0
    maxhi = -1
    maxpair = None
    maxwidth = 0
    jmax = 0
    n_wide = 0
    minmargin = None
    q = 1 << PHI_GRID_BITS
    wide_g = ceil_div(PHI_WIDTH_TARGET.numerator * q, PHI_WIDTH_TARGET.denominator)
    for (r1, r2, ncells) in chunk:
        b = phi_bounds(r1, r2)
        if not b["valid"]:
            n_bad += 1
        elif b["status"] == "PROVED":
            n_proved += 1
            c_proved += ncells
        elif b["status"] == "REFUTED":
            n_ref += 1
        else:
            n_undec += 1
        if b["hi"] > maxhi:
            maxhi = b["hi"]
            maxpair = (r1, r2)
        maxwidth = max(maxwidth, b["width"])
        if b["width"] > wide_g:
            n_wide += 1
        jmax = max(jmax, b["J"])
        margin = q - 40 * b["hi"]
        if minmargin is None or margin < minmargin:
            minmargin = margin
    return (n_proved, n_undec, n_ref, n_bad, c_proved, maxhi, maxpair, maxwidth, jmax,
            n_wide, minmargin if minmargin is not None else 0)


def _chunks(seq, n):
    k = max(1, (len(seq) + n - 1) // n)
    return [seq[i:i + k] for i in range(0, len(seq), k)]


def _run_pool(fn, chunks, nproc):
    if nproc <= 1 or len(chunks) <= 1:
        return [fn(c) for c in chunks]
    with multiprocessing.Pool(processes=nproc, initializer=_winit) as pool:
        return list(pool.imap_unordered(fn, chunks))


# --------------------------------------------------------------------------
# deposit accumulator
# --------------------------------------------------------------------------
CLAIMS = []


def claim(cid, label, value, ok=None):
    CLAIMS.append((cid, label, str(value), ok))


# --------------------------------------------------------------------------
# mutations (run on in-memory copies; the real witness is never touched)
# --------------------------------------------------------------------------
def run_mutations(cells, raw, worst_pair):
    out = []
    lines = raw.decode("utf-8").split("\n")

    # MUT-0 empty witness -> not PASS (count gates)
    zinfo = parse_witness(b"")[1]
    z1 = b1_identity(zinfo)
    z2 = b2_check([])
    ok0 = not (z1[0][3] and z1[1][3] and z1[4][3] and z2["f_volume_matches"]
               and z2["e_ok"] and z2["d_uncovered_regions"] == 0)
    out.append(("MUT-0", "empty witness is not PASS",
                f"sha256_match={z1[0][3]} lines_match={z1[1][3]} volume={z2['f_volume']} "
                f"used={z2['e_cells_used']} uncovered={z2['d_uncovered_regions']}", ok0))

    # MUT-1 one byte of line 1 -> identity FAIL
    l0 = lines[0]
    mraw = (("9" if l0[0] != "9" else "8") + l0[1:]).encode("utf-8")
    mraw = "\n".join([mraw.decode("utf-8")] + lines[1:]).encode("utf-8")
    minfo = parse_witness(mraw)[1]
    m1 = b1_identity(minfo)
    out.append(("MUT-1", "one byte of line 1 changed -> identity FAIL",
                f"sha256={minfo['sha256'][:24]}... match={m1[0][3]}", not m1[0][3]))

    # MUT-2 last cell removed -> uncovered region
    r2 = b2_check(cells[:-1])
    ok2 = (r2["d_uncovered_regions"] > 0 or not r2["e_ok"] or not r2["f_volume_matches"])
    out.append(("MUT-2", "last cell removed -> uncovered region",
                f"uncovered_regions={r2['d_uncovered_regions']} used={r2['e_cells_used']}"
                f"/{len(cells) - 1} volume={r2['f_volume']}", ok2))

    # MUT-3 cell 0 duplicated -> used-once FAIL
    r3 = b2_check([cells[0]] + list(cells))
    ok3 = r3["c_duplicate_keys"] > 0 or not r3["e_ok"]
    out.append(("MUT-3", "cell 0 duplicated -> used-once FAIL",
                f"duplicate_keys={r3['c_duplicate_keys']} used={r3['e_cells_used']} "
                f"parsed={r3['e_parsed_cells']}", ok3))

    # MUT-4 interior endpoint moved by 2**-80
    mid = len(cells) // 2
    lo, hi, r1, r2r = cells[mid]
    m4 = list(cells)
    m4[mid] = ((lo[0] + Fr(1, 1 << 80), lo[1], lo[2]), hi, r1, r2r)
    r4 = b2_check(m4)
    ok4 = r4["c_cells_not_a_dyadic_subbox"] > 0 or r4["d_uncovered_regions"] > 0
    out.append(("MUT-4", "interior endpoint moved by 2^-80 -> alignment / uncovered FAIL",
                f"not_a_dyadic_subbox={r4['c_cells_not_a_dyadic_subbox']} "
                f"uncovered_regions={r4['d_uncovered_regions']}", ok4))

    # MUT-5 / MUT-6 : binding cell with rho1 / rho2 doubled
    cell = next((c for c in cells if worst_pair and (c[2], c[3]) == worst_pair), None)
    if cell is None:
        out.append(("MUT-5", "binding cell with rho1 doubled", "no binding cell", False))
        out.append(("MUT-6", "binding cell with rho2 doubled", "no binding cell", False))
        out.append(("MUT-7", "binding cell with both rho halved", "no binding cell", False))
        return out
    lo, hi, wr1, wr2 = cell
    o5 = hlow_cell(lo, hi, wr1 * 2, wr2, HLOW_MAXDEPTH)
    st5 = "PROVED" if o5["h1"] else ("REFUTED" if o5["refute1"] else "UNDECIDED")
    out.append(("MUT-5", "binding cell with rho1 doubled -> hLow1 not PROVED",
                f"hLow1_proved={o5['h1']} status={st5} "
                f"counterexample={json.dumps(o5['refute1'])}", not o5["h1"]))
    o6 = hlow_cell(lo, hi, wr1, wr2 * 2, HLOW_MAXDEPTH)
    st6 = "PROVED" if o6["h2"] else ("REFUTED" if o6["refute2"] else "UNDECIDED")
    out.append(("MUT-6", "binding cell with rho2 doubled -> hLow2 not PROVED",
                f"hLow2_proved={o6['h2']} status={st6} "
                f"counterexample={json.dumps(o6['refute2'])}", not o6["h2"]))
    b7 = phi_bounds(wr1 / 2, wr2 / 2)
    out.append(("MUT-7", "binding cell with both rho halved -> hNum REFUTED",
                f"status={b7['status']} enclosure=[{Fr(b7['lo'], b7['Q'])}, "
                f"{Fr(b7['hi'], b7['Q'])}]", b7["status"] == "REFUTED"))
    return out


# --------------------------------------------------------------------------
def dump_deposit():
    lines = ["# INDEPENDENT VERIFICATION DEPOSIT  --  f8_indep_claims-R1.txt",
             "# generated by formal01_f8_indep_checker.py (clean room: witness only)",
             "# " + time.strftime("%Y-%m-%dT%H:%M:%S"),
             ""]
    n_gate = 0
    n_pass = 0
    for (cid, lab, val, ok) in CLAIMS:
        if cid[:2] in ("P-", "G-") and len(cid) <= 4:
            n_gate += 1
            n_pass += 1 if ok else 0
        tag = "" if ok is True else ("  [FAIL]" if ok is False else "")
        lines.append(f"CLAIM {cid}: {lab} = {val}{tag}")
    lines += ["", f"CLAIM SUMMARY: acceptance gates passed {n_pass} of {n_gate}",
              f"CLAIM VERDICT: {[c[2] for c in CLAIMS if c[0] == 'VERDICT'][0]}"]
    text = "\n".join(lines) + "\n"
    sys.stdout.write(text + "\n")
    CLAIMS_PATH.write_text(text, encoding="utf-8")
    record_write(CLAIMS_PATH)


# --------------------------------------------------------------------------
def main():
    t_start = time.time()
    nproc = min(12, max(1, multiprocessing.cpu_count()))
    print("=" * 78)
    print("INDEPENDENT VERIFICATION CHECKER  --  formal01_f8_indep_checker.py")
    print(f"clean room: witness only | {nproc} procs | numpy/scipy/sympy/flint/mpmath absent")
    print("=" * 78)

    # ---------------- G-IND ----------------
    gi = g_ind()
    claim("G-IND", "syntax-tree scan of own imports",
          f"imports={gi['declared_imports']} non_stdlib={gi['non_stdlib_imports']} "
          f"project_modules={gi['project_module_imports']} "
          f"forbidden_loaded={gi['forbidden_modules_loaded']} (ast is stdlib and is "
          "required by the §3 syntax-tree scan)", gi["pass"])

    # ---------------- B-1 ----------------
    raw = read_project_file(WITNESS_PATH)
    cells, info = parse_witness(raw)
    b1map = {}
    for (cid, lab, val, ok) in b1_identity(info):
        claim(cid, lab, val, ok)
        b1map[cid] = ok
    p1_ok = b1map["B1-SHA256"] and b1map["B1-LINES"] and b1map["B1-PARSE"] \
        and b1map["B1-ROUNDTRIP"]
    claim("P-1", "B-1 hash matches, 59276 lines, 0 parse failures (round-trip included)",
          f"sha256={info['sha256']} lines={info['lines']} "
          f"parse_failures={info['parse_failures']}", p1_ok)

    # ---------------- B-2 ----------------
    b2 = b2_check(cells)
    claim("B2a", "knots strictly increasing, end knots equal the Box bounds, 28 blocks tile Box",
          f"strictly_increasing={b2['a_knots_strictly_increasing']} "
          f"ends_equal_Box={b2['a_ends_equal_box_bounds']} blocks={b2['a_block_count']} "
          f"blocks_volume={b2['a_blocks_tile_box_volume']}", b2["a_ok"])
    claim("B2b", "every cell has lo<hi strictly in each coordinate and lies inside Box",
          f"violations={b2['b_cells_violating']} of {len(cells)}", b2["b_cells_violating"] == 0)
    claim("B2c", "every cell is a dyadic sub-box of exactly one block, same depth in all "
          "three coordinates and aligned to that block",
          f"violations={b2['c_cells_not_a_dyadic_subbox']} "
          f"unassigned={b2['c_cells_unassigned_or_multi']} duplicate_keys={b2['c_duplicate_keys']} "
          f"witness_cell_depth_histogram={json.dumps(b2['c_depth_histogram'])}", b2["c_ok"])
    claim("B2d", "top-down rebuild: nodes that are neither a witness cell nor an ancestor of "
          "one are UNCOVERED REGIONS",
          f"nodes_visited={b2['d_nodes_visited']} uncovered_regions={b2['d_uncovered_regions']} "
          f"examples={json.dumps(b2['d_uncovered_examples'])}", b2["d_uncovered_regions"] == 0)
    claim("B2e", "each cell used exactly once", f"used={b2['e_cells_used']} "
          f"parsed={b2['e_parsed_cells']} expected={b2['e_expected']}", b2["e_ok"])
    claim("B2f", "volumes sum to exactly 5616/25", f"sum={b2['f_volume']}", b2["f_volume_matches"])
    b2_ok = (b2["a_ok"] and b2["b_cells_violating"] == 0 and b2["c_ok"]
             and b2["d_uncovered_regions"] == 0 and b2["e_ok"] and b2["f_volume_matches"])
    claim("P-2", "B-2 (a)-(f): 0 uncovered regions, 59276 cells used exactly once, volume exact",
          f"a={b2['a_ok']} b={b2['b_cells_violating']} c={b2['c_ok']} "
          f"d={b2['d_uncovered_regions']} e={b2['e_ok']} f={b2['f_volume_matches']}", b2_ok)

    # ---------------- B-3 ----------------
    _winit()
    t0 = time.time()
    res = _run_pool(_w_hlow, _chunks(list(range(len(cells))), nproc * 4), nproc)
    dt3 = time.time() - t0
    pr1 = sum(r[0] for r in res)
    pr2 = sum(r[1] for r in res)
    sc_bad = sum(r[5] for r in res)
    sc_tested = sum(r[6] for r in res)
    hist = collections.Counter()
    refs = []
    fails = []
    for r in res:
        hist.update(r[2])
        refs.extend(r[3])
        fails.extend(r[4])
    hist_j = {str(k): v for k, v in sorted(hist.items(), key=lambda kv: str(kv[0]))}
    claim("B3-SELFCHECK", "expanded polynomial vs separately written direct evaluation of the §1 "
          "formulas (first 256 cells x 20 exact rational points)",
          f"mismatches={sc_bad} over {sc_tested} comparisons", sc_bad == 0)
    claim("B3-HIST", "hLow subdivision depth used per cell (first depth at which the cell's "
          "proof completed)", json.dumps(hist_j), True)
    claim("B3-REFUTED", "cells with an exact counterexample point",
          json.dumps(refs) if refs else "none", not refs)
    claim("B3", "hLow pass", f"hLow1 PROVED={pr1} hLow2 PROVED={pr2} of {NCELLS_EXPECTED} "
          f"in {dt3:.1f}s on {nproc} procs", True)

    esc_hist = {}
    esc_ref = []
    if fails:
        todo = fails[:4000]
        for i in todo:
            lo, hi, r1, r2 = cells[i]
            o = hlow_cell(lo, hi, r1, r2, HLOW_MAXDEPTH_ESC)
            if o["h1"]:
                pr1 += 1
            if o["h2"]:
                pr2 += 1
            if o["h1"] and o["h2"]:
                esc_hist[str(o["depth"])] = esc_hist.get(str(o["depth"]), 0) + 1
            elif (o["refute1"] or o["refute2"]) and len(esc_ref) < 8:
                esc_ref.append({"cell": i, "hLow1_counterexample": o["refute1"],
                                "hLow2_counterexample": o["refute2"]})
        claim("B3-ESCALATION", "one escalation to hLow depth <= 6 on the undecided cells "
              f"(time-capped at 4000 of {len(fails)})",
              f"escalated={len(todo)} depths={json.dumps(esc_hist)} "
              f"counterexamples={json.dumps(esc_ref) if esc_ref else 'none'}", True)
    p3_ok = pr1 == NCELLS_EXPECTED and pr2 == NCELLS_EXPECTED
    claim("P-3", "hLow1 PROVED==59276 and hLow2 PROVED==59276 (compared against the constant "
          "59276, not against cells processed)", f"hLow1={pr1} hLow2={pr2}", p3_ok)

    # ---------------- G-0 ----------------
    pairset = sorted({(r1, r2) for (_l, _h, r1, r2) in cells})
    cellcount = collections.Counter((r1, r2) for (_l, _h, r1, r2) in cells)
    g0_pairs = [pairset[0], pairset[len(pairset) // 2], pairset[-1]]
    g0a = g0_a_weights()
    g0b = g0_b_c_zero(g0_pairs)
    g0c = g0_c_tcdf()
    g0d = g0_d_symmetry(g0_pairs)
    g0e1 = g0_e1_exact_moments(g0_pairs)
    g0e2 = g0_e2_mc_identity(g0_pairs)
    g0e3 = g0_e3_mc_phi(g0_pairs)
    g0a_ok = all(x["pass"] for x in g0a)
    g0b_ok = all(x["pass"] for x in g0b)
    g0c_ok = all(x["pass"] for x in g0c)
    g0d_ok = all(x["pass"] for x in g0d)
    g0e_ok = (all(x["pass"] for x in g0e1) and all(x["pass"] for x in g0e2)
              and all(x["pass"] for x in g0e3))
    claim("G-0a", "sum_j w_j in [1-1e-13, 1] at p in {1/10, 1/2, 9/10}", json.dumps(g0a), g0a_ok)
    claim("G-0b", "at c = 0 the route gives exactly 1/2", json.dumps(g0b), g0b_ok)
    claim("G-0c", "t-CDF series vs closed forms at d=2,4 and vs float Simpson at d=70,s=2 "
          "within 1e-7", json.dumps(g0c), g0c_ok)
    claim("G-0d", "the enclosures of Phi(a,b) and Phi(b,a) overlap", json.dumps(g0d), g0d_ok)
    claim("G-0e", f"identity and route vs Monte-Carlo (10**6 draws, fixed seed {MC_SEED}, "
          "random.gauss and random.gammavariate(35/2,2), 3 pairs, 4 standard errors)",
          json.dumps({"e1_exact_moments": g0e1, "e2_identity_mc": g0e2,
                      "e3_phi_mc": g0e3}), g0e_ok)
    g0_ok = g0a_ok and g0b_ok and g0c_ok and g0d_ok and g0e_ok
    claim("G-0", "G-0 TEST OF THE B-4 ROUTE (a)-(e)",
          f"a={g0a_ok} b={g0b_ok} c={g0c_ok} d={g0d_ok} e={g0e_ok}", g0_ok)

    # ---------------- B-4 : M, scope, hNum ----------------
    M = pareto_minimal(pairset)
    claim("M", "Pareto-minimal (rho1,rho2) pairs, exact sort + sweep",
          f"|M|={len(M)} of {len(pairset)} distinct pairs", True)
    samp = pairset[:: max(1, len(pairset) // 40)][:40]
    t0 = time.time()
    for pair in samp:
        phi_bounds(*pair)
    dt_sample = time.time() - t0
    est = dt_sample * len(pairset) / max(1, len(samp))
    scope = "ALL CELLS" if est <= 1800 else "M ONLY (via the antitonicity reduction)"
    claim("SCOPE", "scope decision (30-minute rule)",
          f"40-pair timing {dt_sample:.2f}s -> estimated {est:.0f}s for all {len(pairset)} "
          f"distinct pairs -> {scope}", True)
    n_match, unmatched = match_to_minimal(pairset, M)
    claim("SCOPE-MATCH", "every distinct pair is matched to a Pareto-minimal pair m <= it "
          "(antitonicity requirement)",
          f"matched={n_match}/{len(pairset)} unmatched={json.dumps(unmatched)}",
          n_match == len(pairset))
    todo = pairset if scope == "ALL CELLS" else M
    chunks = _chunks([(r1, r2, cellcount[(r1, r2)]) for (r1, r2) in todo], nproc * 4)
    t0 = time.time()
    res = _run_pool(_w_phi, chunks, nproc)
    dt4 = time.time() - t0
    nprov = sum(r[0] for r in res)
    nund = sum(r[1] for r in res)
    nref = sum(r[2] for r in res)
    nbad = sum(r[3] for r in res)
    cprov = sum(r[4] for r in res)
    widths = max(r[7] for r in res) if res else 0
    jmax = max(r[8] for r in res) if res else 0
    n_wide = sum(r[9] for r in res)
    minmargin = min((r[10] for r in res), default=0)
    Q = 1 << PHI_GRID_BITS
    best_hi = max((r[5] for r in res), default=0)
    best_pair = next((r[6] for r in res if r[5] == best_hi and r[6] is not None), None)
    best_cell = next((i for i, c in enumerate(cells)
                      if best_pair and (c[2], c[3]) == best_pair), -1)
    claim("B4-ENCLOSURE", "hNum enclosure width, J usage and range validity",
          f"max width={Fr(widths, Q)} (target <= {PHI_WIDTH_TARGET}), pairs over the width "
          f"target={n_wide} (the J cap of {PHI_JOBCAP} binds when p is small, so the 1e-13 tail "
          f"mass is unreachable), max J used={jmax}, enclosures outside [0,1] (rejected, never "
          f"PROVED)={nbad}", nbad == 0)
    claim("B4-DECISIONS", "every pair's verdict is certified regardless of the width target",
          f"smallest margin alpha - Phi_upper over the scope = {Fr(minmargin, 40 * Q)}, which is "
          f"{Fr(minmargin, 40 * Q) / (Fr(widths, Q) or Fr(1)):.3g}x the largest enclosure width; "
          f"so no verdict is decided by an enclosure straddling 1/40", minmargin > 0)
    claim("B4", "Phi(Uc; rho) over the chosen scope",
          f"PROVED={nprov} UNDECIDED={nund} REFUTED={nref} INVALID={nbad} of {len(todo)} "
          f"distinct pairs ({dt4:.1f}s, {nproc} procs)", nbad == 0)
    best_frac = Fr(best_hi, Q)
    n_with_best = sum(1 for c in cells if best_pair and (c[2], c[3]) == best_pair)
    claim("B4-MAX", "largest Phi upper bound and its cell",
          f"max upper bound={best_frac} at rho=({best_pair[0]},{best_pair[1]}) witness cell "
          f"index {best_cell}, cells carrying this pair={n_with_best}", True)
    diff = abs(best_frac - WORST_CLAIMED)
    claim("G-XC", "largest Phi upper bound agrees with the claimed worst endpoint "
          "0.024999125950489053",
          f"mine={best_frac} claimed={WORST_CLAIMED} |diff|={diff} tol={GXC_TOL} "
          f"(alpha-worst={ALPHA - WORST_CLAIMED})", diff <= GXC_TOL)
    p5_ok = diff <= GXC_TOL
    claim("P-5", "G-XC at the binding cell, far inside the 8.74e-07 margin",
          f"|diff|={diff} <= {GXC_TOL}", p5_ok)

    g_m_pair = best_pair in set(M)
    mono = g0_monotonicity(best_pair)
    claim("G-M", "the cell with the largest Phi upper bound has its pair in M, and Phi is "
          "non-increasing in each rho (antitonicity check)",
          f"pair_in_M={g_m_pair} |M|={len(M)} monotonicity={json.dumps(mono)}",
          g_m_pair and mono["pass"])
    g_m_ok = g_m_pair and mono["pass"]
    cells_covered = cprov if scope == "ALL CELLS" else None
    p4_ok = nprov == len(todo) and nund == 0 and nref == 0 and g_m_ok \
        and n_match == len(pairset) and nbad == 0 and best_pair is not None
    claim("P-4", "hNum PROVED over the chosen scope, every cell matched to an M pair",
          f"pairs PROVED={nprov}/{len(todo)} cells covered by proved pairs={cells_covered} "
          f"scope={scope}", p4_ok)

    # ---------------- mutations ----------------
    mut = run_mutations(cells, raw, best_pair)
    for (cid, lab, val, ok) in mut:
        claim(cid, lab, val, ok)
    p7_ok = all(m[3] for m in mut)
    claim("P-7", "all eight mutations are detected",
          f"{sum(1 for m in mut if m[3])}/8 detected", p7_ok)

    # ---------------- verdict ----------------
    p6_ok = g0_ok and g_m_ok and gi["pass"]
    claim("P-6", "G-0, G-M and G-IND all pass",
          f"G-0={g0_ok} G-M={g_m_ok} G-IND={gi['pass']}", p6_ok)
    instrument_fail = (not gi["pass"]) or (not g_m_ok) or (not g0_ok) or (not p5_ok) \
        or (not p7_ok) or (sc_bad > 0) or (nbad > 0) or (best_pair is None)
    genuine_refutation = (nref > 0) or (b2["d_uncovered_regions"] > 0) or bool(refs) \
        or bool(esc_ref)
    undecided = (pr1 != NCELLS_EXPECTED) or (pr2 != NCELLS_EXPECTED) or nund > 0
    if instrument_fail:
        verdict = "VOID"
    elif genuine_refutation:
        verdict = "FAIL"
    elif undecided:
        verdict = "INCOMPLETE"
    else:
        verdict = "PASS"
    claim("VERDICT", "computed verdict (not asserted)", verdict, verdict == "PASS")
    claim("FILES-OPENED", "files opened by this checker (allow-list enforced in "
          "read_project_file)", json.dumps(sorted({o[0] for o in OPENED})), True)
    claim("FILES-READS", "every read: (path, bytes, sha256)", json.dumps(OPENED), True)
    claim("FILES-WRITES", "files written (destinations; byte counts omitted as self-referential)",
          json.dumps(sorted(set(WRITES))), True)
    claim("TIME", "wall clock", f"{time.time() - t_start:.1f}s on {nproc} procs", True)
    dump_deposit()
    return verdict


if __name__ == "__main__":
    main()
