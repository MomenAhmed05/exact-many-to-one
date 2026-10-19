"""P3RA — the Unit 2B falsifier battery F-1..F-5, run against the certified graded cover of
record (deposits/formal01_f8c_claims-R2.txt), and P3RB — the conditional FWER bound witness.

Nothing here re-verifies the enclosure. The cover is REBUILT deterministically by
`formal01_f8c_certified.build_graded_cover()` (the cells, rho's and shaves are a pure
function of the declared parameter vector), its cell count is checked against the deposit,
and then it is ATTACKED.

Rigour basis for every Phi value below:
  * r1(theta), r2(theta) at an exact rational theta are enclosed by `r1_iv` / `r2_iv` on a
    point interval: exact rational arithmetic plus an exact-integer-sqrt enclosure of the
    discriminant (2^-R_PREC grid). No float enters.
  * Phi(c; rho) is an Arb ball from the R-beta route (`phi_r_beta`), radius <= 1e-29 in
    practice, every error term included.
  * Phi is non-increasing in c and in each rho (Lean `Phi_antitone`), so
        Phi(c; r(theta)) <= upper(Phi(c; r_lo(theta)))   and
        Phi(c; r(theta)) >= lower(Phi(c; r_hi(theta))).
    An F-k verdict compares those EXACT rational endpoints with alpha. There is no Monte
    Carlo anywhere, so no SE: the uncertainty statement is the ball radius.

Imports: the checker, the cover builder and the rate engine only. No w2_*, routec_* or
AlignedFactor import. Paper 1's `04-efficient-rule/` and an earlier branch-and-bound run
(`wAtoH_r2/`) are READ as text/JSON only.

Not runnable from this repository as published: it also reads the F-3 tolerance
declaration, the `wAtoH_r2/` files and two cached covers (.pkl, 38 MB), which are not
included. Its results file, deposits/formal01_p3a_claims.txt, is included.

Deposit: deposits/formal01_p3a_claims.txt (P-1: every number in the record comes from here).
"""

from __future__ import annotations

import hashlib
import json
import math
import multiprocessing as mp
import pathlib
import pickle
import random
import re
import sys
import time
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import formal01_f8c_certified as f8c  # noqa: E402
import formal01_phi_engine as pe  # noqa: E402
from formal01_f8a_checker import B_HI, B_LO, LAM_HI, LAM_LO, Iv, r1_iv, r2_iv  # noqa: E402

DEP = ROOT / "deposits"
OUT = DEP / "formal01_p3a_claims.txt"
CACHE_COVER = ROOT / "formal01_p3a_cover.pkl"
CACHE_F1 = ROOT / "formal01_p3a_f1.pkl"
CLAIMS: list[str] = []


def need(m, what, src=""):
    """A regex that must have matched. A missing pattern means the source file changed and
    the clause is OPEN; crashing with a bare AttributeError on None would hide that."""
    if m is None:
        raise RuntimeError(f"could not locate {what} in {src or 'the source'}; this clause is OPEN")
    return m


def read_json(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"could not read {path}: {exc}; this clause is OPEN") from exc


def read_text(path):
    try:
        return path.read_text(encoding="utf-8")
    except OSError as exc:
        raise RuntimeError(f"could not read {path}: {exc}; this clause is OPEN") from exc


def claim(cid: str, label: str, value) -> None:
    line = f"CLAIM {cid}: {label} = {value}"
    CLAIMS.append(line)
    print(line, flush=True)


def sha(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------ constants FROM the Lean file --
LEAN = ROOT / "lean" / "W4lean" / "Unit2BFormalR1.lean"
KAPPA_LEAN = ROOT / "lean" / "W4lean" / "Unit2BKappaStatementR1.lean"
_LEAN_TXT = LEAN.read_text(encoding="utf-8")


def lean_const(name: str) -> F:
    m = re.search(rf"^def {name} : ℝ := ([0-9 /]+)$", _LEAN_TXT, re.M)
    if m is None:
        raise ValueError(f"Lean constant {name} not found")
    parts = [p.strip() for p in m.group(1).split("/")]
    return F(int(parts[0]), int(parts[1])) if len(parts) == 2 else F(int(parts[0]))


LC, UC, ALPHA, TOL, NU = (lean_const(n) for n in ("Lc", "Uc", "alpha", "tol", "nu"))
R_PREC = 100            # 2^-100 grid for the exact sqrt enclosure of r(theta)

# ------------------------------------------------------------ declared probe design --
# Fixed before any run. Every fraction is strictly inside (0,1) and none is 1/2, so no
# probe is a cell corner, a face point, or the cell centre (the locator's representative).
PROBES_DIAG = [(F(r + 1, 9),) * 3 for r in range(8)]          # F-1 rule as declared
PROBES_QUART = [(a, b, c) for a in (F(1, 4), F(3, 4)) for b in (F(1, 4), F(3, 4))
                for c in (F(1, 4), F(3, 4))]
PROBES_NEARV = [(a, b, c) for a in (F(1, 64), F(63, 64)) for b in (F(1, 64), F(63, 64))
                for c in (F(1, 64), F(63, 64))]
PROBES = PROBES_DIAG + PROBES_QUART + PROBES_NEARV
F2_SHAVES = (4, 6, 8, 10, 12, 14, 16, 18, 20, 24)               # perturbation 2^-s, exact
F5_TARGET = (F(1, 3), F(2, 7), F(3, 5))                           # never on a dyadic plane
F5_DEPTH = 8
SEED = 20260928
N_TOP, N_RAND = 16, 48


def point_in(box, fr):
    return tuple(box[i].lo + (box[i].hi - box[i].lo) * fr[i] for i in range(3))


def vertices(box):
    return [tuple(box[i].hi if (n >> i) & 1 else box[i].lo for i in range(3)) for n in range(8)]


def r_enc(theta, prec=R_PREC):
    b = tuple(Iv(x) for x in theta)
    return r1_iv(b, prec), r2_iv(b, prec)


def phi_bounds(c, rho1, rho2):
    v = pe.phi_r_beta(c, rho1, rho2)
    lo, hi = pe.ball_bounds(v)
    return lo, hi, pe.radius(v)


# ------------------------------------------------------------------ F-1 worker -------

def _f1_cell(args):
    """Attack one cell: 24 strictly interior probes + 8 vertices (vertices for r only)."""
    idx, box, rho1, rho2 = args
    klo, khi, krad = phi_bounds(UC, rho1, rho2)
    rec = {"idx": idx, "k_lo": klo, "k_hi": khi, "maxrad": krad, "n": 0, "viol": 0,
           "undec": 0, "dom_fail": 0, "dom_undec": 0, "incons": 0, "not_interior": 0,
           "best_hi": None, "best_lo": None, "best_pt": None,
           "min_r1_hi": None, "min_r2_hi": None, "vdom_fail": 0, "vdom_undec": 0}

    def upd_min(key, v):
        if rec[key] is None or v < rec[key]:
            rec[key] = v

    for fr in PROBES:
        pt = point_in(box, fr)
        if not all(box[i].lo < pt[i] < box[i].hi for i in range(3)):
            rec["not_interior"] += 1
        R1, R2 = r_enc(pt)
        for rho, R in ((rho1, R1), (rho2, R2)):
            if rho > R.hi:
                rec["dom_fail"] += 1
            elif rho > R.lo:
                rec["dom_undec"] += 1
        upd_min("min_r1_hi", R1.hi)
        upd_min("min_r2_hi", R2.hi)
        lo, hi, rad = phi_bounds(UC, R1.lo, R2.lo)
        rec["maxrad"] = max(rec["maxrad"], rad)
        rec["n"] += 1
        if hi > ALPHA:
            if lo > ALPHA:
                rec["viol"] += 1
            else:
                rec["undec"] += 1
        if lo > khi:                     # r_lo >= rho would force Phi(r_lo) <= Phi(rho)
            rec["incons"] += 1
        if rec["best_hi"] is None or hi > rec["best_hi"]:
            rec["best_hi"], rec["best_lo"], rec["best_pt"] = hi, lo, pt
    for v in vertices(box):
        R1, R2 = r_enc(v)
        for rho, R in ((rho1, R1), (rho2, R2)):
            if rho > R.hi:
                rec["vdom_fail"] += 1
            elif rho > R.lo:
                rec["vdom_undec"] += 1
        upd_min("min_r1_hi", R1.hi)
        upd_min("min_r2_hi", R2.hi)
    return rec


def _grid_point(args):
    """A box point (G1 / corner / vertex): Phi(U; r) and dominance vs EVERY containing cell."""
    pt, cells_here = args
    R1, R2 = r_enc(pt)
    lo, hi, rad = phi_bounds(UC, R1.lo, R2.lo)
    dfail = dund = 0
    for rho1, rho2 in cells_here:
        for rho, R in ((rho1, R1), (rho2, R2)):
            if rho > R.hi:
                dfail += 1
            elif rho > R.lo:
                dund += 1
    return {"pt": pt, "lo": lo, "hi": hi, "rad": rad, "ncells": len(cells_here),
            "dom_fail": dfail, "dom_undec": dund}


# ------------------------------------------------------------------ F-3 -----------

def cstar_bracket(theta, c0=2.0, half=F(1, 1 << 44)):
    """Rigorous [c_lo, c_hi] containing cstarPt(theta) = inf{c>0 : Phi(c; r(theta)) <= alpha}.

    Accepted only if  upper(Phi(c_hi; r_lo)) <= alpha  (so cstarPt <= c_hi) and
    lower(Phi(c_lo; r_hi)) > alpha  (so Phi(c) > alpha for all c <= c_lo, cstarPt >= c_lo).
    The float secant only LOCATES the root; the two exact comparisons are the certificate.
    """
    R1, R2 = r_enc(theta)
    # r_enc gives an enclosure of r(theta); its lower endpoint can be non-positive because
    # the interval extension is conservative near the corner, where r2 is small. The engine
    # requires strictly positive rho (Real.sqrt is 0 on negatives, so a negative ratio is not
    # a model point at all). Tighten the enclosure about its own value until the lower
    # endpoints are strictly positive, and only then bracket. If they cannot be made
    # positive, that theta is not certifiable here and must be reported, not crashed on.
    for _ in range(80):
        if R1.lo > 0 and R2.lo > 0:
            break
        mid1, mid2 = (R1.lo + R1.hi) / 2, (R2.lo + R2.hi) / 2
        R1, R2 = Iv(R1.lo + (mid1 - R1.lo) / 2, R1.hi - (mid1 - R1.lo) / 2), \
                 Iv(R2.lo + (mid2 - R2.lo) / 2, R2.hi - (mid2 - R2.lo) / 2)
    if R1.lo <= 0 or R2.lo <= 0:
        raise ArithmeticError(
            f"r(theta) enclosure cannot be made strictly positive at {theta}: "
            f"r1=[{R1.lo}, {R1.hi}] r2=[{R2.lo}, {R2.hi}]")

    def g(c):
        # cstarPt is inf{c > 0 : Phi(c; r) <= alpha}, so the locator may only search c > 0.
        # For c <= 0 the point is outside the domain of the infimum and certainly does not
        # satisfy Phi <= alpha, so it reports "above alpha" rather than evaluating.
        if c <= 0:
            return 1.0
        return float(pe.phi_r_beta(F(c), R1.lo, R2.lo).mid()) - float(ALPHA)

    # Proper bracketing, then bisection. The previous blind secant from c0=2.0 stalled at
    # the same-sign corners, where c* is well below 2.0: the secant cannot tell which side
    # of the root it started on, and a flat Phi makes ga==gb immediately. Instead bracket
    # first (double up or down from c0 until the sign of g changes) and only then bisect.
    a = c0
    ga = g(a)
    lo, hi = a, a
    if ga > 0:
        # root is above c0: double up until g <= 0
        hi = a
        for _ in range(80):
            hi = hi * 2
            if g(hi) <= 0:
                break
        else:
            raise ArithmeticError(f"could not bracket cstar above c0 at {theta}")
        lo = a
    else:
        # root is at or below c0: double down until g > 0 (lo), with hi = c0 where g <= 0
        hi = a
        for _ in range(80):
            lo = hi / 2
            if lo <= 0:
                lo = F(1, 1 << 40)
            if g(lo) > 0:
                break
            hi = lo
        else:
            raise ArithmeticError(f"could not bracket cstar below c0 at {theta}")
    # bisect [lo, hi] to localize, then apply the certified two-sided test
    for _ in range(60):
        mid = (lo + hi) / 2
        if g(mid) > 0:
            lo = mid
        else:
            hi = mid
    cm = hi
    for k in range(40):
        w = half * (1 << k)
        c_lo, c_hi = cm - w, cm + w
        if c_lo <= 0:
            continue
        up = phi_bounds(c_hi, R1.lo, R2.lo)[1]
        dn = phi_bounds(c_lo, R1.hi, R2.hi)[0]
        if up <= ALPHA and dn > ALPHA:
            return c_lo, c_hi
    raise ArithmeticError(f"no certified bracket at {theta}")


def _f3_point(theta):
    lo, hi = cstar_bracket(theta)
    return theta, lo, hi


# ------------------------------------------------------------------ F-5 -----------

def certify_with_shave(p, holds, cert):
    """`f8c.certify_rho`, but also returning WHICH fixed shave was accepted (diagnostic)."""
    if not p > 0.0:
        return None, None
    pf = F(p)
    for i, shave in enumerate(f8c.SHAVES):
        q = f8c._dyadic_down(pf * (1 - shave))
        if q > 0 and holds(q, cert)[0]:
            return q, i
    return None, None


def certify_box(box):
    cert = f8c.CellCert(box)
    p1, p2 = f8c.locate_rhos(cert)
    q1, i1 = certify_with_shave(p1, f8c.hlow1_holds, cert)
    q2, i2 = certify_with_shave(p2, f8c.hlow2_holds, cert)
    return q1, i1, q2, i2


def child_containing(box, t):
    for ch in f8c.children8(box):
        if all(ch[i].lo < t[i] < ch[i].hi for i in range(3)):
            return ch
    raise ValueError("target lies on a split plane")


def _f5_chain(args):
    idx, box, rho1, rho2 = args
    t = point_in(box, F5_TARGET)
    R1, R2 = r_enc(t)
    tlo, thi, _ = phi_bounds(UC, R1.lo, R2.lo)
    levels = []
    b = box
    for lev in range(F5_DEPTH + 1):
        q1, i1, q2, i2 = certify_box(b)
        ph = phi_bounds(UC, q1, q2) if (q1 and q2) else None
        levels.append({"lev": lev, "q1": q1, "i1": i1, "q2": q2, "i2": i2, "phi": ph,
                       "width": max(b[i].hi - b[i].lo for i in range(3))})
        b = child_containing(b, t)
    kids = []
    for ch in f8c.children8(box):
        q1, i1, q2, i2 = certify_box(ch)
        kids.append((q1, i1, q2, i2))
    return {"idx": idx, "rho1": rho1, "rho2": rho2, "t": t, "r1": R1, "r2": R2,
            "t_phi": (tlo, thi), "levels": levels, "kids": kids}


# ------------------------------------------------------------------ helpers --------

def read_claims(path: pathlib.Path) -> dict[str, str]:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"CLAIM ([^:]+): .*? = (.*)$", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def log_grid(lo, hi, m):
    """fcore.log_grid, re-derived (not imported): exp(linspace(log lo, log hi, m))."""
    import numpy as np
    return [float(x) for x in np.exp(np.linspace(np.log(lo), np.log(hi), int(m)))]


def linspace(lo, hi, m):
    import numpy as np
    return [float(x) for x in np.linspace(lo, hi, int(m))]


def fl(x) -> str:
    return repr(float(x))


# The caches are written and read ONLY by this script, in this working tree, and are keyed by
# the sha256 of the sources that produced them. They hold Fractions and the custom Iv class,
# so a plain JSON dump would corrupt them; this codec round-trips those types exactly and
# avoids pickle, which is a deserialisation sink.


def _enc(x: object) -> object:
    if isinstance(x, F):
        return {"#f": [x.numerator, x.denominator]}
    if isinstance(x, Iv):
        return {"#iv": [_enc(x.lo), _enc(x.hi)]}
    if isinstance(x, tuple):
        return {"#t": [_enc(v) for v in x]}
    if isinstance(x, list):
        return [_enc(v) for v in x]
    if isinstance(x, dict):
        return {"#d": [[str(k), _enc(v)] for k, v in x.items()]}
    if isinstance(x, bool) or x is None or isinstance(x, (int, float, str)):
        return x
    raise TypeError(f"cannot cache a {type(x).__name__}")


def _dec(x: object) -> object:
    if isinstance(x, dict):
        if "#f" in x:
            n, d = x["#f"]
            return F(int(n), int(d))
        if "#iv" in x:
            lo, hi = x["#iv"]
            return Iv(_dec(lo), _dec(hi))
        if "#t" in x:
            return tuple(_dec(v) for v in x["#t"])
        if "#d" in x:
            return {k: _dec(v) for k, v in x["#d"]}
        return x
    if isinstance(x, list):
        return [_dec(v) for v in x]
    return x


def save_cache(path: pathlib.Path, obj: dict) -> None:
    path.write_text(json.dumps(_enc(obj), separators=(",", ":")), encoding="utf-8")


def load_cache(path: pathlib.Path) -> dict:
    if not path.exists():
        return {}
    try:
        blob = _dec(json.loads(path.read_text(encoding="utf-8")))
        return blob if isinstance(blob, dict) else {}
    except (OSError, ValueError):
        return {}


def get_cover(workers):
    key = sha(ROOT / "formal01_f8c_certified.py") + sha(ROOT / "formal01_f8a_checker.py")
    if CACHE_COVER.exists():
        blob = load_cache(CACHE_COVER)
        if blob.get("key") == key:
            return blob["cells"], blob["unserved"], True
    cells, unserved = f8c.build_graded_cover(workers)
    save_cache(CACHE_COVER, {"key": key, "cells": cells, "unserved": unserved})
    return cells, unserved, False


def locate(cells_np, cells, pt):
    import numpy as np
    lo, hi = cells_np
    p = np.array([float(x) for x in pt])
    m = np.all((lo <= p + 1e-12) & (p <= hi + 1e-12), axis=1)
    out = []
    for i in np.nonzero(m)[0]:
        b = cells[int(i)]["box"]
        if all(b[j].lo <= pt[j] <= b[j].hi for j in range(3)):
            out.append(int(i))
    return out


# ================================================================== main ============

def main() -> int:
    t_all = time.time()
    workers = max(1, (mp.cpu_count() or 1) - 1)

    # ---------------------------------------------------------------- provenance ----
    r2p = DEP / "formal01_f8c_claims-R2.txt"
    r2sha = sha(r2p)
    R2 = read_claims(r2p)
    F8B = read_claims(DEP / "formal01_f8b_claims.txt")
    claim("P3A-SRC-BATTERY-SHA", "sha256 of formal01_f8a_battery.py", sha(pathlib.Path(__file__)))
    for nm in ("formal01_f8c_certified.py", "formal01_f8a_checker.py", "formal01_phi_engine.py"):
        claim(f"P3A-SRC-{nm}", f"sha256 of {nm}", sha(ROOT / nm))
    claim("P3A-ENCL-DEPOSIT", "enclosure deposit of record and its sha256",
          f"deposits/formal01_f8c_claims-R2.txt {r2sha}")
    claim("P3A-ENCL-SHA-PREFIX-OK", "deposit sha256 begins with the prefix named in advance",
          r2sha.startswith("25fdae5a0ceee53b"))
    claim("P3A-LEAN-SHA", "sha256 of lean/W4lean/Unit2BFormalR1.lean", sha(LEAN))
    claim("P3A-CONST", "Lc, Uc, alpha, tol, nu parsed from the approved Lean file",
          f"Lc={LC} Uc={UC} alpha={ALPHA} tol={TOL} nu={NU}")
    claim("P3A-CONST-UL", "Uc - Lc (exact) and Uc - Lc <= tol", f"{UC - LC} {UC - LC <= TOL}")
    claim("P3A-ENGINE-NU", "engine t-dof equals 2*nu from Lean", pe.NU_T == 2 * NU)
    claim("P3A-DEP-A4-WORST", "deposited A-4 worst upper endpoint (quoted)", R2["F8C-A4-WORST"])
    claim("P3A-DEP-A5-MARGIN", "deposited A-5 hWit margin (quoted)", F8B["A5-MARGIN"])

    # --------------------------------------------------------------------- cover ----
    t0 = time.time()
    cells, unserved, cached = get_cover(workers)
    claim("P3A-COVER-REBUILT", "cover rebuilt by build_graded_cover (or loaded from its keyed cache)",
          f"cells={len(cells)} unserved={len(unserved)} cached={cached} sec={time.time() - t0:.1f}")
    claim("P3A-COVER-MATCHES-DEPOSIT", "rebuilt cell count equals deposited F8C-CELLS, 0 unserved",
          len(cells) == int(R2["F8C-CELLS"]) and not unserved)
    vol = sum((f8c.total_volume(c["box"]) for c in cells), F(0))
    claim("P3A-COVER-VOLUME", "rebuilt leaf volumes sum exactly to the box",
          vol == (LAM_HI - LAM_LO) * (B_HI - B_LO) ** 2)
    import numpy as np
    cells_np = (np.array([[float(c["box"][i].lo) for i in range(3)] for c in cells]),
                np.array([[float(c["box"][i].hi) for i in range(3)] for c in cells]))

    # ======================================================================= F-1 ====
    f1_key = (sha(CACHE_COVER) if CACHE_COVER.exists() else "") + sha(ROOT / "formal01_phi_engine.py") \
        + repr(PROBES) + str(R_PREC) + str(UC)
    recs = None
    if CACHE_F1.exists():
        blob = load_cache(CACHE_F1)
        if blob.get("key") == f1_key:
            recs = blob["recs"]
    if recs is None:
        t0 = time.time()
        one = _f1_cell((0, cells[0]["box"], cells[0]["rho1"], cells[0]["rho2"]))
        claim("F1-TIMING-ONE-CELL", "seconds for one cell (24 probes + 8 vertices + 1 rho eval), serial",
              round(time.time() - t0, 3))
        del one
        t0 = time.time()
        args = [(i, c["box"], c["rho1"], c["rho2"]) for i, c in enumerate(cells)]
        recs = []
        # Worker recycling plus INCREMENTAL checkpointing. The first attempt held every
        # result in memory and wrote the cache only after the whole sweep, so when the
        # Windows spawn pool died at ~1212 s it discarded 20000 completed cells. Now a
        # crash costs at most one chunk, and the resume path above picks the cache up.
        with mp.Pool(workers, maxtasksperchild=64) as pool:
            for k, r in enumerate(pool.imap_unordered(_f1_cell, args, chunksize=32)):
                recs.append(r)
                if (k + 1) % 2000 == 0:
                    tmp = CACHE_F1.with_suffix(".tmp")
                    tmp.write_bytes(pickle.dumps({"key": f1_key, "recs": recs}))
                    tmp.replace(CACHE_F1)
                    print(f"  F-1 {k + 1}/{len(args)} cells, {time.time() - t0:.0f}s", flush=True)
        recs.sort(key=lambda r: r["idx"])
        CACHE_F1.write_bytes(pickle.dumps({"key": f1_key, "recs": recs}))
        claim("F1-SEC", "seconds for the F-1 interior sweep (parallel)", round(time.time() - t0, 1))

    n_pts = sum(r["n"] for r in recs)
    viol = sum(r["viol"] for r in recs)
    undec = sum(r["undec"] for r in recs)
    dfail = sum(r["dom_fail"] for r in recs)
    dund = sum(r["dom_undec"] for r in recs)
    vdfail = sum(r["vdom_fail"] for r in recs)
    vdund = sum(r["vdom_undec"] for r in recs)
    incons = sum(r["incons"] for r in recs)
    notint = sum(r["not_interior"] for r in recs)
    best = max(recs, key=lambda r: r["best_hi"])
    maxrad = max(r["maxrad"] for r in recs)
    binding = max(recs, key=lambda r: r["k_hi"])
    headroom = min(r["k_hi"] - r["best_hi"] for r in recs)

    claim("F1-EXPECTED", "F-1 expected (theory: rho_k <= r(theta) on the cell and Phi antitone)",
          "Phi(U; r(theta)) <= alpha at every sampled interior point; rho_k <= r(theta) at every "
          "point and vertex; upper(Phi(U;r(theta))) never above upper(Phi(U;rho_k)) of its cell")
    claim("F1-TOL", "F-1 tolerance as declared",
          "0 (exact rational comparison of rigorous ball endpoints with alpha = "
          f"{ALPHA}); SE = N/A, exact-arithmetic r(theta) at 2^-{R_PREC} + Arb R-beta balls")
    claim("F1-DESIGN", "probes per cell (8 V43-diagonal (r+1)/9, 8 {1/4,3/4}^3, 8 {1/64,63/64}^3), all strictly interior",
          f"{len(PROBES)} per cell x {len(cells)} cells; the sample is finer than the mesh in every cell")
    claim("F1-N-INTERIOR", "strictly interior points sampled", n_pts)
    claim("F1-NOT-INTERIOR", "probes that were NOT strictly interior (must be 0)", notint)
    claim("F1-MAX-PHI-INTERIOR", "largest rigorous upper(Phi(U; r(theta))) over interior points",
          fl(best["best_hi"]))
    claim("F1-MAX-PHI-INTERIOR-EXACT", "that upper endpoint as an exact dyadic", str(best["best_hi"]))
    claim("F1-MAX-PHI-POINT", "theta = (lam, b0, b1) of that point (float view)",
          tuple(round(float(x), 12) for x in best["best_pt"]))
    claim("F1-MAX-PHI-CELL", "its cell index and that cell's certified upper(Phi(U; rho_k))",
          f"{best['idx']} {fl(recs[best['idx']]['k_hi'])}")
    claim("F1-MARGIN-INTERIOR", "alpha - largest interior upper endpoint", fl(ALPHA - best["best_hi"]))
    claim("F1-VIOLATIONS", "interior points with lower(Phi) > alpha (proven violation)", viol)
    claim("F1-UNDECIDED", "interior points whose ball straddles alpha", undec)
    claim("F1-DOM-FAIL", "interior (point, component) pairs with rho_k > r_hi(theta) (proven A-3 breach)", dfail)
    claim("F1-DOM-UNDECIDED", "interior pairs with r_lo < rho_k <= r_hi (undecided)", dund)
    claim("F1-DOM-VERTEX-FAIL", "vertex pairs with rho_k > r_hi (closed-cell breach)", vdfail)
    claim("F1-DOM-VERTEX-UNDECIDED", "vertex pairs undecided", vdund)
    claim("F1-ENGINE-INCONSISTENT", "points with lower(Phi(U;r_lo)) > upper(Phi(U;rho_k)) (antitone breach)", incons)
    claim("F1-MIN-CELL-HEADROOM", "min over cells of upper(Phi(U;rho_k)) - max interior upper(Phi(U;r))",
          fl(headroom))
    claim("F1-MAX-RADIUS", "largest Phi ball radius in the F-1 sweep", fl(maxrad))
    claim("F1-BINDING-CELL", "cell with the largest upper(Phi(U;rho_k)) and that value (reproduces A-4 worst)",
          f"{binding['idx']} {fl(binding['k_hi'])}")
    claim("F1-A4-REPRODUCED", "rebuilt A-4 worst equals the deposited F8C-A4-WORST (float repr)",
          fl(binding["k_hi"]) == R2["F8C-A4-WORST"])

    # --- box points: the V43 G1 audit grid (1053), the exact corner, the 8 box vertices
    g1 = [(F(lam_), F(a), F(b)) for lam_ in log_grid(0.04, 25.0, 13)
          for a in linspace(-1.5, 1.5, 9) for b in linspace(-1.5, 1.5, 9)]
    corner = (LAM_LO, B_LO, B_HI)
    boxv = [((LAM_LO, LAM_HI)[(n >> 0) & 1], (B_LO, B_HI)[(n >> 1) & 1],
             (B_LO, B_HI)[(n >> 2) & 1]) for n in range(8)]
    in_box = [p for p in g1 if LAM_LO <= p[0] <= LAM_HI]
    claim("F1-G1-IN-BOX", "G1 points inside the declared box (float lambda endpoints)", f"{len(in_box)} of {len(g1)}")
    bpts = in_box + [corner] + boxv
    bargs, uncovered = [], 0
    for p in bpts:
        hit = locate(cells_np, cells, p)
        if not hit:
            uncovered += 1
        bargs.append((p, [(cells[i]["rho1"], cells[i]["rho2"]) for i in hit]))
    with mp.Pool(workers) as pool:
        brec = pool.map(_grid_point, bargs, chunksize=8)
    bbest = max(brec, key=lambda r: r["hi"])
    crec = brec[len(in_box)]
    claim("F1-BOXPTS-N", "box points checked (G1 + exact corner + 8 box vertices)", len(brec))
    claim("F1-BOXPTS-UNCOVERED", "box points contained in NO cell (coverage breach)", uncovered)
    claim("F1-BOXPTS-CONTAINMENTS", "total (point, containing cell) pairs checked (faces shared)",
          sum(r["ncells"] for r in brec))
    claim("F1-BOXPTS-DOM-FAIL", "box-point pairs with rho_k > r_hi", sum(r["dom_fail"] for r in brec))
    claim("F1-BOXPTS-DOM-UNDECIDED", "box-point pairs undecided", sum(r["dom_undec"] for r in brec))
    claim("F1-BOXPTS-VIOLATIONS", "box points with upper(Phi(U;r)) > alpha",
          sum(1 for r in brec if r["hi"] > ALPHA))
    claim("F1-BOXPTS-MAX", "largest upper(Phi(U; r)) over the box points, and where",
          f"{fl(bbest['hi'])} at {tuple(round(float(x), 12) for x in bbest['pt'])}")
    claim("F1-CORNER-PHI", "Phi(U; r(corner)) rigorous ball [lo, hi]", f"[{fl(crec['lo'])}, {fl(crec['hi'])}]")
    claim("F1-CORNER-MARGIN", "alpha - upper(Phi(U; r(corner)))", fl(ALPHA - crec["hi"]))
    f1_pass = (viol == 0 and undec == 0 and dfail == 0 and dund == 0 and vdfail == 0
               and vdund == 0 and incons == 0 and notint == 0 and uncovered == 0
               and all(r["dom_fail"] == 0 and r["dom_undec"] == 0 and r["hi"] <= ALPHA for r in brec))
    claim("F1-VERDICT", "F-1 dominance, decisive, on the real enclosure", "PASS" if f1_pass else "FAIL")

    # ======================================================================= F-2 ====
    order = sorted(range(len(cells)), key=lambda i: recs[i]["k_hi"], reverse=True)
    rng = random.Random(SEED)
    rand = rng.sample(order[N_TOP:], N_RAND)
    sample = order[:N_TOP] + rand
    claim("F2-TOL", "F-2 tolerance as declared",
          "a break counts only if PROVEN: lower ball endpoint > alpha (strict, exact); fixed exact "
          f"perturbations 2^-s, s in {list(F2_SHAVES)}")
    claim("F2-EXPECTED", "F-2 expected",
          "lowering rho_k (or U, or r(theta)) by a small fixed shave pushes Phi above alpha at the "
          "binding cell; RAISING rho cannot raise Phi (Phi_antitone) but breaks A-3 (rho > inf r)")

    def first_break(fn):
        for s in sorted(F2_SHAVES, reverse=True):         # smallest perturbation first
            if fn(F(1, 1 << s)):
                return s
        return None

    rows = []
    for i in sample:
        c = cells[i]
        r1, r2 = c["rho1"], c["rho2"]
        s_dn = first_break(lambda e, r1=r1, r2=r2:
                           phi_bounds(UC, r1 * (1 - e), r2 * (1 - e))[0] > ALPHA)
        s_u = first_break(lambda e, r1=r1, r2=r2:
                          phi_bounds(UC * (1 - e), r1, r2)[0] > ALPHA)
        up_hi = max(phi_bounds(UC, r1 * (1 + F(1, 1 << s)), r2 * (1 + F(1, 1 << s)))[1] for s in F2_SHAVES)
        cert = f8c.CellCert(c["box"])
        a3_rej = first_break(lambda e, r1=r1, r2=r2, cert=cert:
                             not (f8c.hlow1_holds(r1 * (1 + e), cert)[0]
                                  and f8c.hlow2_holds(r2 * (1 + e), cert)[0]))
        gen = first_break(lambda e, r1=r1, r2=r2, i=i, recs=recs:
                          r1 * (1 + e) > recs[i]["min_r1_hi"]
                          or r2 * (1 + e) > recs[i]["min_r2_hi"])
        rows.append({"i": i, "top": i in order[:N_TOP], "k_hi": recs[i]["k_hi"], "s_dn": s_dn,
                     "s_u": s_u, "up_hi": up_hi, "a3": a3_rej, "gen": gen})
    b0 = rows[0]
    claim("F2-BINDING-CELL", "binding cell index and alpha - upper(Phi(U;rho_k))",
          f"{b0['i']} {fl(ALPHA - b0['k_hi'])}")
    claim("F2-BIND-RHO-DOWN", "binding cell: smallest shave s with lower(Phi(U; rho_k(1-2^-s))) > alpha",
          b0["s_dn"])
    claim("F2-BIND-U-DOWN", "binding cell: smallest shave s with lower(Phi(U(1-2^-s); rho_k)) > alpha",
          b0["s_u"])
    claim("F2-BIND-RHO-UP-PHI", "binding cell, LITERAL 'above rho_k': max upper(Phi(U; rho_k(1+2^-s))) over s",
          f"{fl(b0['up_hi'])} (<= alpha: {b0['up_hi'] <= ALPHA})")
    claim("F2-BIND-RHO-UP-A3", "binding cell: smallest s at which A-3 REJECTS rho_k(1+2^-s)", b0["a3"])
    claim("F2-BIND-RHO-UP-GENUINE", "binding cell: smallest s at which rho_k(1+2^-s) exceeds r_hi at a sampled point",
          b0["gen"])
    lphi = phi_bounds(LC, cells[b0["i"]]["rho1"], cells[b0["i"]]["rho2"])
    claim("F2-BIND-AT-L", "binding cell: lower(Phi(Lc; rho_k)) - alpha", fl(lphi[0] - ALPHA))
    top = [r for r in rows if r["top"]]
    rnd = [r for r in rows if not r["top"]]
    for nm, grp in (("TOP", top), ("RAND", rnd)):
        claim(f"F2-{nm}-N", f"{nm} sample size", len(grp))
        claim(f"F2-{nm}-RHO-DOWN-BROKEN", f"{nm}: cells broken by some declared rho shave (down)",
              sum(1 for r in grp if r["s_dn"] is not None))
        claim(f"F2-{nm}-RHO-DOWN-SMALLEST-S", f"{nm}: per-cell smallest breaking s (None = unbroken)",
              [r["s_dn"] for r in grp])
        claim(f"F2-{nm}-UP-PHI-NEVER-ABOVE", f"{nm}: raising rho never gives upper(Phi) > alpha",
              all(r["up_hi"] <= ALPHA for r in grp))
        claim(f"F2-{nm}-UP-A3-REJECTED", f"{nm}: cells where A-3 rejects some rho_k(1+2^-s)",
              sum(1 for r in grp if r["a3"] is not None))
        claim(f"F2-{nm}-UP-GENUINE", f"{nm}: cells where some rho_k(1+2^-s) provably exceeds r at a point",
              sum(1 for r in grp if r["gen"] is not None))
    # certificate slack over ALL cells: how close rho_k is to the sampled infimum of r
    sl1 = sorted((r["min_r1_hi"] - cells[r["idx"]]["rho1"]) / cells[r["idx"]]["rho1"] for r in recs)
    sl2 = sorted((r["min_r2_hi"] - cells[r["idx"]]["rho2"]) / cells[r["idx"]]["rho2"] for r in recs)
    claim("F2-SLACK-R1", "rel. slack (min sampled r1_hi - rho1)/rho1 over all cells: min, median, max",
          f"{fl(sl1[0])} {fl(sl1[len(sl1) // 2])} {fl(sl1[-1])}")
    claim("F2-SLACK-R2", "rel. slack (min sampled r2_hi - rho2)/rho2 over all cells: min, median, max",
          f"{fl(sl2[0])} {fl(sl2[len(sl2) // 2])} {fl(sl2[-1])}")
    margins = [ALPHA - r["k_hi"] for r in recs]
    for thr in (F(1, 10 ** 5), F(1, 10 ** 4), F(1, 10 ** 3)):
        claim(f"F2-CELLS-MARGIN-LT-{thr.denominator}", f"cells with alpha - upper(Phi(U;rho_k)) < {thr}",
              sum(1 for m in margins if m < thr))
    # at the F-1 sampled argmax point
    th = best["best_pt"]
    TR1, TR2 = r_enc(th)
    s_pt = first_break(lambda e, TR1=TR1, TR2=TR2:
                       phi_bounds(UC, TR1.hi * (1 - e), TR2.hi * (1 - e))[0] > ALPHA)
    s_ptu = first_break(lambda e, TR1=TR1, TR2=TR2:
                        phi_bounds(UC * (1 - e), TR1.hi, TR2.hi)[0] > ALPHA)
    # The LITERAL "above r(theta)" probe needs r(theta) > 0 to be meaningful, because
    # Phi(c; rho1, rho2) is the model's rejection rate only for nonnegative ratios and the
    # engine refuses a negative rho. Where the enclosure of r does not certify strict
    # positivity at this point, the literal probe is UNDEFINED, not a failure -- so it is
    # recorded as such instead of raising. Previously this raised
    # "ValueError: c, rho1, rho2 must be >= 0" and killed the whole battery after F-1.
    r_positive = TR1.lo > 0 and TR2.lo > 0
    if r_positive:
        up_pt = max(phi_bounds(UC, TR1.lo * (1 + F(1, 1 << s)), TR2.lo * (1 + F(1, 1 << s)))[1]
                   for s in F2_SHAVES)
    else:
        up_pt = None
    claim("F2-POINT", "sampled point used (the F-1 interior argmax)", tuple(round(float(x), 12) for x in th))
    claim("F2-POINT-R-POSITIVE", "r(theta) enclosure lower endpoints both strictly positive at the point",
          f"r1.lo>0: {TR1.lo > 0} r2.lo>0: {TR2.lo > 0}")
    claim("F2-POINT-R-DOWN", "smallest s with lower(Phi(U; r_hi(theta)(1-2^-s))) > alpha", s_pt)
    claim("F2-POINT-U-DOWN", "smallest s with lower(Phi(U(1-2^-s); r_hi(theta))) > alpha", s_ptu)
    claim("F2-POINT-R-UP", "LITERAL 'above r(theta)': max upper(Phi(U; r_lo(1+2^-s))) over s, and <= alpha",
          "UNDEFINED (r not certified positive at this point)" if up_pt is None
          else f"{fl(up_pt)} {up_pt <= ALPHA}")
    if up_pt is None:
        claim("F2-POINT-R-UP-NOTE",
              "the literal 'above r(theta)' reading of the F-2 clause cannot be evaluated at this point; "
              "it is discharged instead on the binding CELL, where rho_k is strictly positive by "
              "construction, and the direction error in the clause is recorded separately",
              "RECORDED AS A FINDING, not a pass")
    wit = phi_bounds(LC, *[x.hi for x in r_enc(corner)])
    claim("F2-CORNER-AT-L", "constant side: lower(Phi(Lc; r_hi(corner))) - alpha (recomputes hWit)",
          fl(wit[0] - ALPHA))
    f2_pass = (b0["s_dn"] is not None and b0["s_u"] is not None and s_pt is not None
               and b0["a3"] is not None and wit[0] > ALPHA)
    cells_up_all = all(r["up_hi"] <= ALPHA for r in rows)
    pt_up_txt = "UNDEFINED (r not certified positive at that point)" if up_pt is None else (
        f"{up_pt <= ALPHA}")
    claim("F2-LITERAL-UP-VERDICT", "the literal perturbation 'above rho_k / above r(theta)' breaks Phi(U;rho)<=alpha",
          "CANNOT BREAK, by Phi_antitone; confirmed on every sampled cell "
          f"(all cells: {cells_up_all}) and at the sampled point ({pt_up_txt}). This is a direction error in the "
          "clause, recorded as a finding, not a pass")
    claim("F2-VERDICT", "F-2 non-vacuity (corrected direction: rho/U/r lowered; rho raised breaks A-3)",
          "PASS" if f2_pass else "FAIL")

    # ======================================================================= F-3 ====
    w2src = read_text(ROOT.parent / "04-efficient-rule" / "w2_cert.py")
    mb = need(re.search(r"bgrid = np\.linspace\(([-0-9.]+), ([-0-9.]+), (\d+)\)", w2src),
              "the bgrid linspace definition", "w2_cert.py")
    ml = need(re.search(r"lgrid = F\.log_grid\(([0-9.]+), ([0-9.]+), (\d+)\)", w2src),
              "the lgrid log_grid definition", "w2_cert.py")
    bgrid = linspace(float(mb.group(1)), float(mb.group(2)), int(mb.group(3)))
    lgrid = log_grid(float(ml.group(1)), float(ml.group(2)), int(ml.group(3)))
    w2j = read_json(ROOT.parent / "04-efficient-rule" / "deposits" / "w2_cert.json")
    w2d = w2j.get("data", w2j)
    claim("F3-GRID-SOURCE", "w2_cert.py grid read as text (lgrid, bgrid parameters)",
          f"log_grid({ml.group(1)},{ml.group(2)},{ml.group(3)}) x linspace({mb.group(1)},{mb.group(2)},{mb.group(3)})^2")
    claim("F3-GRID-MATCHES-DEPOSIT", "re-derived grid endpoints and size equal w2_cert.json box field",
          [lgrid[0], lgrid[-1]] == w2d["box"]["lam"] and [bgrid[0], bgrid[-1]] == w2d["box"]["beta"]
          and len(lgrid) * len(bgrid) ** 2 == w2d["box"]["points"])
    pts = [(F(lam_), F(a), F(b)) for lam_ in lgrid for a in bgrid for b in bgrid]
    t0 = time.time()
    with mp.Pool(workers) as pool:
        br = pool.map(_f3_point, pts, chunksize=4)
    claim("F3-SEC", "seconds for 175 rigorous c* brackets", round(time.time() - t0, 1))
    wmax = max(hi - lo for _, lo, hi in br)
    gmax_lo = max(lo for _, lo, _ in br)
    gmax_hi = max(hi for _, _, hi in br)
    arg = max(br, key=lambda x: x[2])
    claim("F3-N", "grid points bracketed", len(br))
    claim("F3-BRACKET-WIDTH", "largest bracket width", fl(wmax))
    claim("F3-GRIDMAX", "175-point grid maximum of cstarPt, rigorous [lower, upper]",
          f"[{gmax_lo}, {gmax_hi}] ~ {fl(gmax_lo)}")
    claim("F3-GRIDMAX-AT", "argmax grid point (float view)", tuple(float(x) for x in arg[0]))
    prior = read_claims(DEP / "p3a_battery_claims.txt")["BAT-F3-GRIDMAX"]
    claim("F3-PRIOR", "prior deposit's gridmax (quoted) and |recomputed - prior|",
          f"{prior} {fl(abs(F(prior) - gmax_lo))}")
    claim("F3-CHECK-L", "lower(gridmax) >= L  [exact]", f"{gmax_lo >= LC}  (gap {fl(gmax_lo - LC)})")
    claim("F3-CHECK-U", "upper(gridmax) <= U  [exact]", f"{gmax_hi <= UC}  (gap {fl(UC - gmax_hi)})")
    clo, chi = cstar_bracket(corner)
    corner_ok = gmax_lo >= LC
    uc_val = UC
    upper_ok = uc_val >= gmax_hi
    f3_contains = corner_ok and upper_ok
    claim("F3-CORNER-CSTAR", "cstarPt at the exact corner (1/25,-3/2,3/2), rigorous bracket",
          f"[{fl(clo)}, {fl(chi)}]; lo>=L {corner_ok}; U>=hi {upper_ok}")
    claim("F3-CORNER-VS-GRIDMAX", "cstar(corner) - gridmax (grid lambda is float 0.04 > 1/25)", fl(clo - gmax_hi))
    f3_pass = f3_contains
    claim("F3-TOL", "F-3 tolerance as declared", "0: exact rational containment L <= gridmax <= U")
    claim("F3-VERDICT", "F-3 lower anchor (reviewed reduction): L <= gridmax <= U", "PASS" if f3_pass else "FAIL")
    # V43 declared F-3 row, with its AlignedFactor(MID) anchor
    anc = F(w2d["certified_constant_over_box"])
    v43 = (ROOT / "declarations" / "2B-DECLARATION-V43.md").read_text(encoding="utf-8")
    tolf3 = F(need(re.search(r"TOL_F3 = ([0-9.e-]+)", v43), "TOL_F3", "the prior F-3 record").group(1))
    claim("F3-V43-ANCHOR", "V43 c_grid_anchor = w2_cert.json certified_constant_over_box (AlignedFactor MID)",
          w2d["certified_constant_over_box"])
    claim("F3-V43-LITERAL", "V43 literal: U - anchor >= -1e-10 ; |L - anchor| <= TOL_F3",
          f"{UC - anc >= F(-1, 10 ** 10)} ; {abs(LC - anc) <= tolf3} (U-anchor={fl(UC - anc)}, TOL_F3={tolf3})")
    claim("F3-V43-VERDICT", "V43 literal F-3 row against the AlignedFactor anchor",
          "FAIL — the anchor is an AlignedFactor(MID) value, a resolution the AF-withdrawal packet shows "
          "contradicted at this corner; retiring or replacing it is a decision outside this battery")

    # ======================================================================= F-5 ====
    f5_sample = order[:N_TOP] + rand
    t0 = time.time()
    with mp.Pool(workers) as pool:
        chains = pool.map(_f5_chain, [(i, cells[i]["box"], cells[i]["rho1"], cells[i]["rho2"])
                                      for i in f5_sample], chunksize=1)
    claim("F5-SEC", "seconds for the F-5 chains", round(time.time() - t0, 1))
    rev_rho, rev_phi, dom_bad, repro_bad, uncert, kid_rev, gap_ok = [], [], 0, 0, 0, [], 0
    shave_mix = 0
    max_rev = F(0)
    for ch in chains:
        L0 = ch["levels"]
        if (L0[0]["q1"], L0[0]["q2"]) != (ch["rho1"], ch["rho2"]):
            repro_bad += 1
        for a, b in zip(L0, L0[1:], strict=False):
            if a["q1"] is None or b["q1"] is None or a["q2"] is None or b["q2"] is None:
                uncert += 1
                continue
            for comp in ("q1", "q2"):
                if b[comp] < a[comp]:
                    rel = (a[comp] - b[comp]) / a[comp]
                    max_rev = max(max_rev, rel)
                    rev_rho.append((ch["idx"], b["lev"], comp, float(rel), a["i" + comp[1]], b["i" + comp[1]]))
            if a["phi"] and b["phi"] and b["phi"][0] > a["phi"][1]:
                rev_phi.append((ch["idx"], b["lev"]))
            if a["i1"] != b["i1"] or a["i2"] != b["i2"]:
                shave_mix += 1
        for lv in L0:
            if lv["q1"] is not None and (lv["q1"] > ch["r1"].hi or lv["q2"] > ch["r2"].hi):
                dom_bad += 1
        g0 = ch["r1"].lo - L0[0]["q1"]
        gN = ch["r1"].lo - L0[-1]["q1"]
        if gN < g0:
            gap_ok += 1
        for q1, _i1, q2, _i2 in ch["kids"]:
            if q1 is None or q2 is None or q1 < ch["rho1"] or q2 < ch["rho2"]:
                kid_rev.append(ch["idx"])
    ex = chains[0]
    claim("F5-TOL", "F-5 tolerance as declared", "0: any decrease of a certified rho along a nested chain, "
          "or any child rho below its parent's, is a reversal and an H-4 finding")
    claim("F5-DESIGN", "chains: cells, depth, target fractions",
          f"{len(chains)} cells (top {N_TOP} by A-4 + {N_RAND} seeded random), depth {F5_DEPTH} halvings, "
          f"target {tuple(str(x) for x in F5_TARGET)} of each cell")
    claim("F5-LEVEL0-REPRODUCES", "chains whose level-0 recertification differs from the cover's rho_k", repro_bad)
    claim("F5-UNCERTIFIED", "chain steps where a nested box could not be certified", uncert)
    claim("F5-RHO-REVERSALS", "reversals of certified rho along chains (count)", len(rev_rho))
    claim("F5-RHO-REVERSAL-MAX-REL", "largest relative reversal", fl(max_rev))
    claim("F5-RHO-REVERSAL-DETAIL", "first reversals (cell, level, comp, rel, parent shave idx, child shave idx)",
          rev_rho[:6])
    claim("F5-SHAVE-CHANGES", "chain steps where the accepted fixed shave index changed", shave_mix)
    claim("F5-PHI-REVERSALS", "chain steps with lower(Phi_child) > upper(Phi_parent)", len(rev_phi))
    claim("F5-DOMINANCE", "chain levels with certified rho above r_hi(target) (must be 0)", dom_bad)
    claim("F5-GAP-SHRINKS", "chains where r1_lo(target) - rho1 is smaller at depth 8 than at depth 0",
          f"{gap_ok} of {len(chains)}")
    claim("F5-CHILD-REVERSALS", "cells with some child rho below the parent's rho (or uncertified child)",
          len(set(kid_rev)))
    claim("F5-EXAMPLE-BINDING", "binding cell chain: rho1 gap to r1_lo(target) by level",
          [fl(ex["r1"].lo - lv["q1"]) for lv in ex["levels"]])
    claim("F5-EXAMPLE-BINDING-PHI", "binding cell chain: upper(Phi(U;rho)) by level, then Phi at the target",
          [fl(lv["phi"][1]) for lv in ex["levels"]] + [fl(ex["t_phi"][1])])
    f5_pass = (not rev_rho and not rev_phi and dom_bad == 0 and repro_bad == 0 and uncert == 0
               and not kid_rev and gap_ok == len(chains))
    claim("F5-VERDICT", "F-5 refinement monotone toward the true value",
          "PASS" if f5_pass else "FAIL (H-4 finding — see F5-* detail)")

    # ======================================================================= F-4 ====
    bnb = read_json(ROOT / "wAtoH_r2" / "bnb_final.json")
    bsrc = read_text(ROOT / "wAtoH_r2" / "replica.py")
    claim("F4-ANCHOR-B_df35", "bnb_final.json B_df35 [c_L, c_U], witness lambda (quoted)",
          f"[{bnb['B_df35']['c_L']}, {bnb['B_df35']['c_U']}] {bnb['B_df35']['witness_best']}")
    claim("F4-ANCHOR-A_df3", "bnb_final.json A_df3 [c_L, c_U], witness lambda (quoted)",
          f"[{bnb['A_df3']['c_L']}, {bnb['A_df3']['c_U']}] {bnb['A_df3']['witness_best']}")
    claim("F4-ANCHOR-DIM", "nuisance dimension K of each bnb anchor (len of witness lambda)",
          f"B_df35 K={len(bnb['B_df35']['witness_best'])} A_df3 K={len(bnb['A_df3']['witness_best'])}")
    claim("F4-ANCHOR-TARGET", "bnb target level (replica.py text): ALPHA0 = ABAR - GAMMA",
          need(re.search(r"^ALPHA0 = .*$", bsrc, re.M), "ALPHA0", "replica.py").group(0)
          + " ; " + need(re.search(r"^GAMMA = .*$", bsrc, re.M), "GAMMA", "replica.py").group(0))
    # 1-D specialisation of the 2B problem (b0 = b1 = 0): its critical scale, evidential
    c1 = cstar_bracket((F(1), F(0), F(0)))
    cl = cstar_bracket((LAM_LO, F(0), F(0)))
    claim("F4-2B-1D-CSTAR", "2B at b0=b1=0: cstarPt at lam=1 and lam=1/25 (=lam 25 by symmetry), rigorous",
          f"lam=1 [{fl(c1[0])}, {fl(c1[1])}] ; lam=1/25 [{fl(cl[0])}, {fl(cl[1])}]")
    import mpmath as mpm
    mpm.mp.dps = 50
    tq = mpm.findroot(lambda c: mpm.betainc(35, mpm.mpf(1) / 2, 0, 70 / (70 + c * c), regularized=True) / 2
                      - mpm.mpf(1) / 40, 2)
    claim("F4-CANDIDATE-T70", "candidate replacement anchor: t_70 upper 1/40 quantile by mpmath (independent of flint)",
          mpm.nstr(tq, 25))
    claim("F4-CANDIDATE-CONTAINED", "engine's rigorous cstarPt(1,0,0) bracket contains the mpmath t_70 quantile",
          F(c1[0]) <= F(str(mpm.nstr(tq, 40))) <= F(c1[1]))
    claim("F4-VERDICT", "F-4 precedent recovery at one nuisance dimension",
          "OPEN — the named anchors are not the same object; the t_70 check above is a "
          "CANDIDATE replacement anchor, not F-4, and adopting it is a decision outside this battery")

    # ============================================================== P3RB / F-9 =====
    kap = KAPPA_LEAN.read_text(encoding="utf-8")
    kap_sorry = len(re.findall(r"\bsorry\b", kap))
    claim("P3B-KAPPA-SORRY", "sorry count in Unit2BKappaStatementR1.lean (S-kappa entry point)", kap_sorry)
    claim("P3B-KAPPA-SHA", "sha256 of Unit2BKappaStatementR1.lean", sha(KAPPA_LEAN))
    allmax = max(best["best_hi"], bbest["hi"])
    claim("P3B-WITNESS", "numerical witness: max rigorous upper(Phi(U; r(theta))) over every sampled theta, <= alpha",
          f"{fl(allmax)} {allmax <= ALPHA} over {n_pts + len(brec)} points")
    claim("P3B-BOUND", "conditional FWER bound, as stated",
          "for every theta in Box: Phi(Uc; r1(theta), r2(theta)) <= alpha (cover_sound, premises A-2/A-3/A-4 of "
          "the R2 deposit); = RawRate(Uc; theta) <= alpha ONLY via kappa_reduction, which is sorry")
    claim("P3B-PMISS", "P(miss) status",
          "P(miss) = P(theta outside the declared Box). It is 0 CONDITIONAL on theta in Box (the declared "
          "operating regime) and otherwise equals the prior's mass outside Box; NOT simulated, no prior declared")
    claim("P3B-F9-i", "F-9 (i) finite-grid classifier",
          f"F-3: 175-point gridmax contained in [L,U] = {f3_pass}; EVIDENTIAL tier")
    claim("P3B-F9-ii", "F-9 (ii) continuum statement",
          "cover_sound + comp_thm_2B proved in Lean on [propext, Classical.choice, Quot.sound] as deposited; "
          "their premises hCover/hLow1/hLow2 are discharged by exact-rational checks (A-2, A-3), NOT in Lean")
    claim("P3B-F9-iii", "F-9 (iii) numerical certified bound",
          "hNum (A-4) and hWit (A-5) discharged by Arb ball arithmetic outside Lean; F-1 attacks them with 0 "
          f"violations (verdict {('PASS' if f1_pass else 'FAIL')}); the engine-to-Lean-Phi identity is a numerical edge")
    claim("P3B-F9-iv", "F-9 (iv) conditional FWER theorem + nuisance-region coverage",
          f"CONDITIONAL on S-kappa (kappa_reduction, {kap_sorry} sorry in its file, un-promoted) and on theta in Box")
    claim("P3B-F9-v", "F-9 (v) empirical relevance to observed returns",
          "NOT DISCHARGED: no edge from the Box to observed return data exists in this check")

    claim("P3A-SEC-TOTAL", "total seconds", round(time.time() - t_all, 1))
    claim("P3A-SUMMARY", "F-1 / F-2 / F-3 / F-4 / F-5",
          f"{'PASS' if f1_pass else 'FAIL'} / {'PASS' if f2_pass else 'FAIL'} / "
          f"{'PASS' if f3_pass else 'FAIL'} / OPEN / {'PASS' if f5_pass else 'FAIL'}")
    DEP.mkdir(exist_ok=True)
    OUT.write_text("\n".join(CLAIMS) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
