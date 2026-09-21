"""Rigorous two-sided enclosures of the conditional no-rejection probability.

Reduced form (SPEC section 3; identical to 02-bounded-rule/cert_sup.py so every existing deposit
stays comparable, and cross-checked against 02-bounded-rule/core.py in the self-test):

    F(c, lambda) = E_U[ prod_i Phi(s_i (U + c d_i)) ],   s_i = sqrt(n_i / lambda_i),
    U ~ N(0, 1/n_0)     equivalently     F = int G(z; c, lambda) dPhi(z),
    G(z; c, lambda) = prod_i Phi(s_i (z/sqrt(n_0) + c d_i)),

with `d_i` the drawn standard error of contrast i.

Box enclosure.  For Lambda = prod_i [lam_lo_i, lam_hi_i] put
s_lo_i = sqrt(n_i/lam_hi_i), s_hi_i = sqrt(n_i/lam_lo_i) and x_i(z) = z/sqrt(n_0) + c d_i.
Because Phi(s x) is increasing in s for x >= 0 and decreasing for x < 0, the inner
minimum and maximum over the s-range are attained at a corner:

    Psi_min,i = Phi(s_lo,i x_i) if x_i >= 0 else Phi(s_hi,i x_i)
    Psi_max,i = Phi(s_hi,i x_i) if x_i >= 0 else Phi(s_lo,i x_i)
    G_min = prod_i Psi_min,i,        G_max = prod_i Psi_max,i.

Both are continuous and INCREASING in z, and 0 <= G_min <= G <= G_max <= 1 for every
lambda in the box.  Monotonicity is the whole device.  For ANY partition
z_0 < ... < z_M with cell masses dPhi_j = Phi(z_{j+1}) - Phi(z_j),

    F_lo(c; Lambda) = sum_j G_min(z_j) dPhi_j + G_min(z_M) (1 - Phi(z_M))
    F_hi(c; Lambda) = sum_j G_max(z_{j+1}) dPhi_j + G_max(z_0) Phi(z_0) + (1 - Phi(z_M))

satisfy, rigorously,

    F_lo(c; Lambda) <= min_{lambda in Lambda} F(c; lambda)
    max_{lambda in Lambda} F(c; lambda) <= F_hi(c; Lambda).

(A left Riemann sum lower-bounds the Stieltjes integral of an increasing integrand and a
right sum upper-bounds it; the lower tail contributes at least 0 to the lower bound and at
most G_max(z_0) Phi(z_0) to the upper; the upper tail contributes at least
G_min(z_M)(1 - Phi(z_M)) and at most 1 - Phi(z_M).)  A degenerate box has s_lo = s_hi, so
G_min = G_max = G and the two bounds bracket F(c; lambda).  No derivative bound, Lipschitz
constant, covering radius or vertex-attainment property is used anywhere.

Direction of use (this is the part that must not be confused):

  * `c_upper_box`: the smallest c with F_lo(c; box) >= 1 - alpha.  Because F_lo is a lower
    enclosure of min_lambda F, this c is an UPPER bound on sup_{lambda in box} c*(lambda).
  * `c_lower_s`: the largest c with F_hi(c; {lambda_w}) < 1 - alpha at a WITNESS POINT.
    This is a LOWER bound on c*(lambda_w), hence on the box supremum.
    (An upper bound on max_lambda F does not lower-bound sup c*; the witness point does.)
  * `bracket_point`: both together, giving c*(lambda) in [c_L, c_U] with widths (SPEC 6).

Arithmetic.  `bounds_float` is the numpy search path (NOT certified).  `bounds_arb` is the
certificate: Arb ball arithmetic for every operation, endpoints taken from `.lower()` /
`.upper()`, cells below `mass_floor` dropped with the omitted mass added to the width.
Cell masses are precomputed once per partition (`prepare_partition`) as float endpoints
that bound the true masses from the correct side, and the sums run in Arb; a certified
claim is therefore about exactly the float64 `c`, `alpha`, `n`, `d`, `lambda` reported.
"""

from __future__ import annotations

import contextlib
import os

import numpy as np
import paths
from paths import NAN, _f, _i

Z_MAX = paths._env_float("WA2H_ZMAX", 9.0)
# Upper cap for the critical-value search: an additive widening loop cannot reach the
# O(1e2..1e4) critical values that extreme region vertices require (the P2 rootcert
# battery returned NaN brackets for exactly those witnesses), so the span doubles.
C_MAX = paths._env_float("WA2H_C_MAX", 1e6)

ERRLOG: list[str] = []
MASS_FLOOR = paths._env_float("WA2H_MASS_FLOOR", 1e-18)


# --------------------------------------------------------------------------
# float64 search path (not certified)
# --------------------------------------------------------------------------
def _lo(x, default=None):
    """Lower endpoint of an Arb ball, rounded DOWN (directed).

    Pass-1 verifier finding (MAJOR): `_f(x.lower())` is round-to-nearest, so a mass used as a
    lower bound could exceed the true mass by half an ulp (measured in 24/400 cells), and the
    final endpoint conversion leaked ~5.5e-17.  Directed rounding removes both leaks.
    """
    try:
        v = float(x.lower())
        return float(np.nextafter(v, -np.inf))
    except Exception:  # noqa: BLE001
        return float("nan") if default is None else default


def _hi(x, default=None):
    """Upper endpoint of an Arb ball, rounded UP (directed)."""
    try:
        v = float(x.upper())
        return float(np.nextafter(v, np.inf))
    except Exception:  # noqa: BLE001
        return float("nan") if default is None else default


def s_range(n, lam_lo, lam_hi):
    """Per-coordinate s-range of a box: s_i = sqrt(n_i/lam_i) (lam_hi gives s_lo)."""
    n = np.asarray(n, float)
    lam_lo = np.asarray(lam_lo, float)
    lam_hi = np.asarray(lam_hi, float)
    return np.sqrt(n[1:] / lam_hi), np.sqrt(n[1:] / lam_lo)


def base_partition(z_max=None, nodes=400):
    """Ascending z-nodes spanning [-z_max, z_max]."""
    zm = Z_MAX if z_max is None else _f(z_max, Z_MAX)
    rows = max(3, _i(nodes, 400))
    return np.linspace(-zm, zm, rows + 1)


def g_box(z, c, d_row, s_lo, s_hi, n0, mode="min"):
    """G_min / G_max at the nodes z (vectorised over z)."""
    from scipy.special import ndtr

    zz = np.asarray(z, float)
    cc = _f(c)
    s_lo = np.asarray(s_lo, float)
    s_hi = np.asarray(s_hi, float)
    d_row = np.asarray(d_row, float).reshape(-1)
    u = zz / np.sqrt(_f(n0, 1.0))
    acc = np.ones(zz.size)
    for i in range(len(s_lo)):
        x = u + cc * d_row[i]
        s = (
            np.where(x >= 0.0, s_lo[i], s_hi[i])
            if mode == "min"
            else np.where(x >= 0.0, s_hi[i], s_lo[i])
        )
        acc = acc * ndtr(s * x)
    return acc


def bounds_float(c, z, d_row, s_lo, s_hi, n0, slack=1e-12):
    """(lo, hi) float64 enclosure at the given partition.  NOT certified."""
    from scipy.special import ndtr

    try:
        zz = np.asarray(z, float)
        g_lo = g_box(zz, c, d_row, s_lo, s_hi, n0, "min")
        g_hi = g_box(zz, c, d_row, s_lo, s_hi, n0, "max")
        mass = ndtr(zz[1:]) - ndtr(zz[:-1])
        up_tail = 1.0 - ndtr(zz[-1])
        low_tail = ndtr(zz[0])
        lo = _f(np.dot(g_lo[:-1], mass) + g_lo[-1] * up_tail) - slack
        hi = _f(np.dot(g_hi[1:], mass) + g_hi[0] * low_tail + up_tail) + slack
        return lo, hi
    except Exception:  # noqa: BLE001
        return NAN, NAN


def cell_error(z, c, d_row, s_lo, s_hi, n0):
    """Per-cell contribution to the enclosure width (the refinement criterion)."""
    from scipy.special import ndtr

    zz = np.asarray(z, float)
    mass = ndtr(zz[1:]) - ndtr(zz[:-1])
    g_lo = g_box(zz[:-1], c, d_row, s_lo, s_hi, n0, "min")
    g_hi = g_box(zz[1:], c, d_row, s_lo, s_hi, n0, "max")
    return mass * np.maximum(g_hi - g_lo, 0.0)


def refine_partition(
    c,
    d_row,
    s_lo,
    s_hi,
    n0,
    z=None,
    target=1e-6,
    max_nodes=6000,
    share=0.25,
    rounds=80,
    z_max=None,
):
    """Bisect the worst cells until the float64 width <= target (or the budget is spent).

    Rigour never depends on this: any partition gives a valid enclosure, refinement only
    decides tightness.  Returns the ascending node array.
    """
    zz = np.asarray(z, float) if z is not None else base_partition(z_max)
    cap = max(64, _i(max_nodes, 6000))
    tgt = _f(target, 1e-6)
    for _round in range(max(1, _i(rounds, 80))):
        err = cell_error(zz, c, d_row, s_lo, s_hi, n0)
        total = _f(np.sum(err))
        if not np.isfinite(total) or total <= tgt or zz.size >= cap:
            break
        k = max(2, _i(_i(_f(share, 0.25) * err.size), 2))
        idx = np.argsort(err)[-k:]
        mids = 0.5 * (zz[idx] + zz[idx + 1])
        zz = np.sort(np.concatenate([zz, mids]))
    return zz


def _bisect_float(pred, lo, hi, tol):
    """Float64 bisection on a monotone predicate; returns (c_true_end, c_false_end)."""
    if not pred(hi):
        return None, None
    if pred(lo):
        return lo, lo
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return hi, lo


def c_ub_float(
    d_row,
    n0,
    s_lo,
    s_hi,
    z,
    alpha=paths.ALPHA0,
    c_guess=3.0,
    tol=1e-8,
    widen=0.5,
    grow=40,
):
    """Float64 estimate of the certified upper bound (search path only)."""
    level = 1.0 - _f(alpha, paths.ALPHA0)

    def ok(c):
        lo, _hi = bounds_float(c, z, d_row, s_lo, s_hi, n0)
        return np.isfinite(lo) and lo >= level

    c_lo = max(1e-4, _f(c_guess, 3.0) - _f(widen, 0.5))
    c_hi = max(c_lo + 0.1, _f(c_guess, 3.0) + _f(widen, 0.5))
    n = 0
    span = max(0.5, _f(widen, 0.5))
    while not ok(c_hi):
        span *= 2.0
        c_hi = c_lo + span
        n += 1
        if n > _i(grow, 80) or c_hi > C_MAX:
            return NAN
    while ok(c_lo) and c_lo > 1e-4:
        c_lo = max(1e-4, c_lo - _f(widen, 0.5))
        n += 1
        if n > _i(grow, 40):
            break
    hi, _lo = _bisect_float(ok, c_lo, c_hi, _f(tol, 1e-8))
    return _f(hi, NAN)


def c_lb_float(
    s_pt, d_row, n0, z, alpha=paths.ALPHA0, c_guess=3.0, tol=1e-8, widen=0.5, grow=40
):
    """Float64 estimate of the witness lower bound (search path only)."""
    level = 1.0 - _f(alpha, paths.ALPHA0)

    def below(c):
        _lo, hi = bounds_float(c, z, d_row, s_pt, s_pt, n0)
        return np.isfinite(hi) and hi < level

    c_lo = max(1e-4, _f(c_guess, 3.0) - _f(widen, 0.5))
    c_hi = max(c_lo + 0.1, _f(c_guess, 3.0) + _f(widen, 0.5))
    n = 0
    span = max(0.5, _f(widen, 0.5))
    while below(c_hi):
        span *= 2.0
        c_hi = c_lo + span
        n += 1
        if n > _i(grow, 80) or c_hi > C_MAX:
            return NAN
    while below(c_lo) and c_lo > 1e-4:
        c_lo = max(1e-4, c_lo - _f(widen, 0.5))
        n += 1
        if n > _i(grow, 40):
            break
    # largest c with below(c) True
    lo, hi = c_lo, c_hi
    while hi - lo > _f(tol, 1e-8):
        mid = 0.5 * (lo + hi)
        if below(mid):
            lo = mid
        else:
            hi = mid
    return _f(lo, NAN)


# --------------------------------------------------------------------------
# Arb (certified) path
# --------------------------------------------------------------------------
def _set_prec(prec):
    with contextlib.suppress(Exception):
        from flint import ctx

        ctx.prec = max(20, _i(prec, 64))


def Phi_arb(a, prec=None):
    """Rigorous enclosure of Phi(a) as an Arb ball."""
    from flint import arb

    if prec is not None:
        _set_prec(prec)
    return (1 + (a / arb(2).sqrt()).erf()) / 2


def prepare_partition(z, prec=64, mass_floor=None):
    """Precompute rigorous cell and tail masses for a partition (one Arb pass).

    Masses are stored as float endpoints that bound the true masses from the correct
    side (lower for the lower sum, upper for the upper sum), and the sums themselves run
    in Arb, so the accumulated result is rigorous.
    """
    from flint import arb

    _set_prec(prec)
    zz = np.asarray(z, float)
    mf = MASS_FLOOR if mass_floor is None else _f(mass_floor, MASS_FLOOR)
    za = [arb(_f(v, 0.0)) for v in zz]
    ph = [Phi_arb(v, prec) for v in za]
    mass_lo, mass_up, dropped = [], [], 0.0
    for j in range(len(za) - 1):
        m = ph[j + 1] - ph[j]
        lo_j = max(0.0, _lo(m, 0.0))
        up_j = max(0.0, _hi(m, 0.0))
        if up_j < mf:
            dropped += up_j
            mass_lo.append(0.0)
            mass_up.append(0.0)
        else:
            mass_lo.append(lo_j)
            mass_up.append(up_j)
    return {
        "z": zz,
        "mass_lo": np.asarray(mass_lo, float),
        "mass_up": np.asarray(mass_up, float),
        "up_tail_lo": max(0.0, _lo(1 - ph[-1], 0.0)),
        "up_tail_up": max(0.0, _hi(1 - ph[-1], 0.0)),
        "low_tail_up": max(0.0, _hi(ph[0], 0.0)),
        "dropped_up": dropped,
        "prec": _i(prec, 64),
        "mass_floor": mf,
        "n_nodes": _i(zz.size),
    }


def _g_at_node(z_arb, c_arb, d_a, s_lo_a, s_hi_a, sq_n0, prec):
    """(G_min, G_max) at one node, rigorously."""
    from flint import arb

    u = z_arb / sq_n0
    g_lo = arb(1)
    g_hi = arb(1)
    for i in range(len(s_lo_a)):
        x = u + c_arb * d_a[i]
        xm = s_lo_a[i] * x
        xM = s_hi_a[i] * x
        xl = _f(x.lower(), -1.0)
        xu = _f(x.upper(), 1.0)
        if xl >= 0.0:  # x >= 0 provably
            f_lo, f_hi = Phi_arb(xm, prec), Phi_arb(xM, prec)
        elif xu < 0.0:  # x < 0 provably
            f_lo, f_hi = Phi_arb(xM, prec), Phi_arb(xm, prec)
        else:  # sign not decided: use the corner extremes
            lo = min(_f(xm.lower(), 0.0), _f(xM.lower(), 0.0))
            hi = max(_f(xm.upper(), 0.0), _f(xM.upper(), 0.0))
            f_lo, f_hi = Phi_arb(arb(lo), prec), Phi_arb(arb(hi), prec)
        g_lo = g_lo * f_lo
        g_hi = g_hi * f_hi
    return g_lo, g_hi


def bounds_arb(c, d_row, n0, s_lo, s_hi, prep, need=("lo", "hi")):
    """Rigorous (lo, hi): lo <= min_lambda F(c,lambda), max_lambda F(c,lambda) <= hi.

    `need` selects which direction to compute (the other is returned as nan), because a
    bisection only ever tests one direction.  `s_lo == s_hi` for a point lambda.
    """
    from flint import arb

    try:
        prec = prep.get("prec", 64)
        _set_prec(prec)
        zz = prep["z"]
        c_arb = arb(_f(c))
        n0a = arb(_f(n0, 1.0))
        sq_n0 = n0a.sqrt()
        d_a = [arb(_f(v)) for v in np.asarray(d_row, float).reshape(-1)]
        s_lo_a = [arb(_f(v)) for v in np.asarray(s_lo, float)]
        s_hi_a = [arb(_f(v)) for v in np.asarray(s_hi, float)]
        want_lo = "lo" in need
        want_hi = "hi" in need
        g_min = []
        g_max = []
        for j in range(zz.size):
            za = arb(_f(zz[j], 0.0))
            if want_lo and want_hi:
                glo, ghi = _g_at_node(za, c_arb, d_a, s_lo_a, s_hi_a, sq_n0, prec)
                g_min.append(glo)
                g_max.append(ghi)
            elif want_lo:
                g_min.append(_g_at_node(za, c_arb, d_a, s_lo_a, s_hi_a, sq_n0, prec)[0])
            else:
                g_max.append(_g_at_node(za, c_arb, d_a, s_lo_a, s_hi_a, sq_n0, prec)[1])
        lo_out = NAN
        hi_out = NAN
        if want_lo:
            acc = arb(0)
            for j in range(zz.size - 1):
                m = prep["mass_lo"][j]
                if m > 0.0:
                    acc += g_min[j] * arb(_f(m))
            acc += g_min[-1] * arb(_f(prep["up_tail_lo"]))
            # fold the dropped mass into the Arb sum BEFORE the directed conversion: adjusting the
            # float endpoint afterwards is not directed and leaked ~6.6e-18 (verifier finding 5)
            acc -= arb(_f(prep.get("dropped_up", 0.0), 0.0))
            lo_out = _lo(acc, NAN)
        if want_hi:
            acc = arb(0)
            for j in range(zz.size - 1):
                m = prep["mass_up"][j]
                if m > 0.0:
                    acc += g_max[j + 1] * arb(_f(m))
            acc += g_max[0] * arb(_f(prep["low_tail_up"]))
            acc += arb(_f(prep["up_tail_up"]))
            acc += arb(_f(prep.get("dropped_up", 0.0), 0.0))
            hi_out = _hi(acc, NAN)
        return lo_out, hi_out
    except Exception:  # noqa: BLE001
        return NAN, NAN


def c_upper_box(
    d_row,
    n0,
    s_lo,
    s_hi,
    prep,
    alpha=paths.ALPHA0,
    c_guess=None,
    tol=1e-6,
    widen=0.5,
    grow=40,
):
    """Certified upper bound on sup_{lambda in box} c*(lambda).

    Smallest float64 c with the RIGOROUS bound F_lo(c; box) >= 1 - alpha verified.
    Returns a dict with c_ub, the number of Arb evaluations, and the rigorous F bounds
    at c_ub (so the caller can deposit the enclosure at the answer).
    """
    level = 1.0 - _f(alpha, paths.ALPHA0)
    out = {
        "c_ub": NAN,
        "n_arb": 0,
        "f_lo": NAN,
        "f_hi": NAN,
        "bracket_width": NAN,
        "warm_start": NAN,
    }

    def ok(c):
        lo, _hi = bounds_arb(c, d_row, n0, s_lo, s_hi, prep, need=("lo",))
        out["n_arb"] += 1
        return np.isfinite(lo) and lo >= level

    try:
        if c_guess is None:
            zf = np.linspace(-Z_MAX, Z_MAX, 401)
            c_guess = _f(c_ub_float(d_row, n0, s_lo, s_hi, zf, alpha=alpha), 3.0)
            if not np.isfinite(c_guess):
                c_guess = 3.0
        out["warm_start"] = _f(c_guess)
        c_lo = max(1e-4, _f(c_guess) - _f(widen, 0.5))
        c_hi = max(c_lo + 0.1, _f(c_guess) + _f(widen, 0.5))
        n = 0
        while not ok(c_hi):
            c_hi += _f(widen, 0.5)
            n += 1
            if n > _i(grow, 40):
                return out
        while ok(c_lo) and c_lo > 1e-4:
            c_lo = max(1e-4, c_lo - _f(widen, 0.5))
            n += 1
            if n > _i(grow, 40):
                break
        tgt = _f(tol, 1e-6)
        while c_hi - c_lo > tgt:
            mid = 0.5 * (c_lo + c_hi)
            if ok(mid):
                c_hi = mid
            else:
                c_lo = mid
        out["bracket_width"] = c_hi - c_lo
        flo, fhi = bounds_arb(c_hi, d_row, n0, s_lo, s_hi, prep, need=("lo", "hi"))
        out.update({"c_ub": _f(c_hi), "f_lo": flo, "f_hi": fhi})
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: {exc!r}")
    return out


def c_lower_s(
    s_pt,
    d_row,
    n0,
    prep,
    alpha=paths.ALPHA0,
    c_guess=None,
    tol=1e-6,
    widen=0.5,
    grow=40,
):
    """Certified lower bound on c*(lambda) at a witness point (degenerate box).

    Largest float64 c with F_hi(c; {lambda}) < 1 - alpha verified, hence c*(lambda) > c.
    """
    level = 1.0 - _f(alpha, paths.ALPHA0)
    s_arr = np.asarray(s_pt, float)
    out = {
        "c_lb": NAN,
        "n_arb": 0,
        "f_lo": NAN,
        "f_hi": NAN,
        "bracket_width": NAN,
        "warm_start": NAN,
    }

    def below(c):
        _lo, hi = bounds_arb(c, d_row, n0, s_arr, s_arr, prep, need=("hi",))
        out["n_arb"] += 1
        return np.isfinite(hi) and hi < level

    try:
        if c_guess is None:
            zf = np.linspace(-Z_MAX, Z_MAX, 401)
            c_guess = _f(c_lb_float(s_arr, d_row, n0, zf, alpha=alpha), 3.0)
            if not np.isfinite(c_guess):
                c_guess = 3.0
        out["warm_start"] = _f(c_guess)
        c_lo = max(1e-4, _f(c_guess) - _f(widen, 0.5))
        c_hi = max(c_lo + 0.1, _f(c_guess) + _f(widen, 0.5))
        n = 0
        while below(c_hi):
            c_hi += _f(widen, 0.5)
            n += 1
            if n > _i(grow, 40):
                return out
        while below(c_lo) and c_lo > 1e-4:
            c_lo = max(1e-4, c_lo - _f(widen, 0.5))
            n += 1
            if n > _i(grow, 40):
                break
        tgt = _f(tol, 1e-6)
        while c_hi - c_lo > tgt:
            mid = 0.5 * (c_lo + c_hi)
            if below(mid):
                c_lo = mid
            else:
                c_hi = mid
        out["bracket_width"] = c_hi - c_lo
        flo, fhi = bounds_arb(c_lo, d_row, n0, s_arr, s_arr, prep, need=("lo", "hi"))
        out.update({"c_lb": _f(c_lo), "f_lo": flo, "f_hi": fhi})
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: {exc!r}")
    return out


def bracket_point(
    s_pt, d_row, n0, prep, alpha=paths.ALPHA0, c_guess=None, tol=1e-6, slope_h=1e-4
):
    """SPEC section 6 root certificate at one lambda: c* in [c_L, c_U] with widths.

    Returns c_L, c_U, width, the probability enclosure width Delta p at each end, and a
    local slope bound on dF/dc used to convert Delta p into a c-equivalent.
    """
    res = {
        "c_L": NAN,
        "c_U": NAN,
        "width": NAN,
        "dp": NAN,
        "slope_hi": NAN,
        "c_equiv": NAN,
        "n_arb": 0,
        "warm_start": NAN,
    }
    try:
        up = c_upper_box(
            d_row,
            n0,
            np.asarray(s_pt, float),
            np.asarray(s_pt, float),
            prep,
            alpha=alpha,
            c_guess=c_guess,
            tol=tol,
        )
        lo = c_lower_s(
            s_pt, d_row, n0, prep, alpha=alpha, c_guess=up.get("warm_start"), tol=tol
        )
        res["n_arb"] = _i(up.get("n_arb"), 0) + _i(lo.get("n_arb"), 0)
        res["warm_start"] = up.get("warm_start")
        res.update({"c_U": up.get("c_ub"), "c_L": lo.get("c_lb")})
        res["width"] = _f(res["c_U"]) - _f(res["c_L"])
        dp_up = _f(up.get("f_hi")) - _f(up.get("f_lo"))
        dp_lo = _f(lo.get("f_hi")) - _f(lo.get("f_lo"))
        res["dp"] = _f(np.nanmax([dp_up, dp_lo]))
        c_ref = _f(res["c_U"], 3.0)
        s_arr = np.asarray(s_pt, float)
        f_p, _ = bounds_arb(
            c_ref + slope_h, d_row, n0, s_arr, s_arr, prep, need=("lo",)
        )
        _l, f_m = bounds_arb(
            c_ref - slope_h, d_row, n0, s_arr, s_arr, prep, need=("hi",)
        )
        slope = max(0.0, (_f(f_p, 0.0) - _f(f_m, 0.0)) / (2.0 * slope_h))
        res["slope_hi"] = slope
        res["c_equiv"] = (_f(res["dp"]) / slope) if slope > 1e-12 else NAN
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: {exc!r}")
    return res


# --------------------------------------------------------------------------
# self-test: falsifiers a-g of the plan (P1 step 7)
# --------------------------------------------------------------------------
def _draw(n, lam, seed, D=1):
    from core import draw_var

    s0, si, lh = draw_var(np.asarray(n, float), np.asarray(lam, float), D, seed)
    return s0, si, lh


def self_test(alpha=paths.ALPHA0, prec=64, target=1e-6, seed=20260918, report=True):
    """Falsifiers: containment, collapse, refinement, mass floor, Arb-vs-float."""
    import contextlib as _ctx

    from core import cond_noreject

    n = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    s0, si, _lh = _draw(n, lam, seed, D=1)
    d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    s_pt = np.sqrt(n[1:] / lam)
    out: dict[str, object] = {
        "n": list(n),
        "lam": list(lam),
        "d_row": list(d_row),
        "alpha": _f(alpha),
        "target": _f(target),
        "prec": _i(prec),
    }

    # (a) containment against an independent solver, at a point
    cont = []
    for c in (2.0, 2.5, 2.75, 3.0):
        z = refine_partition(c, d_row, s_pt, s_pt, n[0], target=target)
        prep = prepare_partition(z, prec=prec)
        flo, fhi = bounds_arb(c, d_row, n[0], s_pt, s_pt, prep)
        ref = _f(cond_noreject(np.array([c]), np.array([s0[0]]), si[0:1], n, lam)[0])
        cont.append(
            {
                "c": _f(c),
                "nodes": _i(z.size),
                "F_lo": flo,
                "F_hi": fhi,
                "F_ref": ref,
                "width": _f(fhi) - _f(flo),
                "lo_le_ref": bool(flo <= ref),
                "ref_le_hi": bool(ref <= fhi),
                "ref_inside": bool(flo <= ref <= fhi),
            }
        )
    out["containment"] = cont

    # (b) box containment: the true F at the four corners must sit inside the box bounds
    lo_f, hi_f = lam * 0.75, lam * 1.25
    s_lo, s_hi = s_range(n, lo_f, hi_f)
    c = 2.75
    z = refine_partition(c, d_row, s_lo, s_hi, n[0], target=target)
    prep = prepare_partition(z, prec=prec)
    flo, fhi = bounds_arb(c, d_row, n[0], s_lo, s_hi, prep)
    corners = []
    for a in (lo_f[0], hi_f[0]):
        for b in (lo_f[1], hi_f[1]):
            lc = np.array([a, b])
            corners.append(
                {
                    "lam": [a, b],
                    "F_ref": _f(
                        cond_noreject(np.array([c]), np.array([s0[0]]), si[0:1], n, lc)[
                            0
                        ]
                    ),
                }
            )
    out["box_containment"] = {
        "box": [list(lo_f), list(hi_f)],
        "c": _f(c),
        "F_lo": flo,
        "F_hi": fhi,
        "corners": corners,
        "F_lo_le_min_corner": bool(flo <= min(r["F_ref"] for r in corners)),
        "F_hi_ge_max_corner": bool(fhi >= max(r["F_ref"] for r in corners)),
        "width": _f(fhi) - _f(flo),
    }

    # (c) width vs node count (thinner tails of the refinement criterion)
    widths = []
    for nodes in (100, 400, 1600, 3200):
        zz = base_partition(nodes=nodes)
        l1, h1 = bounds_float(c, zz, d_row, s_pt, s_pt, n[0])
        widths.append({"nodes": _i(nodes), "width": _f(h1) - _f(l1)})
    out["width_vs_nodes"] = widths

    # (d) box collapse: shrinking the box tightens toward the point value
    collapse = []
    for mult in (1.05, 1.01, 1.001, 1.0):
        a, b = lam * (2.0 - mult), lam * mult
        sl, sh = s_range(n, a, b)
        zz = refine_partition(c, d_row, sl, sh, n[0], target=target)
        pr = prepare_partition(zz, prec=prec)
        l2, h2 = bounds_arb(c, d_row, n[0], sl, sh, pr)
        collapse.append(
            {"mult": _f(mult), "width": _f(h2) - _f(l2), "lo": l2, "hi": h2}
        )
    out["collapse"] = collapse

    # (e) mass floor: a larger floor may only widen the enclosure
    tight = prepare_partition(base_partition(nodes=800), prec=prec, mass_floor=0.0)
    coarse = prepare_partition(base_partition(nodes=800), prec=prec, mass_floor=1e-4)
    t1 = bounds_arb(c, d_row, n[0], s_pt, s_pt, tight)
    t2 = bounds_arb(c, d_row, n[0], s_pt, s_pt, coarse)
    out["mass_floor"] = {
        "width_floor_0": _f(t1[1]) - _f(t1[0]),
        "width_floor_1e-4": _f(t2[1]) - _f(t2[0]),
        "wider_or_equal": bool(
            (_f(t2[1]) - _f(t2[0])) >= (_f(t1[1]) - _f(t1[0])) - 1e-15
        ),
        "dropped_up": coarse.get("dropped_up"),
    }

    # (f) Arb vs float on the same partition
    zz = base_partition(nodes=800)
    prep = prepare_partition(zz, prec=prec)
    flo_a, fhi_a = bounds_arb(c, d_row, n[0], s_pt, s_pt, prep)
    flo_f, fhi_f = bounds_float(c, zz, d_row, s_pt, s_pt, n[0], slack=0.0)
    out["arb_vs_float"] = {
        "F_lo_arb": flo_a,
        "F_lo_float": flo_f,
        "F_hi_arb": fhi_a,
        "F_hi_float": fhi_f,
        "abs_diff_lo": abs(_f(flo_a) - _f(flo_f)),
        "abs_diff_hi": abs(_f(fhi_a) - _f(fhi_f)),
        "consistent": bool(
            abs(_f(flo_a) - _f(flo_f)) < 1e-9 and abs(_f(fhi_a) - _f(fhi_f)) < 1e-9
        ),
    }

    # (g) point bracket against the independent oracle solver
    with _ctx.suppress(Exception):
        from core import oracle_c

        c_or = oracle_c(n, lam, _f(alpha))
        br = bracket_point(
            s_pt,
            d_row,
            n[0],
            prepare_partition(
                refine_partition(3.0, d_row, s_pt, s_pt, n[0], target=target), prec=prec
            ),
            alpha=alpha,
        )
        out["bracket"] = {
            "c_oracle_independent": c_or,
            **br,
            "contains_oracle": bool(
                _f(br.get("c_L")) <= c_or <= _f(br.get("c_U")) + 1e-9
            ),
        }
    if report:
        print(f"self_test n={list(n)} lam={list(lam)} alpha={_f(alpha)}")
        for r in cont:
            print(
                f"  c={r['c']:.2f} nodes={r['nodes']:<5} F in [{r['F_lo']:.12f}, "
                f"{r['F_hi']:.12f}] ref {r['F_ref']:.12f} inside={r['ref_inside']} "
                f"width {r['width']:.2e}"
            )
        b = out["box_containment"]
        print(
            f"  box: lo<=min corner {b['F_lo_le_min_corner']}, "
            f"hi>=max corner {b['F_hi_ge_max_corner']}, width {b['width']:.2e}"
        )
        print(
            "  width vs nodes:",
            [(r["nodes"], f"{r['width']:.2e}") for r in out["width_vs_nodes"]],
        )
        print("  collapse widths:", [f"{r['width']:.2e}" for r in out["collapse"]])
        print("  mass floor:", out["mass_floor"])
        print("  arb vs float:", out["arb_vs_float"])
        if "bracket" in out:
            br = out["bracket"]
            print(
                f"  bracket [{br['c_L']:.6f}, {br['c_U']:.6f}] width {br['width']:.2e} "
                f"oracle {br['c_oracle_independent']:.6f} contains="
                f"{br['contains_oracle']} dp={br['dp']:.2e} c_equiv={br['c_equiv']:.2e}"
            )
    return out


if __name__ == "__main__":
    import tiers

    try:
        res = self_test()
        ok = (
            all(r["ref_inside"] for r in res["containment"])
            and res["box_containment"]["F_lo_le_min_corner"]
            and res["box_containment"]["F_hi_ge_max_corner"]
            and res["mass_floor"]["wider_or_equal"]
            and res["arb_vs_float"]["consistent"]
        )
        res["verdict"] = "PASS" if ok else "FAIL"
        print("RESULT:", res["verdict"])
        try:
            tiers.write_deposit(
                res,
                os.path.join(paths.DEPOSITS, "enclose_selftest.json"),
                "region_sup_certified",
                notes="Riemann-path self-test; independent cross-check path "
                "(the certificate path is enclose_gl.py)",
            )
        except Exception as exc:  # noqa: BLE001
            print("deposit failed:", repr(exc))
    except Exception as exc:  # noqa: BLE001
        print("SELF-TEST ERROR:", repr(exc))
