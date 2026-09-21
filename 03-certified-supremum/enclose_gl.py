"""Certificate quadrature: Gauss-Legendre with a rigorous Bernstein-ellipse remainder.

Why not the monotone Riemann rule of `enclose.py`?  It converges like h^2 as a function
of the node count (measured: 160k nodes for a 2.6e-6 enclosure width), which cannot reach
the predeclared `Delta p` target at a feasible Arb cost.  The integrand here is
*entire* (Phi is entire and exp(-z^2/2) is entire), so a high-order rule with an
analyticity remainder converges exponentially instead.

The rule.  For a panel [a, b] with centre m and half-width h, the n-node Gauss-Legendre
rule on [-1, 1], and the panel integral of f = G_min * phi (or G_max * phi):

    int_a^b f dz  =  h * sum_i w_i f(m + h xi_i)  +  E_n ,

and, for f analytic in the Bernstein ellipse E_rho with |f| <= M_rho there (Trefethen,
Approximation Theory and Approximation Practice, Thm 19.3; the constant below is that
theorem's 64/15 multiplied by a safety factor of 2),

    |E_n|  <=  C * M_rho / ( rho^{2n} (1 - rho^{-2}) ) ,     C = 128/15 .

`M_rho` is computed *rigorously without sampling*: the ellipse is inscribed in the
rectangle [m - h cosh t, m + h cosh t] x [-h sinh t, h sinh t] (rho = e^t), and the
integrand is evaluated ONCE on that rectangle in `acb` ball arithmetic.  Interval
arithmetic is inclusion-monotone, so the resulting magnitude bounds every value on the
whole rectangle, hence on the ellipse.  Each factor of G is bounded by
max(|Phi(s_lo x)|, |Phi(s_hi x)|) - valid for both the box-minimum and the box-maximum
form, since each picks one of the two branches.  The float64 node/weight perturbation is NOT
covered by an eps*M term (that claim had no Lipschitz bound behind it and was removed - pass-3
finding #20c): each node and weight carries a tiny Arb radius (NODE_RAD / WEIGHT_RAD), so the
interval arithmetic propagates the perturbation rigorously.

Tails.  Panels cover [-z_max, z_max] with z_max = 9 (mass outside 2.2e-19), and the tails
are added explicitly using monotonicity of G in z:

    lower tail  in [0, G_max(-z_max) Phi(-z_max)]
    upper tail  in [G_min(z_max) (1 - Phi(z_max)), 1 - Phi(z_max)].

Directions (SPEC section 6): `c_upper_box_gl` is the smallest c whose RIGOROUS lower
enclosure of min_lambda F reaches 1 - alpha (an upper bound on the supremum);
`c_lower_s_gl` is the largest c whose RIGOROUS upper enclosure at a witness point is
below 1 - alpha (a lower bound on c* at that point, hence on the supremum).

Arbitrary precision: `prec` (bits) is a parameter; every endpoint reported is an Arb
endpoint, so a certified claim concerns exactly the float64 inputs and the reported c.
"""

from __future__ import annotations

import contextlib
import math

import enclose
import numpy as np
import paths
from enclose import Z_MAX, _g_at_node, _hi, _lo
from paths import NAN, _f, _i

ERRLOG: list[str] = []

GL_N = paths._env_int("WA2H_GL_N", 16)
GL_SAFETY = paths._env_float("WA2H_GL_SAFETY", 2.0)  # factor on the theorem constant
# GL_EPS was the old (unfounded) node/weight rounding constant; it is gone - the perturbation is
# now carried by NODE_RAD / WEIGHT_RAD inside the interval arithmetic (pass-3 finding #20c).
# slack for the float64 summation of ~50 panel errors and tail terms (n * eps * scale)
GL_SUM_SLACK = 1e-14
# Absolute radii that carry the float64 node/weight perturbation INSIDE the interval
# arithmetic (so no derivative/Lipschitz bound is needed - verifier finding: the old
# 2*2*eps*M term was unfounded and breached 13x at a sharp-corner panel).
NODE_RAD = paths._env_float("WA2H_NODE_RAD", 1e-15)
WEIGHT_RAD = paths._env_float("WA2H_WEIGHT_RAD", 4e-15)
T_GRID = (0.25, 0.4, 0.6, 0.8, 1.0, 1.3, 1.7, 2.2, 2.8)

# Panel layout: finer where the transition layer of G lives (z ~ -c d sqrt(n0), i.e.
# roughly [-6, 0] for desk designs), coarser in the tails where f is tiny.
EDGES = (
    -9.0,
    -6.5,
    -5.0,
    -4.0,
    -3.3,
    -2.7,
    -2.2,
    -1.8,
    -1.45,
    -1.15,
    -0.9,
    -0.68,
    -0.48,
    -0.3,
    -0.14,
    0.0,
    0.14,
    0.3,
    0.48,
    0.68,
    0.9,
    1.15,
    1.45,
    1.8,
    2.2,
    2.7,
    3.3,
    4.0,
    5.0,
    6.5,
    9.0,
)


def panel_edges(z_max=None, edges=None):
    """Ascending panel edges covering [-z_max, z_max]."""
    if edges is not None:
        return np.asarray(edges, float)
    zm = Z_MAX if z_max is None else _f(z_max, Z_MAX)
    e = np.asarray(EDGES, float)
    mask = (e > -zm) & (e < zm)
    return np.concatenate([[-zm], e[mask], [zm]])


def layer_edges(
    c_ref,
    d_row,
    n0,
    s_hi,
    z_max=None,
    per_layer=10,
    spread=4.0,
    edges=None,
    max_width=None,
    window=6.0,
    max_edges=800,
):
    """Panel edges that RESOLVE the integrand's transition layer and its kink.

    Two verifier-driven additions to the base layout:

    * KINK EDGES.  For a box the integrand is only piecewise analytic: `Psi_min`/`Psi_max` switch
      branch at `x_i = 0`, i.e. at `z* = -c d_i sqrt(n0)`.  An EVEN `per_layer` places an edge
      exactly at each `z*`, so every panel lies wholly on one side and the Bernstein-ellipse
      remainder applies panel by panel.
    * S-AWARE SUBDIVISION.  The ellipse bound is tight only while `s * (panel half-width) *
      sinh(t)` is moderate; a panel wider than about `8/s` cannot be bounded tightly for ANY `t`
      (large `t` explodes `M`, small `t` kills the decay).  Panels wider than that are subdivided
      where the mass lives.  Measured: at the tiny-lambda corner (`s ~ 109`, lambda ~ 1e-3) the
      bracket width fell from O(10) to the target range, at ~165 extra edges and ~50 ms per
      enclosure evaluation.
    """
    e = list(panel_edges(z_max, edges))
    zm = _f(z_max, Z_MAX)
    sq = math.sqrt(_f(n0, 1.0))
    d_row = np.asarray(d_row, float).reshape(-1)
    per_layer = _i(per_layer, 10) + (_i(per_layer, 10) % 2)  # EVEN: k = per/2 gives z*
    # NOTE (pass-4 finding 4c): the planted edge uses float64 math.sqrt(n0) while the Arb
    # integrand uses arb(n0).sqrt(), so an edge can miss the true kink by up to ~9e-16.  The
    # resulting sliver panel is bounded by its normal mass and the residual is far below
    # GL_SUM_SLACK = 1e-14, but "panel wholly on one side" holds only to that tolerance.
    # iterate the LONGER of the two vectors so a d_row entry beyond len(s_hi) cannot silently
    # skip its kink (pass-4 finding 4b)
    s_hi_arr = np.asarray(s_hi, float)
    for i in range(max(s_hi_arr.size, d_row.size)):
        sh = _f(s_hi_arr[i], 0.0) if i < s_hi_arr.size else 0.0
        if sh <= 0.0:
            continue
        z_star = -_f(c_ref) * _f(d_row[i]) * sq
        w = 1.0 / sh
        for k in range(per_layer + 1):
            zz = z_star + (2.0 * k / per_layer - 1.0) * _f(spread, 4.0) * w
            if -zm < zz < zm:
                e.append(zz)
    s_max = _f(np.max(np.asarray(s_hi, float))) if np.size(s_hi) else 0.0
    if s_max > 0.0:
        mw = _f(max_width) if max_width is not None else 8.0 / s_max
        if mw > 0.0:
            base = sorted(e)
            extra: list = []
            for a, b in zip(base[:-1], base[1:], strict=False):
                if (b - a) <= mw or b < -_f(window, 6.0) or a > _f(window, 6.0):
                    continue
                k = _i(math.ceil((b - a) / mw))
                if k > 1:
                    extra.extend(np.linspace(a, b, k + 1)[1:-1])
                if len(extra) > _i(max_edges, 800):
                    break
            e.extend(extra)
    return np.unique(np.asarray(e, float))


def gl_nodes(n=None):
    """Gauss-Legendre nodes and weights on [-1, 1] (float64)."""
    k = max(2, _i(n, GL_N))
    x, w = np.polynomial.legendre.leggauss(k)
    return x, w


def _sqrt2pi_a():
    from flint import arb

    return (2 * arb.pi()).sqrt()


def Phi_acb(a, prec=None):
    """Rigorous enclosure of Phi(a) for a complex ball argument."""
    from flint import acb

    if prec is not None:
        enclose._set_prec(prec)
    return (1 + (a / acb(2).sqrt()).erf()) / 2


def rect_bound(z_mid, x_half, y_half, c_a, d_a, s_lo_a, s_hi_a, sq_n0_a, prec):
    """Rigorous upper bound on sup |G(z) phi(z)| over the rectangle (both G variants).

    One `acb` evaluation; every factor is bounded by the max of its two branches.
    """
    from flint import acb, arb

    zr = acb(arb(z_mid, max(0.0, _f(x_half, 0.0))), arb(0.0, max(0.0, _f(y_half, 0.0))))
    u = zr / sq_n0_a
    acc = arb(1)
    for i in range(len(s_lo_a)):
        x = u + c_a * d_a[i]
        fa = abs(Phi_acb(s_lo_a[i] * x, prec))
        fb = abs(Phi_acb(s_hi_a[i] * x, prec))
        m = max(_f(fa.upper(), 0.0), _f(fb.upper(), 0.0))
        acc = acc * arb(m)
    phi_mag = abs(((-zr * zr) / 2).exp()).upper() / _sqrt2pi_a()
    return _f(acc.upper(), NAN) * _f(phi_mag, NAN)


def panel_error(m, h, c_a, d_a, s_lo_a, s_hi_a, sq_n0_a, prec, n, edges_ctx=None):
    """Rigorous error bound for one panel: the analyticity (Bernstein-ellipse) term only.

    The analyticity term minimises `C * M_rho / (rho^{2n} (1 - rho^-2))` over the t grid.  The
    float64 node/weight perturbation is NOT a term here: it is carried inside the interval
    arithmetic by the NODE_RAD / WEIGHT_RAD radii applied at the node evaluations (pass-3 finding
    #20c removed the earlier unfounded `2*2*eps*M` claim, whose parenthetical could not yield an
    eps-scaled bound).  The panel's own mass bound `int_panel G phi dz <= Phi(b) - Phi(a)` caps the
    analyticity term, which is what stops far-tail panels - whose ellipse bound explodes because
    `s * (imaginary extent)` is large - from contributing their full width.
    """
    best = NAN
    m_best = 0.0
    for t in T_GRID:
        xh = _f(h) * math.cosh(t)
        yh = _f(h) * math.sinh(t)
        if not (math.isfinite(xh) and math.isfinite(yh)):
            continue
        M = rect_bound(m, xh, yh, c_a, d_a, s_lo_a, s_hi_a, sq_n0_a, prec)
        if M != M:  # NaN -> unusable
            continue
        rho = math.exp(t)
        decay = rho ** (2 * n)
        if decay <= 0.0:
            continue
        err = (128.0 / 15.0) * GL_SAFETY * M / (decay * (1.0 - rho**-2.0))
        val = _f(h) * err
        if best != best or val < best:
            best = val
            m_best = M
    if best != best:
        best = _f(h)  # fallback: the panel's own mass bound
    # The float64 node/weight perturbation is carried INSIDE the interval arithmetic: the node and
    # weight evaluations above use NODE_RAD / WEIGHT_RAD radii, so no separate eps*M term is needed
    # (or available - the earlier "PROVEN" claim had no Lipschitz bound behind it and was removed,
    # pass-3 finding #20c).  WEIGHT_RAD = 4e-15 covers numpy's measured weight error (max |dw| =
    # 6.0e-16 at n=6, 3.3e-15 at n=48); NODE_RAD = 1e-15 covers the node error (<= 7.8e-17).

    # The float64 node/weight perturbation is carried by NODE_RAD/WEIGHT_RAD
    # inside the Arb sums above, so no separate (and unfounded) term is needed.
    round_err = 0.0
    # Rigorous cap on the panel's integral: 0 <= G <= 1, so int_panel G phi <= int_panel phi =
    # Phi(b) - Phi(a), the panel's NORMAL MASS.  Using the panel width h instead (the first
    # version) let far-tail panels - whose ellipse bound explodes because s * (imaginary extent)
    # is large - contribute their full width, which is what kept Delta p at ~1e-2 on the tiny-
    # lambda region corners even after the transition layers were planted.
    cap = _f(h)
    try:
        from flint import arb

        pa = enclose.Phi_arb(arb(_f(m - h)), prec)
        pb = enclose.Phi_arb(arb(_f(m + h)), prec)
        mass_cap = max(0.0, _hi(pb - pa, _f(h)))
        cap = min(cap, mass_cap)
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: panel mass cap unavailable: {exc!r}")
    return min(best, cap) + round_err, m_best


def settings(n=None, prec=64, edges=None):
    """The certificate-affecting constants, recorded in every certificate output.

    Verifier finding N2-b: these constants were environment-overridable and were not
    recorded anywhere, so a certificate could not be reproduced from its own deposit.
    """
    e = panel_edges(edges=edges)
    return {
        "gl_n": _i(n, GL_N),
        "gl_safety": _f(GL_SAFETY),
        "node_rad": _f(NODE_RAD),
        "weight_rad": _f(WEIGHT_RAD),
        "prec": _i(prec, 64),
        "z_max": _f(Z_MAX),
        "n_panels": _i(e.size - 1),
        "edges_first_last": [_f(e[0]), _f(e[-1])],
        "t_grid": [round(_f(t), 3) for t in T_GRID],
        "mass_floor": None,
        "note": "GL certificate path: no mass floor; the ellipse remainder bounds the "
        "quadrature error; the float64 node/weight perturbation is carried by NODE_RAD/WEIGHT_RAD "
        "radii inside the interval arithmetic",
    }


def bounds_gl_arb(
    c,
    d_row,
    n0,
    s_lo,
    s_hi,
    prec=64,
    n=None,
    edges=None,
    need=("lo", "hi"),
    return_diag=False,
    plant_layer=True,
    max_width=None,
):
    """Rigorous (lo, hi): lo <= min_lambda F(c,lambda), max_lambda F(c,lambda) <= hi.

    Gauss-Legendre per panel with the Bernstein-ellipse remainder above; explicit tails.
    """
    from flint import arb

    lo_out = NAN
    hi_out = NAN
    diag = {
        "n_panels": 0,
        "n_gl": _i(n, GL_N),
        "err_total": 0.0,
        "tail_terms": {},
        "panels": [],
        "settings": settings(n=n, prec=prec, edges=edges),
    }
    try:
        enclose._set_prec(prec)
        d_row = np.asarray(d_row, float).reshape(-1)
        s_lo = np.asarray(s_lo, float)
        s_hi = np.asarray(s_hi, float)
        c_a = arb(_f(c))
        n0a = arb(_f(n0, 1.0))
        sq_n0_a = n0a.sqrt()
        d_a = [arb(_f(v)) for v in d_row]
        s_lo_a = [arb(_f(v)) for v in s_lo]
        s_hi_a = [arb(_f(v)) for v in s_hi]
        xs, ws = gl_nodes(n)
        k = xs.size
        e = (
            layer_edges(
                c_ref=c,
                d_row=d_row,
                n0=n0,
                s_hi=s_hi,
                edges=edges,
                max_width=max_width,
            )
            if plant_layer
            else panel_edges(edges=edges)
        )
        want_lo = "lo" in need
        want_hi = "hi" in need
        acc_lo = arb(0)
        acc_hi = arb(0)
        err_total = 0.0
        for j in range(e.size - 1):
            a = _f(e[j])
            b = _f(e[j + 1])
            m = 0.5 * (a + b)
            h = 0.5 * (b - a)
            if h <= 0.0:
                continue
            s_lo_sum = arb(0)
            s_hi_sum = arb(0)
            for i in range(k):
                zi = arb(m) + arb(h) * arb(_f(xs[i], 0.0), NODE_RAD)
                g_lo, g_hi = _g_at_node(zi, c_a, d_a, s_lo_a, s_hi_a, sq_n0_a, prec)
                phi = ((-zi * zi) / 2).exp() / _sqrt2pi_a()
                wi = arb(_f(ws[i], 0.0), WEIGHT_RAD)
                if want_lo:
                    s_lo_sum += wi * g_lo * phi
                if want_hi:
                    s_hi_sum += wi * g_hi * phi
            err, m_b = panel_error(m, h, c_a, d_a, s_lo_a, s_hi_a, sq_n0_a, prec, k)
            err_total += _f(err)
            diag["panels"].append({"a": a, "b": b, "err": _f(err), "M": _f(m_b)})
            if want_lo:
                acc_lo += arb(h) * s_lo_sum - arb(_f(err))
            if want_hi:
                acc_hi += arb(h) * s_hi_sum + arb(_f(err))
        diag["n_panels"] = _i(len(diag["panels"]))
        diag["layer_used"] = bool(plant_layer)
        # tails
        z0 = arb(_f(e[0], -Z_MAX))
        zM = arb(_f(e[-1], Z_MAX))
        g_lo0, g_hi0 = _g_at_node(z0, c_a, d_a, s_lo_a, s_hi_a, sq_n0_a, prec)
        g_loM, g_hiM = _g_at_node(zM, c_a, d_a, s_lo_a, s_hi_a, sq_n0_a, prec)
        phi0 = enclose.Phi_arb(z0, prec)
        phiM = enclose.Phi_arb(zM, prec)
        up_mass = 1 - phiM  # (1 - Phi(z_max))
        low_mass = phi0  # Phi(-z_max)
        if want_lo:
            acc_lo += g_loM * up_mass
        if want_hi:
            acc_hi += g_hi0 * low_mass
            acc_hi += up_mass
        diag["tail_terms"] = {
            "low_mass": _f(low_mass.upper(), 0.0),
            "up_mass": _f(up_mass.upper(), 0.0),
        }
        diag["err_total"] = _f(err_total)
        if want_lo:
            lo_out = _lo(acc_lo, NAN) - GL_SUM_SLACK
        if want_hi:
            hi_out = _hi(acc_hi, NAN) + GL_SUM_SLACK
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: {exc!r}")
    if return_diag:
        return lo_out, hi_out, diag
    return lo_out, hi_out


def c_upper_box_gl(
    d_row,
    n0,
    s_lo,
    s_hi,
    alpha=paths.ALPHA0,
    c_guess=None,
    tol=1e-6,
    prec=64,
    n=None,
    grow=40,
    widen=0.5,
    max_width=None,
):
    """Certified upper bound on sup_{lambda in box} c*(lambda), GL certificate path."""
    level = 1.0 - _f(alpha, paths.ALPHA0)
    out = {
        "c_ub": NAN,
        "n_eval": 0,
        "f_lo": NAN,
        "f_hi": NAN,
        "bracket_width": NAN,
        "warm_start": NAN,
        "err_total": NAN,
        "settings": settings(n=n, prec=prec),
    }

    def ok(c, mw=None):
        lo, _hi = bounds_gl_arb(
            c,
            d_row,
            n0,
            s_lo,
            s_hi,
            prec=prec,
            n=n,
            need=("lo",),
            plant_layer=True,
            max_width=mw,
        )
        out["n_eval"] += 1
        return (lo == lo) and lo >= level

    COARSE = _f(paths._env_float("WA2H_COARSE_WIDTH", 1.0))
    try:
        # If the coarse layout IS the fine layout (no s-aware subdivision needed at this c), the
        # two-phase search would only repeat work - so skip phase 2 in that case.  Measured: this
        # keeps ordinary cells at ~1.8 s while sharp corners drop from ~5 s to ~1.4 s.
        try:
            _ec = layer_edges(
                c_ref=_f(c_guess, 3.0), d_row=d_row, n0=n0, s_hi=s_hi, max_width=COARSE
            )
            _ef = layer_edges(
                c_ref=_f(c_guess, 3.0), d_row=d_row, n0=n0, s_hi=s_hi, max_width=None
            )
            if _ec.size == _ef.size:
                COARSE = None
        except Exception as exc:  # noqa: BLE001
            ERRLOG.append(f"{__name__}: coarse-layout probe failed: {exc!r}")
        if c_guess is None:
            zf = np.linspace(-Z_MAX, Z_MAX, 401)
            c_guess = _f(
                enclose.c_ub_float(d_row, n0, s_lo, s_hi, zf, alpha=alpha), 3.0
            )
            if c_guess != c_guess:
                c_guess = 3.0
        out["warm_start"] = _f(c_guess)
        c_lo = max(1e-4, _f(c_guess) - _f(widen, 0.5))
        c_hi = max(c_lo + 0.1, _f(c_guess) + _f(widen, 0.5))
        it = 0
        span = max(0.5, _f(widen, 0.5))
        while not ok(c_hi):
            span *= 2.0
            c_hi = c_lo + span
            it += 1
            if it > _i(grow, 80) or c_hi > enclose.C_MAX:
                out["error"] = "no bracket within C_MAX"
                return out
        while ok(c_lo) and c_lo > 1e-4:
            c_lo = max(1e-4, c_lo - _f(widen, 0.5))
            it += 1
            if it > _i(grow, 40):
                break
        # PHASE 1: bisect with the COARSE layout (fast; the kink edges are still planted, only the
        # s-aware subdivision is skipped).  A coarse bound is a valid but looser upper bound, so it
        # can only overshoot c* - it never understates it.
        tgt = _f(tol, 1e-6)
        while c_hi - c_lo > tgt:
            mid = 0.5 * (c_lo + c_hi)
            if ok(mid, COARSE):
                c_hi = mid
            else:
                c_lo = mid
        # PHASE 2: certify with the FINE layout.  If the fine check fails, step c_hi up until it
        # holds, then bisect finely - so every returned value is fine-verified (valid AND tight).
        step = 0
        while COARSE is not None and not ok(c_hi, None):
            c_lo = c_hi
            c_hi = c_hi * 1.02 + 1e-3
            step += 1
            if step > _i(grow, 80) or c_hi > enclose.C_MAX:
                out["error"] = "fine certification failed"
                return out
        while COARSE is not None and c_hi - c_lo > tgt:
            mid = 0.5 * (c_lo + c_hi)
            if ok(mid, None):
                c_hi = mid
            else:
                c_lo = mid
        out["bracket_width"] = c_hi - c_lo
        flo, fhi, diag = bounds_gl_arb(
            c_hi,
            d_row,
            n0,
            s_lo,
            s_hi,
            prec=prec,
            n=n,
            need=("lo", "hi"),
            return_diag=True,
            plant_layer=True,
            max_width=None,
        )
        out.update(
            {
                "c_ub": _f(c_hi),
                "f_lo": flo,
                "f_hi": fhi,
                "err_total": diag.get("err_total"),
            }
        )
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: {exc!r}")
    return out


def c_lower_s_gl(
    s_pt,
    d_row,
    n0,
    alpha=paths.ALPHA0,
    c_guess=None,
    tol=1e-6,
    prec=64,
    n=None,
    grow=40,
    widen=0.5,
    max_width=None,
):
    """Certified lower bound on c*(lambda) at a witness point (GL certificate path)."""
    level = 1.0 - _f(alpha, paths.ALPHA0)
    s_arr = np.asarray(s_pt, float)
    out = {
        "c_lb": NAN,
        "n_eval": 0,
        "f_lo": NAN,
        "f_hi": NAN,
        "bracket_width": NAN,
        "warm_start": NAN,
        "err_total": NAN,
        "settings": settings(n=n, prec=prec),
    }

    def below(c, mw=None):
        _lo, hi = bounds_gl_arb(
            c,
            d_row,
            n0,
            s_arr,
            s_arr,
            prec=prec,
            n=n,
            need=("hi",),
            plant_layer=True,
            max_width=mw,
        )
        out["n_eval"] += 1
        return (hi == hi) and hi < level

    COARSE = _f(paths._env_float("WA2H_COARSE_WIDTH", 1.0))
    try:
        # If the coarse layout IS the fine layout (no s-aware subdivision needed at this c), the
        # two-phase search would only repeat work - so skip phase 2 in that case.  Measured: this
        # keeps ordinary cells at ~1.8 s while sharp corners drop from ~5 s to ~1.4 s.
        try:
            _ec = layer_edges(
                c_ref=_f(c_guess, 3.0), d_row=d_row, n0=n0, s_hi=s_arr, max_width=COARSE
            )
            _ef = layer_edges(
                c_ref=_f(c_guess, 3.0), d_row=d_row, n0=n0, s_hi=s_arr, max_width=None
            )
            if _ec.size == _ef.size:
                COARSE = None
        except Exception as exc:  # noqa: BLE001
            ERRLOG.append(f"{__name__}: coarse-layout probe failed: {exc!r}")
        if c_guess is None:
            zf = np.linspace(-Z_MAX, Z_MAX, 401)
            c_guess = _f(enclose.c_lb_float(s_arr, d_row, n0, zf, alpha=alpha), 3.0)
            if c_guess != c_guess:
                c_guess = 3.0
        out["warm_start"] = _f(c_guess)
        c_lo = max(1e-4, _f(c_guess) - _f(widen, 0.5))
        c_hi = max(c_lo + 0.1, _f(c_guess) + _f(widen, 0.5))
        it = 0
        span = max(0.5, _f(widen, 0.5))
        while below(c_hi):
            span *= 2.0
            c_hi = c_lo + span
            it += 1
            if it > _i(grow, 80) or c_hi > enclose.C_MAX:
                out["error"] = "no bracket within C_MAX"
                return out
        while below(c_lo) and c_lo > 1e-4:
            c_lo = max(1e-4, c_lo - _f(widen, 0.5))
            it += 1
            if it > _i(grow, 40):
                break
        # PHASE 1: bisect with the COARSE layout.  For a LOWER bound a looser F_hi makes `below`
        # harder, so the coarse phase undershoots c_lb - conservative, the safe direction.
        tgt = _f(tol, 1e-6)
        while c_hi - c_lo > tgt:
            mid = 0.5 * (c_lo + c_hi)
            if below(mid, COARSE):
                c_lo = mid
            else:
                c_hi = mid
        # PHASE 2: refine with the FINE layout, pushing c_lo up while the fine check holds, so the
        # returned c_lb is fine-verified (valid AND tight).
        hi_fine = c_hi
        step = 0
        while COARSE is not None and below(hi_fine, None) and hi_fine < enclose.C_MAX:
            c_lo = hi_fine
            hi_fine = hi_fine * 1.02 + 1e-3
            step += 1
            if step > _i(grow, 80):
                break
        while COARSE is not None and hi_fine - c_lo > tgt:
            mid = 0.5 * (c_lo + hi_fine)
            if below(mid, None):
                c_lo = mid
            else:
                hi_fine = mid
        out["bracket_width"] = hi_fine - c_lo
        flo, fhi, diag = bounds_gl_arb(
            c_lo,
            d_row,
            n0,
            s_arr,
            s_arr,
            prec=prec,
            n=n,
            need=("lo", "hi"),
            return_diag=True,
            plant_layer=True,
            max_width=None,
        )
        out.update(
            {
                "c_lb": _f(c_lo),
                "f_lo": flo,
                "f_hi": fhi,
                "err_total": diag.get("err_total"),
            }
        )
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: {exc!r}")
    return out


def bracket_point_gl(
    s_pt,
    d_row,
    n0,
    alpha=paths.ALPHA0,
    c_guess=None,
    tol=1e-6,
    prec=64,
    n=None,
    slope_h=0.02,
):
    """SPEC section 6 root certificate at one lambda, GL path, with Delta p and slope."""
    res = {
        "c_L": NAN,
        "c_U": NAN,
        "width": NAN,
        "dp": NAN,
        "slope_hi": NAN,
        "c_equiv": NAN,
        "n_eval": 0,
        "warm_start": NAN,
        "settings": settings(n=n, prec=prec),
    }
    try:
        up = c_upper_box_gl(
            d_row,
            n0,
            np.asarray(s_pt, float),
            np.asarray(s_pt, float),
            alpha=alpha,
            c_guess=c_guess,
            tol=tol,
            prec=prec,
            n=n,
        )
        lo = c_lower_s_gl(
            s_pt,
            d_row,
            n0,
            alpha=alpha,
            c_guess=up.get("warm_start"),
            tol=tol,
            prec=prec,
            n=n,
        )
        res["n_eval"] = _i(up.get("n_eval"), 0) + _i(lo.get("n_eval"), 0)
        res["warm_start"] = up.get("warm_start")
        res.update({"c_U": up.get("c_ub"), "c_L": lo.get("c_lb")})
        res["width"] = _f(res["c_U"]) - _f(res["c_L"])
        dp_up = _f(up.get("f_hi")) - _f(up.get("f_lo"))
        dp_lo = _f(lo.get("f_hi")) - _f(lo.get("f_lo"))
        res["dp"] = _f(np.nanmax([dp_up, dp_lo]))
        s_arr = np.asarray(s_pt, float)
        c_ref = _f(res["c_U"], 3.0)
        f_p, _ = bounds_gl_arb(
            c_ref + slope_h,
            d_row,
            n0,
            s_arr,
            s_arr,
            prec=prec,
            n=n,
            need=("lo",),
            plant_layer=True,
        )
        _l, f_m = bounds_gl_arb(
            c_ref - slope_h,
            d_row,
            n0,
            s_arr,
            s_arr,
            prec=prec,
            n=n,
            need=("hi",),
            plant_layer=True,
        )
        slope = max(0.0, (_f(f_p) - _f(f_m)) / (2.0 * _f(slope_h, 0.02)))
        res["slope_hi"] = slope
        res["c_equiv"] = (_f(res["dp"]) / slope) if slope > 1e-12 else NAN
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: {exc!r}")
    return res


def quickcheck(alpha=paths.ALPHA0, prec=64, seed=20260918):
    """One fast certified containment check against the independent solver.

    Used by `run_modes.quick`; the full battery lives in `test_enclose.py`.
    """
    import importlib

    core = importlib.import_module("core")
    nvec = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    s0v, siv, _lh = core.draw_var(nvec, lam, 1, seed)
    d_row = np.sqrt(siv[0] / nvec[1:] + s0v[0] / nvec[0])
    s_pt = np.sqrt(nvec[1:] / lam)
    c = 2.6
    lo, hi = bounds_gl_arb(c, d_row, nvec[0], s_pt, s_pt, prec=prec)
    ref = _f(
        core.cond_noreject(np.array([c]), np.array([s0v[0]]), siv[0:1], nvec, lam)[0]
    )
    return {
        "contained": bool(lo <= ref <= hi),
        "F_lo": lo,
        "F_hi": hi,
        "F_ref": ref,
        "width": _f(hi) - _f(lo),
        "c": _f(c),
        "alpha0": _f(alpha),
    }


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------
def self_test(alpha=paths.ALPHA0, prec=64, n=None, seed=20260918, n_pts=8, report=True):
    """Containment (vs the independent solver), cross-check vs the Riemann path,
    the population/conditional identity, and the width/cost of the certificate."""
    import importlib

    _core = importlib.import_module("core")
    cond_noreject = _core.cond_noreject
    draw_var = _core.draw_var
    oracle_c = _core.oracle_c

    rng = np.random.default_rng(seed)
    nvec = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    s0v, siv, _lh = draw_var(nvec, lam, 1, seed)
    d_row = np.sqrt(siv[0] / nvec[1:] + s0v[0] / nvec[0])
    s_pt = np.sqrt(nvec[1:] / lam)
    out: dict = {
        "n": list(nvec),
        "lam": list(lam),
        "d_row": list(d_row),
        "alpha": _f(alpha),
        "gl_n": _i(n, GL_N),
    }

    # (a) containment on a point box and on random boxes, against cond_noreject
    rows = []
    for c in (2.0, 2.5, 2.75, 3.0):
        lo, hi, diag = bounds_gl_arb(
            c, d_row, nvec[0], s_pt, s_pt, prec=prec, n=n, return_diag=True
        )
        ref = _f(
            cond_noreject(np.array([c]), np.array([s0v[0]]), siv[0:1], nvec, lam)[0]
        )
        rows.append(
            {
                "c": _f(c),
                "F_lo": lo,
                "F_hi": hi,
                "F_ref": ref,
                "width": _f(hi) - _f(lo),
                "inside": bool(lo <= ref <= hi),
                "err_total": diag.get("err_total"),
            }
        )
    out["point_containment"] = rows
    for _ in range(max(1, _i(n_pts, 8))):
        fac = 1.0 + 0.4 * rng.random(2)
        lo_f, hi_f = lam / fac, lam * fac
        s_lo, s_hi = enclose.s_range(nvec, lo_f, hi_f)
        c = 2.0 + rng.random()
        lo, hi = bounds_gl_arb(c, d_row, nvec[0], s_lo, s_hi, prec=prec, n=n)
        fmin, fmax = 1.0, 0.0
        for a in (lo_f[0], hi_f[0]):
            for b in (lo_f[1], hi_f[1]):
                v = _f(
                    cond_noreject(
                        np.array([c]),
                        np.array([s0v[0]]),
                        siv[0:1],
                        nvec,
                        np.array([a, b]),
                    )[0]
                )
                fmin, fmax = min(fmin, v), max(fmax, v)
        out.setdefault("box_containment", []).append(
            {
                "c": _f(c),
                "box": [list(lo_f), list(hi_f)],
                "F_lo": lo,
                "F_hi": hi,
                "min_corner": fmin,
                "max_corner": fmax,
                "width": _f(hi) - _f(lo),
                "lo_ok": bool(lo <= fmin),
                "hi_ok": bool(hi >= fmax),
                "inside": bool(lo <= fmin and hi >= fmax),
            }
        )

    # (b) cross-check against the Riemann-Arb path (agreement within combined widths)
    c = 2.75
    prep = enclose.prepare_partition(enclose.base_partition(nodes=20000), prec=prec)
    r_lo, r_hi = enclose.bounds_arb(c, d_row, nvec[0], s_pt, s_pt, prep)
    g_lo, g_hi = bounds_gl_arb(c, d_row, nvec[0], s_pt, s_pt, prec=prec, n=n)
    out["cross_check_riemann"] = {
        "riemann": [r_lo, r_hi],
        "riemann_width": _f(r_hi) - _f(r_lo),
        "gl": [g_lo, g_hi],
        "gl_width": _f(g_hi) - _f(g_lo),
        "agree": bool(r_lo <= g_hi and g_lo <= r_hi),
        "width_ratio": (
            (_f(r_hi) - _f(r_lo)) / (_f(g_hi) - _f(g_lo))
            if (_f(g_hi) - _f(g_lo)) > 0
            else NAN
        ),
    }

    # (c) population/conditional identity: the draw-average at the oracle c is 1 - alpha
    with contextlib.suppress(Exception):
        c_or = oracle_c(nvec, lam, _f(alpha))
        D = 4000
        s0b, sib, _ = draw_var(nvec, lam, D, seed + 7)
        p = cond_noreject(np.full(D, c_or), s0b, sib, nvec, lam)
        out["population_identity"] = {
            "c_oracle": c_or,
            "mean_cond_noreject": _f(np.mean(p)),
            "target": 1.0 - _f(alpha),
            "mcse": _f(np.std(p, ddof=1) / np.sqrt(D)),
            "z": _f(
                (np.mean(p) - (1.0 - _f(alpha))) / (np.std(p, ddof=1) / np.sqrt(D))
            ),
            "draws": _i(D),
        }

    # (d) the bracket at the point (and its cost)
    br = bracket_point_gl(s_pt, d_row, nvec[0], alpha=alpha, prec=prec, n=n)
    out["bracket"] = br

    if report:
        print(
            f"GL self_test n={list(nvec)} lam={list(lam)} alpha={_f(alpha)} "
            f"gl_n={_i(n, GL_N)} panels={len(panel_edges()) - 1}"
        )
        for r in out["point_containment"]:
            print(
                f"  c={r['c']:.2f} F in [{r['F_lo']:.12f}, {r['F_hi']:.12f}] "
                f"ref {r['F_ref']:.12f} inside={r['inside']} width {r['width']:.2e} "
                f"err {_f(r['err_total']):.1e}"
            )
        bc = out.get("box_containment", [])
        print(
            f"  random boxes: {sum(1 for r in bc if r['inside'])}/{len(bc)} contained, "
            f"widths {[f'{r[chr(119) + chr(105) + chr(100) + chr(116) + chr(104)]:.1e}' for r in bc[:4]]}"
        )
        cc = out["cross_check_riemann"]
        print(
            f"  vs Riemann: agree={cc['agree']} widths {cc['riemann_width']:.2e} "
            f"(riemann) / {cc['gl_width']:.2e} (gl) ratio {_f(cc['width_ratio']):.1f}x"
        )
        if "population_identity" in out:
            pi = out["population_identity"]
            print(
                f"  population identity: mean {pi['mean_cond_noreject']:.6f} vs "
                f"{pi['target']:.6f} (z={pi['z']:.2f}, {pi['draws']} draws)"
            )
        print(
            f"  bracket [{_f(br['c_L']):.6f}, {_f(br['c_U']):.6f}] width "
            f"{_f(br['width']):.2e} dp {_f(br['dp']):.2e} slope {_f(br['slope_hi']):.4f} "
            f"c_equiv {_f(br['c_equiv']):.2e} evals {br['n_eval']}"
        )
    return out


if __name__ == "__main__":
    import tiers

    try:
        res = self_test()
        okp = all(r["inside"] for r in res["point_containment"])
        okb = all(r["inside"] for r in res.get("box_containment", []))
        okc = bool(res["cross_check_riemann"]["agree"])
        res["verdict"] = "PASS" if (okp and okb and okc) else "FAIL"
        print(
            "RESULT:",
            res["verdict"],
            f"(point {okp}, box {okb}, cross {okc})",
        )
        try:
            tiers.write_deposit(
                res,
                paths.ensure_dirs() + "/enclose_gl_selftest.json",
                "region_sup_certified",
                notes="GL certificate path self-test (Bernstein-ellipse "
                "remainder); cross-checked against the Riemann path",
            )
        except Exception as exc:  # noqa: BLE001
            print("deposit failed:", repr(exc))
    except Exception as exc:  # noqa: BLE001
        print("SELF-TEST ERROR:", repr(exc))
