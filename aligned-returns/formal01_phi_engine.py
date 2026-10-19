"""F-8 CHECKER — the rigorous rate engine for Lean's `Phi`.

    Phi(c; rho1, rho2) = P( N > c * sqrt((rho1*X1 + rho2*X2) / nu) )
    N ~ N(0,1), X1, X2 ~ chi2_35 independent, nu = 35, rho1, rho2 >= 0.

Two INDEPENDENT rigorous routes, each returning an Arb ball whose radius INCLUDES every
error term (quadrature enclosure from acb.integral, truncation, input rounding).

  R-beta.  X1 = W*B, X2 = W*(1-B), W ~ chi2_70, B ~ Beta(35/2, 35/2) independent, so
           rho1*X1 + rho2*X2 = W*m(B), m(B) = rho2 + (rho1-rho2)*B. Conditioning on B,
               P(N > c*sqrt(W*m/35)) = (1/2) P(N^2/(W/70) > 2*m*c^2)
                                     = (1/2) P(F_{1,70} > 2*m*c^2)
                                     = (1/2) P(|T_70| > c*sqrt(2m)) = t70sf(c*sqrt(2m)).
           Hence   Phi = int_0^1 t70sf(c*sqrt(2*m(b))) Beta(db)   — the ISSUED SPEC FORM,
           confirmed by both FK-8 anchors. The alternative E_B[(1/2) t70sf(2*m*c^2)]
           confuses the F_{1,70} tail with the t_70 tail and misses anchor A by ~0.0249.
           (W is NOT 70*T_70; W/70 is the chi2_70/70 denominator of F_{1,70}.)

           Integrated in theta, b = sin^2(theta): Beta density * db = 2 sin^34 cos^34 / B(a,a)
           is ENTIRE (no log, no fractional power), and t70sf is the exact even-nu finite
           sum (holomorphic, sqrt taken with the analytic flag). The old integrand used
           b.log(), (1-b).log(): on the integrator's bounding balls touching theta = 0 or
           pi/2 those are nan, and bisection can never shed an endpoint -> nan always.

  R-Imhof. Phi = (1/2) P(W > 0), W = 35*Z^2 - c^2*rho1*X1 - c^2*rho2*X2, and
               Phi = 1/4 + (1/(2*pi)) int_0^inf Im(phi(t))/t dt,
               phi(t) = (1-70it)^(-1/2) (1+2i c^2 rho1 t)^(-35/2) (1+2i c^2 rho2 t)^(-35/2).
           Im(phi)/t is not holomorphic as written (Im, and a 0/0 at t=0). Its holomorphic
           extension is sinc(theta(t)) * (theta(t)/t) / rho(t) with
               theta(t)/t = sum h_j a_j 2F1(1/2,1;3/2;-4 a_j^2 t^2)      (= atan(2at)/t),
               1/rho(t)   = prod [(1-2i a_j t)(1+2i a_j t)]^(-h_j/4),
           analytic on the strip |Im t| < 1/(2 max|a_j|); the integrand checks that strip
           itself when `analytic` is set. |phi(t)| <= prod (2|a_j| t)^(-h_j/2), so the
           integrand decays like t^(-1-p), p = sum_{a_j != 0} h_j/2 (35.5, or 18 when one
           rho is 0) — NOT t^(-71.5). The tail K*T^(-p)/(p*2*pi) is added to the radius.

Exact cases: c = 0 or rho1 = rho2 = 0 gives W >= 0 a.s., so Phi = P(N > 0) = 1/2 exactly.
"""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F
from math import comb

from flint import acb, arb, ctx, fmpq

ROOT = pathlib.Path(__file__).resolve().parent
DEP = ROOT / "deposits"
CLAIMS: list[str] = []
BITS = 192
GOAL = 2.0 ** -100          # integrator goal only; the certificate is the returned ball
TAIL_TARGET = F(1, 10 ** 20)


def claim(cid: str, label: str, value) -> None:
    line = f"CLAIM {cid}: {label} = {value}"
    CLAIMS.append(line)
    print(line)


def q_arb(q) -> arb:
    """Exact rational -> arb ball at the CURRENT precision (no float ever enters).

    An existing `arb` is passed straight through. A-5 evaluates Phi at r(corner), which is
    irrational, so rho arrives as a BALL that must be propagated rather than
    point-evaluated; forcing such an input through Fraction raises "argument should be a
    string or a Rational instance". An `acb` is NOT passed through here: every call site in
    this file requires a real ball.
    """
    if isinstance(q, arb):
        return q
    q = F(q)
    return arb(fmpq(q.numerator, q.denominator))


def dyadic(x) -> F:
    m, e = x.man_exp()
    m, e = int(m), int(e)
    return F(m * (1 << e)) if e >= 0 else F(m, 1 << -e)


def real_part(x) -> arb:
    return x.real if isinstance(x, acb) else x


def ball_bounds(x) -> tuple[F, F]:
    """Exact rational lower/upper endpoints of the real part of a ball."""
    r = real_part(x)
    return dyadic(r.lower()), dyadic(r.upper())


def radius(x) -> F:
    return dyadic(real_part(x).rad())


def imag_is_zero_ball(x) -> bool:
    return not isinstance(x, acb) or x.imag.contains(0)


# ------------------------------------------------------------ t_70 survival ---------

NU_T = 70
# P(|T_v| < x) = sin(phi) * sum_{j<v/2} C(2j,j)/4^j cos^{2j}(phi), tan(phi) = x/sqrt(v), v even
R_COEF = [F(comb(2 * j, j), 4 ** j) for j in range(NU_T // 2)]


def t70sf_holo(x, x2, analytic):
    """P(T_70 > x) by the exact even-nu finite sum; x2 must equal x*x.

    Holomorphic wherever 70 + x2 avoids (-inf, 0]; the sqrt carries the analytic flag.
    """
    u = x2 + NU_T
    w = NU_T / u
    s = 0
    for r in reversed(R_COEF):
        s = s * w + q_arb(r)
    return (1 - x / u.sqrt(analytic=analytic) * s) / 2


def t_sf_ref(x2: arb, v: int = NU_T) -> arb:
    """Reference P(T_v > sqrt(x2)), x2 >= 0 real, via the built-in regularized incomplete
    beta. Independent of t70sf_holo; used for the rho1 = rho2 closed form."""
    V = q_arb(v)
    half, hv = q_arb(F(1, 2)), q_arb(F(v, 2))
    if x2 < V:
        z = x2 / (V + x2)
        return (1 - z.beta_lower(half, hv, regularized=True)) / 2
    z = V / (V + x2)
    return z.beta_lower(hv, half, regularized=True) / 2


def _check_inputs(c, rho1, rho2):
    """Accept exact rationals OR Arb/Acb balls.

    A-5 evaluates Phi at r(corner), which is irrational, so rho arrives as a BALL that must
    be propagated rather than point-evaluated. Forcing such an input through Fraction
    raised "argument should be a string or a Rational instance". A ball is passed straight
    through; the sign check is only applied to exact inputs, and a ball whose lower
    endpoint is negative is rejected rather than silently clipped.
    """
    out = []
    for v in (c, rho1, rho2):
        if isinstance(v, arb):
            if v.lower() < 0:
                raise ValueError("a ball input has a negative lower endpoint")
            out.append(v)
        elif isinstance(v, acb):
            # A complex ball has no sign; its validity is the caller's business.
            out.append(v)
        else:
            v = F(v)
            if v < 0:
                raise ValueError("c, rho1, rho2 must be >= 0")
            out.append(v)
    return out[0], out[1], out[2]


def _exact_half(c, rho1, rho2) -> bool:
    """True exactly when Phi is 1/2 without integrating: c = 0, or rho1 = rho2 = 0.

    Only exact zeros count. A ball is NOT tested for equality with zero: a ball is an
    enclosure, and treating a non-degenerate ball as 'zero' would silently return the
    exact-case value for a general input.
    """
    def is_zero(v):
        return (not isinstance(v, (arb, acb))) and F(v) == 0
    return is_zero(c) or (is_zero(rho1) and is_zero(rho2))


# ---------------------------------------------------------------- R-beta route -----

def phi_r_beta(c, rho1, rho2, bits=BITS, rejected_form=False) -> arb:
    """R-beta. `rejected_form=True` integrates E_B[(1/2) t70sf(2*m*c^2)] instead, ONLY to
    deposit how far that candidate formula lands from the anchors; never a certified Phi."""
    c, rho1, rho2 = _check_inputs(c, rho1, rho2)
    with ctx.workprec(bits):
        if _exact_half(c, rho1, rho2):
            return q_arb(F(1, 2))
        C, R1, R2 = q_arb(c), q_arb(rho1), q_arb(rho2)
        c2 = C * C
        csq2 = (2 * c2).sqrt()                     # c*sqrt(2), c > 0
        a = q_arb(F(35, 2))
        inv_B = (2 * a).gamma() / (a.gamma() ** 2)
        sR1, sR2 = R1.sqrt(), R2.sqrt()

        def integrand(th, analytic):
            s, co = th.sin(), th.cos()
            s2, co2 = s * s, co * co
            if rho2 == 0:
                sqrt_m = sR1 * s                  # m = rho1 sin^2: entire square root
                m = R1 * s2
            elif rho1 == 0:
                sqrt_m = sR2 * co
                m = R2 * co2
            else:
                m = R2 * co2 + R1 * s2
                sqrt_m = m.sqrt(analytic=analytic)
            if rejected_form:
                x = 2 * c2 * m
                sf = t70sf_holo(x, x * x, analytic) / 2
            else:
                x = csq2 * sqrt_m
                sf = t70sf_holo(x, 2 * c2 * m, analytic)
            return sf * (2 * inv_B) * (s2 * co2) ** 17

        v = acb.integral(integrand, 0, acb.pi() / 2, rel_tol=GOAL, abs_tol=GOAL)
        if not v.is_finite() or not imag_is_zero_ball(v):
            raise ArithmeticError(f"R-beta: non-finite or non-real enclosure {v}")
        return v.real


# ---------------------------------------------------------------- R-Imhof route -----

def _imhof_terms(c, rho1, rho2):
    """(a_j, h_j) for W = sum a_j chi2_{h_j}; a == 0 terms dropped (they are exactly 1)."""
    c2 = c * c
    terms = [(F(35), 1), (-c2 * rho1, 35), (-c2 * rho2, 35)]
    return [(a, h) for a, h in terms if a != 0]


def imhof_truncation(terms, target=TAIL_TARGET):
    """Rigorous bound on |Phi - truncated Imhof integral| and the dyadic T achieving it.

    |Im phi(t)/t| <= |phi(t)|/t <= K t^(-1-p), K = prod (2|a|)^(-h/2), p = sum h/2.
    So the Phi-tail beyond T is <= K T^(-p) / (2 pi p). Returned as an arb upper bound.
    """
    p = sum(F(h, 2) for _, h in terms)
    logK = arb(0)
    for a, h in terms:
        logK -= q_arb(F(h, 2)) * (2 * q_arb(abs(a))).log()
    T = F(1, 8)
    for _ in range(200):
        tail = (logK - q_arb(p) * q_arb(T).log()).exp() / (2 * arb.pi() * q_arb(p))
        if dyadic(tail.upper()) <= target:
            return T, arb(tail.upper())
        T *= 2
    raise ArithmeticError("Imhof truncation not found")


def phi_r_imhof(c, rho1, rho2, bits=BITS) -> arb:
    c, rho1, rho2 = _check_inputs(c, rho1, rho2)
    with ctx.workprec(bits):
        if _exact_half(c, rho1, rho2):
            return q_arb(F(1, 2))
        terms = _imhof_terms(c, rho1, rho2)
        A = [(q_arb(a), h) for a, h in terms]
        amax = max(abs(a) for a, _ in terms)
        strip = q_arb(1 / (2 * amax))            # branch points of all factors at |Im t| = 1/(2|a|)
        half, one, three_half = q_arb(F(1, 2)), arb(1), q_arb(F(3, 2))
        T, tail = imhof_truncation(terms)

        def integrand(t, analytic):
            if analytic and not (abs(t.imag) < strip):
                return acb("nan")
            th_over_t = 0
            inv_rho = 1
            for a, h in A:
                at2 = 2 * a * t
                th_over_t += h * a * (-(at2 * at2)).hypgeom_2f1(half, one, three_half)
                u = ((1 - 1j * at2).sqrt(analytic=analytic)
                     * (1 + 1j * at2).sqrt(analytic=analytic)).sqrt(analytic=analytic)
                inv_rho *= u ** (-h)
            return (t * th_over_t).sinc() * th_over_t * inv_rho

        v = acb.integral(integrand, 0, q_arb(T), rel_tol=GOAL, abs_tol=GOAL)
        if not v.is_finite() or not imag_is_zero_ball(v):
            raise ArithmeticError(f"R-Imhof: non-finite or non-real enclosure {v}")
        val = v.real / (2 * arb.pi()) + q_arb(F(1, 4))
        return val + arb(0, tail)                # truncation error folded into the ball


# ---------------------------------------------------------------- driver -----------

C0 = F("2.7173398035")
ANCHORS = {
    "A": ((C0, F("0.50504286"), F("0.04674997")), F("0.025")),
    "B": ((C0, F("5.54") / F("10.04"), F(0)), F("0.02562587081")),
}
ANCHOR_TOL = F(1, 10 ** 8)
RAD_BUDGET = F(1, 10 ** 9)


def fmt_rad(r: F) -> str:
    n, d = r.numerator, r.denominator
    return f"{n}/2^{d.bit_length() - 1}" if d & (d - 1) == 0 else f"{n}/{d}"


def within(ball, centre: F, tol: F) -> bool:
    lo, hi = ball_bounds(ball)
    return centre - tol <= lo and hi <= centre + tol


def overlap(x, y) -> bool:
    (a, b), (c, d) = ball_bounds(x), ball_bounds(y)
    return a <= d and c <= b


def report(cid, label, ball):
    r = radius(ball)
    ok = r <= RAD_BUDGET
    claim(cid, label, f"{ball.str(30, radius=True)}  radius={fmt_rad(r)} "
                      f"(~{float(r):.3e}) radius<=1e-9:{ok}")
    return ok


def main() -> int:
    ok = True
    with ctx.workprec(BITS):
        # t70sf finite sum vs built-in incomplete beta (independent implementations)
        for i, x2 in enumerate([F(1, 3), F(29, 4), F(70), F(211)]):
            X2 = q_arb(x2)
            h = t70sf_holo(acb(X2.sqrt()), acb(X2), False).real
            ref = t_sf_ref(X2)
            agree = overlap(h, ref)
            ok &= agree
            claim(f"TSF{i}", f"t70sf(sqrt({x2})) finite-sum vs beta_lower overlap", agree)

    for k, ((c, r1, r2), target) in ANCHORS.items():
        rb = phi_r_beta(c, r1, r2)
        ri = phi_r_imhof(c, r1, r2)
        ok &= report(f"{k}.beta", f"Phi_Rbeta(c0;{r1},{r2})", rb)
        ok &= report(f"{k}.imhof", f"Phi_RImhof(c0;{r1},{r2})", ri)
        for nm, b in (("beta", rb), ("imhof", ri)):
            w = within(b, target, ANCHOR_TOL)
            ok &= w
            claim(f"{k}.{nm}.anchor", f"ball within {target} +/- 1e-8", w)
        o = overlap(rb, ri)
        ok &= o
        claim(f"{k}.overlap", "R-beta and R-Imhof balls overlap", o)
        rj = phi_r_beta(c, r1, r2, rejected_form=True)
        claim(f"{k}.rejected", "E_B[0.5*t70sf(2*m*c0^2)] [rejected form]",
              f"{rj.str(20, radius=True)} within anchor+/-1e-8:{within(rj, target, ANCHOR_TOL)}")
        with ctx.workprec(BITS):
            gap = rb - q_arb(target)
        claim(f"{k}.gap", f"Phi_Rbeta - {target}", gap.str(12, radius=True))

    # exact cases
    for cid, args in (("X1", (C0, 0, 0)), ("X2", (0, F("0.3"), F("0.7")))):
        for nm, fn in (("beta", phi_r_beta), ("imhof", phi_r_imhof)):
            v = fn(*args)
            e = ball_bounds(v) == (F(1, 2), F(1, 2))
            ok &= e
            claim(f"{cid}.{nm}", f"Phi{args} == 1/2 exactly", e)

    # rho1 = rho2 = rho self-test against the closed form t70sf(c sqrt(2 rho))
    for i, rho in enumerate([F("0.3"), F("0.27589641")]):
        with ctx.workprec(BITS):
            closed = t_sf_ref(2 * q_arb(C0) ** 2 * q_arb(rho))
        rb = phi_r_beta(C0, rho, rho)
        ri = phi_r_imhof(C0, rho, rho)
        report(f"S{i}.closed", f"t70sf(c0*sqrt(2*{rho})) [beta_lower]", closed)
        ok &= report(f"S{i}.beta", f"Phi_Rbeta(c0;{rho},{rho})", rb)
        ok &= report(f"S{i}.imhof", f"Phi_RImhof(c0;{rho},{rho})", ri)
        for nm, b in (("beta", rb), ("imhof", ri)):
            o = overlap(b, closed)
            ok &= o
            claim(f"S{i}.{nm}.closed", "overlaps closed form", o)
        with ctx.workprec(BITS):
            wrong = t_sf_ref((2 * q_arb(rho) * q_arb(C0) ** 2) ** 2) / 2
        claim(f"S{i}.altform", "0.5*t70sf(2*rho*c0^2) [rejected form] overlaps beta route",
              overlap(rb, wrong))

    claim("ALL", "every check passed", ok)
    DEP.mkdir(exist_ok=True)
    (DEP / "formal01_phi_claims.txt").write_text("\n".join(CLAIMS) + "\n", encoding="utf-8")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
