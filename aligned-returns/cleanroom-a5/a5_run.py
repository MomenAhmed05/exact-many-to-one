"""A5 clean-room driver: certificate + validations + cross-checks + mutations -> a5_claims.txt.
Decision path: a5core (Fraction + isqrt only).  Floats appear ONLY in validations/cross-checks,
which are labelled as such."""
import math
import os
import random
import sys
import time
from fractions import Fraction as Q

import a5core
from a5core import check, corner_exact, ratio_intervals, phi_bounds, sqrt_lo, sqrt_hi, rd, ru

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "a5_claims.txt")
WITNESS = os.path.join(HERE, "..", "deposits", "f8c_cover_witness-R1.txt")

NU = 35
ALPHA = Q(1, 40)
LC = Q(271731, 100000)
CORNER = (Q(1, 25), Q(-3, 2), Q(3, 2))

claims = []


def C(cid, label, value):
    claims.append(f"CLAIM {cid}: {label} = {value}")
    print(claims[-1])
    sys.stdout.flush()


def fr(q):
    return f"{q.numerator}/{q.denominator}" if q.denominator != 1 else f"{q.numerator}"


def iv(t):
    return f"[{fr(t[0])}, {fr(t[1])}]  (~[{float(t[0]):.17g}, {float(t[1]):.17g}])"


all_valid = True


def V(cid, label, ok, detail):
    global all_valid
    all_valid = all_valid and ok
    C(cid, label, f"{'PASS' if ok else 'FAIL'} ; {detail}")


# ======================= 1. exact corner quantities =======================
ex = corner_exact(*CORNER)
for key in ["aF", "g00", "g11", "g01", "trG", "detG", "discG"]:
    C(f"A5-EX-{key}", f"{key}(corner) exact", fr(ex[key]))
d = ex["discG"]
C("A5-EX-discG-square", "discG is a rational square", str(sqrt_lo(d) == sqrt_hi(d)))
C("A5-EX-sqrt_discG", "sqrt(discG) = sqrt(52929)/50 ; 52929 = 230^2 + 29 (not a square)",
  f"{230**2 + 29 == 52929 and d == Q(52929, 2500)}")
ivs = ratio_intervals(ex)
# rigorous sanity on the sqrt enclosure
sd = ivs["sqrt_discG"]
C("A5-EX-sqrt_encl", "sqrt(discG) enclosure (lo^2 <= discG <= hi^2 checked exactly)",
  f"{iv(sd)} ; check={sd[0]**2 <= d <= sd[1]**2}")
C("A5-EX-mu1", "mu1(corner) = (277/50 + sqrt(52929)/50)/2 = (277+sqrt(52929))/100, enclosure", iv(ivs["mu1"]))
C("A5-EX-mu2", "mu2(corner) = (277-sqrt(52929))/100, enclosure", iv(ivs["mu2"]))
C("A5-EX-r1", "r1(corner) = (277+sqrt(52929))/1004, enclosure", iv(ivs["r1"]))
C("A5-EX-r2", "r2(corner) = (277-sqrt(52929))/1004, enclosure", iv(ivs["r2"]))
C("A5-EX-r1r2-sum", "r1+r2 = trG/aF exact", fr(ex["trG"] / ex["aF"]))
C("A5-EX-r1r2-prod", "r1*r2 = detG/aF^2 exact", fr(ex["detG"] / ex["aF"] ** 2))
C("A5-EX-width", "enclosure widths r1, r2 (upper-lower)",
  f"{float(ivs['r1'][1]-ivs['r1'][0]):.3e}, {float(ivs['r2'][1]-ivs['r2'][0]):.3e}")

# ======================= 2. the certificate =======================
t0 = time.time()
ok, lo, hi, info = check(ALPHA, LC, *CORNER)
C("A5-ROUTE", "route",
  "NegBin(r=35/2,p=rho_lo/rho_hi) mixture of even-d t tails, rewritten via x*sum_k a_k y^k = 1 as "
  "Phi = (x/2) sum_{k>=35} a_k y^k P(J<=k-35), all terms >= 0; exact Fraction + isqrt, outward rounding to 2^-320")
C("A5-PARAM-y", "y = 1/(1+Lc^2 rho_lo/35) enclosure", iv(info["y"]))
C("A5-PARAM-x", "x = sqrt(1-y) enclosure", iv(info["x"]))
C("A5-PARAM-p", "p = rho_lo/rho_hi enclosure", iv(info["p"]))
C("A5-SERIES-terms", "series terms used (k = 35..K)", f"{info['terms']} (K = {info['K_last']}); "
  f"analytic tail bound a_(K+1) y_hi^(K+1)/(1-y_hi) <= {float(info['tail_hi']):.3e}; "
  f"1 - F_lo(K-35) = {float(1-info['F_lo_final']):.3e}")
C("A5-PHI-LO", "certified LOWER bound on Phi(Lc; r(corner)) (exact rational)", fr(lo))
C("A5-PHI-LO-dec", "certified LOWER bound, decimal (floor, 30 digits)",
  f"0.{(lo.numerator * 10**30 // lo.denominator):030d}")
C("A5-PHI-HI", "certified UPPER bound on Phi(Lc; r(corner)) (exact rational)", fr(hi))
C("A5-PHI-HI-dec", "certified UPPER bound, decimal (ceil, 30 digits)",
  f"0.{(-((-hi.numerator * 10**30) // hi.denominator)):030d}")
C("A5-PHI-WIDTH", "width hi - lo", f"{fr(hi-lo)}  (~{float(hi-lo):.6e})")
C("A5-CMP", "exact rational inequality lower_bound > 1/40",
  f"{lo > ALPHA}  ; lower_bound - 1/40 = {fr(lo-ALPHA)}  (~{float(lo-ALPHA):.10e}) ; "
  f"cross-multiplied: 40*num(lo) = {40*lo.numerator} vs den(lo) = {lo.denominator} ; "
  f"40*num > den : {40*lo.numerator > lo.denominator}")
C("A5-CMP-rel", "relative margin (lo - 1/40)/(1/40)", f"~{float((lo-ALPHA)/ALPHA):.6e}")
C("A5-TIME", "certificate wall time (s)", f"{time.time()-t0:.2f}")

# precision re-run: P = 200 bits must give a (wider) enclosure that contains the P = 320 one
_saveP = a5core.P_BITS
a5core.P_BITS = 200
iv160 = ratio_intervals(ex)
lo160, hi160, _ = phi_bounds(LC, iv160["r1"], iv160["r2"], tol_exp=30)
a5core.P_BITS = _saveP
V("A5-VAL-precision", "independent-precision re-run at P=200 bits (tol 1e-30) encloses the P=320 enclosure",
  lo160 <= lo and hi <= hi160, f"P200 = [{float(lo160):.17g}, {float(hi160):.17g}], width {float(hi160-lo160):.3e}")


# ======================= 3. validations =======================
# ---- float helpers for validation only ----
def t_series_tail(d, s):
    """P(t_d > s), even d, via 1/2 - (x/2) sum_{k<d/2} a_k y^k  (float)"""
    y = d / (d + s * s)
    x = s / math.sqrt(d + s * s)
    a, S, yk = 1.0, 0.0, 1.0
    for k in range(d // 2):
        if k > 0:
            a *= (2 * k - 1) / (2 * k)
        S += a * yk
        yk *= y
    return 0.5 - 0.5 * x * S


def t_tail_literal_plus(d, s):
    """The formula exactly as written in the specification: 1/2 + (x/2) sum ...  (tested, not used)"""
    return 1.0 - t_series_tail(d, s)


def t_pdf(d, t):
    lc = math.lgamma((d + 1) / 2) - math.lgamma(d / 2) - 0.5 * math.log(d * math.pi)
    return math.exp(lc - (d + 1) / 2 * math.log1p(t * t / d))


def simpson_tail(d, s, upper=80.0, n=400000):
    h = (upper - s) / n
    tot = t_pdf(d, s) + t_pdf(d, upper)
    for i in range(1, n):
        tot += (4 if i % 2 else 2) * t_pdf(d, s + i * h)
    return tot * h / 3


# ---- (a) NegBin weights sum to 1 at several p (exact rational enclosures) ----
def weight_sum_encl(p_iv, M):
    r = Q(35, 2)
    omp = (rd(1 - p_iv[1]), ru(1 - p_iv[0]))
    pr = a5core.imul(a5core.ipow(p_iv, 17), a5core.isqrt_iv(p_iv))
    w = pr
    s_lo, s_hi = w[0], w[1]
    for j in range(1, M + 1):
        g = (r + j - 1) / j
        w = (rd(w[0] * g * omp[0]), ru(w[1] * g * omp[1]))
        s_lo += w[0]
        s_hi += w[1]
    # tail after M: ratio_j = (r+j)/(j+1)*(1-p) decreasing in j; bound by geometric series
    rho = (r + M) / (M + 1) * omp[1]
    assert rho < 1
    tail = ru(w[1] * rho / (1 - rho))
    return s_lo, s_hi + tail


for name, p_iv, M in [("1/2", (Q(1, 2), Q(1, 2)), 400), ("1/10", (Q(1, 10), Q(1, 10)), 3000),
                      ("p_corner", info["p"], 3000), ("1/37", (Q(1, 37), Q(1, 37)), 12000),
                      ("9/10", (Q(9, 10), Q(9, 10)), 200)]:
    slo, shi = weight_sum_encl(p_iv, M)
    okw = slo <= 1 <= shi and (shi - slo) < Q(1, 10**40)
    V(f"A5-VAL-a-{name}", f"(a) NegBin weights: sum_(j<=M) w_j + geometric tail bound encloses 1, p={name}, M={M}",
      okw, f"enclosure [1 - {float(1-slo):.3e}, 1 + {float(shi-1):.3e}]")

# ---- (b) c = 0 gives exactly 1/2 ----
# direct form Phi = 1/2 - (x/2) sum_j w_j S_{35+j}; at c=0: s=0 -> x = 0 exactly, and each t tail is
# exactly 1/2, independent of the weights (which sum to 1).  Evaluate exactly:
x0 = sqrt_lo(Q(0)), sqrt_hi(Q(0))
# the actual direct-form sum at c=0 (y = 1 exactly, so S_m = sum_{k<m} a_k), weights at p_corner, j <= 400
def _direct_sum_c0(p_iv, M=400):
    r = Q(35, 2)
    omp = (rd(1 - p_iv[1]), ru(1 - p_iv[0]))
    w = a5core.imul(a5core.ipow(p_iv, 17), a5core.isqrt_iv(p_iv))[1]
    a, S = Q(1), Q(0)
    for k in range(35):
        if k:
            a *= Q(2 * k - 1, 2 * k)
        S += a  # y^k = 1
    tot = w * S
    for j in range(1, M + 1):
        w = ru(w * (r + j - 1) / j * omp[1])
        k = 34 + j
        a *= Q(2 * k - 1, 2 * k)
        S += a
        tot += w * S
    return tot
_c0sum = _direct_sum_c0(info["p"])
y0 = 1 / (1 + Q(0) ** 2 * info["rho_lo"][1] / 35)
x0 = sqrt_lo(1 - y0), sqrt_hi(1 - y0)
phi_c0 = Q(1, 2) - x0[1] / 2 * _c0sum
V("A5-VAL-b", "(b) c=0: direct form 1/2 - (x/2)*sum, x = sqrt(1-y) with y = 1 -> x = 0 exactly",
  phi_c0 == Q(1, 2) and x0 == (0, 0), f"y={fr(y0)}, x=[{fr(x0[0])},{fr(x0[1])}], direct sum (j<=400) ~{float(_c0sum):.6g}; Phi(0) = {fr(phi_c0)} exactly")
# also: the c -> 0+ limit of the NONNEGATIVE series form (float), vs 1/2 - Q-dist; small c
def pos_series_float(c, rlo, rhi, tol=1e-15):
    z = c * c * rlo / 35
    y = 1 / (1 + z)
    x = math.sqrt(z / (1 + z))
    p = rlo / rhi
    r = 17.5
    w = p ** r
    F = w
    a = 1.0
    for k in range(1, 36):
        a *= (2 * k - 1) / (2 * k)
    yk = y ** 35
    S = 0.0
    k, j = 35, 0
    while True:
        S += a * yk * F
        k += 1
        a *= (2 * k - 1) / (2 * k)
        yk *= y
        j += 1
        w *= (r + j - 1) / j * (1 - p)
        F += w
        if (a * yk / (1 - y) < tol and 1 - F < 1e-13) or k > 2_000_000:  # float F cannot reach 1-1e-15
            break
    return x * S / 2


def direct_form_float(c, rlo, rhi, J=4000):
    """1/2 - (x/2) sum_j w_j S_{35+j}  (float), independent code path from the positive series."""
    z = c * c * rlo / 35
    y = 1 / (1 + z)
    x = math.sqrt(z / (1 + z))
    p = rlo / rhi
    r = 17.5
    # S_m = sum_{k<m} a_k y^k, incrementally
    a, yk, S = 1.0, 1.0, 0.0
    for k in range(35):
        if k:
            a *= (2 * k - 1) / (2 * k)
        S += a * yk
        yk *= y
    tot = 0.0
    w = p ** r
    wsum = 0.0
    for j in range(J):
        if j:
            w *= (r + j - 1) / j * (1 - p)
            k = 34 + j
            a *= (2 * k - 1) / (2 * k)
            S += a * yk
            yk *= y
        tot += w * S
        wsum += w
    return 0.5 - 0.5 * x * tot, wsum


r1f, r2f = float(ivs["r1"][0]), float(ivs["r2"][0])
dirv, wsum = direct_form_float(float(LC), r2f, r1f)
posv = pos_series_float(float(LC), r2f, r1f)
V("A5-VAL-b2", "(b') direct form 1/2-(x/2)sum_j w_j S_(35+j) vs nonneg. series form, float, at Lc",
  abs(dirv - posv) < 1e-12 and abs(dirv - float(lo)) < 1e-12,
  f"direct={dirv:.15g} (sum w={wsum:.15g}), positive={posv:.15g}, cert_lo={float(lo):.15g}")
smallc = [(0.5, pos_series_float(0.5, r2f, r1f, 1e-13))]
# a c-> 0 trend check against the normal-free limit: Phi(c) -> 1/2
V("A5-VAL-b3", "(b'') nonneg. series at small c=1/2 is below and near 1/2 (limit c->0+ is 1/2)",
  0.35 < smallc[0][1] < 0.5, f"Phi(0.5) = {smallc[0][1]:.12g}")

# ---- (c) t-series vs closed forms at d=2, d=4 (EXACT at rational points), Simpson at d=2,4,70 ----
def rat_sqrt_exact(q):
    """exact sqrt of a rational square (None if not a square) -- via isqrt of num and den"""
    from math import isqrt as _is
    q = Q(q)
    a, b = _is(q.numerator), _is(q.denominator)
    return Q(a, b) if a * a == q.numerator and b * b == q.denominator else None


def t_series_tail_exact(d, s):
    """exact when d+s^2 is a rational square"""
    s = Q(s)
    D = d + s * s
    rt = rat_sqrt_exact(D)
    assert rt is not None and rt * rt == D, "choose s with d+s^2 a rational square"
    y = Q(d) / D
    x = s / rt
    a, S = Q(1), Q(0)
    for k in range(d // 2):
        if k:
            a *= Q(2 * k - 1, 2 * k)
        S += a * y ** k
    return Q(1, 2) - x / 2 * S, rt


# d=2: textbook F(t) = 1/2 + t/(2 sqrt(2+t^2))  -> tail 1/2 - t/(2 sqrt(2+t^2)); at t=1/2: sqrt(9/4)=3/2
d2_pts = [Q(1, 2), Q(7, 4), Q(23, 2), Q(7, 1)]  # 9/4, 81/16, (537/4 no), 51 no; + 1/12? below  # 2+1/4=9/4, 2+49/16=81/16, 2+529/4=537/4? -> filtered below
d2_pts = [t for t in d2_pts if rat_sqrt_exact(2 + t * t) is not None]
ok_c2 = True
det = []
for t in d2_pts:
    ser, rt = t_series_tail_exact(2, t)
    closed = Q(1, 2) - t / (2 * rt)
    ok_c2 &= (ser == closed)
    det.append(f"t={fr(t)}: series={fr(ser)} closed={fr(closed)}")
V("A5-VAL-c-d2", "(c) d=2 series == textbook closed form 1/2 - t/(2 sqrt(2+t^2)), EXACT rationals",
  ok_c2 and len(d2_pts) >= 2, "; ".join(det))
# d=4: textbook (e.g. standard tables) F(t) = 1/2 + (3/8) (t/sqrt(1+t^2/4)) [1 - t^2/(12(1+t^2/4))]
d4_pts = [Q(3, 2), Q(16, 3), Q(5, 6), Q(9, 20)]  # 25/4, (292/9 no), 169/36, 1681/400  # 4+9/4=25/4 ; 4+256/9=292/9? ; 4+25/36=169/36
d4_pts = [t for t in d4_pts if rat_sqrt_exact(4 + t * t) is not None]
ok_c4 = True
det = []
for t in d4_pts:
    ser, rt = t_series_tail_exact(4, t)
    u = 1 + t * t / 4
    su = rt / 2  # sqrt(1+t^2/4) = sqrt(4+t^2)/2
    closed = 1 - (Q(1, 2) + Q(3, 8) * (t / su) * (1 - t * t / (12 * u)))
    ok_c4 &= (ser == closed)
    det.append(f"t={fr(t)}: series={fr(ser)} closed={fr(closed)}")
V("A5-VAL-c-d4", "(c) d=4 series == textbook closed form (3/8 form), EXACT rationals",
  ok_c4 and len(d4_pts) >= 2, "; ".join(det))
# Simpson (float) at d = 2, 4, 70
for dd, svals, tol in [(2, [0.5, 1.7], 1e-7), (4, [1.5, 2.3], 1e-7),
                       (70, [1.0, 1.9, 2.0186, 2.5], 1e-7)]:
    worst = 0.0
    det = []
    for s in svals:
        up = 80.0 if dd == 70 else 1e4
        n = 400000 if dd == 70 else 2000000
        if dd in (2, 4):
            # heavy tails: integrate to U then add the closed tail of the density beyond U analytically is
            # circular; instead use the substitution t = s + u/(1-u) on u in [0,1)
            def g(u, dd=dd, s=s):
                if u >= 1.0:
                    return 0.0
                t = s + u / (1 - u)
                return t_pdf(dd, t) / (1 - u) ** 2
            nn = 200000
            h = 1.0 / nn
            tot = g(0.0) + g(1.0)
            for i in range(1, nn):
                tot += (4 if i % 2 else 2) * g(i * h)
            simp = tot * h / 3
        else:
            simp = simpson_tail(dd, s, up, n)
        ser = t_series_tail(dd, s)
        worst = max(worst, abs(simp - ser))
        det.append(f"s={s}: series={ser:.12g} simpson={simp:.12g}")
    V(f"A5-VAL-c-simpson-d{dd}", f"(c) d={dd} t-series vs float Simpson of the t density, tol {tol:g}",
      worst < tol, f"max|diff|={worst:.3e}; " + "; ".join(det))
# the literal "1/2 + (x/2)sum" of the specification, tested as-is (NOT used): it is the CDF, not the tail
lit = t_tail_literal_plus(2, 0.5)
C("A5-NOTE-spec-sign", "the specification's literal tail formula '1/2 + (x/2)sum' at d=2,t=1/2 (textbook tail = 1/3)",
  f"{lit:.12g} -> it is the CDF F(t)=2/3, not the tail; I derived and use 1/2 - (x/2)sum (tail); "
  f"validation (c) confirms the minus form, the plus form fails (c) outright")

# ---- (d) Monte Carlo, 10^6 draws, fixed seed, stdlib random ----
SEED = 20260930
rng = random.Random(SEED)
nmc = 10 ** 6
cf = float(LC)
hits = 0
for _ in range(nmc):
    X1 = 2.0 * rng.gammavariate(17.5, 1.0)
    X2 = 2.0 * rng.gammavariate(17.5, 1.0)
    Nn = rng.gauss(0.0, 1.0)
    if Nn > cf * math.sqrt((r1f * X1 + r2f * X2) / 35.0):
        hits += 1
phat = hits / nmc
se = math.sqrt(phat * (1 - phat) / nmc)
zmc = (phat - float(lo)) / se
V("A5-VAL-d", f"(d) Monte Carlo n=10^6, seed={SEED}, stdlib random (chi2_35 = 2*Gamma(17.5,1))",
  abs(zmc) < 4, f"hits={hits}, phat={phat:.6f}, SE={se:.3e}, z vs cert={zmc:+.3f} (|z|<4). "
  f"NOTE: SE ({se:.1e}) >> margin ({float(lo-ALPHA):.1e}); MC validates the value, it cannot decide the claim")

# ---- (e) independent cross-checks (NOT decision path) ----
try:
    import mpmath as mp
    mp.mp.dps = 40
    r1m = mp.mpf(277) / 1004 + mp.sqrt(52929) / 1004
    r2m = mp.mpf(277) / 1004 - mp.sqrt(52929) / 1004
    cm = mp.mpf(271731) / 100000
    lam3 = mp.mpf(35) / cm ** 2
    # Phi = 1/2 * P(Q < 0), Q = r1 X1 + r2 X2 - (35/c^2) Z^2 ; Imhof: P(Q>0) = 1/2 + (1/pi) int sin(th)/(u rho) du
    def integrand(u):
        th = (35 * mp.atan(r1m * u) + 35 * mp.atan(r2m * u) - mp.atan(lam3 * u)) / 2
        rh = (1 + (r1m * u) ** 2) ** (mp.mpf(35) / 4) * (1 + (r2m * u) ** 2) ** (mp.mpf(35) / 4) \
            * (1 + (lam3 * u) ** 2) ** (mp.mpf(1) / 4)
        return mp.sin(th) / (u * rh)
    I = mp.quad(integrand, [0, 0.5, 1, 2, 4, 8, 16, 32, 64, mp.inf])
    PQpos = mp.mpf(1) / 2 + I / mp.pi
    phi_imhof = (1 - PQpos) / 2
    lo_m = mp.mpf(lo.numerator) / lo.denominator
    hi_m = mp.mpf(hi.numerator) / hi.denominator
    okI = lo_m - mp.mpf(10) ** -25 <= phi_imhof <= hi_m + mp.mpf(10) ** -25
    V("A5-XCHK-imhof", "mpmath 1.3.0 (40 dps) Imhof/Gil-Pelaez integral for Phi = P(Q<0)/2, Q = r1X1+r2X2-(35/Lc^2)Z^2 "
      "(independent route; checks the certificate value)", bool(okI),
      f"Phi_imhof = {mp.nstr(phi_imhof, 25)} ; cert=[{float(lo):.17g},{float(hi):.17g}] ; "
      f"Phi_imhof - 1/40 = {mp.nstr(phi_imhof - mp.mpf(1)/40, 10)}")
except Exception as e:  # pragma: no cover
    V("A5-XCHK-imhof", "mpmath Imhof cross-check", False, f"error {e!r}")

try:
    from scipy import integrate, stats
    # 2-D float quadrature: Phi = E_{X1,X2}[ Qnorm(c sqrt((r1 X1 + r2 X2)/35)) ]  (another independent route)
    f35 = stats.chi2(35).pdf
    def inner(v):
        return integrate.quad(lambda u: f35(u) * stats.norm.sf(cf * math.sqrt((r1f * u + r2f * v) / 35.0)),
                              0, 200, epsabs=1e-14, epsrel=1e-12, limit=400)[0]
    phi2d = integrate.quad(lambda v: f35(v) * inner(v), 0, 200, epsabs=1e-14, epsrel=1e-12, limit=400)[0]
    V("A5-XCHK-scipy2d", "scipy 1.17.1 2-D quadrature of E[Qnorm(c sqrt((r1X1+r2X2)/35))] (float; independent route)",
      abs(phi2d - float(lo)) < 1e-9, f"Phi_2d = {phi2d:.15g} ; |diff vs cert| = {abs(phi2d-float(lo)):.2e}")
except Exception as e:  # pragma: no cover
    V("A5-XCHK-scipy2d", "scipy 2-D cross-check", False, f"error {e!r}")

# ======================= 4. mutation tests (checker must FAIL) =======================
def mut(cid, label, alpha=ALPHA, c=LC, corner=CORNER, expect_fail=True):
    okm, lom, him, infom = check(alpha, c, *corner)
    behaved = (not okm) if expect_fail else okm
    C(cid, label, f"checker says {'PASS' if okm else 'FAIL'} (expected {'FAIL' if expect_fail else 'PASS'}) "
      f"-> {'OK' if behaved else 'MUTATION NOT DETECTED'} ; bounds [{float(lom):.12g}, {float(him):.12g}]")
    return behaved


mres = []
mres.append(mut("A5-MUT-1", "alpha = 1/4", alpha=Q(1, 4)))
mres.append(mut("A5-MUT-2", "b0 sign flipped, corner (1/25, +3/2, 3/2) (flips sign of b0*b1; changes aF to 26/25)",
                corner=(Q(1, 25), Q(3, 2), Q(3, 2))))
mres.append(mut("A5-MUT-3", "b1 sign flipped, corner (1/25, -3/2, -3/2) (flips sign of b0*b1)",
                corner=(Q(1, 25), Q(-3, 2), Q(-3, 2))))
mres.append(mut("A5-MUT-4", "alpha = 1/40 + 2e-6 (just above the true value)", alpha=ALPHA + Q(2, 10**6)))
mres.append(mut("A5-MUT-5", "alpha = certified lower bound itself (strict inequality must fail)", alpha=lo))
mres.append(mut("A5-MUT-6", "Lc = 2.7175 (+1.9e-4)", c=Q(27175, 10000)))
mres.append(mut("A5-MUT-7", "lambda = 1/20 instead of 1/25", corner=(Q(1, 20), Q(-3, 2), Q(3, 2))))
mres.append(mut("A5-MUT-ctrl", "CONTROL (should pass): alpha = 1/40 - 1e-6", alpha=ALPHA - Q(1, 10**6),
                expect_fail=False))
# inert mutation: sign of g01 alone (g01 enters only via g01^2)
exm = dict(ex)
exm["g01"] = -ex["g01"]
exm["detG"] = exm["g00"] * exm["g11"] - exm["g01"] ** 2
C("A5-MUT-inert", "flip sign of g01 ONLY (b0*b1 -> -b0*b1 inside g01)",
  f"detG unchanged ({exm['detG'] == ex['detG']}) -> r1,r2 and Phi unchanged; this mutation is INERT by algebra "
  f"(g01 enters squared), so the checker correctly does not flag it. The detectable b0b1-sign mutations are MUT-2/MUT-3.")

# rounding-direction mutation at coarse precision: outward must enclose truth, inward may not
save_rd, save_ru, saveP = a5core.rd, a5core.ru, a5core.P_BITS
try:
    truth = float(phi_imhof)
except Exception:
    truth = float(lo)
a5core.P_BITS = 64
lo_o, hi_o, _ = phi_bounds(LC, ratio_intervals(ex)["r1"], ratio_intervals(ex)["r2"], tol_exp=12, K_max=40000)
a5core.rd, a5core.ru = save_ru, save_rd  # deliberately round INWARD
try:
    lo_i, hi_i, _ = phi_bounds(LC, ratio_intervals(ex)["r1"], ratio_intervals(ex)["r2"], tol_exp=12, K_max=40000)
except Exception as e:
    lo_i, hi_i = None, None
a5core.rd, a5core.ru, a5core.P_BITS = save_rd, save_ru, saveP
enc_o = float(lo_o) <= truth <= float(hi_o)
enc_i = lo_i is not None and float(lo_i) <= truth <= float(hi_i)
C("A5-MUT-round", "rounding-direction mutation at P=64 bits, tol 1e-12, K_max 40000 (outward vs deliberately inward)",
  f"outward [{float(lo_o):.16g},{float(hi_o):.16g}] (truth ~{truth:.16g}) encloses truth: {enc_o}; "
  f"inward [{float(lo_i) if lo_i is not None else 'err'},{float(hi_i) if hi_i is not None else 'err'}] "
  f"encloses truth: {enc_i} -> {'OK: inward rounding is detected as non-enclosing' if (enc_o and not enc_i) else 'inward rounding NOT distinguished at this precision'}")
C("A5-MUT-summary", "all mutations behaved as expected", str(all(mres)))

# ======================= 5. witness record (the one allowed project file) =======================
with open(WITNESS, "r", encoding="utf-8") as fh:
    wl = [ln.rstrip("\n") for ln in fh]
hitl = [(i + 1, ln) for i, ln in enumerate(wl)
        if ln.startswith("[1/25,") and "|[-3/2," in ln and ", 3/2]#" in ln]
C("A5-WIT-lines", "witness file line count; lines whose sub-box contains the corner (1/25,-3/2,3/2)",
  f"{len(wl)} lines; {len(hitl)} match: " + " | ".join(f"line {i}: {ln}" for i, ln in hitl))
for i, ln in hitl:
    parts = ln.split("#")
    wa, wb = Q(parts[1]), Q(parts[2])
    C("A5-WIT-cmp", f"witness line {i} numbers vs my corner r1, r2 (informational only; format undocumented)",
      f"{float(wa):.17g} vs r1~{float(ivs['r1'][0]):.17g} (w <= r1_lo: {wa <= ivs['r1'][0]}); "
      f"{float(wb):.17g} vs r2~{float(ivs['r2'][0]):.17g} (w <= r2_lo: {wb <= ivs['r2'][0]})")

C("A5-NOTE-harness-fixes", "harness bugs fixed after run 1 (the Phi formula, series, bounds logic and decision rule are unchanged; a5core.py changed ONLY in item (v), how the precision parameter is resolved)",
  "run 1 had: (i) float-only helper loop could not terminate (1-F<1e-15 unreachable in doubles) -> threshold 1e-13 + cap; "
  "(ii) (c)-d4 FAILED only its '>=2 test points' clause because my square detector (sqrt_lo==sqrt_hi) recognises only "
  "dyadic squares and dropped t=5/6; the one point evaluated (t=3/2) matched exactly 13/125; detector replaced by exact "
  "isqrt(num),isqrt(den); (iii) mpmath cross-check raised TypeError(mpf from Fraction) -> explicit num/den conversion; "
  "(iv) rounding mutation re-parameterised and (v) a5core rd/ru/sqrt_lo/sqrt_hi had default arg P=P_BITS bound at definition time, so run-1's P=160 re-run and P=24 rounding mutation silently ran at 320 bits; now P resolves at call time (default 320 unchanged, main certificate unaffected); precision re-run now P=200, rounding mutation now P=64")

# ======================= 6. provenance + verdict =======================
C("A5-FILES-OPENED", "files opened",
  "deposits/f8c_cover_witness-R1.txt (read-only) ; this folder's own files only. "
  "No other file of the project was opened, searched or imported.")
C("A5-LIBS", "libraries", "decision path: Python 3.12 stdlib (fractions, math.isqrt). Cross-checks only: "
  "mpmath 1.3.0 (Imhof integral, checks the Phi value), scipy 1.17.1 (2-D quadrature, checks the Phi value). "
  "numpy not used directly.")
C("A5-VALID-ALL", "all validations passed", str(all_valid))
verdict = "CONFIRMED" if (ok and all_valid and all(mres)) else "REFUTED"
claims.append(f"CLAIM A5-VERDICT: {verdict}")
with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(claims) + "\n")
lines_out = [c for c in claims]
print("verdict line: CLAIM A5-VERDICT:", verdict)
