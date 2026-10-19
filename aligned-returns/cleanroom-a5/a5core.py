"""A5 clean-room certificate core.  Standard library only.  No floats on the decision path.

Route: negative-binomial mixture of t tails, rewritten as a series of NONNEGATIVE terms.

  rho_lo*X1 + rho_hi*X2  ~  rho_lo * chi2_{70+2J},  J ~ NegBin(r=35/2, p=rho_lo/rho_hi),
  w_j = C(r+j-1, j) p^r (1-p)^j.
  Given J=j (d = 70+2j):  P(N > c sqrt(rho_lo chi2_d / 35)) = P(t_d > s_j),  s_j^2 = c^2 rho_lo d / 35.
  Even-d t tail:  P(t_d > s) = 1/2 - (x/2) * sum_{k<d/2} a_k y^k,
      a_0 = 1, a_k = a_{k-1}(2k-1)/(2k),  y = d/(d+s^2) = 1/(1+c^2 rho_lo/35),  x = sqrt(1-y)  (>= 0).
  (NB: this is 1/2 MINUS; "1/2 + ..." is the CDF, not the tail.)
  Since sum_{k>=0} a_k y^k = (1-y)^(-1/2) = 1/x for 0<=y<1, the tail equals
      (x/2) * sum_{k >= d/2} a_k y^k,
  so, exchanging the (nonnegative) double sum,
      Phi = (x/2) * sum_{k>=35} a_k y^k F(k-35),   F(m) = P(J <= m) = sum_{j<=m} w_j.
  Every term is >= 0, so any partial sum built from LOWER bounds of x, y, a_k, w_j is a rigorous
  lower bound.  Upper bound: upper bounds of everything, F <= min(1, sum w_hi), plus the tail
      sum_{k>K} a_k y^k F <= a_{K+1} y^{K+1} / (1-y)   (a_k is nonincreasing).
"""
from fractions import Fraction as Q
from math import isqrt

P_BITS = 320  # working precision of outward rounding (denominator 2^P_BITS)


# ---------- outward rounding ----------
def rd(q, P=None):
    """largest m/2^P <= q"""
    P = P_BITS if P is None else P  # resolved at CALL time
    q = Q(q)
    return Q((q.numerator << P) // q.denominator, 1 << P)


def ru(q, P=None):
    """smallest m/2^P >= q"""
    P = P_BITS if P is None else P  # resolved at CALL time
    q = Q(q)
    return Q(-((-(q.numerator << P)) // q.denominator), 1 << P)


def sqrt_lo(q, P=None):
    """rational <= sqrt(q), q >= 0"""
    P = P_BITS if P is None else P  # resolved at CALL time
    q = Q(q)
    assert q >= 0
    n = (q.numerator << (2 * P)) // q.denominator  # floor(q*4^P)
    return Q(isqrt(n), 1 << P)  # isqrt(floor(q 4^P)) <= sqrt(q 4^P)


def sqrt_hi(q, P=None):
    """rational >= sqrt(q), q >= 0"""
    P = P_BITS if P is None else P  # resolved at CALL time
    q = Q(q)
    assert q >= 0
    n = -((-(q.numerator << (2 * P))) // q.denominator)  # ceil(q*4^P)
    u = isqrt(n)
    if u * u < n:
        u += 1
    return Q(u, 1 << P)  # u >= sqrt(ceil(q 4^P)) >= sqrt(q 4^P)


# ---------- intervals of POSITIVE rationals (lo, hi) ----------
def imul(a, b):
    assert a[0] >= 0 and b[0] >= 0
    return (rd(a[0] * b[0]), ru(a[1] * b[1]))


def idiv(a, b):
    assert a[0] >= 0 and b[0] > 0
    return (rd(a[0] / b[1]), ru(a[1] / b[0]))


def isqrt_iv(a):
    return (sqrt_lo(a[0]), sqrt_hi(a[1]))


def ipow(a, n):
    r = (Q(1), Q(1))
    for _ in range(n):
        r = imul(r, a)
    return r


# ---------- exact corner quantities ----------
def corner_exact(lam, b0, b1):
    lam, b0, b1 = Q(lam), Q(b0), Q(b1)
    aF = (b1 - b0) ** 2 + 1 + lam
    g00 = b0 ** 2 + 1
    g11 = b1 ** 2 + lam
    g01 = b0 * b1
    trG = g00 + g11
    detG = g00 * g11 - g01 ** 2
    discG = trG ** 2 - 4 * detG
    return dict(aF=aF, g00=g00, g11=g11, g01=g01, trG=trG, detG=detG, discG=discG)


def ratio_intervals(ex):
    """intervals for mu1, mu2, r1, r2 (outward).  sqrt(discG) is the only irrational."""
    sd = (sqrt_lo(ex['discG']), sqrt_hi(ex['discG']))
    tr, aF = ex['trG'], ex['aF']
    mu1 = (rd((tr + sd[0]) / 2), ru((tr + sd[1]) / 2))
    mu2 = (rd((tr - sd[1]) / 2), ru((tr - sd[0]) / 2))  # outward: subtract the UPPER sqrt for lo
    assert mu2[0] > 0, "mu2 lower bound not positive"
    r1 = (rd(mu1[0] / aF), ru(mu1[1] / aF))
    r2 = (rd(mu2[0] / aF), ru(mu2[1] / aF))
    return dict(sqrt_discG=sd, mu1=mu1, mu2=mu2, r1=r1, r2=r2)


# ---------- the certificate ----------
def phi_bounds(c, rho1, rho2, nu=35, tol_exp=60, K_max=200000):
    """Certified [lo, hi] for Phi(c; rho1, rho2) with rho* given as positive intervals.
    nu must be such that 2*nu is the total d.o.f. with each X ~ chi2_nu; here nu=35 (r = 35/2)."""
    c = Q(c)
    assert c > 0
    # order the weights rigorously
    if rho2[1] < rho1[0]:
        rlo, rhi = rho2, rho1
    elif rho1[1] < rho2[0]:
        rlo, rhi = rho1, rho2
    else:
        raise ValueError("cannot order rho1, rho2 rigorously")
    r = Q(nu, 2)
    half = nu  # d/2 at j=0 is (2*nu)/2 = nu = 35
    c2 = c * c
    # z = c^2 rho_lo / nu ; y = 1/(1+z) ; x = sqrt(1-y) = sqrt(z/(1+z))
    z = (rd(c2 * rlo[0] / nu), ru(c2 * rlo[1] / nu))
    y = (rd(1 / (1 + z[1])), ru(1 / (1 + z[0])))
    xx = (rd(z[0] / (1 + z[0])), ru(z[1] / (1 + z[1])))  # z/(1+z) increasing in z
    x = isqrt_iv(xx)
    one_minus_y_lo = rd(1 - y[1])
    assert one_minus_y_lo > 0
    # p = rho_lo/rho_hi ; 1-p
    p = idiv(rlo, rhi)
    assert p[1] < 1
    omp = (rd(1 - p[1]), ru(1 - p[0]))
    # p^r = p^(nu//2) * sqrt(p)   (nu odd -> r half-integer)
    assert nu % 2 == 1
    pr = imul(ipow(p, nu // 2), isqrt_iv(p))
    # a_k and y^k for k < 35 : advance to k = 35
    a = (Q(1), Q(1))
    yk = (Q(1), Q(1))
    for k in range(1, half + 1):
        f = Q(2 * k - 1, 2 * k)
        a = (rd(a[0] * f), ru(a[1] * f))
        yk = imul(yk, y)
    # now a = a_35, yk = y^35
    w = pr  # w_0
    F_lo, F_hi = w[0], w[1]
    S_lo, S_hi = Q(0), Q(0)
    tol = Q(1, 10 ** tol_exp)
    k = half
    j = 0  # m = k - 35
    while True:
        # term k uses F(k-35) = sum_{j<=k-35} w_j
        S_lo += rd(a[0] * yk[0] * F_lo)
        S_hi += ru(a[1] * yk[1] * min(Q(1), F_hi))
        # next k
        k += 1
        f = Q(2 * k - 1, 2 * k)
        a = (rd(a[0] * f), ru(a[1] * f))
        yk = imul(yk, y)
        # next weight w_{j+1} = w_j * (r+j)/(j+1) * (1-p)
        j += 1
        g = (r + j - 1) / j
        w = (rd(w[0] * g * omp[0]), ru(w[1] * g * omp[1]))
        F_lo += w[0]
        F_hi += w[1]
        tail_hi = ru(a[1] * yk[1] / one_minus_y_lo)  # bound on sum_{k' >= k} a y^k'
        if (tail_hi < tol and 1 - F_lo < tol) or k > K_max:
            break
    lo = rd(x[0] * S_lo / 2)
    hi = ru(x[1] * (S_hi + tail_hi) / 2)
    info = dict(K_last=k - 1, terms=k - half, y=y, x=x, z=z, p=p, pr=pr, F_lo_final=F_lo,
                tail_hi=tail_hi, rho_lo=rlo, rho_hi=rhi)
    return lo, hi, info


def check(alpha, c, lam, b0, b1):
    """The decision: returns (passed, lo, hi, info). passed iff certified lo > alpha."""
    ex = corner_exact(lam, b0, b1)
    iv = ratio_intervals(ex)
    lo, hi, info = phi_bounds(c, iv['r1'], iv['r2'])
    return (lo > Q(alpha)), lo, hi, dict(ex=ex, iv=iv, **info)
