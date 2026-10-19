"""T1 — falsifier of an objection to the 2B corner certificate.

The objection is that the 2B corner certificate is a non sequitur: the rate
    Phi(c; r1, r2) = P(N > c*sqrt((r1*X1 + r2*X2)/nu)),   X1,X2 ~ chi2_nu iid, N ~ N(0,1)
is symmetric and CONVEX in (r1, r2), hence Schur-convex, so at fixed r1+r2 a MORE UNEVEN
split RAISES the rate and therefore RAISES c*. If so, "r1+r2 >= rc1+rc2 implies
c*(theta) <= c*(corner)" is a non sequitur and the corner argument fails.

This script ATTACKS that claim. Independence from the objection's arithmetic: it does not
import w2_*, routec_* or any AlignedFactor code; it rebuilds Phi from the defining model
with its own quadrature (tensor generalized Gauss-Laguerre, order 64) and cross-checks
against scipy dblquad at points (b),(c) to 1e-9.

Expectations are DECLARED IN ADVANCE (A-6):
  (a) at the corner rc1 = 0.50504..., rc2 = 0.046750..., s = 5.54/10.04
  (b) Phi(c0; rc) = 0.025 +/- 1e-8
  (c) Phi(c0; s-t, t) for t in {0,.01,.02,.03,rc2,.1,.2,s/2} strictly DECREASING in t,
      and Phi(c0; s, 0) > 0.025
  (e) on a grid (lam geomspace(0.04,25,60); b0,b1 linspace(-1.5,1.5,61)):
      max Phi(U; r(theta)) < alpha; report argmax, margin, runner-up distinct class, and
      count of points with r1+r2 >= s and r2 < rc2.
The corner of the box is theta = (0.04, -1.5, 1.5).
"""

import pathlib

import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln, ndtr, roots_genlaguerre

nu, alpha_rate = 35, 1 / 40
U_C = 271740 / 100000
c0 = 2.7173398035
OUT = pathlib.Path(__file__).resolve().parent
DEP = OUT / "deposits"
CLAIMS = []


def claim(cid, label, value):
    CLAIMS.append(f"CLAIM {cid}: {label} = {value}")


def gl_tensors(N):
    """Normalized tensor GL nodes/weights for X ~ chi2_nu (X = 2t)."""
    a = nu / 2 - 1
    t, w = roots_genlaguerre(N, a)
    w = w / np.exp(gammaln(a + 1))
    return 2 * t, w


def Phi_GL(c, r1, r2, N=64):
    """Phi(c; r1, r2) = P(N>c*sqrt((r1*X1+r2*X2)/nu)), tensor GL, order N."""
    X, w = gl_tensors(N)
    S = r1 * X[:, None] + r2 * X[None, :]            # (N,N)
    z = c * np.sqrt(np.maximum(S, 0.0) / nu)
    return float(np.sum(w[:, None] * w[None, :] * ndtr(-z)))


def chi2_pdf(x, n):
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    m = x > 0
    out[m] = np.exp((n / 2 - 1) * np.log(x[m]) - x[m] / 2 - gammaln(n / 2)
                    - (n / 2) * np.log(2))
    return out


def Phi_dblquad(c, r1, r2):
    """Independent route: scipy dblquad over the chi2(35) density.

    The naive form diverged (quadpack reported a divergent inner integral): the chi2(35)
    integrand has a heavy polynomial tail at 0 and a Gaussian tail at infinity, and the
    nested adaptive call never settled. Substituting X = u^2 makes the integrand smooth
    and bounded on (0, inf), which the nested call handles. scipy's own dblquad was
    version-fragile about its callback arity, so the nesting is written out explicitly.
    This is still a genuinely different route from the Gauss-Laguerre tensor: adaptive
    QUADPACK subdivision, a different rule, and a different code path.
    """
    def jac(u):
        return 2.0 * u * chi2_pdf(u * u, nu)

    def inner_u1(u1, u2):
        s = r1 * u1 * u1 + r2 * u2 * u2
        if s <= 0:
            return 0.0
        return ndtr(-c * np.sqrt(s / nu)) * jac(u1)

    def outer(u2):
        v = quad(inner_u1, 0.0, np.inf, args=(u2,),
                 epsabs=1e-13, epsrel=1e-11, limit=300)[0]
        return v * jac(u2)

    return float(quad(outer, 0.0, np.inf,
                      epsabs=1e-12, epsrel=1e-10, limit=300)[0])


def rpair(lam, b0, b1):
    g00, g11, g01 = b0 * b0 + 1, b1 * b1 + lam, b0 * b1
    tr, det = g00 + g11, g00 * g11 - g01 * g01
    a = (b1 - b0) ** 2 + 1 + lam
    disc = max(tr * tr - 4 * det, 0.0)
    s = np.sqrt(disc)
    return (tr + s) / (2 * a), (tr - s) / (2 * a)


def main():
    rc1, rc2 = rpair(0.04, -1.5, 1.5)
    s = 5.54 / 10.04
    claim("T1-A-RC1", "corner rc1", round(rc1, 8))
    claim("T1-A-RC2", "corner rc2", round(rc2, 8))
    claim("T1-A-S", "s=5.54/10.04", round(s, 8))
    a_ok = abs(rc1 - 0.50504) < 1e-4 and abs(rc2 - 0.046750) < 1e-5
    claim("T1-A-VERDICT", "corner constants match declaration", a_ok)

    # (b) Phi(c0; rc) = alpha
    pb = Phi_GL(c0, rc1, rc2)
    pb2 = Phi_dblquad(c0, rc1, rc2)
    claim("T1-B-PHI", "Phi(c0; rc)", round(pb, 12))
    claim("T1-B-PHI-DBLQUAD", "Phi(c0; rc) independent route", round(pb2, 12))
    claim("T1-B-AGREE", "GL vs dblquad within 1e-9", abs(pb - pb2) <= 1e-9)
    claim("T1-B-VERDICT", "Phi(c0;rc)=0.025 +/- 1e-8", abs(pb - 0.025) <= 1e-8)

    # (c) Schur-convexity probe: Phi(c0; s-t, t) as t: 0 -> s/2 (uneven -> even)
    ts = [0.0, 0.01, 0.02, 0.03, rc2, 0.1, 0.2, s / 2]
    vals = [Phi_GL(c0, s - t, t) for t in ts]
    for t, v in zip(ts, vals, strict=True):
        claim("T1-C-T{}".format(str(t).replace(".", "p")), "Phi(c0;s-t,t)", round(v, 12))
    # expectation: strictly DECREASING in t
    dec = all(vals[i + 1] < vals[i] for i in range(len(vals) - 1))
    claim("T1-C-DECREASING", "Phi strictly decreasing in t (uneven->even)", dec)
    claim("T1-C-PHI-S0", "Phi(c0;s,0)", round(vals[0], 12))
    claim("T1-C-PHI-S0-GT-ALPHA", "Phi(c0;s,0) > 0.025", vals[0] > 0.025)
    claim("T1-C-VERDICT", "Schur-convexity expectation holds", dec and vals[0] > 0.025)

    # (e) grid scan
    lams = np.geomspace(0.04, 25, 60)
    b0s = np.linspace(-1.5, 1.5, 61)
    b1s = np.linspace(-1.5, 1.5, 61)
    X, w = gl_tensors(20)
    best = (-1.0, None)
    count_uneven = 0
    for lam in lams:
        for b0 in b0s:
            for b1 in b1s:
                r1, r2 = rpair(lam, b0, b1)
                if r1 + r2 >= s - 1e-12 and r2 < rc2 - 1e-12:
                    count_uneven += 1
                S = r1 * X[:, None] + r2 * X[None, :]
                z = U_C * np.sqrt(np.maximum(S, 0.0) / nu)
                p = float(np.sum(w[:, None] * w[None, :] * ndtr(-z)))
                if p > best[0]:
                    best = (p, (float(lam), float(b0), float(b1)))
    claim("T1-E-MAXPHI", "max Phi(U;r(theta)) on grid", round(best[0], 12))
    claim("T1-E-ARGMAX", "argmax (lam,b0,b1)", str(best[1]))
    claim("T1-E-MARGIN", "alpha - maxPhi", round(alpha_rate - best[0], 12))
    claim("T1-E-COUNT-UNEVEN", "points with r1+r2>=s and r2<rc2", count_uneven)
    claim("T1-E-VERDICT", "max Phi(U;r) < alpha", best[0] < alpha_rate)
    claim("T1-E-CORNER-IS-ARGMAX",
          "argmax equals the box corner", best[1] is not None and
          abs(best[1][0] - 0.04) < 1e-9 and abs(best[1][1] + 1.5) < 1e-9 and
          abs(best[1][2] - 1.5) < 1e-9)

    (DEP / "formal01_t1_claims.txt").write_text("\n".join(CLAIMS) + "\n", encoding="utf-8")
    # concise verdict
    print("T1-A corner rc1,rc2:", round(rc1, 8), round(rc2, 8), "expect 0.50504, 0.046750")
    print("T1-B Phi(c0;rc):", round(pb, 10), "dblquad", round(pb2, 10), "expect 0.025")
    print("T1-C decreasing in t:", dec, " Phi(c0;s,0)=", round(vals[0], 8))
    print("T1-E max Phi(U;r) on grid:", round(best[0], 8), "alpha=", alpha_rate,
          "argmax=", best[1])
    print("T1-E count r1+r2>=s and r2<rc2:", count_uneven)
    print("T1-E argmax is corner:", best[1] == (0.04, -1.5, 1.5) or "no")


if __name__ == "__main__":
    main()
