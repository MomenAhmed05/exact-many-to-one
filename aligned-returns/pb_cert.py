# -*- coding: utf-8 -*-
"""pb_cert.py — option B: certified bounds for the calibrated table rule at K = 2.

Every transcendental value comes from an Arb ball (python-flint) and is stored as a float rounded
outward. Floating point is used only for +, -, *, / and sqrt, which IEEE 754 rounds correctly; every
float argument passed to a table is padded outward by PAD * (1 + |y|), far above the accumulated
rounding of the few operations that produce it, and every float sum of positive terms is inflated
by SUMPAD.

Coordinates (Bartlett): x = sqrt(u), u ~ chi2_nu; e ~ N(0,1); q' = sqrt(w), w ~ chi2_{nu-1}; all
independent and free of rho. m = rho x + s e, q = s q', s = sqrt(1 - rho^2). Then a11 = x^2,
a22 = m^2 + q^2, r_hat = m / sqrt(m^2 + q^2), and t_hat = m / q = r_hat / sqrt(1 - r_hat^2) is
increasing in r_hat, so the table's bins are t_hat-intervals with edges tau_e.

One candidate's rejection rate. With g_j(x) = P(Z > c_j x / sqrt(nu)) nondecreasing in j (the
deployed constants are nonincreasing) and T_e(y) = P(t_hat > tau_e | x) = E_{q'} Phi(y - tau_e q'),
y = rho x / s:

    P1(rho) = E_x[ g_0(x) + sum_{j>=1} (g_j(x) - g_{j-1}(x)) T_{j-1}(rho x / s) ].

T_e is increasing in y and rho / s is increasing in rho, so P1 is nondecreasing in rho: on a
rho-interval its supremum is its value at the upper end.
"""
import math

import numpy as np
from flint import arb, ctx, fmpq

ctx.prec = 96
PAD = 1e-10
SUMPAD = 1e-9


def flo(a):
    """Float lower bound of an arb ball."""
    return math.nextafter(float(a.lower()), -math.inf)


def fhi(a):
    return math.nextafter(float(a.upper()), math.inf)


def q(x):
    """Exact rational of a float."""
    n, d = float(x).as_integer_ratio()
    return arb(fmpq(n, d))


# ---------------------------------------------------------------- Phi table --------------
class PhiTable:
    """Phi(y) = P(N(0,1) <= y) on the dyadic grid y_i = Y0 + i D, outward-rounded."""

    def __init__(self, y0=-16.0, y1=16.0, logD=-15):
        self.Y0, self.D = y0, 2.0 ** logD
        n = int(round((y1 - y0) / self.D)) + 1
        lo, hi = np.empty(n), np.empty(n)
        r2 = arb(2).sqrt()
        for i in range(n):
            y = q(y0 + i * self.D)
            v = (-y / r2).erfc() / 2
            lo[i], hi[i] = max(flo(v), 0.0), min(fhi(v), 1.0)
        self.lo, self.hi, self.N = lo, hi, n

    def up(self, y):
        """Upper bound of Phi(y) for float y (array)."""
        y = np.asarray(y, float)
        yy = y + PAD * (1 + np.abs(y))
        i = np.ceil((yy - self.Y0) / self.D + 1e-9)
        out = np.ones_like(y)
        ok = i <= self.N - 1
        i = np.clip(i, 0, self.N - 1).astype(np.int64)
        out[ok] = self.hi[i[ok]]
        return out

    def lo_(self, y):
        """Lower bound of Phi(y)."""
        y = np.asarray(y, float)
        yy = y - PAD * (1 + np.abs(y))
        i = np.floor((yy - self.Y0) / self.D - 1e-9)
        out = np.zeros_like(y)
        ok = i >= 0
        i = np.clip(i, 0, self.N - 1).astype(np.int64)
        out[ok] = self.lo[i[ok]]
        return out


# ---------------------------------------------------------------- chi cells --------------
def chi_cells(df, ncell, tail=1e-16):
    """Uniform cells in v = sqrt(chi2_df) on [0, vmax]; certified masses (upper and lower).

    Returns edges (floats, dyadic multiples), mass_hi, mass_lo, tail_hi (P(v > vmax), upper)."""
    from scipy.stats import chi2
    vmax = math.sqrt(chi2.isf(tail, df))
    step = 2.0 ** math.floor(math.log2(vmax / ncell))
    k = int(math.ceil(vmax / step))
    edges = np.arange(k + 1) * step
    s = arb(df) / 2
    cdf = [ (q(v) * q(v) / 2).gamma_lower(s, regularized=1) if v > 0 else arb(0) for v in edges]
    mhi = np.array([fhi(cdf[i + 1] - cdf[i]) for i in range(k)])
    mlo = np.array([max(flo(cdf[i + 1] - cdf[i]), 0.0) for i in range(k)])
    tail_hi = fhi(1 - cdf[-1])
    return edges, mhi, mlo, tail_hi


# ---------------------------------------------------------------- the deployed table ----
class Table:
    """Deployed rule: rational r-edges and rational constants, constants nonincreasing."""

    def __init__(self, r_edges, consts):
        self.r = [fmpq(*float(x).as_integer_ratio()) if not isinstance(x, fmpq) else x for x in r_edges]
        self.c = [fmpq(*float(x).as_integer_ratio()) if not isinstance(x, fmpq) else x for x in consts]
        assert len(self.c) == len(self.r) + 1
        assert all(self.c[j] <= self.c[j - 1] for j in range(1, len(self.c))), "constants must be nonincreasing"
        tau = [arb(r) / (1 - arb(r) ** 2).sqrt() for r in self.r]
        self.tau_lo = np.array([flo(t) for t in tau])
        self.tau_hi = np.array([fhi(t) for t in tau])
        self.c_lo = np.array([flo(arb(c)) for c in self.c])
        self.c_hi = np.array([fhi(arb(c)) for c in self.c])


# ---------------------------------------------------------------- P1 ---------------------
class P1Bound:
    """Certified upper bound of one candidate's rejection rate, for one n and one table."""

    def __init__(self, n, table, phi, nx=4000, nq=2000, y_max=60.0, dy=0.005):
        self.n, self.nu, self.T, self.phi = n, n - 1, table, phi
        nu = self.nu
        self.xe, self.xm_hi, _, self.xtail = chi_cells(nu, nx)
        self.qe, qm_hi, _, qtail = chi_cells(nu - 1, nq)
        # T_e(y) upper bounds on the grid y_k = -y_max + k dy
        self.ys = -y_max + dy * np.arange(int(round(2 * y_max / dy)) + 1)
        self.dy, self.ymin = dy, -y_max
        ne = len(table.tau_lo)
        Tup = np.empty((ne, len(self.ys)))
        qa, qb = self.qe[:-1], self.qe[1:]
        for e in range(ne):
            tl = table.tau_lo[e]
            qw = qa if tl >= 0 else qb                      # q' that maximises y - tau q'
            for k0 in range(0, len(self.ys), 2000):
                y = self.ys[k0:k0 + 2000]
                arg = y[:, None] - tl * qw[None, :]
                Tup[e, k0:k0 + 2000] = (phi.up(arg) @ qm_hi) * (1 + SUMPAD) + qtail
        self.Tup = np.minimum(Tup, 1.0)

    def T_at(self, y):
        """Upper bound of T_e(y) for every edge e, at float y (array) -> (ne, len(y))."""
        y = np.asarray(y, float)
        yy = y + PAD * (1 + np.abs(y))
        k = np.ceil((yy - self.ymin) / self.dy + 1e-9)
        out = np.ones((self.Tup.shape[0], len(y)))
        ok = k <= len(self.ys) - 1
        k = np.clip(k, 0, len(self.ys) - 1).astype(np.int64)
        out[:, ok] = self.Tup[:, k[ok]]
        return out

    def upper(self, rho):
        """Certified upper bound of P1(rho), rho a float in (-1, 1)."""
        nu = self.nu
        xa, xb = self.xe[:-1], self.xe[1:]
        s = math.sqrt(1.0 - rho * rho)
        if rho > 0:                                          # y = rho x / s, maximised over the cell
            y = rho * xb / (s * (1 - 1e-15)) * (1 + 1e-14)
        elif rho < 0:
            y = rho * xa / (s * (1 + 1e-15)) * (1 - 1e-14)
        else:
            y = np.zeros_like(xa)
        Tm = self.T_at(y)                                    # (ne, ncell)
        rnu = math.sqrt(nu) * (1 + 1e-15)
        # g_j(xa) upper bounds: Phi(-c_j xa / sqrt(nu)), c at its lower end
        g = self.phi.up(-(self.T.c_lo[:, None] * xa[None, :]) / rnu)   # (nbin, ncell)
        g = np.maximum.accumulate(g, axis=0)                 # nondecreasing in j, still upper bounds
        dg = np.diff(g, axis=0)                              # >= 0
        inner = g[0] + (dg * Tm).sum(axis=0)
        inner = np.minimum(inner * (1 + SUMPAD), 1.0)
        return float((inner * self.xm_hi).sum() * (1 + SUMPAD) + self.xtail)


# ---------------------------------------------------------------- the joint term -----------
def normal_cells(step_log2=-4, emax=9.0):
    """Uniform cells in e ~ N(0,1) on [-emax, emax]; certified masses; tail mass (both sides)."""
    step = 2.0 ** step_log2
    k = int(round(2 * emax / step))
    edges = -emax + step * np.arange(k + 1)
    r2 = arb(2).sqrt()
    cdf = [(-q(v) / r2).erfc() / 2 for v in edges]
    mhi = np.array([fhi(cdf[i + 1] - cdf[i]) for i in range(k)])
    mlo = np.array([max(flo(cdf[i + 1] - cdf[i]), 0.0) for i in range(k)])
    tail_hi = fhi(cdf[0] + (1 - cdf[-1]))
    return edges, mhi, mlo, tail_hi


def _tail12_lo(phi, h, k, rho, K=64, L=6.0):
    """Lower bound of P(Z1 > h, Z2 > k) for corr rho (float), arrays h, k (same shape)."""
    s = math.sqrt(1.0 - rho * rho)
    z = h[:, None] + (L / K) * np.arange(K + 1)[None, :]
    dm = np.maximum(phi.lo_(z[:, 1:]) - phi.up(z[:, :-1]), 0.0)
    zz = z[:, :-1] if rho >= 0 else z[:, 1:]
    cond = phi.lo_((rho * zz - k[:, None]) / s)
    return (dm * cond).sum(axis=1) * (1 - SUMPAD)


def _q_up(phi, h, k, rho, K=64):
    """Upper bound of P(Z1 <= h, Z2 > k) for corr rho (float), arrays h, k."""
    s = math.sqrt(1.0 - rho * rho)
    if rho > 0:
        zl = np.maximum((k - 7.0 * s) / rho, -9.0)
    else:
        zl = np.full_like(h, -9.0)
    zl = np.minimum(zl, h)
    z = zl[:, None] + ((h - zl) / K)[:, None] * np.arange(K + 1)[None, :]
    z[:, -1] = h
    dm = np.maximum(phi.up(z[:, 1:]) - phi.lo_(z[:, :-1]), 0.0)
    zz = z[:, 1:] if rho >= 0 else z[:, :-1]
    cond = phi.up((rho * zz - k[:, None]) / s)
    body = (dm * cond).sum(axis=1)
    if rho > 0:
        tail = phi.up(zl) * phi.up((rho * zl - k) / s)
    else:
        tail = phi.up(zl)
    return (body + tail) * (1 + SUMPAD)


class JointBound:
    """Certified lower bound of P12 = P(both reject) and upper bound of Q = P(t1 <= c < t2),
    each over a whole rho-interval, by enclosure on 3-D Bartlett cells."""

    def __init__(self, n, table, phi, nx=128, nq=128, e_log2=-4, K=64):
        self.n, self.nu, self.T, self.phi, self.K = n, n - 1, table, phi, K
        self.xe, self.xmh, self.xml, self.xt = chi_cells(n - 1, nx)
        self.qe, self.qmh, self.qml, self.qt = chi_cells(n - 2, nq)
        self.ee, self.emh, self.eml, self.et = normal_cells(e_log2)
        self.tail = (self.xt + self.qt + self.et) * (1 + 1e-12)

    def bounds(self, ra, rb, want="both", cut=1e-13):
        T, phi = self.T, self.phi
        rnu_lo, rnu_hi = math.sqrt(self.nu) * (1 - 1e-15), math.sqrt(self.nu) * (1 + 1e-15)
        if ra < 0 < rb:
            s_hi = 1.0
        else:
            s_hi = math.sqrt(1.0 - min(ra * ra, rb * rb)) * (1 + 1e-15)
        s_lo = math.sqrt(1.0 - max(ra * ra, rb * rb)) * (1 - 1e-15)
        qa, qb = self.qe[:-1][None, :, None], self.qe[1:][None, :, None]
        ea, eb = self.ee[:-1][None, None, :], self.ee[1:][None, None, :]
        P12, Q = 0.0, 0.0
        for i in range(len(self.xe) - 1):
            xa, xb = self.xe[i], self.xe[i + 1]
            rx_lo, rx_hi = min(ra * xa, ra * xb), max(rb * xa, rb * xb)
            se_lo = np.minimum(s_lo * ea, s_hi * ea)
            se_hi = np.maximum(s_lo * eb, s_hi * eb)
            m_lo = (rx_lo + se_lo) * np.ones_like(qa)
            m_hi = (rx_hi + se_hi) * np.ones_like(qa)
            m_lo = m_lo - 1e-14 * (1 + np.abs(m_lo))
            m_hi = m_hi + 1e-14 * (1 + np.abs(m_hi))
            q_lo, q_hi = s_lo * qa * (1 - 1e-15), s_hi * qb * (1 + 1e-15)
            q_lo, q_hi = np.broadcast_to(q_lo, m_lo.shape), np.broadcast_to(q_hi, m_lo.shape)
            straddle = (m_lo <= 0) & (m_hi >= 0)
            m2min = np.where(straddle, 0.0, np.minimum(m_lo * m_lo, m_hi * m_hi))
            m2max = np.maximum(m_lo * m_lo, m_hi * m_hi)
            R_lo = np.sqrt(m2min + q_lo * q_lo) * (1 - 1e-14)
            R_hi = np.sqrt(m2max + q_hi * q_hi) * (1 + 1e-14)
            with np.errstate(divide="ignore", invalid="ignore"):
                tl = np.where(m_lo >= 0, m_lo / q_hi, np.where(q_lo > 0, m_lo / q_lo, -np.inf))
                th = np.where(m_hi >= 0, np.where(q_lo > 0, m_hi / q_lo, np.inf), m_hi / q_hi)
            tl = tl - 1e-13 * (1 + np.abs(tl))
            th = th + 1e-13 * (1 + np.abs(th))
            Jlo = np.searchsorted(T.tau_hi, tl.ravel(), side="left")
            Jhi = np.searchsorted(T.tau_lo, th.ravel(), side="left")
            cmax, cmin = T.c_hi[Jlo], T.c_lo[Jhi]
            Mx_h, Mx_l = self.xmh[i], self.xml[i]
            Mh = (Mx_h * self.qmh[:, None] * self.emh[None, :]).ravel() * (1 + 1e-14)
            Ml = (Mx_l * self.qml[:, None] * self.eml[None, :]).ravel() * (1 - 1e-14)
            h_up = cmax * xb / rnu_lo
            if want in ("both", "p12"):
                k_up = cmax * R_hi.ravel() / rnu_lo
                crude = phi.up(-np.maximum(h_up, k_up))          # P12 <= min tail
                sel = (Ml * crude > cut)
                if sel.any():
                    P12 += float((Ml[sel] * _tail12_lo(phi, h_up[sel], k_up[sel], ra, K=self.K)).sum())
            if want in ("both", "q"):
                k_lo = cmin * R_lo.ravel() / rnu_hi
                crude = np.minimum(phi.up(h_up), phi.up(-k_lo))
                sel = (Mh * crude > cut)
                Q += float((Mh[~sel] * crude[~sel]).sum())
                if sel.any():
                    Q += float((Mh[sel] * _q_up(phi, h_up[sel], k_lo[sel], ra, K=self.K)).sum())
        P12 = P12 * (1 - SUMPAD)
        Q = Q * (1 + SUMPAD) + self.tail
        return P12, Q
