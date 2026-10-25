# -*- coding: utf-8 -*-
"""pb_joint.py — option B: adaptive certified enclosure of the joint term over a rho-interval.

P12 = P(t1 > c(r_hat), t2 > c(r_hat)) is bounded below, or Q = P(t1 <= c(r_hat) < t2) above.
Cells are index boxes on fine dyadic grids of the Bartlett coordinates x, q', e (pb_cert), so a
cell's mass is an exact difference of certified CDF values. A cell's bound uses the whole
rho-interval: the geometry by interval arithmetic, the normal correlation by monotonicity
(P12 is increasing and Q decreasing in the correlation of (Z1, Z2) at fixed h, k). The bound is
therefore valid at every rho in the interval.
"""
import math

import numpy as np
from flint import arb

from pb_cert import SUMPAD, fhi, flo, q


class FineCDF:
    """Certified CDF of one Bartlett coordinate on a dyadic grid v_i = v0 + i 2^L."""

    def __init__(self, kind, df=None, L=-9):
        from scipy.stats import chi2
        self.step = 2.0 ** L
        if kind == "chi":
            vmax = math.sqrt(chi2.isf(1e-16, df))
            self.v0 = 0.0
            nn = int(math.ceil(vmax / self.step))
            v = np.arange(nn + 1) * self.step
            s = arb(df) / 2
            vals = [(q(x) * q(x) / 2).gamma_lower(s, regularized=1) if x > 0 else arb(0) for x in v]
            self.tail_hi = fhi(1 - vals[-1])
        else:
            self.v0 = -9.0
            nn = int(round(18.0 / self.step))
            v = self.v0 + np.arange(nn + 1) * self.step
            r2 = arb(2).sqrt()
            vals = [(-q(x) / r2).erfc() / 2 for x in v]
            self.tail_hi = fhi(1 - vals[-1]) + fhi(vals[0])
        self.v = v
        self.lo = np.array([max(flo(a), 0.0) for a in vals])
        self.hi = np.array([min(fhi(a), 1.0) for a in vals])

    def mass(self, ia, ib):
        return np.maximum(self.lo[ib] - self.hi[ia], 0.0), self.hi[ib] - self.lo[ia]


def tail12(phi, h, k, rho, dz, L, upper):
    """Bound of P(Z1 > h, Z2 > k) at correlation rho: an upper bound if `upper`, else a lower."""
    s = math.sqrt(1.0 - rho * rho)
    K = int(math.ceil(L / dz))
    z = h[:, None] + dz * np.arange(K + 1)[None, :]
    if upper:
        dm = np.maximum(phi.up(z[:, 1:]) - phi.lo_(z[:, :-1]), 0.0)
        zz = z[:, 1:] if rho >= 0 else z[:, :-1]
        val = (dm * phi.up((rho * zz - k[:, None]) / s)).sum(axis=1) + phi.up(-z[:, -1])
        return val * (1 + SUMPAD)
    dm = np.maximum(phi.lo_(z[:, 1:]) - phi.up(z[:, :-1]), 0.0)
    zz = z[:, :-1] if rho >= 0 else z[:, 1:]
    return (dm * phi.lo_((rho * zz - k[:, None]) / s)).sum(axis=1) * (1 - SUMPAD)


def qterm(phi, h, k, rho, dz, upper):
    """Bound of P(Z1 <= h, Z2 > k) at correlation rho >= 0."""
    assert rho >= 0
    s = math.sqrt(1.0 - rho * rho)
    zl = np.maximum((k - 7.0 * s) / rho, -9.0) if rho > 0 else np.full_like(h, -9.0)
    zl = np.minimum(zl, h)
    K = max(int(math.ceil(float(np.max(h - zl)) / dz)), 1)
    w = (h - zl) / K
    z = zl[:, None] + w[:, None] * np.arange(K + 1)[None, :]
    z[:, -1] = h
    if upper:
        dm = np.maximum(phi.up(z[:, 1:]) - phi.lo_(z[:, :-1]), 0.0)
        val = (dm * phi.up((rho * z[:, 1:] - k[:, None]) / s)).sum(axis=1)
        val += phi.up(zl) * (phi.up((rho * zl - k) / s) if rho > 0 else 1.0)
        return val * (1 + SUMPAD)
    dm = np.maximum(phi.lo_(z[:, 1:]) - phi.up(z[:, :-1]), 0.0)
    return (dm * phi.lo_((rho * z[:, :-1] - k[:, None]) / s)).sum(axis=1) * (1 - SUMPAD)


class JointAdaptive:
    def __init__(self, n, table, phi, L=-9, start=(-3, -3, -2)):
        self.n, self.nu, self.T, self.phi = n, n - 1, table, phi
        self.X = FineCDF("chi", n - 1, L)
        self.Qd = FineCDF("chi", n - 2, L)
        self.E = FineCDF("normal", None, L)
        edges = []
        for F, sl in zip((self.X, self.Qd, self.E), start):
            w = int(round(2.0 ** sl / F.step))
            e = np.arange(0, len(F.v) - 1, w)
            edges.append(np.append(e, len(F.v) - 1))
        g = np.meshgrid(*[np.arange(len(e) - 1) for e in edges], indexing="ij")
        self.cells = np.stack([edges[0][g[0].ravel()], edges[0][g[0].ravel() + 1],
                               edges[1][g[1].ravel()], edges[1][g[1].ravel() + 1],
                               edges[2][g[2].ravel()], edges[2][g[2].ravel() + 1]], axis=1)
        self.tail = (self.X.tail_hi + self.Qd.tail_hi + self.E.tail_hi) * (1 + 1e-12)

    def hk(self, C, ra, rb):
        X, Qd, E, T = self.X, self.Qd, self.E, self.T
        xa, xb = X.v[C[:, 0]], X.v[C[:, 1]]
        qa, qb = Qd.v[C[:, 2]], Qd.v[C[:, 3]]
        ea, eb = E.v[C[:, 4]], E.v[C[:, 5]]
        s_hi = 1.0 if ra < 0 < rb else math.sqrt(1.0 - min(ra * ra, rb * rb)) * (1 + 1e-15)
        s_lo = math.sqrt(1.0 - max(ra * ra, rb * rb)) * (1 - 1e-15)
        m_lo = np.minimum(ra * xa, ra * xb) + np.minimum(s_lo * ea, s_hi * ea)
        m_hi = np.maximum(rb * xa, rb * xb) + np.maximum(s_lo * eb, s_hi * eb)
        m_lo = m_lo - 1e-14 * (1 + np.abs(m_lo))
        m_hi = m_hi + 1e-14 * (1 + np.abs(m_hi))
        q_lo, q_hi = s_lo * qa * (1 - 1e-15), s_hi * qb * (1 + 1e-15)
        strad = (m_lo <= 0) & (m_hi >= 0)
        m2min = np.where(strad, 0.0, np.minimum(m_lo * m_lo, m_hi * m_hi))
        m2max = np.maximum(m_lo * m_lo, m_hi * m_hi)
        R_lo = np.sqrt(m2min + q_lo * q_lo) * (1 - 1e-14)
        R_hi = np.sqrt(m2max + q_hi * q_hi) * (1 + 1e-14)
        with np.errstate(divide="ignore", invalid="ignore"):
            tl = np.where(m_lo >= 0, m_lo / q_hi, np.where(q_lo > 0, m_lo / q_lo, -np.inf))
            th = np.where(m_hi >= 0, np.where(q_lo > 0, m_hi / q_lo, np.inf), m_hi / q_hi)
        tl = tl - 1e-13 * (1 + np.abs(tl))
        th = th + 1e-13 * (1 + np.abs(th))
        cmax = T.c_hi[np.searchsorted(T.tau_hi, tl, side="left")]
        cmin = T.c_lo[np.searchsorted(T.tau_lo, th, side="left")]
        rl, rh = math.sqrt(self.nu) * (1 - 1e-15), math.sqrt(self.nu) * (1 + 1e-15)
        h_lo, h_hi = cmin * xa / rh, cmax * xb / rl
        k_lo, k_hi = cmin * R_lo / rh, cmax * R_hi / rl
        mx, mq, me = X.mass(C[:, 0], C[:, 1]), Qd.mass(C[:, 2], C[:, 3]), E.mass(C[:, 4], C[:, 5])
        mlo = mx[0] * mq[0] * me[0] * (1 - 1e-14)
        mhi = mx[1] * mq[1] * me[1] * (1 + 1e-14)
        return h_lo, h_hi, k_lo, k_hi, mlo, mhi

    def cell_bounds(self, C, ra, rb, mode, lat_lo, lat_up):
        phi = self.phi
        h_lo, h_hi, k_lo, k_hi, mlo, mhi = self.hk(C, ra, rb)
        if mode == "p12":
            lo = mlo * lat_lo.lower(h_hi, k_hi)
            up = mhi * lat_up.upper(h_lo, k_lo)
        else:                                   # Q = P(Z2 > k) - P(Z1 > h, Z2 > k)
            up = mhi * np.maximum(phi.up(-k_lo) - lat_lo.lower(h_hi, k_lo), 0.0) * (1 + SUMPAD)
            lo = mlo * np.maximum(phi.lo_(-k_hi) - lat_up.upper(h_lo, k_hi), 0.0) * (1 - SUMPAD)
        return lo, up

    def run(self, ra, rb, mode, need, max_rounds=14, max_cells=3_000_000, log=None):
        """Refine until the bound meets `need` (P12 lower >= need, or Q upper <= need) or the
        budget ends. Returns (bound, cells, rounds); the final cells seed the next interval."""
        lat_lo, lat_up = BivLattice(self.phi, ra, False), BivLattice(self.phi, rb, True)
        C = self.cells
        lo, up = self.cell_bounds(C, ra, rb, mode, lat_lo, lat_up)
        for rnd in range(max_rounds + 1):
            if mode == "p12":
                val = float(lo.sum()) * (1 - SUMPAD)
                ok = val >= need
            else:
                val = float(up.sum()) * (1 + SUMPAD) + self.tail
                ok = val <= need
            if log:
                log(f"   round {rnd}: cells {len(C)}, bound {val:.7f}, need {need:.7f}, gap {float((up - lo).sum()):.2e}")
            if ok or rnd == max_rounds or len(C) > max_cells:
                break
            gap = up - lo
            order = np.argsort(-gap)
            cum = np.cumsum(gap[order])
            sel = order[:int(np.searchsorted(cum, 0.6 * cum[-1])) + 1]
            sel = sel[((C[sel, 1] - C[sel, 0]) > 1) | ((C[sel, 3] - C[sel, 2]) > 1) | ((C[sel, 5] - C[sel, 4]) > 1)]
            if len(sel) == 0:
                break
            P = C[sel]
            for ax in range(3):
                a, b = P[:, 2 * ax], P[:, 2 * ax + 1]
                mid = (a + b) // 2
                can = (b - a) > 1
                left, right = P.copy(), P.copy()
                left[can, 2 * ax + 1] = mid[can]
                right[can, 2 * ax] = mid[can]
                P = np.concatenate([left, right[can]])
            keep = np.ones(len(C), bool)
            keep[sel] = False
            nlo, nup = self.cell_bounds(P, ra, rb, mode, lat_lo, lat_up)
            C = np.concatenate([C[keep], P])
            lo = np.concatenate([lo[keep], nlo])
            up = np.concatenate([up[keep], nup])
        self.cells = C
        return val, len(C), rnd


class BivLattice:
    """Certified bound of G(h, k) = P(Z1 > h, Z2 > k), corr rho, on the lattice h, k = i d, d = 2^L.

    Segment i covers Z1 in [h_i, h_{i+1}]; on it the conditional P(Z2 > k | Z1 = z) = Phi((rho z - k)/s)
    is monotone in z, so its value at one end bounds it. G(h_i, k) is the reverse cumulative sum of the
    segment bounds (plus, for the upper bound, the mass beyond the last lattice point)."""

    def __init__(self, phi, rho, upper, L=-8, H=8.5):
        self.phi, self.upper_, self.d = phi, upper, 2.0 ** L
        N = int(round(H / self.d))
        h = self.d * np.arange(N + 1)
        s = math.sqrt(1.0 - rho * rho)
        if upper:
            mass = np.maximum(phi.up(h[1:]) - phi.lo_(h[:-1]), 0.0)
            zend = h[1:] if rho >= 0 else h[:-1]
            cond = phi.up((rho * zend[:, None] - h[None, :]) / s)
        else:
            mass = np.maximum(phi.lo_(h[1:]) - phi.up(h[:-1]), 0.0)
            zend = h[:-1] if rho >= 0 else h[1:]
            cond = phi.lo_((rho * zend[:, None] - h[None, :]) / s)
        seg = mass[:, None] * cond                                   # (N, N+1): segment i, k index l
        G = np.zeros((N + 1, N + 1))
        G[:N] = np.cumsum(seg[::-1], axis=0)[::-1]
        if upper:
            G = G * (1 + SUMPAD) + phi.up(-h[-1:])[0]
        else:
            G = G * (1 - SUMPAD)
        self.G, self.N = G, N

    def lower(self, h, k):
        """Lower bound of G at arrays h, k (G decreasing in both: round the indices up)."""
        i = np.ceil(h * (1 + 1e-12) / self.d + 1e-9)
        j = np.ceil(k * (1 + 1e-12) / self.d + 1e-9)
        out = np.zeros(len(h))
        ok = (i <= self.N) & (j <= self.N)
        out[ok] = self.G[i[ok].astype(np.int64), j[ok].astype(np.int64)]
        return out

    def upper(self, h, k):
        """Upper bound of G (round the indices down; beyond the lattice use its last point)."""
        i = np.clip(np.floor(h * (1 - 1e-12) / self.d - 1e-9), 0, self.N).astype(np.int64)
        j = np.clip(np.floor(k * (1 - 1e-12) / self.d - 1e-9), 0, self.N).astype(np.int64)
        return self.G[i, j]
