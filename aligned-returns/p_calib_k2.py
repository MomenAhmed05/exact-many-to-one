# -*- coding: utf-8 -*-
"""p_calib_k2.py — option B, Phase 0: the calibrated own-denominator rule at K = 2 (measured).

The desk estimates the correlation r_hat of the two difference series and applies the constant a
deployed table gives for r_hat's bin. The table holds the exact constant c*(r; a, n) of p_law at an
internal level a <= alpha, and a is calibrated so that the rule's family-wise error rate is at most
alpha at every correlation on a declared grid.

Bartlett in (x, m, q): x = sqrt(u), u ~ chi2_nu; m | x ~ N(rho x, s^2), s^2 = 1 - rho^2;
q = s sqrt(w), w ~ chi2_{nu-1}. Then a11 = x^2, a22 = m^2 + q^2 and r_hat = m / sqrt(m^2 + q^2),
so a table bin is an exact interval of m given q, and

    P(no rejection) = E_{x,q} sum_j  int_{m in bin j}  Phi2(c_j x/sqrt(nu) - d, c_j sqrt(m^2+q^2)/sqrt(nu); rho) dP(m | x)

with d = sqrt(n) * edge of candidate 1 (d = 0 under the null).

    python p_calib_k2.py gates
    python p_calib_k2.py sweep N
"""
import math
import os
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri
from scipy.stats import chi2, t as tdist

import p_law as P

ALPHA = 0.025
ZMAX = math.atanh(0.999)
B = 80
HERE = os.path.dirname(os.path.abspath(__file__))


def edges(B=B, zmax=ZMAX):
    return np.tanh(np.linspace(-zmax, zmax, B + 1))


def cstar(n, r, a):
    L = P.Law2(n, r, panels=8, order=8)
    return brentq(lambda c: L.F(c) - (1 - a), 0.3, 40.0, xtol=1e-11)


def table(n, a, B=B):
    """Constants for the B + 2 bins: Bonferroni at a on (-1, r_0], then c*(left edge)."""
    r = edges(B)
    cs = [float(tdist.isf(a / 2, n - 1))] + [cstar(n, x, a) for x in r]
    return r, np.array(cs)


class Eval:
    """Non-rejection probability of a table rule at one (n, rho)."""

    def __init__(self, n, rho, r_edges, panels=8, order=8, inner=10):
        self.n, self.nu, self.rho = n, n - 1, rho
        nu = self.nu
        s = math.sqrt(1 - rho * rho)
        x, w = P.score_rule(panels, order)
        w = w / w.sum()
        X = np.sqrt(chi2.ppf(ndtr(x), nu))
        Q = s * np.sqrt(chi2.ppf(ndtr(x), nu - 1))
        Xg, Qg = np.meshgrid(X, Q, indexing="ij")
        self.W = np.outer(w, w)                              # (Nx, Nq)
        # m-boundaries of the bins, as probability levels of m | x
        ext = np.concatenate([[-np.inf], r_edges, [np.inf]])
        with np.errstate(divide="ignore", invalid="ignore"):
            mb = Qg[..., None] * ext / np.sqrt(1 - ext ** 2)  # (Nx, Nq, B+3)
        mb[..., 0], mb[..., -1] = -np.inf, np.inf
        pb = ndtr((mb - rho * Xg[..., None]) / s)            # levels of m | x
        g, wg = leggauss(inner)
        lo, hi = pb[..., :-1], pb[..., 1:]                   # (Nx, Nq, Bins)
        pn = lo[..., None] + (hi - lo)[..., None] * (g + 1) / 2   # (Nx, Nq, Bins, inner)
        pn = np.clip(pn, 1e-300, 1 - 1e-16)
        m = rho * Xg[..., None, None] + s * ndtri(pn)
        self.wm = (hi - lo)[..., None] * wg / 2              # weights sum to the bin's mass
        self.hx = Xg / math.sqrt(nu)                          # x / sqrt(nu)
        self.hr = np.sqrt(m * m + Qg[..., None, None] ** 2) / math.sqrt(nu)
        self.mass = (hi - lo)                                 # (Nx, Nq, Bins)

    def nonreject(self, cs, d=0.0):
        cs = np.asarray(cs, float)
        h = cs[None, None, :, None] * self.hx[..., None, None] - d
        k = cs[None, None, :, None] * self.hr
        h = np.broadcast_to(h, k.shape)
        v = P.phi2(h.ravel(), k.ravel(), self.rho).reshape(k.shape)
        return float(((v * self.wm).sum(axis=(2, 3)) * self.W).sum())

    def fwer(self, cs):
        return 1.0 - self.nonreject(cs)

    def bin_mass(self):
        return (self.mass * self.W[..., None]).sum(axis=(0, 1))


def bin_index(rh, r_edges):
    return np.searchsorted(r_edges, rh, side="left")       # 0 = (-1, r_0], j = (r_{j-1}, r_j]


def brute_table(n, rho, r_edges, cs, N, seed, d=0.0, batch=100_000):
    rng = np.random.default_rng(seed)
    L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]))
    hits = done = 0
    while done < N:
        m = min(batch, N - done)
        D = rng.standard_normal((m, n, 2)) @ L.T
        Dc = D - D.mean(1, keepdims=True)
        a = (Dc ** 2).sum(1)
        rh = (Dc[:, :, 0] * Dc[:, :, 1]).sum(1) / np.sqrt(a[:, 0] * a[:, 1])
        c = cs[bin_index(rh, r_edges)]
        t = math.sqrt(n) * D.mean(1) / np.sqrt(a / (n - 1))
        t[:, 0] += d / np.sqrt(a[:, 0] / (n - 1))
        hits += int((t.max(1) > c).sum())
        done += m
    p = hits / N
    return p, math.sqrt(p * (1 - p) / N)


RHO_GRID = np.tanh(np.linspace(-ZMAX, ZMAX, 41))


def out(name):
    return os.path.join(HERE, "deposits", name)


def gates():
    lines = []

    def claim(cid, label, value):
        lines.append(f"CLAIM {cid}: {label} = {value}")
        print(lines[-1], flush=True)

    n = 12
    r = edges()
    for rho in (-0.5, 0.0, 0.9):
        E = Eval(n, rho, r)
        Lw = P.Law2(n, rho)
        dev = max(abs(E.nonreject(np.full(B + 2, c)) - Lw.F(c)) for c in (1.5, 2.5, 3.5))
        claim(f"G0-1-{rho:+.1f}", f"constant table vs Law2.F, n=12, rho={rho}, max |diff| over c in 1.5/2.5/3.5", f"{dev:.2e} ({'PASS' if dev < 1e-8 else 'FAIL'})")
    rr, cs = table(n, ALPHA)
    worst = 0.0
    for rho in (-0.9, 0.0, 0.5, 0.95):
        f1 = Eval(n, rho, rr).fwer(cs)
        f2 = Eval(n, rho, rr, panels=12, inner=16).fwer(cs)
        worst = max(worst, abs(f1 - f2))
    claim("G0-2", "refinement move (panels 8->12, inner 10->16), n=12, plug-in at alpha, max over 4 rho", f"{worst:.2e} ({'PASS' if worst < 2e-6 else 'FAIL'})")
    for k, rho in enumerate((-0.5, 0.5, 0.95)):
        f = Eval(n, rho, rr).fwer(cs)
        p, se = brute_table(n, rho, rr, cs, 1_000_000, 20261025 + k)
        z = (p - f) / se
        claim(f"G0-3-{rho:+.2f}", f"brute 1e6 vs evaluator, plug-in at alpha, n=12, rho={rho}", f"brute {p:.5f}, eval {f:.6f}, z {z:+.2f} ({'PASS' if abs(z) < 3 else 'FAIL'})")
    with open(out("p_calib_k2_gates_claims.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def sweep(n):
    lines = []

    def claim(cid, label, value):
        lines.append(f"CLAIM {cid}: {label} = {value}")
        print(lines[-1], flush=True)

    r = edges()
    def worst(cs):
        f = {float(rho): Eval(n, float(rho), r).fwer(cs) for rho in RHO_GRID}
        rho = max(f, key=f.get)
        return f[rho], rho, f

    _, c_a = table(n, ALPHA)
    w0, at0, f0 = worst(c_a)
    claim(f"PB0-{n}-PLUGIN", f"uncalibrated plug-in table at a = alpha, n={n}: worst FWER over rho grid (41 pts in [-0.999, 0.999])", f"{w0:.6f} at rho={at0:+.4f} ({w0 / ALPHA:.4f} alpha)")
    claim(f"PB0-{n}-PLUGIN-PROFILE", "FWER at rho = -0.999, -0.9, -0.5, 0, 0.5, 0.9, 0.999", ", ".join(f"{f0[min(f0, key=lambda x: abs(x - t))]:.6f}" for t in (-0.999, -0.9, -0.5, 0, 0.5, 0.9, 0.999)))
    lo, hi = 0.8 * ALPHA, ALPHA
    best = None
    for _ in range(9):
        a = 0.5 * (lo + hi)
        _, cs = table(n, a)
        w, at, _ = worst(cs)
        if w <= ALPHA:
            lo, best = a, (a, cs, w, at)
        else:
            hi = a
    a, cs, w, at = best
    claim(f"PB0-{n}-LEVEL", f"calibrated internal level (largest feasible after 9 bisection steps), n={n}", f"{a:.6f} ({a / ALPHA:.4f} alpha); worst FWER {w:.6f} at rho={at:+.4f}")
    cB = float(tdist.isf(ALPHA / 2, n - 1))
    for rho in (0.0, 0.45, 0.75, 0.9, 0.95, 0.99):
        E = Eval(n, rho, r)
        co = cstar(n, rho, ALPHA)
        powr = lambda c_arr, d: 1.0 - E.nonreject(c_arr, d)
        d80 = brentq(lambda d: powr(np.full(B + 2, co), d) - 0.8, 0.0, 10.0, xtol=1e-7)
        pc, pb = powr(cs, d80), powr(np.full(B + 2, cB), d80)
        claim(f"PB0-{n}-POW-{rho:.2f}", f"power at oracle-80% edge, n={n}, rho={rho}: oracle c*, Bonferroni c; calibrated / Bonferroni power; gain over Bonferroni (pp)",
              f"{co:.6f}, {cB:.6f}; {pc:.4f} / {pb:.4f}; {100 * (pc - pb):+.2f}")
    with open(out(f"p_calib_k2_sweep_n{n}_claims.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "gates":
        gates()
    else:
        sweep(int(sys.argv[2]))
