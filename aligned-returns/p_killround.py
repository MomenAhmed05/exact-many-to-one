# -*- coding: utf-8 -*-
"""p_killround.py — the paired kill round at K = 2 (Section 8).

Executes the eight kill tests declared before the run,
with the evaluator p_law.Law2 (gates E-0..E-3 passed: deposits/p_law_gates_claims.txt).
Writes deposits/p_killround_claims.txt ALWAYS.
    python p_killround.py
"""
import math
import sys
import time
import traceback

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import chi2, t as tdist

import p_law as P

ALPHA = 0.025
rows = []


def claim(cid, label, value):
    rows.append("CLAIM %s: %s = %s" % (cid, label, value))
    print(rows[-1], flush=True)


def cstar(n, rho, alpha=ALPHA):
    L = P.Law2(n, rho, panels=8, order=8)
    return brentq(lambda c: L.F(c) - (1 - alpha), 0.5, 30.0, xtol=1e-9)


RHO_GRID = np.tanh(np.linspace(np.arctanh(-0.995), np.arctanh(0.995), 81))


def table(n, alpha=ALPHA):
    vals = np.array([cstar(n, r, alpha) for r in RHO_GRID])
    sp = CubicSpline(RHO_GRID, vals)
    mids = 0.5 * (RHO_GRID[::8][:-1] + RHO_GRID[::8][1:])
    err = max(abs(sp(m) - cstar(n, m, alpha)) for m in mids)
    return sp, vals, err


def simulate_plugin(n, rho, sp, N, seed):
    """Realised FWER of the plug-in rule c*(rho_hat) under true correlation rho."""
    rng = np.random.default_rng(seed)
    L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]))
    hits = done = 0
    while done < N:
        m = min(100_000, N - done)
        D = rng.standard_normal((m, n, 2)) @ L.T
        Dc = D - D.mean(1, keepdims=True)
        s = np.sqrt((Dc ** 2).sum(1) / (n - 1))
        rh = np.clip((Dc[:, :, 0] * Dc[:, :, 1]).sum(1) / (n - 1) / (s[:, 0] * s[:, 1]), -0.995, 0.995)
        t = math.sqrt(n) * D.mean(1) / s
        hits += int((t.max(1) > sp(rh)).sum())
        done += m
    p = hits / N
    return p, math.sqrt(ALPHA * (1 - ALPHA) / N)


def c_mvt(n, rho):
    """Common-denominator bivariate t constant: P(T1<=c, T2<=c) = 1 - alpha, nu = n - 1."""
    nu = n - 1
    x, w = P.score_rule(64, 8)
    s = np.sqrt(chi2.ppf(ndtr(x), nu) / nu)
    f = lambda c: float(P.phi2(c * s, c * s, rho) @ w / w.sum()) - (1 - ALPHA)
    return brentq(f, 0.5, 30.0, xtol=1e-10)


def main():
    t0 = time.time()
    # tables of c*(rho) for the plug-in tests
    tabs = {}
    for n in (12, 36, 120):
        sp, vals, err = table(n)
        tabs[n] = (sp, vals)
        claim("PKR-TAB-n%d" % n, "c*(rho) table on 81 tanh-spaced rho; max interpolation error at 10 midpoints",
              "%.2e" % err)
    sp9, _, err9 = table(36, 0.9 * ALPHA)
    claim("PKR-TAB-n36-a0", "c*(rho) table at alpha0 = 0.9 alpha, n = 36; interpolation error", "%.2e" % err9)

    # 1. plug-in at n = 12
    res1 = {}
    for j, rho in enumerate((-0.5, 0.0, 0.5, 0.9)):
        p, se = simulate_plugin(12, rho, tabs[12][0], 400_000, 81100 + j)
        res1[rho] = p
        claim("PKR-1-r%+.1f" % rho, "plug-in realised FWER, n=12 (KILL unless >1.05a and >a+3SE)",
              "%.5f (z %+.1f)" % (p, (p - ALPHA) / se))
    surv1 = any(p > 1.05 * ALPHA and p > ALPHA + 3 * math.sqrt(ALPHA * (1 - ALPHA) / 400_000)
                for p in res1.values())
    claim("PKR-1-VERDICT", "conjecture 1 (plug-in liberal by >5% at n=12)", "SURVIVES" if surv1 else "KILLED")

    # 2. 1/n shrinkage at the worst rho of test 1
    rw = max(res1, key=res1.get)
    ex = {12: res1[rw] - ALPHA}
    for j, n in enumerate((36, 120)):
        p, se = simulate_plugin(n, rw, tabs[n][0], 400_000, 81200 + j)
        ex[n] = p - ALPHA
        claim("PKR-2-n%d" % n, "plug-in realised FWER at worst rho %.1f" % rw, "%.5f (z %+.1f)" % (p, (p - ALPHA) / se))
    r1 = ex[12] / ex[36] if ex[36] != 0 else float("inf")
    r2 = ex[36] / ex[120] if ex[120] != 0 else float("inf")
    surv2 = 1.5 <= r1 <= 6 and 1.5 <= r2 <= 6
    claim("PKR-2-VERDICT", "conjecture 2 (excess ~ 1/n): ratios 12/36, 36/120 in [1.5, 6]",
          "%.2f, %.2f -> %s" % (r1, r2, "SURVIVES" if surv2 else "KILLED"))

    # 3. common-denominator multivariate t at the TRUE R, n = 12
    worst3 = 0.0
    for rho in (-0.5, 0.0, 0.5, 0.9):
        c = c_mvt(12, rho)
        rate = 1 - P.Law2(12, rho).F(c)
        worst3 = max(worst3, rate)
        claim("PKR-3-r%+.1f" % rho, "multivariate-t constant and its TRUE rate, n=12",
              "c_mvt %.4f, rate %.5f (%.2fx)" % (c, rate, rate / ALPHA))
    claim("PKR-3-VERDICT", "conjecture 3 (mvt liberal at true R, n=12)",
          "SURVIVES" if worst3 > 1.05 * ALPHA else "KILLED")

    # 4 / 5. where c*(rho) is largest
    grid = (-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 0.95)
    arg_is_zero, neg_beats = True, False
    for n in (12, 36):
        cs = {r: cstar(n, r) for r in grid}
        claim("PKR-45-n%d" % n, "c*(rho) on the declared grid", "; ".join("%+.2f:%.6f" % kv for kv in cs.items()))
        top = max(cs, key=cs.get)
        if top != 0.0 and cs[top] - cs[0.0] > 1e-6:
            arg_is_zero = False
        if any(cs[r] - cs[0.0] > 1e-6 for r in grid if r < 0):
            neg_beats = True
    claim("PKR-4-VERDICT", "conjecture 4 (argmax at rho=0)", "SURVIVES" if arg_is_zero else "KILLED")
    claim("PKR-5-VERDICT", "conjecture 5 (a negative rho beats rho=0)", "SURVIVES" if neg_beats else "KILLED")

    # 6. Sidak product, two-sided, own denominators, K = 2
    worst6, at6 = 1.0, None
    for n in (2, 3, 5, 12, 36):
        for rho in (-0.99, -0.9, -0.5, -0.1, 0.1, 0.5, 0.9, 0.99):
            L = P.Law2(n, rho)
            for c in (0.5, 1.0, 2.0, 3.0, 5.0):
                d = L.rect(c) - (2 * tdist.cdf(c, n - 1) - 1) ** 2
                if d < worst6:
                    worst6, at6 = d, (n, rho, c)
    claim("PKR-6-NEARMISS", "min over grid of P(|t1|,|t2|<=c) - P(|t|<=c)^2 (200 configs)",
          "%.3e at n=%d rho=%+.2f c=%.1f" % ((worst6,) + at6))
    claim("PKR-6-VERDICT", "conjecture 6 (Sidak product holds, own denominators, K=2)",
          "SURVIVES" if worst6 >= -1e-6 else "KILLED")

    # 7. Bonferroni vs oracle power, n = 36
    n = 36
    cB = tdist.isf(ALPHA / 2, n - 1)
    surv7 = True
    for j, rho in enumerate((0.0, 0.45)):
        co = cstar(n, rho)
        rng = np.random.default_rng(81700 + j)
        L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]))
        N = 200_000
        D = rng.standard_normal((N, n, 2)) @ L.T
        Dc = D - D.mean(1, keepdims=True)
        s = np.sqrt((Dc ** 2).sum(1) / (n - 1))
        base = math.sqrt(n) * D.mean(1)
        def powers(delta):
            t = (base + math.sqrt(n) * np.array([delta, 0.0])) / s
            m = t.max(1)
            return (m > co).mean(), (m > cB).mean()
        dl = brentq(lambda d: powers(d)[0] - 0.8, 0.0, 2.0, xtol=1e-6)
        po, pb = powers(dl)
        loss = po - pb
        se = math.sqrt(2 * 0.8 * 0.2 / N)
        if loss >= 0.02 + 3 * se:
            surv7 = False
        claim("PKR-7-r%.2f" % rho, "oracle c*, Bonferroni c, effect, power oracle/Bonferroni, loss",
              "%.4f, %.4f, %.4f, %.4f/%.4f, %.4f" % (co, cB, dl, po, pb, loss))
    claim("PKR-7-VERDICT", "conjecture 7 (Bonferroni loss < 2pp)", "SURVIVES" if surv7 else "KILLED")

    # 8. calibrated internal level 0.9 alpha, n = 36
    surv8 = True
    for j, rho in enumerate((-0.5, 0.0, 0.5, 0.9, 0.95)):
        p, se = simulate_plugin(36, rho, sp9, 400_000, 81800 + j)
        if p > ALPHA + 3 * se:
            surv8 = False
        claim("PKR-8-r%+.2f" % rho, "plug-in at alpha0 = 0.9a, n=36, realised FWER", "%.5f (z %+.1f)" % (p, (p - ALPHA) / se))
    claim("PKR-8-VERDICT", "conjecture 8 (internal level >= 0.9a feasible at n=36)", "SURVIVES" if surv8 else "KILLED")
    claim("PKR-WALL", "seconds", "%.0f" % (time.time() - t0))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        claim("PKR-ERROR", "the run failed", traceback.format_exc().replace("\n", " | "))
    finally:
        open("deposits/p_killround_claims.txt", "w", encoding="utf-8").write("\n".join(rows) + "\n")
