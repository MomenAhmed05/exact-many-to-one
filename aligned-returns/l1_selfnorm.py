# -*- coding: utf-8 -*-
"""l1_selfnorm.py — L-1, Shao's self-normalised test under AR(1) returns.

The self-normalised test on Table 31's twelve cells (paper 1, Section 11.5).
Declarations D-1..D-6 (seeds, anchors, tolerances) were fixed before the run.

Stages (each writes deposits/l1_selfnorm_<stage>.json and prints CLAIM lines):
  uq       F-1  the upper-0.05 point q of U1 by two routes
  large    F-2  size at n = 4000, phi in {0, 0.6}, equal lengths
  anchors  F-3  Wang-Zhang Table 1 (S15) and Table 3 (WSNS) cells
  table         the twelve Table-31 cells

    python l1_selfnorm.py <stage>
"""
import hashlib
import json
import math
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
DEP = ROOT / "deposits"
sys.path.insert(0, str(ROOT.parent / "04-efficient-rule"))
ALPHA = 0.025
P_UPPER = 2 * ALPHA            # D-3: one-sided level alpha uses the upper-2alpha point


def claim(rows, cid, label, value):
    rows.append("CLAIM %s: %s = %s" % (cid, label, value))
    print(rows[-1], flush=True)


def save(stage, obj, rows):
    obj["script_sha256"] = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    (DEP / ("l1_selfnorm_%s.json" % stage)).write_text(json.dumps(obj, indent=1), encoding="utf-8")
    (DEP / ("l1_selfnorm_%s_claims.txt" % stage)).write_text("\n".join(rows) + "\n", encoding="utf-8")


# ----------------------------------------------------------------- statistics

def sn_stat(d):
    """Shao's T = n Dbar^2 / W_n, W_n = n^-2 sum_t (S_t - t Dbar)^2, per row."""
    n = d.shape[1]
    db = d.mean(1)
    s = np.cumsum(d - db[:, None], axis=1)
    w = (s * s).sum(1) / n ** 2
    return n * db * db / w, db


def q_from_uq():
    return json.loads((DEP / "l1_selfnorm_uq.json").read_text())["q_of_record"]


class Warped:
    """Wang-Zhang warped self-normalised statistic (their (3.1)), vectorised over rows.

    X is observed at times k = 1..m (array index s-1), Y at times l..n (array index s-l).
    """

    def __init__(self, m, l, n):
        self.m, self.l, self.n = m, l, n
        self.Nx, self.Ny = m, n - l + 1
        self.N = self.Nx + self.Ny

    def sets(self, i, j):
        """Index sets I_{x,i,j}, I_{y,i,j} as (lo, hi) times; None if empty (D-1)."""
        l, N = self.l, self.N
        xlo = l + i - 1 - (j * (l - 1)) // N
        xhi = l - (j * (l - 1)) // N + (j * self.Nx) // N - 1 + i - 1
        ylo = l + i - 1
        yhi = l + (j * self.Ny) // N - 1 + i - 1
        if xhi < xlo or yhi < ylo:
            return None
        return xlo, xhi, ylo, yhi

    def MN(self, B):
        l, N = self.l, self.N
        return min(self.m - l + 2 + (B * (l - 1)) // N - (B * self.Nx) // N,
                   self.n - l + 2 - (B * self.Ny) // N)

    def _wmeans(self, CX, CY, xlo, xhi, ylo, yhi):
        # X time s -> cumsum index s; Y time s -> cumsum index s - l + 1
        mx = (CX[:, xhi] - CX[:, xlo - 1]) / (xhi - xlo + 1)
        my = (CY[:, yhi - self.l + 1] - CY[:, ylo - self.l]) / (yhi - ylo + 1)
        return mx - my

    def stat(self, CX, CY, i, B, delta):
        """T_{i,B}(delta) per row; delta is an array (per row) or a scalar."""
        full = self.sets(i, B)
        assert full is not None
        dB = self._wmeans(CX, CY, *full)
        den = np.zeros(CX.shape[0])
        for j in range(1, B + 1):
            s = self.sets(i, j)
            if s is None:
                continue
            dj = self._wmeans(CX, CY, *s)
            den += j * j * (dj - dB) ** 2
        den /= B * B
        return B * (dB - delta) ** 2 / den, dB

    @staticmethod
    def cums(X, Y):
        z = np.zeros((X.shape[0], 1))
        return np.hstack([z, np.cumsum(X, 1)]), np.hstack([z, np.cumsum(Y, 1)])

    def subsample_q(self, CX, CY, B, dhat, p):
        vals = np.column_stack([self.stat(CX, CY, i, B, dhat)[0]
                                for i in range(1, self.MN(B) + 1)])
        return np.quantile(vals, p, axis=1, method="linear")       # D-2


# ----------------------------------------------------------------- data

def ar1_stationary(R, n, rho, rng):
    """AR(1) with standard normal innovations, stationary start (Wang-Zhang 4.1)."""
    x = np.empty((R, n))
    x[:, 0] = rng.standard_normal(R) / math.sqrt(1 - rho * rho)
    e = rng.standard_normal((R, n))
    for t in range(1, n):
        x[:, t] = rho * x[:, t - 1] + e[:, t]
    return x


# ----------------------------------------------------------------- stages

def stage_uq():
    rows, out = [], {}
    N, NB = 10_000_000, 10
    res = {}
    for tag, m, seed in (("brownian_m1000", 1000, 61001), ("brownian_m4000", 4000, 61002)):
        rng = np.random.default_rng(seed)
        qs, allv = [], []
        per = N // NB
        for b in range(NB):
            v = []
            done = 0
            while done < per:
                r = min(2500 if m == 4000 else 10000, per - done)
                t, _ = sn_stat(rng.standard_normal((r, m)))
                v.append(t)
                done += r
            v = np.concatenate(v)
            qs.append(np.quantile(v, 1 - P_UPPER))
            allv.append(v)
        allv = np.concatenate(allv)
        res[tag] = {"q": float(np.quantile(allv, 1 - P_UPPER)),
                    "se": float(np.std(qs, ddof=1) / math.sqrt(NB)), "batches": [float(x) for x in qs]}
        claim(rows, "L1UQ-%s" % tag, "upper-0.05 point of U1, N=1e7, batch SE",
              "%.4f (se %.4f)" % (res[tag]["q"], res[tag]["se"]))
    # Karhunen-Loeve: U1 = Z0^2 / (sum_{k<=K} Zk^2/(k pi)^2 + tail mean)
    K = 2000
    lam = 1.0 / (np.arange(1, K + 1) * math.pi) ** 2
    tail = 1.0 / 6.0 - lam.sum()                     # sum over k > K of 1/(k pi)^2, exactly
    rng = np.random.default_rng(61003)
    qs, allv = [], []
    per = N // NB
    for b in range(NB):
        v = []
        done = 0
        while done < per:
            r = min(5000, per - done)
            z = rng.standard_normal((r, K))
            den = (z * z) @ lam + tail
            v.append(rng.standard_normal(r) ** 2 / den)
            done += r
        v = np.concatenate(v)
        qs.append(np.quantile(v, 1 - P_UPPER))
        allv.append(v)
    allv = np.concatenate(allv)
    res["kl_K2000"] = {"q": float(np.quantile(allv, 1 - P_UPPER)),
                       "se": float(np.std(qs, ddof=1) / math.sqrt(NB)), "tail_mean": tail}
    claim(rows, "L1UQ-kl_K2000", "upper-0.05 point of U1 by Karhunen-Loeve, N=1e7, batch SE",
          "%.4f (se %.4f)" % (res["kl_K2000"]["q"], res["kl_K2000"]["se"]))
    a, b4, k = res["brownian_m1000"], res["brownian_m4000"], res["kl_K2000"]
    agree = abs(b4["q"] - k["q"]) <= 3 * math.hypot(b4["se"], k["se"])
    refine = abs(b4["q"] - a["q"]) <= b4["se"]
    claim(rows, "L1UQ-F1-AGREE", "Brownian(m=4000) vs KL within 3 combined SE",
          "diff %.4f, 3SE %.4f, pass %s" % (b4["q"] - k["q"], 3 * math.hypot(b4["se"], k["se"]), agree))
    claim(rows, "L1UQ-F1-REFINE", "Brownian m=1000 -> 4000 moves q by less than one SE",
          "move %.4f, SE %.4f, pass %s" % (b4["q"] - a["q"], b4["se"], refine))
    q = float((b4["q"] / b4["se"] ** 2 + k["q"] / k["se"] ** 2) / (1 / b4["se"] ** 2 + 1 / k["se"] ** 2))
    claim(rows, "L1UQ-Q", "q of record (inverse-variance mean of the two routes)", "%.4f" % q)
    claim(rows, "L1UQ-Q-T", "equivalent one-sided t-type critical value sqrt(q)", "%.4f" % math.sqrt(q))
    out.update(routes=res, q_of_record=q, f1_agree=agree, f1_refine=refine)
    save("uq", out, rows)


def stage_large():
    rows = []
    from w3_asymptotic import ar1
    q = q_from_uq()
    n, R = 4000, 400_000
    rng = np.random.default_rng(63001)
    res = {}
    for phi in (0.0, 0.6):
        hits = done = 0
        while done < R:
            r = min(5000, R - done)
            b = ar1(r, n, phi, 1.0, rng)
            c = ar1(r, n, phi, 1.0, rng)
            t, db = sn_stat(c - b)
            hits += int(((db > 0) & (t > q)).sum())
            done += r
        p = hits / R
        se = math.sqrt(ALPHA * (1 - ALPHA) / R)
        z = (p - ALPHA) / se
        res[str(phi)] = {"size": p, "z": z, "pass": abs(z) <= 4}
        claim(rows, "L1F2-phi%.1f" % phi, "size at n=4000, N=4e5 (pass |z|<=4)",
              "%.5f (z %+.2f, pass %s)" % (p, z, abs(z) <= 4))
    save("large", {"q": q, "cells": res}, rows)


def stage_anchors():
    rows = []
    q = q_from_uq()
    R = 4000
    res = {}
    # S15, Table 1, Scenario 1: k = l = 1, m = n = 1500, r = 0
    rng = np.random.default_rng(62001)
    cov = 0
    for _ in range(R // 500):
        X = ar1_stationary(500, 1500, 0.3, rng)
        Y = ar1_stationary(500, 1500, 0.3, rng)
        t, _ = sn_stat(X - Y)
        cov += int((t <= q).sum())
    p = cov / R
    tol = 3 * math.sqrt(0.951 * 0.049 * (1 / 1000 + 1 / R))
    res["S15_T1"] = {"published": 0.951, "ours": p, "tol": tol, "pass": abs(p - 0.951) <= tol}
    claim(rows, "L1F3-S15", "S15 coverage, WZ Table 1 Scen.1 r=0 mean 95% (published 0.951)",
          "%.4f, |diff| %.4f, tol %.4f, pass %s" % (p, abs(p - 0.951), tol, abs(p - 0.951) <= tol))
    # WSNS, Table 3, Scenario 3: k = 1, l = 751, m = n = 2250, B_N = 100, r = 0
    W = Warped(2250, 751, 2250)
    rng = np.random.default_rng(62002)
    cov = 0
    for _ in range(R // 100):
        X = ar1_stationary(100, 2250, 0.3, rng)
        Y = ar1_stationary(100, 1500, 0.3, rng)
        CX, CY = W.cums(X, Y)
        T, dhat = W.stat(CX, CY, 1, W.N, 0.0)
        qh = W.subsample_q(CX, CY, 100, dhat, 0.95)
        cov += int((T <= qh).sum())
    p = cov / R
    tol = 3 * math.sqrt(0.945 * 0.055 * (1 / 1000 + 1 / R))
    res["WSNS_T3"] = {"published": 0.945, "ours": p, "tol": tol, "pass": abs(p - 0.945) <= tol,
                      "MN": W.MN(100)}
    claim(rows, "L1F3-WSNS", "WSNS coverage, WZ Table 3 Scen.3 r=0 mean 95%, B_N=100 (published 0.945)",
          "%.4f, |diff| %.4f, tol %.4f, MN %d, pass %s"
          % (p, abs(p - 0.945), tol, W.MN(100), abs(p - 0.945) <= tol))
    save("anchors", {"q": q, "anchors": res}, rows)


def stage_table():
    rows = []
    from w3_asymptotic import ar1
    q = q_from_uq()
    R = 400_000
    se = math.sqrt(ALPHA * (1 - ALPHA) / R)
    cells = []
    k = 0
    W = Warped(36, 24, 36)                                    # D-4
    grid = [B for B in range(4, 25)
            if W.sets(1, B) is not None and W.MN(B) >= 8]
    claim(rows, "L1T-BGRID", "feasible B_N grid for (36,13) (both windows nonempty, MN>=8)",
          "%s (MN %s)" % (grid, [W.MN(B) for B in grid]))
    for dtag, n0, n1 in (("n36", 36, 36), ("n120", 120, 120), ("n36_13", 36, 13)):
        for phi in (0.0, 0.2, 0.4, 0.6):
            rng = np.random.default_rng(64001 + k)
            k += 1
            if n0 == n1:
                hits = done = 0
                while done < R:
                    r = min(50000, R - done)
                    b = ar1(r, n0, phi, 1.0, rng) + 2.0
                    c = ar1(r, n1, phi, 1.0, rng) + 2.0
                    t, db = sn_stat(c - b)
                    hits += int(((db > 0) & (t > q)).sum())
                    done += r
                p = hits / R
                cells.append({"design": dtag, "phi": phi, "method": "self-normalised (Shao)",
                              "size": p, "mcse": se, "z": (p - ALPHA) / se})
                claim(rows, "L1T-%s-%.1f" % (dtag, phi), "size, Shao SN", "%.5f (z %+.1f)" % (p, (p - ALPHA) / se))
            else:
                fixed = (6, 10, 16)
                hit_mv = 0
                hit_fix = {B: 0 for B in fixed}
                done = 0
                while done < R:
                    r = min(20000, R - done)
                    X = ar1(r, n0, phi, 1.0, rng) + 2.0       # benchmark, dates 1..36
                    Y = ar1(r, n1, phi, 1.0, rng) + 2.0       # candidate, dates 24..36
                    CX, CY = W.cums(X, Y)
                    T, dxy = W.stat(CX, CY, 1, W.N, 0.0)       # dxy = Xbar - Ybar
                    up = -dxy > 0                              # candidate above benchmark
                    qB = {B: W.subsample_q(CX, CY, B, dxy, 1 - P_UPPER) for B in grid}
                    Q = np.column_stack([qB[B] for B in grid])
                    # minimum volatility: sd of q over a +-2 neighbourhood, interior B only
                    vol = np.column_stack([Q[:, j - 2:j + 3].std(1) for j in range(2, len(grid) - 2)])
                    pick = vol.argmin(1) + 2
                    qmv = Q[np.arange(r), pick]
                    hit_mv += int((up & (T > qmv)).sum())
                    for B in fixed:
                        hit_fix[B] += int((up & (T > qB[B])).sum())
                    done += r
                p = hit_mv / R
                cells.append({"design": dtag, "phi": phi, "method": "warped self-normalised (Wang-Zhang)",
                              "size": p, "mcse": se, "z": (p - ALPHA) / se,
                              "fixed_B": {str(B): hit_fix[B] / R for B in fixed}})
                claim(rows, "L1T-%s-%.1f" % (dtag, phi), "size, WSNS min-volatility B (fixed B 6/10/16)",
                      "%.5f (z %+.1f) (%s)" % (p, (p - ALPHA) / se,
                                              " / ".join("%.5f" % (hit_fix[B] / R) for B in fixed)))
    worst = max(cells, key=lambda c: c["size"])
    nbreach = sum(1 for c in cells if c["z"] > 3)
    claim(rows, "L1T-WORST", "worst size, times target, cells above target",
          "%.4f, %.2f, %d of %d" % (worst["size"], worst["size"] / ALPHA, nbreach, len(cells)))
    save("table", {"q": q, "draws": R, "cells": cells, "b_grid": grid,
                   "worst_size": worst["size"], "worst_ratio": worst["size"] / ALPHA,
                   "cells_above": nbreach}, rows)


if __name__ == "__main__":
    t0 = time.time()
    {"uq": stage_uq, "large": stage_large, "anchors": stage_anchors, "table": stage_table}[sys.argv[1]]()
    print("seconds %.0f" % (time.time() - t0))
