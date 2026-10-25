# -*- coding: utf-8 -*-
"""pb_certify.py — option B, Phase 1: certify the calibrated table rule at K = 2 over rho in D.

Claim: for the deployed table at internal level a, FWER(rho) <= 1/40 for every rho in
D = [-0.99, 0.99]. FWER = 2 P1 - P12 = P1 + Q; P1 is nondecreasing in rho (pb_cert), so on an interval
[ra, rb] it is bounded by its certified upper bound at rb; the joint term is enclosed over the whole
interval (pb_joint).

    python pb_certify.py N LEVEL            e.g.  python pb_certify.py 12 0.0247
    python pb_certify.py falsifier          (C2: table at 1.03 alpha must fail near rho = 0.19)
"""
import json
import math
import os
import sys
import time

import numpy as np
from flint import fmpq

import p_calib_k2 as PC2
import pb_cert as PC
import pb_joint as PJ

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = math.nextafter(0.025, 0.0)            # the largest double below 1/40
DZ = 0.003
RMAX = 0.99


def dec6(x, up=False):
    """Exact 6-decimal rational, rounded to nearest (edges) or up (constants)."""
    v = math.ceil(x * 1e6 - 1e-9) if up else round(x * 1e6)
    return fmpq(int(v), 10 ** 6)


def deployed(n, a):
    r = PC2.edges()
    r6 = [dec6(x) for x in r]
    rf = np.array([int(e.p) / int(e.q) for e in r6])
    cs = [PC2.cstar(n, x, a) for x in rf]        # c*(left edge) at the rounded edges
    from scipy.stats import t as tdist
    cs = [float(tdist.isf(a / 2, n - 1))] + cs
    c6 = [dec6(c, up=True) for c in cs]
    for j in range(1, len(c6)):                  # nonincreasing
        if c6[j] > c6[j - 1]:
            c6[j] = c6[j - 1]
    return r6, c6, rf, np.array([int(c.p) / int(c.q) for c in c6])


def intervals():
    zs = np.arange(-math.atanh(RMAX), math.atanh(RMAX) + DZ / 2, DZ)
    rs = np.tanh(zs)
    rs[0], rs[-1] = -RMAX, RMAX
    rs = np.concatenate([rs[rs < RMAX], [RMAX]])
    return list(zip(rs[:-1], rs[1:]))


class Certifier:
    def __init__(self, n, r6, c6, log, rounds=12, max_cells=3_000_000, depth=3):
        self.n, self.log = n, log
        self.rounds, self.max_cells, self.depth = rounds, max_cells, depth
        self.phi = PC.PhiTable()
        self.T = PC.Table(r6, c6)
        self.P1 = PC.P1Bound(n, self.T, self.phi)
        self.J = PJ.JointAdaptive(n, self.T, self.phi)

    def one(self, ra, rb, depth=0):
        p1u = self.P1.upper(rb)
        if 2 * p1u <= ALPHA:
            return [dict(ra=ra, rb=rb, route="P1", U=2 * p1u, p1u=p1u, depth=depth)]
        modes = ("p12", "q") if rb <= 0.8 else ("q", "p12")
        if ra < 0:
            modes = ("p12",)
        best = None
        for mode in modes:
            need = 2 * p1u - ALPHA if mode == "p12" else ALPHA - p1u
            val, nc, rnd = self.J.run(ra, rb, mode, need, max_rounds=self.rounds, max_cells=self.max_cells)
            U = 2 * p1u - val if mode == "p12" else p1u + val
            rec = dict(ra=ra, rb=rb, route=mode, U=U, p1u=p1u, joint=val, cells=nc, rounds=rnd, depth=depth)
            if best is None or U < best["U"]:
                best = rec
            if U <= ALPHA:
                return [rec]
        if depth < self.depth:
            zm = 0.5 * (math.atanh(ra) + math.atanh(rb))
            rm = math.tanh(zm)
            self.log(f"  bisect [{ra:.6f}, {rb:.6f}] (U {best['U']:.7f}) at depth {depth}")
            return self.one(ra, rm, depth + 1) + self.one(rm, rb, depth + 1)
        best["route"] = "FAIL"
        return [best]


def c1_check(n, rf, cf, rec):
    for rho in (rec["ra"], rec["rb"]):
        f = PC2.Eval(n, float(rho), rf).fwer(cf)
        if f > rec["U"] + 1e-15:
            return f"DEFECT: float FWER {f:.8f} > certified {rec['U']:.8f} at rho={rho}"
    return None


def certify(n, a, out_tag=None, only=None, rounds=12, max_cells=3_000_000, depth=3):
    t0 = time.time()
    lines = []
    logf = open(os.path.join(HERE, "deposits", f"pb_certify_n{n}{out_tag or ''}.log"), "w", encoding="utf-8")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    def claim(cid, label, value):
        lines.append(f"CLAIM {cid}: {label} = {value}")
        log(lines[-1])

    r6, c6, rf, cf = deployed(n, a)
    tag = f"PB1-{n}{out_tag or ''}"
    claim(f"{tag}-RULE", f"deployed table, n={n}, internal level a={a}: 81 edges r_j (6 dp) | 82 constants (6 dp, rounded up, nonincreasing)",
          " ".join(str(x) for x in r6) + " | " + " ".join(str(x) for x in c6))
    C = Certifier(n, r6, c6, log, rounds, max_cells, depth)
    log(f"budget: rounds {rounds}, max cells {max_cells}, bisection depth {depth}")
    log(f"setup {time.time() - t0:.0f}s")
    ivs = intervals() if only is None else only
    recs = []
    for i, (ra, rb) in enumerate(ivs):
        rr = C.one(ra, rb)
        for rec in rr:
            if (i % 10 == 0 and rec is rr[0]) or rec["depth"] > 0:
                msg = c1_check(n, rf, cf, rec)
                rec["c1"] = "ok" if msg is None else msg
                if msg:
                    log(msg)
                    raise SystemExit(msg)
        recs += rr
        worst = max(rr, key=lambda d: d["U"])
        if i % 25 == 0 or worst["route"] == "FAIL" or worst["U"] > 0.02495:
            log(f"[{i}/{len(ivs)}] [{ra:+.6f}, {rb:+.6f}] {worst['route']} U {worst['U']:.7f} ({time.time() - t0:.0f}s)")
    fails = [d for d in recs if d["route"] == "FAIL"]
    W = max(recs, key=lambda d: d["U"])
    routes = {k: sum(1 for d in recs if d["route"] == k) for k in ("P1", "p12", "q", "FAIL")}
    claim(f"{tag}-INTERVALS", "intervals (after bisection) by route P1-only / P12 / Q / FAIL",
          f"{len(recs)}: {routes['P1']} / {routes['p12']} / {routes['q']} / {routes['FAIL']}; bisected {sum(1 for d in recs if d['depth'] > 0)}")
    claim(f"{tag}-MAX", "largest certified FWER upper bound over D = [-0.99, 0.99] and where",
          f"{W['U']:.8f} on [{W['ra']:+.6f}, {W['rb']:+.6f}] via {W['route']}")
    claim(f"{tag}-C1", "consistency checks (float FWER <= certified bound) run / defects",
          f"{sum(1 for d in recs if 'c1' in d)} / 0")
    claim(f"{tag}-VERDICT", f"FWER <= 1/40 on all of D for the table at a={a}, n={n}",
          "CERTIFIED" if not fails else f"NOT CERTIFIED ({len(fails)} intervals fail)")
    claim(f"{tag}-TIME", "run time (s)", f"{time.time() - t0:.0f}")
    with open(os.path.join(HERE, "deposits", f"pb_certify_n{n}{out_tag or ''}_claims.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(HERE, "deposits", f"pb_certify_n{n}{out_tag or ''}_intervals.json"), "w", encoding="utf-8") as fh:
        json.dump(recs, fh)
    return not fails


if __name__ == "__main__":
    if sys.argv[1] == "falsifier":
        ivs = [iv for iv in intervals() if iv[0] >= 0.17][:5]
        ok = certify(12, 0.025 * 1.03, out_tag="_falsifier", only=ivs)
        print("C2:", "FAIL OF THE CHECK (an invalid rule certified)" if ok else "PASS (the invalid rule was not certified)")
    else:
        a = sys.argv
        certify(int(a[1]), float(a[2]), out_tag=a[3] if len(a) > 3 else None,
                rounds=int(a[4]) if len(a) > 4 else 12, max_cells=int(a[5]) if len(a) > 5 else 3_000_000,
                depth=int(a[6]) if len(a) > 6 else 3)
