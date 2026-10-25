# -*- coding: utf-8 -*-
"""pb_checks.py — option B, Phase 1 checks on a certified rule: C3 (brute force) and power.

Reads the deployed table from deposits/pb_certify_n{N}_claims.txt (the RULE claim), so the checks run
on exactly the rule that was certified.

    python pb_checks.py N
"""
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.stats import t as tdist

import p_calib_k2 as PC2

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.025


def load_rule(n, tag=""):
    path = os.path.join(HERE, "deposits", f"pb_certify_n{n}{tag}_claims.txt")
    for line in open(path, encoding="utf-8"):
        if f"PB1-{n}{tag}-RULE" in line:
            body = line.split(" = ", 1)[1]
            er, ec = body.split(" | ")
            frac = lambda s: int(s.split("/")[0]) / int(s.split("/")[1]) if "/" in s else float(s)
            return np.array([frac(x) for x in er.split()]), np.array([frac(x) for x in ec.split()])
    raise SystemExit("rule not found")


def main(n, tag=""):
    r, c = load_rule(n, tag)
    lines = []

    def claim(cid, label, value):
        lines.append(f"CLAIM {cid}: {label} = {value}")
        print(lines[-1], flush=True)

    for k, rho in enumerate((0.19, 0.6, 0.95)):
        f = PC2.Eval(n, rho, r).fwer(c)
        p, se = PC2.brute_table(n, rho, r, c, 1_000_000, 20261026 + 10 * n + k)
        z = (p - f) / se
        claim(f"PB1C3-{n}-{rho:.2f}", f"certified rule, n={n}, rho={rho}: brute 1e6 vs evaluator",
              f"brute {p:.5f}, eval {f:.6f}, z {z:+.2f} ({'PASS' if abs(z) < 3 else 'FAIL'})")
    cB = float(tdist.isf(ALPHA / 2, n - 1))
    for rho in (0.0, 0.45, 0.75, 0.9, 0.95, 0.99):
        E = PC2.Eval(n, rho, r)
        co = PC2.cstar(n, rho, ALPHA)
        nb = len(c)
        powr = lambda arr, d: 1.0 - E.nonreject(arr, d)
        d80 = brentq(lambda d: powr(np.full(nb, co), d) - 0.8, 0.0, 10.0, xtol=1e-7)
        pc, pb = powr(c, d80), powr(np.full(nb, cB), d80)
        claim(f"PB1POW-{n}-{rho:.2f}", f"power at oracle-80% edge, n={n}, rho={rho}: certified / Bonferroni / oracle; gain over Bonferroni (pp); loss to oracle (pp)",
              f"{pc:.4f} / {pb:.4f} / 0.8000; {100 * (pc - pb):+.2f}; {100 * (0.8 - pc):+.2f}")
    with open(os.path.join(HERE, "deposits", f"pb_checks_n{n}{tag}_claims.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else "")
