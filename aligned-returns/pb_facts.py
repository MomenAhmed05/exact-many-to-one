# -*- coding: utf-8 -*-
"""pb_facts.py — option B: evaluator values quoted in the text of Section 8.9 (Tested evaluator).

    python pb_facts.py
"""
import os

import numpy as np
from scipy.special import ndtr

import p_calib_k2 as PC2
import p_law as P
import pb_certify as PB
import pb_checks as CK

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    lines = []

    def claim(cid, label, value):
        lines.append(f"CLAIM {cid}: {label} = {value}")
        print(lines[-1], flush=True)

    # the falsifier's table (1.03 alpha, n = 12): its true FWER where it was run
    _, _, rf, cf = PB.deployed(12, 0.025 * 1.03)
    for rho in (0.172, 0.18, 0.186):
        claim(f"PBF-FALS-{rho}", f"float FWER of the 1.03-alpha table, n=12, rho={rho}", f"{PC2.Eval(12, rho, rf).fwer(cf):.6f}")
    # the n = 36 tables at the levels that failed, and the certified one, at high rho
    for tag, a in (("", 0.02488), ("_r2a02468", 0.02468)):
        r, c = CK.load_rule(36, tag) if tag else PB.deployed(36, a)[2:]
        for rho in (0.96, 0.98, 0.99):
            f = PC2.Eval(36, rho, r).fwer(c)
            claim(f"PBF-36-{a}-{rho}", f"float FWER, n=36 table at a={a}, rho={rho}; margin below 1/40", f"{f:.6f}; {0.025 - f:.2e}")
    # P12 of the certified n = 12 rule near where its calibration binds
    r, c = CK.load_rule(12, "")
    for rho in (0.0, 0.2, 0.5, 0.9):
        E = PC2.Eval(12, rho, r)
        h = np.broadcast_to(c[None, None, :, None] * E.hx[..., None, None], E.hr.shape)
        k = c[None, None, :, None] * E.hr
        W = E.wm * E.W[..., None, None]
        p1 = float((ndtr(-h) * W).sum())
        p12 = float((P.phi2(-h.ravel(), -k.ravel(), rho).reshape(k.shape) * W).sum())
        claim(f"PBF-12-P-{rho}", f"certified n=12 rule at rho={rho}: FWER, P1, P12", f"{2 * p1 - p12:.6f}, {p1:.6f}, {p12:.6f}")
    with open(os.path.join(HERE, "deposits", "pb_facts_claims.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
