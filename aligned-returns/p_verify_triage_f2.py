# -*- coding: utf-8 -*-
"""Triage of falsifier F-2 at rho = +/-0.999 (Section 8.3).

Compares the K = 2 evaluator p_law.Law2 with 8,000,000-draw brute force and with the
rho -> +/-1 limit, at n = 12 and c = 1, 2.  Writes deposits/p_verify_triage_F2_claims.txt.
"""
import pathlib

from scipy.stats import t as tdist

import p_law as P

rows = []
for n in (12,):
    for rho in (-0.999, 0.999):
        L = P.Law2(n, rho, panels=32)
        for c in (1.0, 2.0):
            f = L.F(c)
            lim = (2*tdist.cdf(c, n-1)-1) if rho < 0 else tdist.cdf(c, n-1)
            pb, se = P.brute(n, rho, c, 8_000_000, 110000 + int(10*c) + (rho>0))
            rows.append("CLAIM TRIAGE-F2-n%d-r%+.3f-c%.0f: evaluator / brute (8e6) / limit = %.6f / %.6f (z %+.2f) / %.6f; evaluator-limit %.2e"
                        % (n, rho, c, f, pb, (pb-f)/se, lim, f-lim))
            print(rows[-1])
out = pathlib.Path(__file__).resolve().parent / "deposits" / "p_verify_triage_F2_claims.txt"
out.write_text("\n".join(rows) + "\n", encoding="utf-8")
