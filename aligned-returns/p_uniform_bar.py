# -*- coding: utf-8 -*-
"""p_uniform_bar.py — Unit P: the fixed critical value that is valid at every correlation (K = 2).

Section 8.8 proves that a fixed c >= 0 keeps the two-candidate family-wise error rate at most alpha
for every correlation in [-1, 1] exactly when P(T_nu > c) <= alpha / 2, i.e. at Bonferroni's
critical value t_{nu, 1 - alpha/2}. This script supplies the numbers that illustrate it:

  * Bonferroni's critical value at n = 12 and n = 36;
  * the exact own-denominator critical value c*(rho) of p_law.Law2 as rho approaches -1, which
    must stay at or below Bonferroni's and approach it;
  * an anchor: c*(-0.5) reproduces Table 6;
  * a brute-force check of F at c*(-0.99), n = 12.

    python p_uniform_bar.py
"""
import os

from scipy.stats import t as tdist

from p_law import Law2, brute

ALPHA = 0.025
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "p_uniform_bar_claims.txt")
TABLE6 = {12: 2.593012, 36: 2.341953}          # Table 6, rho = -0.5


def main():
    lines = []

    def claim(cid, label, value):
        lines.append(f"CLAIM {cid}: {label} = {value}")
        print(lines[-1])

    for n in (12, 36):
        nu = n - 1
        cb = float(tdist.ppf(1 - ALPHA / 2, nu))
        claim(f"UB-{n}-BONF", f"Bonferroni critical value t_(nu={nu}, 1-alpha/2), alpha = {ALPHA}", f"{cb:.6f}")
        c05 = Law2(n, -0.5).cstar(ALPHA)
        claim(f"UB-{n}-ANCHOR", f"c*(rho=-0.5), n={n}, reproduces Table 6 ({TABLE6[n]})",
              f"{c05:.6f} (|diff| {abs(c05 - TABLE6[n]):.1e})")
        for rho in (-0.9, -0.99, -0.999):
            c = Law2(n, rho).cstar(ALPHA)
            claim(f"UB-{n}-{abs(rho)}", f"c*(rho={rho}), n={n}; gap to Bonferroni",
                  f"{c:.6f} (gap {cb - c:.2e}, below Bonferroni: {c <= cb + 1e-9})")
    n, rho = 12, -0.99
    c = Law2(n, rho).cstar(ALPHA)
    p, se = brute(n, rho, c, 4_000_000, seed=20261024)
    claim("UB-BRUTE", f"brute-force P(max t <= c*(-0.99)), n=12, 4,000,000 draws; target {1 - ALPHA}",
          f"{p:.5f} (se {se:.5f}, z {(p - (1 - ALPHA)) / se:+.2f})")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
