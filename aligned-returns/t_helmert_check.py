"""Independent numerical check of the Helmert identity claims for H-kappa.

Three claims, each checked against simulation:
  1. Nstd = sqrt(n)*(ybar1 - ybar0)/sqrt(aF)  has mean 0 and variance 1;
  2. Nstd is INDEPENDENT of Sres = nu*(s0^2 + s1^2)  (the reduction's load-bearing step);
  3. Sres has the law mu1*X1 + mu2*X2, X_i ~ chi2_nu independent, where mu1, mu2 are the
     eigenvalues of G = M M^T and r_j = mu_j / aF.
"""
import numpy as np

rng = np.random.default_rng(7)
n, nu = 36, 35
lam, b0, b1 = 1 / 25, -3 / 2, 3 / 2
aF = (b1 - b0) ** 2 + 1 + lam
M = np.array([[b0, 1, 0], [b1, 0, np.sqrt(lam)]], float)
G = M @ M.T
ev = np.linalg.eigvalsh(G)
mu2, mu1 = ev[0], ev[1]
print("aF =", aF)
print("eigenvalues of G =", mu2, mu1)
print("r1 = mu1/aF =", mu1 / aF)
print("r2 = mu2/aF =", mu2 / aF)

N = 400_000
f = rng.standard_normal((N, n))
e0 = rng.standard_normal((N, n))
e1 = np.sqrt(lam) * rng.standard_normal((N, n))
y0 = b0 * f + e0
y1 = b1 * f + e1
s0 = y0.var(axis=1, ddof=1)
s1 = y1.var(axis=1, ddof=1)
Nstd = np.sqrt(n) * (y1.mean(axis=1) - y0.mean(axis=1)) / np.sqrt(aF)
Sres = nu * (s0 + s1)

print()
print("claim 1: Nstd mean =", round(float(Nstd.mean()), 5), " var =", round(float(Nstd.var()), 5),
      " (target 0, 1)")
corr = float(np.corrcoef(Nstd, Sres)[0, 1])
print("claim 2: corr(Nstd, Sres) =", round(corr, 5), " (target 0, independent)")

X1 = rng.chisquare(nu, N)
X2 = rng.chisquare(nu, N)
ref = mu1 * X1 + mu2 * X2
print()
print("claim 3: quantile of Sres vs mu1*X1 + mu2*X2")
for q in (0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99):
    a = float(np.quantile(Sres, q))
    b = float(np.quantile(ref, q))
    print(f"   q={q:5.2f}  Sres={a:12.4f}  ref={b:12.4f}  rel diff={abs(a-b)/max(abs(b),1e-12):.2e}")
