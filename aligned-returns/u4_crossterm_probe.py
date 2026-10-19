# -*- coding: utf-8 -*-
"""u4_crossterm_probe.py — check of the cross term behind Section 7.3.

Question: with m = 2 latent factors, does a candidate's variance estimate depend
on the factor paths only through each factor's own mean and sum of squares (the
unit sheet's "2m + 2" outer rule), or also through the cross-product
sum_t (f1_t - f1bar)(f2_t - f2bar)?

Theory (to be falsified, not trusted): given the factor paths, nu * s_i^2 / sigma_i^2
is noncentral chi-square with nu = n - 1 degrees of freedom and noncentrality
b_i' S_FF b_i / sigma_i^2, where S_FF is the full 2 x 2 centred sum-of-squares
matrix. So two factor paths with the SAME diagonal of S_FF but DIFFERENT
off-diagonal must give DIFFERENT laws for s_i^2 whenever b_i1 * b_i2 != 0.

Test, declared before running: fix two factor paths A and B with identical means
and identical diagonal of S_FF but off-diagonals of opposite sign. Draw the
idiosyncratic noise D times under each and compare the mean of s_i^2.
  * Sanity anchor (P-4): with b_i = (1, 0) the two laws must agree (|z| < 4).
  * Test: with b_i = (1, 1) the means must differ by the predicted amount
    2 * b1 * b2 * (S12_A - S12_B) / nu, and |z| against that prediction < 4.
Tolerance: z-gate 4 (two-sided normal tail 6e-5 per check, two checks).
"""
import numpy as np

n, D, sigma2 = 36, 400_000, 1.0
nu = n - 1
rng = np.random.default_rng(20261006)

# build two factor paths with equal means and equal diagonal S_FF, opposite S12
g = rng.standard_normal((2, n))
g -= g.mean(1, keepdims=True)
# orthonormalise the centred rows, then set the target Gram matrices exactly
q, _ = np.linalg.qr(np.column_stack([np.ones(n) / np.sqrt(n), g.T]))
u1, u2 = q[:, 1], q[:, 2]              # centred, orthonormal
s11, s22, rho = 30.0, 40.0, 0.6
def path(r):
    f1 = np.sqrt(s11) * u1
    f2 = np.sqrt(s22) * (r * u1 + np.sqrt(1 - r * r) * u2)
    return np.vstack([f1, f2]) + 0.3   # common mean 0.3 for both factors
FA, FB = path(rho), path(-rho)
SA = (FA - FA.mean(1, keepdims=True)) @ (FA - FA.mean(1, keepdims=True)).T
SB = (FB - FB.mean(1, keepdims=True)) @ (FB - FB.mean(1, keepdims=True)).T
print('CLAIM U4P-DIAG: S_FF diagonals A, B = %s, %s' % (np.round(np.diag(SA), 10), np.round(np.diag(SB), 10)))
print('CLAIM U4P-OFFDIAG: S12 A, B = %.6f, %.6f' % (SA[0, 1], SB[0, 1]))

def s2_draws(F, b):
    eps = rng.standard_normal((D, n)) * np.sqrt(sigma2)
    r = (b @ F)[None, :] + eps
    return r.var(1, ddof=1)

for label, b in (('anchor b=(1,0)', np.array([1.0, 0.0])), ('test b=(1,1)', np.array([1.0, 1.0]))):
    a, bb = s2_draws(FA, b), s2_draws(FB, b)
    diff = a.mean() - bb.mean()
    se = np.sqrt(a.var() / D + bb.var() / D)
    pred = (b @ SA @ b - b @ SB @ b) / nu
    print('CLAIM U4P-%s: mean s2 under A minus under B = %.5f, predicted %.5f, se %.5f, z vs prediction %+.2f'
          % (label.split()[0].upper(), diff, pred, se, (diff - pred) / se))
