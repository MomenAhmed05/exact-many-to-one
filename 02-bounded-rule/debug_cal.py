"""Localise the calibrate mismatch: diff G/H/cgrid/s2 between implementations."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import cond_noreject  # noqa: E402
from h12_calibrated import BOOT, NGRID, SPAN, plugin_sigma  # noqa: E402
from h12_hetero_exact import HeteroMaxT  # noqa: E402
from h12_quad import normal_nodes  # noqa: E402
from h12_tamhane_headtohead import cond_reject  # noqa: E402

n = np.array([4.0, 13.0, 13.0])
sig = np.array([10.0, 10.0, 50.0])
sw = sig / sig[0]
nu = n - 1.0
alpha = 0.025
un = normal_nodes(16, 10)
u_nodes = (un[0] * sw[0] / np.sqrt(n[0]), un[1])

rng_p = np.random.default_rng(31337)
c_ref_p = HeteroMaxT(n, sw, grid=BOOT).critical_value(alpha)
cgrid_p = c_ref_p + np.linspace(SPAN[0], SPAN[1], NGRID)
s2_p = (rng_p.chisquare(nu, size=n.size) / nu) * sw**2
m_p = HeteroMaxT(n, plugin_sigma(s2_p), grid=BOOT)
G_p = np.array([1.0 - m_p.cdf(float(c)) for c in cgrid_p])
H_p = np.array([cond_reject(n, sw, s2_p, float(c), u_nodes) for c in cgrid_p])

rng_m = np.random.default_rng(31337)
c_ref_m = HeteroMaxT(n, sw, grid=BOOT).critical_value(alpha)
cgrid_m = c_ref_m + np.linspace(SPAN[0], SPAN[1], NGRID)
s2_m = (rng_m.chisquare(nu, size=n.size) / nu) * sw**2
m_m = HeteroMaxT(n, plugin_sigma(s2_m), grid=BOOT)
G_m = np.array([1.0 - m_m.cdf(float(c)) for c in cgrid_m])
H_m = 1.0 - cond_noreject(
  np.asarray(cgrid_m, float), np.array([s2_m[0]]), s2_m[None, 1:], n, sw**2
)

print("c_ref  ", c_ref_p, c_ref_m)
print("cgrid equal:", np.allclose(cgrid_p, cgrid_m))
print("s2_p", s2_p, "\ns2_m", s2_m)
print("G_p", G_p)
print("G_m", G_m)
print("H_p", H_p)
print("H_m", H_m)
print("max|dG|", np.abs(G_p - G_m).max(), "max|dH|", np.abs(H_p - H_m).max())
