"""Well-conditioned quadrature primitives for the heteroscedastic max-t reduction.

The natural rules for these integrals are Gauss-Hermite for the normal and
generalized Gauss-Laguerre for the chi-square. Both were tried first and both
fail here, in two different ways:

  * generalized Gauss-Laguerre OVERFLOWS for large degrees of freedom, because
    its weights sum to Gamma(nu/2), which passes the double-precision ceiling at
    nu around 340;
  * at high node counts BOTH rules lose conditioning, so refining the grid makes
    the answer worse rather than better. Measured on a nu = 5 control, the
    critical value moved 2.603140 -> 2.603162 -> 2.603445 -> 2.603209 as the
    chi-square rule went 24 -> 48 -> 96 -> 192 nodes, and 384 nodes returned NaN.
    A rule that cannot be refined cannot certify its own answer.

So both are replaced by the same construction: map the variable to a probability
with its own quantile function, map the probability to a z-score, and integrate
against the standard normal density in z with COMPOSITE Gauss-Legendre over
[-zmax, zmax]. Every panel is a low-order rule, so conditioning never degrades;
accuracy is refined by adding panels, not by raising the order. The integrands
here are bounded by 1, so truncating at zmax costs at most 2*(1 - Phi(zmax)),
which is 1.2e-15 at zmax = 8.

Both branches are cross-checked against the classical rules in the range where
those are trustworthy.
"""

import numpy as np
from scipy.special import ndtr
from scipy.stats import chi2, norm


def zpanel_nodes(n_panel=16, order=12, zmax=8.0):
    """Composite Gauss-Legendre nodes and normal-density weights on [-zmax, zmax].

    Returns (z, w) with sum(w) = Phi(zmax) - Phi(-zmax).
    """
    if n_panel < 1 or order < 2:
        raise RuntimeError("need n_panel >= 1 and order >= 2")
    x, wl = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(-zmax, zmax, n_panel + 1)
    half = 0.5 * (edges[1:] - edges[:-1])
    mid = 0.5 * (edges[1:] + edges[:-1])
    z = (mid[:, None] + half[:, None] * x[None, :]).ravel()
    w = (half[:, None] * wl[None, :]).ravel() * norm.pdf(z)
    return z, w


def normal_nodes(n_panel=16, order=12, zmax=8.0):
    """Nodes and probability weights for a standard normal."""
    return zpanel_nodes(n_panel, order, zmax)


def chi2_over_df_nodes(nu, n_panel=16, order=12, zmax=8.0):
    """Nodes and probability weights for S^2/sigma^2 = chi2_nu / nu.

    Stable at every nu: the quantile function carries the shape, and the
    quadrature runs on the z-scale where the rule is fixed.
    """
    if nu <= 0:
        raise RuntimeError("nu must be positive, got %r" % (nu,))
    z, w = zpanel_nodes(n_panel, order, zmax)
    r = chi2.ppf(ndtr(z), nu) / nu
    if not np.all(np.isfinite(r)):
        raise RuntimeError("non-finite chi2 quantile at nu=%r" % (nu,))
    return r, w
