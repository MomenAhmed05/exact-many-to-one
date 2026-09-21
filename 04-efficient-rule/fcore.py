"""Shared machinery for the efficiency-frontier programme (04-efficient-rule).

The question this folder answers is NOT validity.  A finite-sample valid rule
already exists (`03-certified-supremum`).  The question is efficiency: how much power does hard
uniform FWER control cost when the variance ratios are unknown, and how much of
that cost is removable.

Three objects must never be conflated, and every function here is named for the
one it computes:

  * `cond_reject`      P(reject | S), the CONDITIONAL rejection probability at a
                       given threshold and a given true ratio vector.  Exact.
  * `fwer`             E_S[cond_reject], the UNCONDITIONAL familywise error rate
                       of a rule whose threshold is a function of the realised
                       variance estimates.  Exact by quadrature at K = 1, and
                       Rao-Blackwellised (exact inner, sampled outer) above it.
  * `power`            the same expectation at a declared alternative.

THE K=1 CLOSED FORM.  With one candidate the inner expectation collapses:
E_z[Phi(a + b z)] = Phi(a / sqrt(1 + b^2)), so

    P(no reject | S) = Phi( c * d / sigma_diff ),
    d^2 = S_1^2/n_1 + S_0^2/n_0 ,  sigma_diff^2 = sigma_1^2/n_1 + sigma_0^2/n_0 .

Everything at K = 1 is therefore deterministic to quadrature precision, with no
draws anywhere, which is what makes the lower-bound work in `minimax.py`
possible.  The identity is checked against the general engine in
`tests_frontier.py` (falsifier X1).

SCALE INVARIANCE.  sigma_0 = 1 throughout; lambda_i = sigma_i^2 / sigma_0^2.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from math import sqrt

import numpy as np
from scipy.stats import norm

_HERE = os.path.dirname(os.path.abspath(__file__))
_AUDIT = os.path.dirname(_HERE)
for _p in (_AUDIT, os.path.join(_AUDIT, "01-exact-law"),
           os.path.join(_AUDIT, "02-bounded-rule")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from h12_hetero_exact import HeteroMaxT  # noqa: E402
from h12_quad import chi2_over_df_nodes, normal_nodes  # noqa: E402

DEP = os.path.join(_HERE, "deposits")
os.makedirs(DEP, exist_ok=True)

ABAR = 0.025          # the hard unconditional ceiling
GRID = {"u": (14, 14), "v": (14, 14), "s": (4, 14)}
FAST = {"u": (8, 10), "v": (8, 10), "s": (3, 10)}

# tiers, in the 03-certified-supremum sense; every deposit carries one
TIERS = (
    "oracle_exact_law",
    "plugin_measured",
    "envelope_measured",
    "region_sup_certified",
    "unconditional_candidate_measured",
    "unconditional_valid_certified",
    "frontier_measured",
    "frontier_certified",
    # the beyond-Welch programme (W1-W4)
    "exact_in_model",
    "certified_bound",
    "measured",
    "asymptotic",
    "diagnostic",
)


# ----------------------------------------------------------------- deposits
def _san(o):
    if isinstance(o, dict):
        return {k: _san(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_san(v) for v in o]
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, (np.integer, int)):
        return int(o)
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _san(o.tolist())
    return o


def write_deposit(name, obj, tier):
    """Atomic, strict-JSON deposit carrying a tier label."""
    if tier not in TIERS:
        raise ValueError("unknown tier %r" % (tier,))
    obj = dict(obj)
    obj["tier"] = tier
    path = os.path.join(DEP, name)
    fd, tmp = tempfile.mkstemp(dir=DEP, suffix=".tmp")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(_san(obj), fh, indent=1, allow_nan=False)
    os.replace(tmp, path)
    return path


def read_deposit(name):
    with open(os.path.join(DEP, name), encoding="utf-8") as fh:
        return json.load(fh)


# ------------------------------------------------------------ the state law
def state_nodes(nu, npanel=24, order=12):
    """Nodes and probability weights for S^2/sigma^2 = chi2_nu / nu."""
    return chi2_over_df_nodes(float(nu), npanel, order)


def state_grid(n, npanel=24, order=12):
    """Tensor grid over the K+1 scaled chi-squares, as (points, weights).

    `points[:, i]` is u_i = S_i^2 / sigma_i^2.  Only feasible for small K; above
    K = 2 use `draw_states`.
    """
    n = np.asarray(n, float)
    axes = [state_nodes(ni - 1.0, npanel, order) for ni in n]
    qs = [a[0] for a in axes]
    ws = [a[1] for a in axes]
    mesh = np.meshgrid(*qs, indexing="ij")
    wmesh = np.meshgrid(*ws, indexing="ij")
    pts = np.stack([m.ravel() for m in mesh], axis=1)
    wt = np.ones(pts.shape[0])
    for m in wmesh:
        wt = wt * m.ravel()
    return pts, wt


def draw_states(n, D, seed):
    """D draws of u_i = S_i^2/sigma_i^2 ~ chi2_{n_i-1}/(n_i-1)."""
    n = np.asarray(n, float)
    rng = np.random.default_rng(int(seed))
    nu = n - 1.0
    return rng.chisquare(nu[None, :], size=(int(D), n.size)) / nu[None, :]


# ------------------------------------------------- conditional rejection law
def _u_nodes(npanels=64, order=12):
    return normal_nodes(npanels, order)


def cond_noreject(c, u, n, lam, quad=None):
    """P(max_i T_i <= c | S) with S_i^2 = lam_i * u_i (sigma_0 = 1, lam_0 = 1).

    `c` is a scalar or a length-K vector; `u` is a length-(K+1) array of the
    scaled chi-squares.  Exact: the only approximation is the z-quadrature.
    """
    n = np.asarray(n, float)
    u = np.asarray(u, float)
    lam = np.concatenate([[1.0], np.asarray(lam, float)])
    s2 = lam * u                      # the realised sample variances
    K = n.size - 1
    d = np.sqrt(s2[1:] / n[1:] + s2[0] / n[0])
    cv = np.broadcast_to(np.asarray(c, float), (K,))
    if K == 1:
        sig = sqrt(lam[1] / n[1] + 1.0 / n[0])
        return float(norm.cdf(cv[0] * d[0] / sig))
    z, w = quad if quad is not None else _u_nodes()
    sfac = np.sqrt(n[1:] / lam[1:])
    arg = (z[:, None] / sqrt(n[0]) + cv[None, :] * d[None, :]) * sfac[None, :]
    return float(w @ np.prod(norm.cdf(arg), axis=1))


def cond_reject(c, u, n, lam, quad=None):
    return 1.0 - cond_noreject(c, u, n, lam, quad)


def cond_power(c, u, n, lam, delta, quad=None):
    """P(at least one rejection | S) when candidate means are shifted by delta.

    Under the shift every candidate's numerator gains `delta`, so the event
    becomes Z_i - U <= c d_i - delta and the same factorisation applies.
    """
    n = np.asarray(n, float)
    u = np.asarray(u, float)
    lam = np.concatenate([[1.0], np.asarray(lam, float)])
    s2 = lam * u
    K = n.size - 1
    d = np.sqrt(s2[1:] / n[1:] + s2[0] / n[0])
    cv = np.broadcast_to(np.asarray(c, float), (K,))
    dl = np.broadcast_to(np.asarray(delta, float), (K,))
    if K == 1:
        sig = sqrt(lam[1] / n[1] + 1.0 / n[0])
        return float(norm.cdf((dl[0] - cv[0] * d[0]) / sig))
    z, w = quad if quad is not None else _u_nodes()
    sfac = np.sqrt(n[1:] / lam[1:])
    arg = (z[:, None] / sqrt(n[0]) + cv[None, :] * d[None, :] - dl[None, :]) * sfac[None, :]
    return 1.0 - float(w @ np.prod(norm.cdf(arg), axis=1))


# ------------------------------------------------------- oracle references
def oracle_c_pop(n, lam, alpha=ABAR, grid=None):
    """The population fixed-ratio critical value (the paper's exact law)."""
    n = np.asarray(n, float)
    sig = np.sqrt(np.concatenate([[1.0], np.asarray(lam, float)]))
    return float(HeteroMaxT(n, sig, grid=grid or GRID).critical_value(alpha))


def oracle_c_cond(u, n, lam, alpha=ABAR, quad=None, lo=1e-4, hi=60.0):
    """The CONDITIONAL oracle: the c with P(max T_i > c | S) = alpha at the true
    ratios.  This is the right like-for-like reference for any rule whose
    threshold is computed from the realised variance estimates."""
    tgt = 1.0 - alpha
    a, b = lo, hi
    for _ in range(80):
        m = 0.5 * (a + b)
        if cond_noreject(m, u, n, lam, quad) < tgt:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def oracle_power_k1(n, lam, delta, alpha=ABAR):
    """At K = 1 the known-lambda oracle test is X1 - X0 > z_{1-a} sigma_diff, so
    its power is free of the realised variances."""
    n = np.asarray(n, float)
    sig = sqrt(float(lam[0] if np.ndim(lam) else lam) / n[1] + 1.0 / n[0])
    return float(norm.cdf(delta / sig - norm.ppf(1.0 - alpha)))


def delta_for_oracle_power_k1(n, lam, target, alpha=ABAR):
    """The effect at which the known-lambda oracle attains `target` power."""
    n = np.asarray(n, float)
    sig = sqrt(float(lam[0] if np.ndim(lam) else lam) / n[1] + 1.0 / n[0])
    return float(sig * (norm.ppf(target) + norm.ppf(1.0 - alpha)))


def delta_for_oracle_power(n, lam, target, alpha=ABAR, quad=None):
    """General-K analogue: the common shift at which the population oracle rule
    promotes at least one candidate with probability `target`.

    The oracle here is the exact fixed-ratio constant at the true ratios, which
    is the reference the paper uses; its power is averaged over the state law.
    """
    n = np.asarray(n, float)
    K = n.size - 1
    if K == 1:
        return delta_for_oracle_power_k1(n, lam, target, alpha)
    c = oracle_c_pop(n, lam, alpha)
    pts, wt = state_grid(n, 10, 8)

    def pw(dl):
        return sum(
            w * cond_power(c, u, n, lam, dl, quad) for u, w in zip(pts, wt)
        ) / wt.sum()

    lo, hi = 0.0, 1.0
    while pw(hi) < target and hi < 64.0:
        hi *= 2.0
    for _ in range(40):
        m = 0.5 * (lo + hi)
        if pw(m) < target:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


# -------------------------------------------------------------- rule shells
class Rule:
    """A threshold function of the realised variance estimates.

    `c(u, n)` must depend on the data ONLY through quantities a desk can
    compute: the sample variances and the sample sizes.  It may not see lambda.
    """

    name = "rule"
    tier = "unconditional_candidate_measured"

    def c(self, u, n):  # pragma: no cover - interface
        raise NotImplementedError


class OracleCond(Rule):
    """Not deployable; the like-for-like conditional reference."""

    name = "oracle_conditional"
    tier = "oracle_exact_law"

    def __init__(self, lam, alpha=ABAR):
        self.lam, self.alpha = np.asarray(lam, float), alpha

    def c(self, u, n):
        return oracle_c_cond(u, n, self.lam, self.alpha)


class Plugin(Rule):
    """The exact law at the estimated ratios, at an internal level."""

    name = "plugin"
    tier = "plugin_measured"

    def __init__(self, alpha=ABAR, grid=None, pop=True):
        self.alpha, self.grid, self.pop = alpha, grid or FAST, pop

    def lam_hat(self, u, n):
        s2 = np.asarray(u, float)     # sigma_i = 1 scaling cancels in the ratio
        return s2[1:] / s2[0]

    def c(self, u, n):
        lh = self.lam_hat(u, n)
        if self.pop:
            return oracle_c_pop(n, lh, self.alpha, self.grid)
        return oracle_c_cond(u, n, lh, self.alpha)


def lam_hat_of(u, lam):
    """The estimate a desk actually sees: S_i^2/S_0^2 = lam_i u_i / u_0."""
    u = np.asarray(u, float)
    lam = np.asarray(lam, float)
    return lam * u[1:] / u[0]


# ------------------------------------------------------- FWER and power maps
def fwer_exact_k1(cfun, n, lam, npanel=24, order=12, alpha_shift=0.0):
    """Exact (quadrature) unconditional FWER at K = 1 for a threshold function.

    `cfun(lam_hat)` returns the threshold.  No draws: the expectation over the
    two scaled chi-squares is a tensor Gauss-Legendre rule on the probability
    scale, so the number is deterministic and reproducible.
    """
    n = np.asarray(n, float)
    lam = float(np.asarray(lam, float).ravel()[0])
    q0, w0 = state_nodes(n[0] - 1.0, npanel, order)
    q1, w1 = state_nodes(n[1] - 1.0, npanel, order)
    U0, U1 = np.meshgrid(q0, q1, indexing="ij")
    W = np.outer(w0, w1)
    lh = lam * U1 / U0
    cv = cfun(lh.ravel()).reshape(lh.shape)
    d = np.sqrt(lam * U1 / n[1] + U0 / n[0])
    sig = sqrt(lam / n[1] + 1.0 / n[0])
    rej = norm.cdf(-(cv * d) / sig + alpha_shift)
    return float((W * rej).sum() / W.sum())


def power_exact_k1(cfun, n, lam, delta, npanel=24, order=12):
    n = np.asarray(n, float)
    lam = float(np.asarray(lam, float).ravel()[0])
    q0, w0 = state_nodes(n[0] - 1.0, npanel, order)
    q1, w1 = state_nodes(n[1] - 1.0, npanel, order)
    U0, U1 = np.meshgrid(q0, q1, indexing="ij")
    W = np.outer(w0, w1)
    lh = lam * U1 / U0
    cv = cfun(lh.ravel()).reshape(lh.shape)
    d = np.sqrt(lam * U1 / n[1] + U0 / n[0])
    sig = sqrt(lam / n[1] + 1.0 / n[0])
    pw = norm.cdf((delta - cv * d) / sig)
    return float((W * pw).sum() / W.sum())


def sup_fwer_k1(cfun, n, lam_grid, npanel=24, order=12):
    """max over a declared lambda grid of the exact unconditional FWER."""
    vals = [fwer_exact_k1(cfun, n, l, npanel, order) for l in lam_grid]
    i = int(np.argmax(vals))
    return float(vals[i]), float(lam_grid[i]), [float(v) for v in vals]


def log_grid(lo, hi, m):
    return np.exp(np.linspace(np.log(lo), np.log(hi), int(m)))
