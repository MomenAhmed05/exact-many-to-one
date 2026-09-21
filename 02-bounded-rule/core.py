"""Core library for the W2-W6 programme on the exact many-to-one null.

Everything here is built on the paper's verified HeteroMaxT (base grid certified
to ~1e-13 by refinement; the FAST grid is certified against the base grid in
w0_discrepancy.py).  The Rao-Blackwellised FWER harness draws only the sample
variances, computes each rule's critical value from those draws exactly as on
real data, and evaluates the conditional rejection probability in closed form
using the TRUE variances.  The oracle constant - the exact critical value at the
true ratios - is the control variate, so every estimate is nominal plus a paired
difference.
"""

import os
import sys
from math import sqrt

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "01-exact-law"))  # tools/audit
sys.path.insert(0, _HERE)
from h12_hetero_exact import HeteroMaxT  # noqa: E402
from h12_quad import normal_nodes  # noqa: E402
from io_util import f64  # noqa: E402

FAST = {"u": (8, 10), "v": (8, 10), "s": (3, 10)}
GRID = {"u": (14, 14), "v": (14, 14), "s": (4, 14)}
BOOT = {"u": (4, 8), "v": (4, 8), "s": (3, 8)}
GRIDS = {"fast": FAST, "base": GRID, "boot": BOOT}
ALPHA = 0.025


# ----------------------------------------------------------------------
# critical values
# ----------------------------------------------------------------------
def _solve(args):
    n, lam, alpha, grid = args
    try:
        c = HeteroMaxT(
            np.asarray(n, float), np.sqrt([1.0] + list(lam)), grid=GRIDS[grid]
        ).critical_value(alpha)
        return f64(c)
    except Exception:  # noqa: BLE001
        return float("nan") if False else np.nan


def batch_c(n, lam_mat, alpha=ALPHA, grid="fast", workers=10, pool=None):
    """Critical values at each row of lam_mat (shape D x K), true sigma0 = 1.

    Pass an existing multiprocessing pool via `pool` to avoid creating one per
    batch (Windows worker churn is the main failure mode of long runs).
    """
    lam_mat = np.atleast_2d(np.asarray(lam_mat, float))
    args = [(tuple(np.asarray(n, float)), tuple(row), alpha, grid) for row in lam_mat]
    if pool is not None and len(args) > 1:
        out = pool.map(_solve, args)
    elif workers and workers > 1 and len(args) > 4:
        import multiprocessing as mp

        with mp.Pool(min(workers, len(args))) as p:
            out = p.map(_solve, args)
    else:
        out = [_solve(a) for a in args]
    return np.asarray(out, float)


def oracle_c(n, lam, alpha=ALPHA, grid="base"):
    return f64(
        HeteroMaxT(
            np.asarray(n, float),
            np.sqrt([1.0] + list(np.asarray(lam, float))),
            grid=GRIDS[grid],
        ).critical_value(alpha)
    )


# ----------------------------------------------------------------------
# Rao-Blackwellised evaluation
# ----------------------------------------------------------------------
def draw_var(n, lam, D, seed):
    """Draw sample variances under the true model (sigma_0 = 1)."""
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    s0 = rng.chisquare(nu[0], D) / nu[0]
    si = np.column_stack(
        [lam[i] * rng.chisquare(nu[i + 1], D) / nu[i + 1] for i in range(len(lam))]
    )
    return s0, si, si / s0[:, None]


def cond_noreject(cvec, s0, si, n, lam_true, npanels=40, order=10):
    """P(no rejection | drawn sample variances), TRUE ratios; 1-D in U."""
    n0 = n[0]
    ni = np.asarray(n[1:], float)
    z, a = normal_nodes(npanels, order)
    d = len(cvec)
    prod = np.ones((d, len(z)))
    for i in range(len(ni)):
        b = cvec[:, None] * np.sqrt(si[:, i][:, None] / ni[i] + s0[:, None] / n0)
        arg = np.sqrt(ni[i] / lam_true[i]) * (z[None, :] / np.sqrt(n0) + b)
        prod = prod * np.exp(_log_ndtr(arg))
    return (prod * a[None, :]).sum(axis=1)


def _log_ndtr(x):
    from scipy.special import log_ndtr

    return log_ndtr(x)


def rb_fwer(cvec, s0, si, n, lam_true, c_or, alpha=ALPHA):
    """Paired FWER estimate and its MC standard error."""
    d = cond_noreject(np.full(len(cvec), c_or), s0, si, n, lam_true) - cond_noreject(
        cvec, s0, si, n, lam_true
    )
    return {
        "fwer": f64(alpha + d.mean()),
        "mcse": f64(d.std(ddof=1) / sqrt(len(d))),
        "mean_c": f64(np.mean(cvec)),
        "mean_c_norej": f64(d.mean()),
        "c_oracle": f64(c_or),
    }


def rule_eval(n, lam_true, cvec, s0, si, alpha=ALPHA):
    c_or = oracle_c(n, lam_true, alpha)
    out = rb_fwer(cvec, s0, si, n, lam_true, c_or, alpha)
    out["rel_err_pct"] = 100 * (out["fwer"] - alpha) / alpha
    out["power_cost_pct"] = 100 * (out["mean_c"] - c_or) / c_or
    return out


# ----------------------------------------------------------------------
# rules that are a transform of lam_hat
# ----------------------------------------------------------------------
def f_quantile(nu_i, nu_0, gamma):
    from scipy.stats import f as fdist

    return fdist.ppf(gamma, nu_i, nu_0)


def lam_lower(n, lh, gamma):
    """Per-coordinate lower confidence bound, coverage gamma."""
    nu = np.asarray(n, float) - 1.0
    q = np.array([f_quantile(nu[i + 1], nu[0], gamma) for i in range(lh.shape[1])])
    return lh / q[None, :]


def lam_upper(n, lh, gamma):
    nu = np.asarray(n, float) - 1.0
    q = np.array([f_quantile(nu[i + 1], nu[0], gamma) for i in range(lh.shape[1])])
    return lh * q[None, :]


def rule_plugin(n, lh, s0, si):
    return lh


def rule_cb(n, gamma):
    def fn(lh, s0, si):
        return lam_lower(n, lh, gamma)

    return fn


def rule_dir_cb(n, gamma, cfun):
    """Direction-aware CB: shrink whichever way raises c.

    cfun(lam_matrix) -> critical values.  Per draw, evaluate both the lower and
    the upper bound and keep the larger critical value (conservative).
    """

    def fn(lh, s0, si):
        lo = lam_lower(n, lh, gamma)
        hi = lam_upper(n, lh, gamma)
        return ("both", lo, hi)

    return fn
