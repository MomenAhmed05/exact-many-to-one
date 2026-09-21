"""S4.1: low-df frontier -- which admission REGION shape is cheapest?

The desk rule admits every lambda in a rectangular region built from
chi-square quantiles: lambda_i in [lam_hat_i * w_lo/y_hi_i, lam_hat_i * w_hi/y_lo_i].
Three choices are free, at the SAME certified coverage 1-gamma:

  split    how the level budget is split between the shared W term and each
           candidate group (the desk uses the uniform split a_w = a_i = gamma/2/K)
  shape    equal-tail quantiles vs the SHORTEST (HPD) chi-square interval --
           for a skewed chi-square the HPD is strictly narrower at equal mass
  scale    the largest total level mass whose EXACT quadrature coverage is still
           >= 1-gamma (coverage falls as the mass grows)

We enumerate: uniform split (current), width-optimal split, each in equal-tail and
HPD shape; all at their maximal feasible total mass.  Metric = the mean log-width
of the admitted region (per candidate group) and the resulting mean c_env over the
same draws, plus the paired realised FWER against the oracle.  The
HPD-with-optimal-split value is the FLOOR: no region of coverage 1-gamma can cost
less in c, because the width is already minimal and the budget is already fully
spent (the minimax floor attempt -- a best-possible statement, not a theorem).

Outputs: deposits/frontier/shapes.json
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, batch_c, draw_var  # noqa: E402
from env_rule import _vertex_lams, fwer_of  # noqa: E402

# Desk designs: the two families and the six benchmark dfs of the desk map.
DESK_FAMS = (("A", (0.04, 0.04)), ("B", (1.0, 25.0)))
GAMMA_DF = (3, 5, 7, 11, 19, 35)
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "frontier")
GAMMA = 0.0005
D_DRAWS = 40
RHO_GRID = tuple(np.round(np.arange(0.05, 0.96, 0.05), 3))
try:
    NAN = float("nan")
    INF = float("inf")
except Exception:  # noqa: BLE001
    NAN = 0.0
    INF = 1e308


def _f(x, default=NAN):
    """Guarded float conversion."""
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return float(default)


def _i(x, default=0):
    """Guarded int conversion."""
    try:
        return int(x)
    except Exception:  # noqa: BLE001
        return int(default)


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


POOL_N = _env_int("W26_POOL", 3)


def hpd_interval(nu, a):
    """Shortest chi-square interval with mass 1-a (equal density at the ends)."""
    nu = _f(nu)
    a = _f(a)
    if not np.isfinite(nu) or not np.isfinite(a) or a <= 0.0:
        return 0.0, INF
    try:
        if nu <= 2.0:
            return 0.0, _f(chi2.ppf(1.0 - a, nu))
        b_lo = _f(chi2.ppf(1.0 - a, nu))
        b_hi = _f(chi2.ppf(1.0 - a / 2.0, nu))

        def eq(x, b):
            return (nu / 2.0 - 1.0) * (np.log(x) - np.log(b)) - (x - b) / 2.0

        def mass_gap(b):
            try:
                alo = _f(brentq(lambda x: eq(x, b), 1e-12, b * 0.999999, maxiter=200))
            except Exception:  # noqa: BLE001
                alo = 0.0
            return _f(chi2.cdf(b, nu) - chi2.cdf(alo, nu)) - (1.0 - a)

        b = _f(brentq(mass_gap, b_lo, b_hi, maxiter=200))
        try:
            a_pt = _f(brentq(lambda x: eq(x, b), 1e-12, b * 0.999999, maxiter=200))
        except Exception:  # noqa: BLE001
            a_pt = 0.0
        return a_pt, b
    except Exception:  # noqa: BLE001
        return _f(chi2.ppf(a / 2.0, nu)), _f(chi2.ppf(1.0 - a / 2.0, nu))


def levels(n, a_w, a_i, shape="eqtail"):
    """Raw levels (w_lo, w_hi, y_lo, y_hi) as /nu ratios."""
    n = np.asarray(n, float)
    nu = n - 1.0
    K = n.size - 1
    if shape == "hpd":
        w_lo, w_hi = hpd_interval(nu[0], a_w)
        ys = [hpd_interval(nu[i + 1], a_i) for i in range(K)]
    else:
        w_lo = _f(chi2.ppf(a_w / 2.0, nu[0]))
        w_hi = _f(chi2.ppf(1.0 - a_w / 2.0, nu[0]))
        ys = [
            (
                _f(chi2.ppf(a_i / 2.0, nu[i + 1])),
                _f(chi2.ppf(1.0 - a_i / 2.0, nu[i + 1])),
            )
            for i in range(K)
        ]
    y_lo = np.array([ys[i][0] / nu[i + 1] for i in range(K)])
    y_hi = np.array([ys[i][1] / nu[i + 1] for i in range(K)])
    return _f(w_lo / nu[0]), _f(w_hi / nu[0]), y_lo, y_hi


def coverage(n, w_lo, w_hi, y_lo, y_hi, nquad=200):
    """Exact P(lambda in region) by 1-D quadrature (same integral as env_rule)."""
    n = np.asarray(n, float)
    nu = n - 1.0
    K = n.size - 1
    alo = y_lo / w_hi
    ahi = y_hi / w_lo

    def integrand(w):
        p = 1.0
        for i in range(K):
            lo = _f(chi2.cdf(nu[i + 1] * w * alo[i], nu[i + 1]))
            hi = _f(chi2.cdf(nu[i + 1] * w * ahi[i], nu[i + 1]))
            p *= max(0.0, hi - lo)
        return p * _f(chi2.pdf(nu[0] * w, nu[0])) * nu[0]

    try:
        cov, _err = quad(integrand, 0.0, np.inf, limit=nquad)
        return _f(cov)
    except Exception:  # noqa: BLE001
        return NAN


def max_total(n, gamma, rho, shape, margin=1e-7):
    """Largest total level mass t at split ratio rho with coverage >= 1-gamma."""
    n = np.asarray(n, float)
    K = n.size - 1

    def cov_at(t):
        w_lo, w_hi, y_lo, y_hi = levels(n, t * rho, t * (1.0 - rho) / K, shape)
        return coverage(n, w_lo, w_hi, y_lo, y_hi)

    lo, hi = 1e-6, 1e-1
    try:
        if cov_at(hi) >= 1.0 - gamma + margin:
            return hi
        for _ in range(40):
            mid = _f(np.sqrt(lo * hi))
            if cov_at(mid) >= 1.0 - gamma + margin:
                lo = mid
            else:
                hi = mid
        return _f(lo)
    except Exception:  # noqa: BLE001
        return NAN


def width(n, w_lo, w_hi, y_lo, y_hi):
    """Mean log-width per candidate group (lam_hi/lam_lo = (w_hi/y_lo)/(w_lo/y_hi))."""
    try:
        lw = _f(np.log(w_hi / w_lo))
        return lw + _f(np.mean(np.log(y_hi / y_lo)))
    except Exception:  # noqa: BLE001
        return NAN


def c_of_region(n, lam_lo, lam_hi, lam_hat, pool):
    """Replicates env_rule.c_env's final step for an explicitly supplied region."""
    rows = _vertex_lams(n, lam_lo, lam_hi, lam_hat, 64)
    lam_mat = np.vstack(rows)
    D = np.atleast_2d(lam_hat).shape[0]
    cs = batch_c(n, lam_mat, alpha=ALPHA, pool=pool, workers=10)
    return cs.reshape(len(rows), D).max(axis=0)


def shape_spec(n, gamma, rho, shape):
    t = max_total(n, gamma, rho, shape)
    if not np.isfinite(_f(t)):
        return None
    K = n.size - 1
    w_lo, w_hi, y_lo, y_hi = levels(n, t * rho, t * (1.0 - rho) / K, shape)
    return {
        "rho": _f(rho),
        "shape": shape,
        "total_mass": _f(t),
        "w_lo": w_lo,
        "w_hi": w_hi,
        "y_lo": y_lo.tolist(),
        "y_hi": y_hi.tolist(),
        "coverage": coverage(n, w_lo, w_hi, y_lo, y_hi),
        "width": width(n, w_lo, w_hi, y_lo, y_hi),
    }


def main():
    t0 = time.time()
    rows = []
    with mp.Pool(POOL_N) as pool:
        for fam_name, lam_vec in DESK_FAMS:
            for df0 in GAMMA_DF:
                n = np.array([df0 + 1, 13, 13], float)
                lam_true = np.asarray(lam_vec, float)
                s0, si, lam_hat = draw_var(n, lam_true, D_DRAWS, 4242 + df0)
                cands = []
                for shape in ("eqtail", "hpd"):
                    for rho in RHO_GRID:
                        sp = shape_spec(n, GAMMA, rho, shape)
                        if sp is not None:
                            cands.append(sp)
                ref = {"uniform": shape_spec(n, GAMMA, 0.5, "eqtail")}
                for shape in ("eqtail", "hpd"):
                    sh = [
                        sp
                        for sp in cands
                        if sp["shape"] == shape and np.isfinite(_f(sp["width"]))
                    ]
                    if sh:
                        ref[f"best_{shape}"] = min(sh, key=lambda s: _f(s["width"]))
                rec = {
                    "family": fam_name,
                    "df0": i64(df0),
                    "n": [f64(x) for x in n],
                    "draws": i64(D_DRAWS),
                }
                uni_c = None
                for tag, sp in ref.items():
                    if sp is None or not np.isfinite(_f(sp["width"])):
                        continue
                    y_lo = np.asarray(sp["y_lo"], float)
                    y_hi = np.asarray(sp["y_hi"], float)
                    lam_lo = lam_hat * (_f(sp["w_lo"]) / y_hi)[None, :]
                    lam_hi = lam_hat * (_f(sp["w_hi"]) / y_lo)[None, :]
                    try:
                        cvec = c_of_region(n, lam_lo, lam_hi, lam_hat, pool)
                        fw = fwer_of(n, lam_true, cvec, s0, si)
                    except Exception as exc:  # noqa: BLE001
                        rec[f"{tag}_error"] = repr(exc)
                        continue
                    entry = {
                        "rho": f64(sp["rho"]),
                        "shape": sp["shape"],
                        "total_mass": f64(sp["total_mass"]),
                        "coverage": f64(sp["coverage"]),
                        "width": f64(sp["width"]),
                        "mean_c": fw["mean_c"],
                        "fwer": fw["fwer"],
                        "mcse": fw["mcse"],
                        "c_oracle": fw["c_oracle"],
                        "power_cost_pct": fw["power_cost_pct"],
                    }
                    if tag == "uniform":
                        uni_c = _f(entry["mean_c"])
                    if uni_c:
                        entry["price_vs_uniform_pct"] = f64(
                            100.0 * (_f(entry["mean_c"]) - uni_c) / uni_c
                        )
                    rec[tag] = entry
                rows.append(rec)
                seen = {
                    k: (rec[k]["width"] if k in rec else None)
                    for k in ("uniform", "best_eqtail", "best_hpd")
                }
                print(
                    f"{fam_name} df{df0}: widths "
                    + " | ".join(f"{k} {_f(v):.4f}" for k, v in seen.items()),
                    flush=True,
                )
    dump_json(
        {
            "rows": rows,
            "gamma": GAMMA,
            "alpha": ALPHA,
            "rho_grid": list(RHO_GRID),
            "seconds": f64(time.time() - t0),
            "note": (
                "width = mean log-width per candidate group at certified "
                "coverage >= 1-gamma; best_hpd is the floor (widest budget "
                "spend, narrowest interval shape) -- a best-possible "
                "statement among the shapes tried, not a theorem."
            ),
        },
        os.path.join(OUT, "shapes.json"),
    )
    print(f"shapes done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
