"""P3 step 18: the three implementation classes and the analytic limits, with widths recorded
(plan section 3.5 and section 6.1).

Comparisons, each like-for-like (the population law and the conditional law are different objects,
so they are never compared directly):

  (i)  **K = 1 closed form.**  For one candidate the conditional statistic is a scaled normal:
       conditional on the draw, `T_1 = N(0, sigma^2)/D` with `sigma = sqrt(lambda_1/n_1 + 1/n_0)` and
       `D = sqrt(S_1^2/n_1 + S_0^2/n_0)`, so the conditional critical value is EXACTLY
       `c* = z_{1-alpha} * sigma / D`.  The enclosure's bracket must contain it.
  (ii) **Population law, independent route.**  `core.oracle_c` (the production `HeteroMaxT` solver)
       against the draw-average identity: solve for the c whose draw-average conditional no-rejection
       probability equals `1 - alpha`.  Agreement is a cross-implementation check of the population
       law against the conditional one.
  (iii) **Vanishing-ratio limit.**  As `lambda_1 -> 0` the candidate factor becomes a step, and the
       conditional no-rejection probability is exactly `Phi(c d_1 sqrt(n_0))`, so
       `c* = z_{1-alpha} / (d_1 sqrt(n_0))`.  The enclosure's bracket must contain it - this exercises
       the kink handling directly.

Widths and tolerances are recorded with every comparison.

Run:  python impl_compare.py        (~2-3 min)
Out:  deposits/impl_compare.json    (tier: region_sup_certified)
"""

from __future__ import annotations

import importlib
import os
import time

import enclose_gl
import numpy as np
import paths
import tiers
from paths import _f, _i

OUT = os.path.join(paths.ensure_dirs(), "impl_compare.json")
CORE = importlib.import_module("core")
cond_noreject = CORE.cond_noreject
draw_var = CORE.draw_var
oracle_c = CORE.oracle_c


def _z(alpha):
    from scipy.stats import norm

    return _f(norm.ppf(1.0 - _f(alpha)))


def k1_closed_form(lam1, n0, n1, s0, s1, alpha=paths.ALPHA0):
    """Exact conditional critical value for K = 1 (scaled normal)."""
    sigma = np.sqrt(_f(lam1) / _f(n1) + 1.0 / _f(n0))
    D = np.sqrt(_f(s1) / _f(n1) + _f(s0) / _f(n0))
    return _z(alpha) * sigma / D


def check_i(alpha=paths.ALPHA0, seed=20260918):
    rows = []
    for n0, n1, lam1 in (
        (13.0, 13.0, 1.0),
        (4.0, 7.0, 1.0),
        (30.0, 4.0, 25.0),
        (13.0, 13.0, 0.04),
    ):
        n = np.array([n0, n1])
        lam = np.array([lam1])
        s0, si, _lh = draw_var(n, lam, 1, seed)
        d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
        s_pt = np.sqrt(n[1:] / lam)
        exact = k1_closed_form(lam1, n0, n1, s0[0], si[0][0], alpha)
        br = enclose_gl.bracket_point_gl(s_pt, d_row, n0, alpha=alpha)
        rows.append(
            {
                "n": [n0, n1],
                "lam": [lam1],
                "c_exact": _f(exact),
                "c_L": br["c_L"],
                "c_U": br["c_U"],
                "width": br["width"],
                "contains": bool(_f(br["c_L"]) <= exact <= _f(br["c_U"]) + 1e-12),
                "abs_diff_to_mid": _f(
                    abs(exact - 0.5 * (_f(br["c_L"]) + _f(br["c_U"])))
                ),
            }
        )
    return {
        "rows": rows,
        "all_contain": bool(all(r["contains"] for r in rows)),
        "max_width": _f(max(_f(r["width"]) for r in rows)),
    }


def check_ii(alpha=paths.ALPHA0, seed=20260919, draws=3000):
    """Population critical value by the draw-average identity vs the production solver."""
    n = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    c_prod = _f(oracle_c(n, lam, alpha))
    s0, si, _lh = draw_var(n, lam, _i(draws, 3000), seed)
    lo, hi = c_prod - 0.05, c_prod + 0.05
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        p = _f(
            np.mean(
                cond_noreject(
                    np.full(_i(draws, 3000), mid),
                    s0,
                    si,
                    n,
                    lam,
                    npanels=paths.REF_PANELS,
                )
            )
        )
        if p >= 1.0 - _f(alpha):
            hi = mid
        else:
            lo = mid
    c_avg = 0.5 * (lo + hi)
    return {
        "c_production": c_prod,
        "c_draw_average": _f(c_avg),
        "abs_diff": _f(abs(c_avg - c_prod)),
        "bracket_width": _f(hi - lo),
        "draws": _i(draws),
        "agree": bool(abs(c_avg - c_prod) <= 5e-3),
    }


def check_iii(alpha=paths.ALPHA0, seed=20260920):
    """Vanishing-ratio limit: c* = z_{1-alpha} / (d_1 sqrt(n_0))."""
    rows = []
    for n0, n1 in ((13.0, 13.0), (4.0, 13.0)):
        n = np.array([n0, n1])
        lam = np.array([1e-8])
        s0, si, _lh = draw_var(n, lam, 1, seed)
        d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
        exact = _z(alpha) / (_f(d_row[0]) * np.sqrt(n0))
        s_pt = np.sqrt(n[1:] / lam)
        br = enclose_gl.bracket_point_gl(s_pt, d_row, n0, alpha=alpha)
        rows.append(
            {
                "n": [n0, n1],
                "lam": [1e-8],
                "c_exact": _f(exact),
                "c_L": br["c_L"],
                "c_U": br["c_U"],
                "width": br["width"],
                "contains": bool(_f(br["c_L"]) <= exact <= _f(br["c_U"]) + 1e-9),
                "n_eval": br["n_eval"],
            }
        )
    return {"rows": rows, "all_contain": bool(all(r["contains"] for r in rows))}


def main():
    t0 = time.time()
    out: dict = {"alpha0": _f(paths.ALPHA0), "tol": _f(paths.TOL)}
    out["i_k1_closed_form"] = check_i()
    print(
        "(i) K=1 closed form:",
        out["i_k1_closed_form"]["all_contain"],
        f"max width {_f(out['i_k1_closed_form']['max_width']):.2e}",
        flush=True,
    )
    out["ii_population_route"] = check_ii()
    print(
        "(ii) population law vs draw-average route:",
        out["ii_population_route"],
        flush=True,
    )
    out["iii_vanishing_ratio"] = check_iii()
    print(
        "(iii) vanishing-ratio limit:",
        out["iii_vanishing_ratio"]["all_contain"],
        flush=True,
    )
    out["seconds"] = _f(time.time() - t0)
    ok = (
        bool(out["i_k1_closed_form"]["all_contain"])
        and bool(out["ii_population_route"]["agree"])
        and bool(out["iii_vanishing_ratio"]["all_contain"])
    )
    out["verdict"] = "PASS" if ok else "FAIL"
    out["note"] = (
        "three implementation classes compared like-for-like: an exact K=1 closed form, "
        "an independent population route (draw-average identity), and a vanishing-ratio "
        "limit; widths are recorded with every comparison.  A homoscedastic comparison "
        "against published Dunnett tables is NOT included: scipy has no Dunnett routine "
        "and the 02-bounded-rule kernel regression already covers the population solver against an "
        "independent transformation/panel architecture."
    )
    print(f"\nRESULT: {out['verdict']} ({_f(out['seconds']):.0f}s)")
    tiers.write_deposit(
        out,
        OUT,
        "region_sup_certified",
        notes="analytic limits and cross-implementation comparisons with widths",
    )
    return out


if __name__ == "__main__":
    import sys

    try:
        res = main()
        sys.exit(0 if res.get("verdict") == "PASS" else 1)
    except Exception as exc:  # noqa: BLE001
        print("IMPL COMPARE ERROR:", repr(exc))
        sys.exit(2)
