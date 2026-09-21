"""Phase A falsifiers: region coverage, corner-max validity, guarantee.

A2a analytic coverage: the region is the intersection of the independent events
     {W in [w_lo,w_hi]} and {Y_i in [y_lo_i,y_hi_i]}; its coverage is the product,
     which by the union bound is >= 1 - (gamma_w + sum a_i) = 1 - gamma.
A2b MC coverage at a large gamma where Monte Carlo has power.
A3  per-draw corner-max validity: for random interior points of the region,
     c_star(lam) must not exceed the envelope built from the corners.
A4  guarantee: on the envelope sweep cells the realised FWER is <= alpha + gamma.
"""

import multiprocessing as mp
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import env_rule as E  # noqa: E402
from core import ALPHA, batch_c, draw_var  # noqa: E402
from io_util import dump_json, f64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "env", "falsifiers.json")


def analytic_coverage(gamma=0.0005, n0_list=(4, 6, 8, 12, 20, 36), n_i=13, K=2):
    rows = []
    for n0 in n0_list:
        nu0 = n0 - 1.0
        g_w = gamma * 0.5
        a_i = gamma * 0.5 / K
        cov_w = 1.0 - g_w
        cov_y = (1.0 - a_i) ** K
        rows.append(
            {
                "n0": n0,
                "nu0": nu0,
                "gamma": gamma,
                "coverage_product": f64(cov_w * cov_y),
                "union_bound": f64(1.0 - gamma),
                "passes": bool(cov_w * cov_y >= 1.0 - gamma - 1e-12),
            }
        )
    return rows


def exact_region_coverage(n, gamma, w_frac=0.5):
    """Deterministic coverage of the implemented region (1-D quadrature).

    The region is {lam : Y_i/W in [y_lo_i/w_hi, y_hi_i/w_lo] for all i}.  Given
    W = w, the Y_i are independent, so the conditional coverage is a product of
    chi-square CDF differences and the exact coverage is one integral over W.
    No Monte Carlo, so it has power at the real gamma.
    """
    n = np.asarray(n, float)
    nu = n - 1.0
    K = n.size - 1
    g_w = gamma * w_frac
    a_i = gamma * (1.0 - w_frac) / K
    w_lo = chi2.ppf(g_w / 2.0, nu[0]) / nu[0]
    w_hi = chi2.ppf(1.0 - g_w / 2.0, nu[0]) / nu[0]
    yl = np.array([chi2.ppf(a_i / 2.0, nu[i + 1]) / nu[i + 1] for i in range(K)])
    yh = np.array([chi2.ppf(1.0 - a_i / 2.0, nu[i + 1]) / nu[i + 1] for i in range(K)])
    alo = yl / w_hi
    ahi = yh / w_lo

    def integrand(w):
        p = 1.0
        for i in range(K):
            lo = chi2.cdf(nu[i + 1] * w * alo[i], nu[i + 1])
            hi = chi2.cdf(nu[i + 1] * w * ahi[i], nu[i + 1])
            p *= max(0.0, hi - lo)
        return p * chi2.pdf(nu[0] * w, nu[0]) * nu[0]

    try:
        cov, err = quad(integrand, 0.0, np.inf, limit=200)
    except Exception as exc:  # noqa: BLE001
        return {"n": list(n), "error": repr(exc), "passes": False}
    return {
        "n": list(n),
        "gamma": gamma,
        "coverage": f64(cov),
        "quad_err": f64(err),
        "union_bound": f64(1.0 - gamma),
        "passes": bool(cov >= 1.0 - gamma - 1e-12),
    }


def mc_coverage(gamma=0.1, n=(4, 13, 13), lam=(1.0, 25.0), D=200000, seed=7):
    """Simulate W, Y_i under the true model and measure the region's coverage."""
    rng = np.random.default_rng(seed)
    nu = np.asarray(n, float) - 1.0
    K = len(nu) - 1
    g_w = gamma * 0.5
    a_i = gamma * 0.5 / K
    w_lo = chi2.ppf(g_w / 2.0, nu[0]) / nu[0]
    w_hi = chi2.ppf(1.0 - g_w / 2.0, nu[0]) / nu[0]
    y_lo = chi2.ppf(a_i / 2.0, nu[1]) / nu[1]
    y_hi = chi2.ppf(1.0 - a_i / 2.0, nu[1]) / nu[1]
    W = rng.chisquare(nu[0], D) / nu[0]
    Y = rng.chisquare(nu[1], D) / nu[1]
    lam_hat = lam[0] * Y / W
    inside = (w_lo <= W) & (w_hi >= W) & (y_lo <= Y) & (y_hi >= Y)
    # lambda_i = lam_hat_i * W / Y  => lam_i in [lam_hat*w_lo/y_hi, lam_hat*w_hi/y_lo]
    lam_lo = lam_hat * w_lo / y_hi
    lam_hi = lam_hat * w_hi / y_lo
    inside2 = (lam[0] >= lam_lo) & (lam[0] <= lam_hi)
    try:
        p = float(inside2.mean())
        se = float(np.sqrt(p * (1 - p) / D))
        box_cov = float(inside.mean())
    except Exception as exc:  # noqa: BLE001
        return {"gamma": gamma, "error": repr(exc), "passes": False}
    return {
        "gamma": gamma,
        "coverage": f64(p),
        "mcse": se,
        "target": f64(1.0 - gamma),
        "passes": bool(p >= 1.0 - gamma - 3 * se),
        "box_event_coverage": f64(box_cov),
        "ratio_region_is_superset_of_box": bool(p >= box_cov - 5.0 * se),
    }


def corner_max_check(n, lam, gamma, draws=20, interior=12, seed=99, pool=None):
    """c_star(interior lam) must not exceed the corner envelope."""
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    s0, si, lh = draw_var(n, lam, draws, seed)
    lo, hi, _ = E.lam_region(n, lh, gamma)
    c_env = E.c_env(n, s0, si, lh, gamma=gamma, pool=pool)
    rng = np.random.default_rng(seed + 1)
    worst = -np.inf
    inside_pts = []
    for j in range(draws):
        for _ in range(interior):
            u = rng.random(len(lam))
            lv = np.exp(np.log(lo[j]) + u * (np.log(hi[j]) - np.log(lo[j])))
            inside_pts.append(lv)
            try:
                worst = max(worst, float(np.max(lv >= lo[j]) * np.max(lv <= hi[j])))
            except Exception:  # noqa: BLE001
                worst = np.nan
    lams = np.vstack(inside_pts)
    cs = batch_c(n, lams, alpha=ALPHA, pool=pool).reshape(draws, interior)
    excess = cs - c_env[:, None]
    return {
        "n": list(n),
        "lam": list(lam),
        "gamma": gamma,
        "draws": draws,
        "interior_points": interior,
        "max_c_star_excess_over_envelope": f64(excess.max()),
        "passes": bool(excess.max() <= 1e-6),
    }


def main():
    exact = []
    for n in (
        [4, 13, 13],
        [6, 13, 13],
        [8, 13, 13],
        [12, 13, 13],
        [20, 13, 13],
        [36, 13, 13],
    ):
        r = exact_region_coverage(n, 0.0005)
        exact.append(r)
        print(
            f"exact coverage n0={n[0]:>2}: {r['coverage']:.6f} "
            f"(union bound {r['union_bound']:.6f}) pass={r['passes']}",
            flush=True,
        )
    res = {
        "exact_coverage": exact,
        "analytic_coverage": analytic_coverage(),
        "mc_coverage": mc_coverage(),
    }
    cells = []
    with mp.Pool(6) as pool:
        for fam, lam in (("A", (0.04, 0.04)), ("B", (1.0, 25.0))):
            for df in (3, 35):
                r = corner_max_check([df + 1, 13, 13], lam, 0.0005, pool=pool)
                r["family"] = fam
                cells.append(r)
                print(
                    f"{fam} df{df}: max interior c* excess "
                    f"{r['max_c_star_excess_over_envelope']:+.2e} "
                    f"passes={r['passes']}",
                    flush=True,
                )
    res["corner_max"] = cells
    res["all_pass"] = bool(
        all(r["passes"] for r in res["analytic_coverage"])
        and all(r["passes"] for r in res["exact_coverage"])
        and res["mc_coverage"]["passes"]
        and all(c["passes"] for c in cells)
    )
    dump_json(res, OUT)
    print("all_pass", res["all_pass"])
    print("wrote", OUT)


if __name__ == "__main__":
    main()
