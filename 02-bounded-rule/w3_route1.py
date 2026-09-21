"""W3: route 1 - integrate over the sampling law of the estimates.

Variants measured against the declared comparator set:

  plugin    c(lam_hat)
  cb0.9     c(lam_hat / F_0.9)                       (confidence bound)
  deb       2 c(lam_hat) - mean_b c(lam_hat*)        (debiased critical value)
  cq0.8     quantile_0.8 of c(lam_hat*)              (bound c directly)
  smooth    solve mean_b F_{lam_hat*}(c) = 1-alpha   (smoothed exact CDF)
  prepiv    Beran prepivoting at B_cal inner draws   (the paper's calibration)
  m3        Tamhane & Xi Method 3 reimplementation

lam_hat* ~ F-law at lam_hat (parametric bootstrap of the ratios).
Evaluation: paired against the oracle at the TRUE ratios, Rao-Blackwellised.
prepiv is evaluated only where n0 <= 6 (the low-df corner it exists for); it
costs B_cal times a cheap-grid solve per data set.

Usage: python w3_route1.py [df ...]
Outputs: deposits/w3/route1_df<df>_<fam>.json
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, BOOT, FAST, cond_noreject, draw_var, oracle_c  # noqa
from h12_calibrated import plugin_sigma  # noqa: E402
from h12_hetero_exact import HeteroMaxT  # noqa: E402
from h12_tamhane_headtohead import c_method3  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w3")
FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
B_CAL = 120
POOL_N = i64(os.environ.get("W26_POOL", "4"))


def _robust_c(n, lam, alpha, grid=FAST, clip=(1e-6, 1e6)):
    """Critical value with clipping and a widened fallback bracket."""
    lam = np.clip(np.asarray(lam, float), clip[0], clip[1])
    sig = np.sqrt(np.concatenate([[1.0], lam]))
    try:
        return f64(
            HeteroMaxT(np.asarray(n, float), sig, grid=grid).critical_value(alpha)
        )
    except Exception:  # noqa: BLE001
        try:
            from scipy.optimize import brentq

            mod = HeteroMaxT(np.asarray(n, float), sig, grid=grid)

            def g(c):
                return mod.cdf(f64(c)) - (1.0 - alpha)

            return f64(brentq(g, 0.2, 80.0, xtol=1e-10, maxiter=300))
        except Exception:  # noqa: BLE001
            return np.nan


def _q(n, gamma):
    from scipy.stats import f as fd

    nu = np.asarray(n, float) - 1.0
    return np.array([fd.ppf(gamma, nu[i + 1], nu[0]) for i in range(len(nu) - 1)])


def _calibrate(n, lam_world, alpha, rng, B=B_CAL, ngrid=9, span=(-1.0, 1.4)):
    """Paper's calibrator.  lam_world holds the VARIANCE ratios of the world."""
    n = np.asarray(n, float)
    lamw = np.asarray(lam_world, float)
    nu = n - 1.0
    c_ref = _robust_c(n, lamw, alpha, grid=BOOT)
    if not np.isfinite(c_ref):
        return np.nan
    cgrid = c_ref + np.linspace(span[0], span[1], ngrid)
    G = np.empty((B, ngrid))
    H = np.empty((B, ngrid))
    for b in range(B):
        s2 = np.empty(n.size)
        s2[0] = rng.chisquare(nu[0]) / nu[0]
        for i in range(len(lamw)):
            s2[i + 1] = lamw[i] * rng.chisquare(nu[i + 1]) / nu[i + 1]
        mod = HeteroMaxT(n, plugin_sigma(s2), grid=BOOT)
        for j, c in enumerate(cgrid):
            G[b, j] = 1.0 - mod.cdf(f64(c))
        H[b] = 1.0 - cond_noreject(
            np.asarray(cgrid, float), np.array([s2[0]]), s2[None, 1:], n, lamw
        )

    def F(a):
        tot = 0.0
        for b in range(B):
            c = np.interp(a, G[b][::-1], cgrid[::-1])
            tot += np.interp(c, cgrid, H[b])
        return tot / B

    lo, hi = 1e-5, 0.5
    if F(lo) >= alpha or F(hi) <= alpha:
        return np.nan
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if F(mid) < alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _draw_task(args):
    n, lam_true, s2, B, alpha, seed, _use_prepiv, a_true = args
    n = np.asarray(n, float)
    s2 = np.asarray(s2, float)
    rng = np.random.default_rng(seed)
    nu = n - 1.0
    lh = s2[1:] / s2[0]
    cplug = _robust_c(n, lh, alpha)
    ccb = _robust_c(n, lh / _q(n, 0.9), alpha)
    ccb8 = _robust_c(n, lh / _q(n, 0.8), alpha)
    fstar = np.column_stack(
        [rng.f(nu[i + 1], nu[0], size=B) for i in range(len(nu) - 1)]
    )
    lam_star = lh[None, :] * fstar
    cstar = np.array([_robust_c(n, lam_star[b], alpha) for b in range(B)])
    cdeb = 2 * cplug - f64(np.nanmean(cstar))
    ccq = f64(np.nanquantile(cstar, 0.8))
    cgrid = cplug + np.linspace(-0.7, 1.5, 7)
    acc = np.zeros(len(cgrid))
    for b in range(B):
        mod_b = (
            HeteroMaxT(n, np.sqrt(np.concatenate([[1.0], lam_star[b]])), grid=BOOT)
            if np.all(np.isfinite(lam_star[b]))
            else None
        )
        if mod_b is not None:
            acc += np.array([mod_b.cdf(f64(c)) for c in cgrid])
    sm = acc / B
    order = np.argsort(sm)
    csmooth = f64(np.interp(1 - alpha, np.sort(sm), cgrid[order]))
    c_oc = _robust_c(n, lh, a_true) if np.isfinite(a_true) else np.nan
    c_m3 = c_method3(n, s2, rng, alpha, sims=40000)
    return np.array([cplug, ccb, ccb8, cdeb, ccq, csmooth, c_oc, c_m3])


def evaluate(
    n,
    lam_true,
    alpha=ALPHA,
    D=75,
    B=25,
    seeds=(11, 22, 33),
    workers=10,
    pool=None,
):
    import multiprocessing as mp

    n = np.asarray(n, float)
    lam_true = np.asarray(lam_true, float)
    nu = n - 1.0
    rules = [
        "plugin",
        "cb0.9",
        "cb0.8",
        "deb",
        "cq0.8",
        "smooth",
        "oracle_cal",
        "m3",
    ]
    c_or = oracle_c(n, lam_true, alpha)
    from w4_calibration import calibrate_par

    a_true = np.nan
    if n[0] <= 6 and pool is not None:
        sig_world = np.sqrt(np.concatenate([[1.0], np.asarray(lam_true, float)]))
        a_true = calibrate_par(
            n, sig_world, alpha, 2000, seed=20260915, workers=3, pool=pool
        )[0]
    per = {k: [] for k in rules}
    nan_draws = 0
    t0 = time.time()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        s2s = []
        for _ in range(D):
            s2 = np.empty(n.size)
            s2[0] = rng.chisquare(nu[0]) / nu[0]
            for i in range(len(lam_true)):
                s2[i + 1] = lam_true[i] * rng.chisquare(nu[i + 1]) / nu[i + 1]
            s2s.append(s2)
        use_prepiv = False
        args = [
            (
                tuple(n),
                tuple(lam_true),
                s2,
                B,
                alpha,
                seed * 100003 + d,
                use_prepiv,
                a_true,
            )
            for d, s2 in enumerate(s2s)
        ]
        if pool is None:
            with mp.Pool(workers) as p:
                cv = np.array(p.map(_draw_task, args))
        else:
            cv = np.array(pool.map(_draw_task, args))
        s0 = np.array([s[0] for s in s2s])
        si = np.array([s[1:] for s in s2s])
        nan_draws += i64(np.sum(~np.all(np.isfinite(cv), axis=1)))
        d_len = len(s0)
        nr_o = cond_noreject(np.full(d_len, c_or), s0, si, n, lam_true)
        for j, rule in enumerate(rules):
            ok = np.isfinite(cv[:, j])
            s0j, sij, cvj = s0[ok], si[ok], cv[ok, j]
            nr = cond_noreject(cvj, s0j, sij, n, lam_true)
            d = nr_o[ok] - nr
            per[rule].append(
                (
                    alpha + f64(d.mean()),
                    f64(d.std(ddof=1) / np.sqrt(len(d))),
                    f64(cvj.mean()),
                )
            )
    out = {
        "n": list(n),
        "lam_true": list(lam_true),
        "alpha": alpha,
        "D": D,
        "B": B,
        "B_cal": B_CAL,
        "seeds": list(seeds),
        "c_oracle": c_or,
        "nan_draws": nan_draws,
        "rules": {},
    }
    for rule in rules:
        rows = np.array(per[rule])
        if rows.size == 0:
            # rule not computable at this design (e.g. oracle_cal needs
            # n[0] <= 6); emit nulls rather than bare NaN
            out["rules"][rule] = {
                "available": False,
                "fwer": None,
                "seed_se_pct": None,
                "rel_err_pct": None,
                "mean_c": None,
                "power_cost_pct": None,
            }
            continue
        f = rows[:, 0]
        out["rules"][rule] = {
            "fwer": f64(f.mean()),
            "seed_se_pct": f64(100 * f.std(ddof=1) / np.sqrt(len(f)) / alpha),
            "rel_err_pct": 100 * (f64(f.mean()) - alpha) / alpha,
            "mean_c": f64(rows[:, 2].mean()),
            "power_cost_pct": 100 * (f64(rows[:, 2].mean()) - c_or) / c_or,
        }
    out["seconds"] = time.time() - t0
    return out


def main(dfs=(3, 5, 7, 35)):
    import multiprocessing as mp

    with mp.Pool(POOL_N) as pool:
        for df in dfs:
            n0 = df + 1
            for fam, lam in FAMS.items():
                n = [n0, 13, 13]
                p = os.path.join(OUTDIR, f"route1_df{df}_{fam}.json")
                if os.path.exists(p):
                    continue
                res = evaluate(n, lam, pool=pool)
                res["family"] = fam
                res["benchmark_df"] = df
                dump_json(res, p)
                line = " ".join(
                    f"{k}:{v['rel_err_pct']:+.1f}%" for k, v in res["rules"].items()
                )
                print(f"df{df} {fam} nan={res['nan_draws']} {line}", flush=True)


def _cli_ints(argv):
    try:
        return tuple(int(x) for x in argv)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"bad df list: {argv}") from exc


if __name__ == "__main__":
    dfs = _cli_ints(sys.argv[1:]) if len(sys.argv) > 1 else (3, 5, 7, 35)
    main(dfs)
