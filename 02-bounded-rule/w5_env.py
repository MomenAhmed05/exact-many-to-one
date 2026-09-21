"""Phase A3: envelope on the W5 grid, K <= 3 (falsifier A3 + power map).

Per (K, shape, vol, alpha) cell and per benchmark m, deposit:
  - the deterministic guarantee check: max over draws of the conditional
    rejection at the TRUE lambda under c_env must be <= alpha, restricted to
    draws where lambda_true lies inside the region (that is the event on which
    the envelope argument applies);
  - the empirical region-miss rate (frac_outside), which enters the bound
    FWER <= alpha + P(lambda not in Lambda);
  - the power cost of the envelope vs the oracle and vs the plug-in rule.

Cheap by construction: no oracle bisection per draw beyond the vertex solves.
"""

import multiprocessing as mp
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import env_rule as _er  # noqa: E402
from core import batch_c, cond_noreject, draw_var, oracle_c  # noqa: E402
from env_rule import c_env, lam_region  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402
from w5_boundary import MS_BY_K, lamvec, sizes  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "env")
_SHAPES = ("bench_small", "short_cand", "balanced")
_VOLS = ("noisy_cand", "noisy_bench", "graded")


def _env_floats(name, default):
    try:
        raw = os.environ.get(name) or default
        return tuple(float(x) for x in str(raw).split(","))
    except Exception:  # noqa: BLE001
        return tuple(float(x) for x in str(default).split(","))


_ALPHAS = _env_floats("W26_ENV_ALPHAS", "0.025,0.05")


def _env_tuple(name, default):
    try:
        raw = os.environ.get(name) or default
        return tuple(int(x) for x in str(raw).split(","))
    except Exception:  # noqa: BLE001
        return tuple(int(x) for x in str(default).split(","))


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


def _env_float(name, default):
    try:
        return float(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return float(default)


KS = _env_tuple("W26_ENV_KS", "1,2,3")
MS = _env_tuple("W26_ENV_MS", "4,8,36")
DRAWS = _env_int("W26_ENV_D", 100)
GAMMA = _env_float("W26_ENV_GAMMA", 0.0005)
TIGHT = os.environ.get("W26_ENV_TIGHT", "1") == "1"
POOL_N = _env_int("W26_POOL", 4)
MAXV = _env_int("W26_ENV_MAXV", 64)
SEEDS = _env_tuple("W26_ENV_SEEDS", "101")


def _tight_scale(n, gamma):
    """Exact-coverage scale if env_rule provides it (post-A5b), else None."""
    fn = getattr(_er, "tight_scale", None)
    if fn is None:
        return None
    try:
        return fn(n, gamma)
    except Exception:  # noqa: BLE001
        return None


def _region(n, lam_hat, gamma, tight):
    """Region bounds; `tight` selects the exact-coverage scaled region."""
    scale = _tight_scale(n, gamma) if tight else None
    if scale is None:
        return lam_region(n, lam_hat, gamma)
    try:
        return lam_region(n, lam_hat, gamma, scale=scale)
    except TypeError:
        return lam_region(n, lam_hat, gamma)


def _cenv(n, s0, si, lh, gamma, tight, pool, alpha):
    if tight and getattr(_er, "tight_scale", None) is not None:
        try:
            return c_env(
                n,
                s0,
                si,
                lh,
                alpha=alpha,
                gamma=gamma,
                pool=pool,
                tight=True,
                max_vertices=MAXV,
            )
        except TypeError:
            pass
    return c_env(n, s0, si, lh, alpha=alpha, gamma=gamma, pool=pool)


def _draw_multi(n, lam, draws, seeds):
    """Paired multi-seed draws: one draw_var batch per seed, stacked.

    The conditional-rejection estimator stays Rao-Blackwellised (it conditions on
    s^2), and stacking seeds widens the draw window without changing the per-row
    arithmetic; the seed list is recorded in every row.
    """
    s0s, sis, lhs = [], [], []
    for seed in seeds:
        try:
            s0, si, lh = draw_var(n, lam, draws, int(seed))
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(f"draw_var failed for seed {seed}") from exc
        s0s.append(s0)
        sis.append(si)
        lhs.append(lh)
    return (np.concatenate(s0s, axis=0), np.concatenate(sis, axis=0),
            np.concatenate(lhs, axis=0))


def run_cell(K, shape, vol, alpha, pool):
    lam = np.asarray(lamvec(vol, K), float)
    rows = []
    for m in [x for x in MS_BY_K[K] if x in MS]:
        n = np.asarray(sizes(shape, m, K), float)
        s0, si, lh = _draw_multi(n, lam, DRAWS, SEEDS)
        lo, hi, reg = _region(n, lh, GAMMA, TIGHT)
        c_e = _cenv(n, s0, si, lh, GAMMA, TIGHT, pool, alpha)
        c_p = batch_c(n, lh, alpha=alpha, pool=pool)
        c_o = oracle_c(n, list(lam), alpha)
        rej_e = 1.0 - cond_noreject(c_e, s0, si, n, lam)
        rej_o = 1.0 - cond_noreject(np.full(len(s0), c_o), s0, si, n, lam)
        rej_p = 1.0 - cond_noreject(c_p, s0, si, n, lam)
        inside = np.all((lam[None, :] >= lo) & (lam[None, :] <= hi), axis=1)
        try:
            n_out = int(np.count_nonzero(~inside))
            excess = rej_e - rej_o
            worst = (
                float(np.max(excess[inside])) if bool(np.any(inside)) else float("nan")
            )
            mean_excess = float(np.mean(excess))
            mean_rej = float(np.mean(rej_e))
            ndraw = int(len(rej_e))
            mcse = float(np.std(rej_e, ddof=1) / np.sqrt(ndraw))
            mean_c = float(np.mean(c_e))
            mean_cp = float(np.mean(c_p))
            mean_rejp = float(np.mean(rej_p))
            cov = reg.get("coverage_exact")
            scale = reg.get("scale")
            covv = float(cov) if cov is not None else float("nan")
            scalev = float(scale) if scale is not None else 1.0
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(
                f"draw aggregation failed for K={K} {shape}/{vol} a={alpha} m={m}"
            ) from exc
        rows.append(
            {
                "m": m,
                "df": m - 1,
                "n": list(n),
                "lam": list(lam),
                "alpha": alpha,
                "gamma": GAMMA,
                "draws": DRAWS,
                "draws_total": i64(ndraw),
                "seeds": [i64(s) for s in SEEDS],
                "c_oracle": f64(c_o),
                "mean_c_env": f64(mean_c),
                "c_plugin_mean": f64(mean_cp),
                "power_cost_pct": f64(100 * (mean_c - c_o) / c_o),
                "plugin_cost_pct": f64(100 * (mean_cp - c_o) / c_o),
                "mean_rej_env": f64(mean_rej),
                "mcse_rej_env": f64(mcse),
                "max_excess_vs_oracle_inside": f64(worst),
                "mean_excess_vs_oracle": f64(mean_excess),
                "n_outside": i64(n_out),
                "frac_outside": f64(n_out / ndraw),
                "fwer_est": f64(alpha + mean_excess),
                "fwer_bound": f64(alpha + mean_excess + n_out / ndraw),
                "plugin_fwer": f64(mean_rejp),
                "region_scale": f64(scalev),
                "region_coverage": f64(covv),
                "guarantee": f64(alpha + GAMMA),
                "passes_bound": bool(
                    alpha + mean_excess + n_out / ndraw <= alpha + GAMMA
                ),
                "passes_pointwise": bool(np.isnan(worst) or worst <= 1e-9),
            }
        )
    return {
        "K": K,
        "shape": shape,
        "vol": vol,
        "alpha": alpha,
        "gamma": GAMMA,
        "tight": bool(TIGHT),
        "rows": rows,
    }


def main():
    tag = "tight" if TIGHT else "ub"
    summary = []
    with mp.Pool(POOL_N) as pool:
        for K in KS:
            for shape in _SHAPES:
                for vol in _VOLS:
                    if K == 1 and vol == "noisy_cand":
                        continue
                    for alpha in _ALPHAS:
                        res = run_cell(K, shape, vol, alpha, pool)
                        fn = (
                            f"w5env_{tag}_K{K}_{res['shape']}_{res['vol']}"
                            f"_a{str(alpha).replace('.', '')}.json"
                        )
                        dump_json(res, os.path.join(OUT, fn))
                        for r in res["rows"]:
                            print(
                                f"K{K} {res['shape']:<11s} "
                                f"{res['vol']:<11s} a{alpha} m{r['m']:>2}: "
                                f"cost {r['power_cost_pct']:+6.1f}% | "
                                f"fwer {r['fwer_est']:.4f} | exc_oracle "
                                f"{r['max_excess_vs_oracle_inside']:+.1e} | out "
                                f"{r['n_outside']:>3d} | bound "
                                f"{r['fwer_bound']:.4f} | plug "
                                f"{r['plugin_fwer']:.4f}",
                                flush=True,
                            )
                            summary.append(
                                {
                                    "K": K,
                                    "shape": res["shape"],
                                    "vol": res["vol"],
                                    "alpha": alpha,
                                    "m": r["m"],
                                    **{
                                        k: r[k]
                                        for k in (
                                            "power_cost_pct",
                                            "fwer_est",
                                            "max_excess_vs_oracle_inside",
                                            "n_outside",
                                            "fwer_bound",
                                            "plugin_fwer",
                                            "passes_bound",
                                            "passes_pointwise",
                                        )
                                    },
                                }
                            )
    viol = [s for s in summary if not s["passes_bound"]]
    cond_viol = [s for s in summary if not s["passes_pointwise"]]
    dump_json(
        {
            "tag": tag,
            "n_rows": len(summary),
            "gamma": GAMMA,
            "tight": bool(TIGHT),
            "violations": viol,
            "conditional_violations": cond_viol,
            "all_pass": bool(not viol and not cond_viol),
        },
        os.path.join(OUT, f"w5env_{tag}_summary.json"),
    )
    print(
        f"rows {len(summary)} | bound violations {len(viol)} | "
        f"pointwise violations {len(cond_viol)}"
    )


if __name__ == "__main__":
    main()
