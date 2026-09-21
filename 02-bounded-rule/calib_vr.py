"""Track B: variance-reduced calibration of the plug-in rule (alpha*).

Problem (W4): the calibrated level alpha* is itself estimated from B inner worlds,
and its noise dominates the end-to-end falsifier -- a single-estimate |z| of
-5.86 / +2.70 / -3.58 at B = 600 / 600 / 2400, and a pass needs B of order 1e5.

Two variance reductions, both reusing W4's inner-world computation unchanged:
  1. ANTITHETIC inner worlds: half the draws are generated from uniforms u and
     half from 1-u through the chi-square ppf, so pairs are negatively correlated.
  2. CONTROL VARIATES with EXACTLY known means: for each inner world,
     z = s2_0/sigma_0^2 - 1 has mean 0 exactly (E[chi2_nu/nu] = 1) and
     w = log(s2_0/sigma_0^2) - (digamma(nu0/2) + log(2/nu0)) has mean 0 exactly.
     The benchmark variance drives every comparison through the shared
     denominator, so regressing the per-world term on (z, w) at each level a and
     subtracting the fitted control part leaves the mean unchanged and removes
     the correlated noise.
The estimator remains unbiased: every control is centered at its exact mean.
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np
from scipy.special import digamma
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w4_calibration as W  # noqa: E402
from core import ALPHA, batch_c, cond_noreject, draw_var, oracle_c  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "calib_vr")
CASES = W.CASES


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


B_SWEEP = _env_tuple("W26_VR_BS", "600,1500,3000")
SEEDS = _env_tuple("W26_VR_SEEDS", "11,22,33,44")
D_OUTER = _env_int("W26_VR_D", 1200)
POOL_N = _env_int("W26_POOL", 4)
TOL = 1e-12
try:
    NAN = float("nan")
except Exception:  # noqa: BLE001
    NAN = 0.0


def _finite_or_none(x):
    """JSON-safe value: non-finite becomes null (bare NaN is invalid strict JSON)."""
    try:
        v = float(x)
        return v if v == v and abs(v) != float("inf") else None
    except Exception:  # noqa: BLE001
        return None


def inner_worlds(n, sw, B, seed, antithetic=True):
    """B inner variance worlds plus exactly-centered controls.

    Returns (s2s, ctrl, names).  With antithetic=True, draws come in pairs from
    uniforms u and 1-u so the pair average has lower variance.

    NOTE: the antithetic branch must sit under `if antithetic:`.  An earlier
    lint-driven edit left it inside the `except` handler, which made the flag
    dead code -- plain and antithetic draws were identical, so the reading
    "antithetic pairing does nothing" was a tautology rather than a measurement.
    """
    n = np.asarray(n, float)
    sw = np.asarray(sw, float)
    nu = n - 1.0
    K = n.size
    rng = np.random.default_rng(seed)
    try:
        if antithetic:
            half = int(np.ceil(B / 2.0))
            u = np.clip(rng.random((half, K)), TOL, 1.0 - TOL)
            first = chi2.ppf(u, nu) / nu * sw**2
            second = chi2.ppf(1.0 - u, nu) / nu * sw**2
            s2s = np.vstack([first, second])[:B]
        else:
            s2s = rng.chisquare(nu, size=(B, K)) / nu * sw**2
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("inner-world draw failed") from exc
    try:
        z0 = s2s[:, 0] / sw[0] ** 2 - 1.0
        w0 = np.log(s2s[:, 0] / sw[0] ** 2) - (
            digamma(nu[0] / 2.0) + np.log(2.0 / nu[0])
        )
        zbar = (s2s[:, 1:] / sw[1:] ** 2 - 1.0).mean(axis=1)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("control construction failed") from exc
    ctrl = np.column_stack([z0, w0, zbar])
    return s2s, ctrl, ["z0", "w0", "zbar"]


def _curves(n, sw, cgrid, s2s, pool):
    args = [(tuple(n), tuple(sw), tuple(s2), tuple(cgrid)) for s2 in s2s]
    if pool is None:
        with mp.Pool(POOL_N) as p:
            out = p.map(W._boot_task, args)
    else:
        out = pool.map(W._boot_task, args)
    G = np.array([o[0] for o in out])
    H = np.array([o[1] for o in out])
    return G, H


def _term_at(a, G, H, cgrid):
    """Per-world term H_b(c_b(a)) with c_b(a) from the world's own size curve."""
    cs = np.array([np.interp(a, G[b][::-1], cgrid[::-1]) for b in range(G.shape[0])])
    return np.array([np.interp(cs[b], cgrid, H[b]) for b in range(G.shape[0])]), cs


def calibrate(n, sigma_world, alpha, B, seed, variant="vr", pool=None, antithetic=True):
    """Calibrated level a* with variant in {plain, anti, vr}.

    plain : W4's estimator (no antithetic, no controls)
    anti  : antithetic worlds only
    vr    : antithetic worlds + regression control variates
    """
    n = np.asarray(n, float)
    sw = np.asarray(sigma_world, float) / np.asarray(sigma_world, float)[0]
    try:
        raw_ref = W.HeteroMaxT(np.asarray(n, float), sw, grid=W.BOOT).critical_value(
            alpha
        )
        c_ref = float(np.asarray(raw_ref, dtype=object).reshape(-1)[0])
        cgrid = c_ref + np.linspace(W.SPAN[0], W.SPAN[1], W.NGRID)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("reference critical value failed") from exc
    use_anti = variant in ("anti", "vr")
    use_ctrl = variant == "vr"
    s2s, ctrl, names = inner_worlds(n, sw, B, seed, antithetic=use_anti)
    G, H = _curves(n, sw, cgrid, s2s, pool)
    diag = {
        "B": i64(B),
        "variant": variant,
        "seed": i64(seed),
        "controls": names if use_ctrl else [],
    }

    def F(a):
        try:
            y, _cs = _term_at(a, G, H, cgrid)
            y_mean = float(y.mean())
            y_se = float(y.std(ddof=1) / np.sqrt(y.size))
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError("term evaluation failed") from exc
        if not use_ctrl:
            return y_mean, y_se, NAN
        X = np.column_stack([np.ones(y.size), ctrl])
        try:
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            fitted = X @ beta
            resid = y - fitted
            ss_tot = float(((y - y.mean()) ** 2).sum())
            r2 = 1.0 - float((resid**2).sum()) / ss_tot if ss_tot > 0 else float("nan")
            # controls are exactly centered: E[ctrl] = 0, so this stays unbiased
            adj = float(resid.mean()) + float(beta[0])
            return adj, float(resid.std(ddof=1) / np.sqrt(y.size)), r2
        except Exception:  # noqa: BLE001
            return y_mean, y_se, NAN

    lo, hi = 1e-5, 0.5
    f_lo = F(lo)[0]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if F(mid)[0] < alpha:
            lo = mid
        else:
            hi = mid
    a_star = 0.5 * (lo + hi)
    f_mid, se_mid, r2 = F(a_star)
    diag.update(
        {
            "a_star": f64(a_star),
            "F_at_root": f64(f_mid),
            "F_se": f64(se_mid),
            "control_r2": _finite_or_none(r2),
            "F_lo_at_zero": f64(f_lo),
            "B_used": i64(s2s.shape[0]),
        }
    )
    return a_star, diag


def fwer_case(n, sig, a_star, D, seed, pool=None):
    """End-to-end falsifier: calibrate at the truth, apply the level to the
    plug-in rule, and compare with the oracle (W4's oracle_cal_fwer)."""
    n = np.asarray(n, float)
    sig = np.asarray(sig, float)
    lam = list((sig / sig[0]) ** 2)[1:]
    s0, si, lh = draw_var(n, lam, D, seed)
    c_or = oracle_c(n, lam, ALPHA)
    c_rule = batch_c(n, lh, alpha=a_star, pool=pool)
    d = cond_noreject(np.full(D, c_or), s0, si, n, lam) - cond_noreject(
        c_rule, s0, si, n, lam
    )
    try:
        f = ALPHA + float(d.mean())
        se = float(d.std(ddof=1) / np.sqrt(D))
        z = (f - ALPHA) / se if se > 0 else float("nan")
        payload = {
            "fwer": f64(f),
            "mcse": f64(se),
            "residual": f64(f - ALPHA),
            "z": f64(z),
        }
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("end-to-end falsifier aggregation failed") from exc
    return payload


def main():
    t0 = time.time()
    try:
        os.makedirs(OUT, exist_ok=True)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"cannot create {OUT}") from exc
    rows = []
    with mp.Pool(POOL_N) as pool:
        for case in CASES[:2]:
            name, n, sig = case
            for B in B_SWEEP:
                for seed in SEEDS:
                    diags = {}
                    for variant in ("plain", "anti", "vr"):
                        a, d = calibrate(n, sig, ALPHA, B, seed, variant, pool=pool)
                        diags[variant] = d
                    rows.append(
                        {
                            "case": name,
                            "n": list(n),
                            "sig": list(sig),
                            "B": i64(B),
                            "seed": i64(seed),
                            "diag": diags,
                        }
                    )
                    print(
                        f"{name[:22]:<22s} B={B:<5d} s={seed:<3d} "
                        f"plain {diags['plain']['a_star']:.5f} | "
                        f"anti {diags['anti']['a_star']:.5f} | "
                        f"vr {diags['vr']['a_star']:.5f} "
                        f"(R2 {diags['vr']['control_r2']:.3f})",
                        flush=True,
                    )
    try:
        secs = f64(time.time() - t0)
        out_path = os.path.join(OUT, "variance.json")
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("summary plumbing failed") from exc
    summary = {
        "B_sweep": [i64(b) for b in B_SWEEP],
        "seeds": [i64(s) for s in SEEDS],
        "rows": rows,
        "seconds": secs,
    }
    # sd of a* across seeds per (case, B, variant)
    sds = []
    for case in CASES[:2]:
        name, _n, _sig = case
        for B in B_SWEEP:
            sub = [r for r in rows if r["case"] == name and r["B"] == B]
            for v in ("plain", "anti", "vr"):
                vals = [r["diag"][v]["a_star"] for r in sub]
                try:
                    sd = float(np.std(vals, ddof=1))
                except Exception:  # noqa: BLE001
                    sd = float("nan")
                sds.append(
                    {
                        "case": name,
                        "B": i64(B),
                        "variant": v,
                        "sd_a_star": f64(sd),
                        "n_seeds": i64(len(vals)),
                    }
                )
    summary["sd_across_seeds"] = sds
    # Paired noise ratios: with 4 seeds the across-seed sd is too noisy to resolve
    # a factor-2 gain, but the within-run F_se is a direct paired measure of the
    # estimator noise at the same draws and B.
    ratios = []
    for case in CASES[:2]:
        name, _n2, _s2 = case
        for B in B_SWEEP:
            sub = [r for r in rows if r["case"] == name and r["B"] == B]
            try:
                ratio = float(
                    np.median(
                        [
                            r["diag"]["plain"]["F_se"] / r["diag"]["vr"]["F_se"]
                            for r in sub
                        ]
                    )
                )
                ratio_anti = float(
                    np.median(
                        [
                            r["diag"]["plain"]["F_se"] / r["diag"]["anti"]["F_se"]
                            for r in sub
                        ]
                    )
                )
                r2 = float(np.median([r["diag"]["vr"]["control_r2"] for r in sub]))
            except Exception:  # noqa: BLE001
                ratio = ratio_anti = r2 = float("nan")
            ratios.append(
                {
                    "case": name,
                    "B": i64(B),
                    "median_se_ratio_vr": f64(ratio),
                    "median_se_ratio_anti": f64(ratio_anti),
                    "median_r2_vr": f64(r2),
                    "n_seeds": i64(len(sub)),
                }
            )
    summary["paired_noise_ratios"] = ratios
    dump_json(summary, out_path)
    for s in sds:
        print(
            f"sd(a*) {s['case'][:22]:<22s} B={s['B']:<5d} {s['variant']:<6s} "
            f"{s['sd_a_star']:.6f}",
            flush=True,
        )
    try:
        print(f"elapsed {time.time() - t0:.0f}s")
    except Exception:  # noqa: BLE001
        print("elapsed: unknown")


if __name__ == "__main__":
    main()
