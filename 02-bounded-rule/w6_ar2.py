"""Track D1: AR(1) candidate dependence -- block-bootstrap rule vs plug-in / envelope.

W6 measured the failure: with AR(1) series the naive chi-square plug-in rule
rejects at FWER 0.297 (n=36, phi=0.6) and even the HAC-refitted rule is liberal
(0.164), so a serial-correlation correction that *keeps* the chi-square law cannot
work.  This module evaluates three data-computable rules against the true
long-run variances on the same simulated stream:

  plugin    : exact chi-square critical value at the estimated ratios (wrong law)
  envelope  : the Phase A tight-region envelope (same wrong law, conservative)
  bootstrap : moving-block bootstrap of the CENTRED series, critical value = the
              (1-alpha) quantile of the bootstrap max-t statistic -- a rule that
              does not assume a chi-square denominator law

Acceptance is stated against the plan: a rule within 5% relative for phi <= 0.3
at n = 60, or the boundary where it fails.
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, batch_c  # noqa: E402
from env_rule import c_env as c_env_rule  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402
from w6_model import ar1_batch, t_stats  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "ar2")
PHIS = (0.0, 0.2, 0.4, 0.6)
NS = (36, 60)
LAM = (1.0, 4.0)


def _env_floats(name, default):
    try:
        raw = os.environ.get(name) or default
        return tuple(float(x) for x in str(raw).split(","))
    except Exception:  # noqa: BLE001
        return tuple(float(x) for x in str(default).split(","))


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


R = _env_int("W26_AR_R", 200)
R_BOOT = _env_int("W26_AR_RB", 200)
POOL_N = _env_int("W26_POOL", 4)
PHIS = _env_floats("W26_AR_PHIS", "0,0.2,0.4,0.6")
NS = _env_tuple("W26_AR_NS", "36,60")


def block_indices(rng, n, L, m):
    """(m, n) moving-block bootstrap indices for series of length n."""
    try:
        nblk = int(np.ceil(n / L))
        starts = rng.integers(0, max(1, n - L + 1), size=(m, nblk))
        offs = np.arange(L)[None, None, :]
        idx = (starts[:, :, None] + offs).reshape(m, nblk * L)[:, :n]
        return np.clip(idx, 0, n - 1)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("block index construction failed") from exc


def bootstrap_c(series_centred, n_vec, alpha, m, L, rng):
    """Block-bootstrap critical value for one replicate.

    `series_centred` must already be centred by its own group mean (that is how
    the null of equal means is imposed).  The resampled series is NOT re-centred:
    the resampling-induced variation of the resampled group means IS the bootstrap
    analogue of the sampling variability.  Re-centring the resample forces every
    bootstrap group mean to zero, which drives the critical value to ~0 and makes
    the rule reject almost always (measured FWER 0.66-0.75 before this fix).
    """
    nb = series_centred.shape[0]
    idx = block_indices(rng, series_centred.shape[1], L, m)
    xs = np.stack([series_centred[i][idx] for i in range(nb)])  # (nb, m, n)
    xb = xs.mean(axis=2)  # (nb, m)
    s2 = xs.var(axis=2, ddof=1)  # (nb, m)
    t = t_stats(xb.T, s2.T, n_vec)  # (m, K)
    try:
        vals = np.sort(t.max(axis=1))
        k = int(np.ceil((1.0 - alpha) * vals.size)) - 1
        return float(vals[max(0, min(k, vals.size - 1))])
    except Exception:  # noqa: BLE001
        return float("nan")


def cell(n_i, phi, pool):
    n_vec = np.array([n_i] * 3, float)
    try:
        rng = np.random.default_rng(9000 + int(1000 * phi) + n_i)
        L = max(2, int(round(n_i ** (1.0 / 3.0))))
        sds = [1.0] + [float(np.sqrt(v)) for v in LAM]
        series = [ar1_batch(rng, phi, sds[i], n_i, R) for i in range(3)]
        xb = np.empty((R, 3))
        s2n = np.empty((R, 3))
        boot_c = np.empty(R)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"setup failed for n={n_i} phi={phi}") from exc
    for r in range(R):
        one = np.stack([series[i][r] for i in range(3)])
        xb[r] = one.mean(axis=1)
        s2n[r] = one.var(axis=1, ddof=1)
        boot_c[r] = bootstrap_c(
            one - one.mean(axis=1, keepdims=True), n_vec, ALPHA, R_BOOT, L, rng
        )
    lh = s2n[:, 1:] / s2n[:, 0][:, None]
    t_obs = t_stats(xb, s2n, n_vec)
    tmax = t_obs.max(axis=1)
    c_plug = batch_c(n_vec, lh, alpha=ALPHA, pool=pool)
    c_env = c_env_rule(
        n_vec,
        s2n[:, 0],
        s2n[:, 1:],
        lh,
        alpha=ALPHA,
        gamma=0.0005,
        pool=pool,
        tight=True,
    )
    rules = {
        "plugin": c_plug,
        "envelope": c_env,
        "bootstrap": boot_c,
    }
    try:
        mc_p = float(np.mean(c_plug))
        mc_e = float(np.mean(c_env))
        mc_b = float(np.mean(boot_c))
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("mean critical value aggregation failed") from exc
    out = {
        "n": i64(n_i),
        "phi": f64(phi),
        "R": i64(R),
        "R_boot": i64(R_BOOT),
        "block_len": i64(L),
        "alpha": ALPHA,
        "lam": list(LAM),
        "c_mean_plugin": f64(mc_p),
        "c_mean_envelope": f64(mc_e),
        "c_mean_bootstrap": f64(mc_b),
    }
    for name, cv in rules.items():
        try:
            rej = tmax > cv
            p = float(np.mean(rej))
            se = float(np.sqrt(max(1e-300, p * (1 - p)) / R))
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(f"rejection aggregation failed for {name}") from exc
        out[f"fwer_{name}"] = f64(p)
        out[f"mcse_{name}"] = f64(se)
        out[f"rel_err_pct_{name}"] = f64(100 * (p - ALPHA) / ALPHA)
        out[f"within_5pct_{name}"] = bool(abs(p - ALPHA) <= 0.05 * ALPHA)
        out[f"within_10pct_{name}"] = bool(abs(p - ALPHA) <= 0.10 * ALPHA)
    return out


def main():
    t0 = time.time()
    rows = []
    with mp.Pool(POOL_N) as pool:
        for n_i in NS:
            for phi in PHIS:
                r = cell(n_i, phi, pool)
                rows.append(r)
                print(
                    f"n={n_i} phi={phi:.1f}: plugin FWER {r['fwer_plugin']:.3f} "
                    f"({r['rel_err_pct_plugin']:+.1f}%) | env "
                    f"{r['fwer_envelope']:.3f} ({r['rel_err_pct_envelope']:+.1f}%) | "
                    f"boot {r['fwer_bootstrap']:.3f} "
                    f"({r['rel_err_pct_bootstrap']:+.1f}%)",
                    flush=True,
                )
    boundary = []
    for r in rows:
        for name in ("plugin", "envelope", "bootstrap"):
            if not r[f"within_5pct_{name}"]:
                boundary.append(
                    {
                        "n": r["n"],
                        "phi": r["phi"],
                        "rule": name,
                        "fwer": r[f"fwer_{name}"],
                        "rel_err_pct": r[f"rel_err_pct_{name}"],
                    }
                )
    try:
        phi_ok = [r["phi"] for r in rows if r["within_5pct_bootstrap"]]
        phi_max = f64(max(phi_ok)) if phi_ok else None
    except Exception:  # noqa: BLE001
        phi_max = None
    out = {
        "rows": rows,
        "alpha": ALPHA,
        "phi_list": list(PHIS),
        "n_list": list(NS),
        "out_of_5pct": boundary,
        "bootstrap_boundary_phi_max": phi_max,
        "note": (
            "block length = round(n^(1/3)); the bootstrap rule resamples "
            "the centred series, so it does not assume a chi-square "
            "denominator law"
        ),
        "seconds": f64(time.time() - t0),
    }
    dump_json(out, os.path.join(OUT, "ar1.json"))
    print(
        f"out-of-5% cells: {len(boundary)} | bootstrap max phi within 5%: "
        f"{out['bootstrap_boundary_phi_max']}"
    )


if __name__ == "__main__":
    main()
