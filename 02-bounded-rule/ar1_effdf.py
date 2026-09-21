"""S2.2: effective-degrees-of-freedom plug-in for weak autocorrelation.

The exact law assumes i.i.d. normal data, so under AR(1) its chi-square
denominators are wrong.  The cheapest desk-ready rescue is to keep the law but
replace the degrees of freedom by an *effective* value read off the data:

  (a) sample-size route   nu_eff = (n - 1) * mean_i (1 - phi_i)/(1 + phi_i)
  (b) Satterthwaite route nu_eff = (tr(AV))^2 / tr((AV)^2),  A = (I - 11'/n)/(n-1),
      V the AR(1) correlation matrix at phi_hat (exact for Var(s^2) under AR(1):
      Var(y'Ay) = 2 sigma^4 tr((AV)^2))

Both are plugged into the CRITICAL VALUE only: the observed t-statistics keep the
real sample sizes, and the critical value comes from `core.batch_c(n_eff, lam_hat)`.
Every rule sees the same replicates, so differences resolve far better than
absolute FWERs.

Outputs: deposits/ar2/effdf.json
"""

import contextlib
import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, batch_c  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402
from w6_model import ar1_batch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "ar2")
PHIS = (0.0, 0.1, 0.2, 0.3, 0.4)
NS = (36, 60)
LAM = (1.0, 4.0)
NU_MIN = 1.0
try:
    NAN = float("nan")
except Exception:  # noqa: BLE001
    NAN = 0.0


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


R = _env_int("W26_EFF_R", 2000)
POOL_N = _env_int("W26_POOL", 3)
BIN_W = _f(os.environ.get("W26_EFF_BIN") or 0.25, 0.25)  # nu bin width for grouping


def _binned_c(n_i, nu, lh, pool):
    """Critical value per row, grouping rows into bins of width BIN_W in nu.

    `batch_c` takes one sample-size vector for all rows, so per-row nu must be
    grouped; the rounding error is returned alongside the values and reported.
    """
    nb = np.round(nu / BIN_W).astype(int)
    out = np.empty(len(nu))
    worst = 0.0
    for b in np.unique(nb):
        m = nb == b
        nu_b = _f(np.clip(_i(b) * BIN_W, NU_MIN, n_i - 1.0), NU_MIN)
        out[m] = batch_c(np.full(3, nu_b + 1.0), lh[m], alpha=ALPHA, pool=pool)
        dev = np.abs(nu[m] - nu_b)
        with contextlib.suppress(Exception):
            worst = max(worst, _f(np.max(dev), 0.0))
    return out, worst


def phi_hat(x):
    """Lag-1 autocorrelation of a demeaned series (method of moments)."""
    try:
        xc = np.asarray(x, float) - float(np.mean(x))
        den = float(np.dot(xc, xc))
        if den <= 0.0:
            return 0.0
        return float(np.dot(xc[1:], xc[:-1]) / den)
    except Exception:  # noqa: BLE001
        return 0.0


def nu_eff_satterthwaite(n_i, phi):
    """(tr(AV))^2 / tr((AV)^2) for A = (I - 11'/n)/(n-1), V = AR(1) at phi."""
    try:
        n_i = int(n_i)
        V = np.zeros((n_i, n_i))
        for k in range(n_i):
            V[k, k:] = phi ** np.arange(n_i - k)
            V[k:, k] = phi ** np.arange(n_i - k)
        A = (np.eye(n_i) - np.ones((n_i, n_i)) / n_i) / (n_i - 1.0)
        AV = A @ V
        tr = float(np.trace(AV))
        tr2 = float(np.trace(AV @ AV))
        if tr2 <= 0.0:
            return float(n_i - 1)
        return float((tr * tr) / tr2)
    except Exception:  # noqa: BLE001
        return float("nan")


def _draws(n_i, phi):
    try:
        rng = np.random.default_rng(9000 + int(1000 * phi) + int(n_i))
        sds = [1.0] + [float(np.sqrt(v)) for v in LAM]
        return [ar1_batch(rng, phi, sds[i], int(n_i), R) for i in range(3)]
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"AR(1) draw failed for n={n_i} phi={phi}") from exc


def cell(n_i, phi, pool):
    x = _draws(n_i, phi)
    t_max = np.empty(R)
    lh = np.empty((R, 2))
    nu_a = np.empty(R)
    nu_b = np.empty(R)
    try:
        for r in range(R):
            one = np.stack([x[i][r] for i in range(3)])
            xb = one.mean(axis=1)
            s2 = one.var(axis=1, ddof=1)
            t = (xb[1:] - xb[0]) / np.sqrt(s2[1:] / n_i + s2[0] / n_i)
            t_max[r] = float(np.max(t))
            lh[r] = s2[1:] / s2[0]
            ph = [min(max(phi_hat(one[i]), -0.99), 0.99) for i in range(3)]
            nu_a[r] = (n_i - 1.0) * float(
                np.mean([(1.0 - ph[i]) / (1.0 + ph[i]) for i in range(3)])
            )
            nu_b[r] = float(
                np.mean([nu_eff_satterthwaite(n_i, ph[i]) for i in range(3)])
            )
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"replicate loop failed for n={n_i} phi={phi}") from exc
    nu_a = np.clip(nu_a, NU_MIN, n_i - 1.0)
    nu_b = np.clip(np.nan_to_num(nu_b, nan=_f(n_i - 1.0)), 1.0, _f(n_i - 1.0))
    n_vec = np.array([n_i] * 3, float)
    try:
        c_naive = batch_c(n_vec, lh, alpha=ALPHA, pool=pool)
        c_a, err_a = _binned_c(n_i, nu_a, lh, pool)
        c_b, err_b = _binned_c(n_i, nu_b, lh, pool)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"critical values failed for n={n_i} phi={phi}") from exc
    out = {
        "n": i64(n_i),
        "phi": f64(phi),
        "R": i64(R),
        "alpha": ALPHA,
        "lam": list(LAM),
        "nu_min": f64(NU_MIN),
        "bin_w": f64(BIN_W),
        "max_nu_bin_err_sample": f64(err_a),
        "max_nu_bin_err_satt": f64(err_b),
        "mean_nu_eff_sample": f64(_f(np.mean(nu_a))),
        "mean_nu_eff_satt": f64(_f(np.mean(nu_b))),
    }
    for name, cv in (("naive", c_naive), ("effdf_sample", c_a), ("effdf_satt", c_b)):
        try:
            p = float(np.mean(t_max > cv))
            se = float(np.sqrt(max(1e-300, p * (1 - p)) / R))
            mc = float(np.mean(cv))
        except Exception:  # noqa: BLE001
            p = se = mc = NAN
        out[f"fwer_{name}"] = f64(p)
        out[f"mcse_{name}"] = f64(se)
        out[f"rel_err_pct_{name}"] = f64(100 * (p - ALPHA) / ALPHA)
        out[f"within_5pct_{name}"] = bool(abs(p - ALPHA) <= 0.05 * ALPHA)
        out[f"mean_c_{name}"] = f64(mc)
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
                    f"n={n_i} phi={phi:.1f}: naive {r['fwer_naive']:.4f} "
                    f"({r['rel_err_pct_naive']:+.1f}%) | effdf_sample "
                    f"{r['fwer_effdf_sample']:.4f} "
                    f"({r['rel_err_pct_effdf_sample']:+.1f}%) | effdf_satt "
                    f"{r['fwer_effdf_satt']:.4f} "
                    f"({r['rel_err_pct_effdf_satt']:+.1f}%)",
                    flush=True,
                )
    win = [r for r in rows if r["within_5pct_effdf_satt"]]
    boundary = []
    for r in rows:
        for name in ("naive", "effdf_sample", "effdf_satt"):
            if not r[f"within_5pct_{name}"]:
                boundary.append(
                    {
                        "n": r["n"],
                        "phi": r["phi"],
                        "rule": name,
                        "fwer": r[f"fwer_{name}"],
                    }
                )
    try:
        phi_max = f64(max(r["phi"] for r in win)) if win else f64(float("nan"))
    except Exception:  # noqa: BLE001
        phi_max = f64(float("nan"))
    summary = {
        "rows": rows,
        "alpha": ALPHA,
        "R": i64(R),
        "effdf_satt_within_5pct": i64(len(win)),
        "out_of_5pct": boundary,
        "phi_max_within_5pct_satt": phi_max,
        "seconds": f64(time.time() - t0),
        "note": (
            "nu_eff replaces the degrees of freedom used to BUILD the "
            "critical value only; the observed t-statistics keep the real "
            "sample sizes, and every rule sees the same replicates."
        ),
    }
    dump_json(summary, os.path.join(OUT, "effdf.json"))
    print(
        f"effdf_satt within 5% in {len(win)}/{len(rows)} cells | max phi "
        f"{phi_max:.2f} | {time.time() - t0:.0f}s"
    )


if __name__ == "__main__":
    main()
