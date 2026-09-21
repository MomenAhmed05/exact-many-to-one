"""W6: the model layer - autocorrelation and correlated candidates.

Part 1 (autocorrelation).  Under AR(1) errors the sample mean's variance is the
long-run variance, not sigma^2/n, so the Welch denominator understates the noise
and the exact-law plug-in rule - exact under the model it assumes - is applied
to data that violate it.  Measured directly by simulating data sets and counting
rejections of the plug-in rule, with and without a Newey-West correction.

Part 2 (correlated candidates).  Candidates sharing a factor break the
conditional independence the reduction rests on.  Measured directly; plus a
check of the algebraic claim that conditional on the factor path the candidate
sample variance is noncentral chi-square (so the product formula fails).

Outputs: deposits/w6/ar1.json, deposits/w6/factor.json
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import batch_c  # noqa: E402
from io_util import dump_json, f64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w6")


def _int(x):
    try:
        return int(x)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("not an int") from exc


def nw_lag(n):
    try:
        return max(1, int(np.floor(4.0 * (n / 100.0) ** (2.0 / 9.0))))
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("bad n") from exc


def ar1_batch(rng, phi, sig, n, B):
    """B AR(1) series of length n with marginal sd sig."""
    eta = rng.standard_normal((B, n)) * sig * np.sqrt(max(1e-12, 1 - phi**2))
    x = np.empty((B, n))
    x[:, 0] = eta[:, 0] / np.sqrt(max(1e-12, 1 - phi**2))
    for t in range(1, n):
        x[:, t] = phi * x[:, t - 1] + eta[:, t]
    return x


def long_run_var(x, lag):
    """Newey-West (Bartlett) long-run variance of one series' mean, x (B,n)."""
    B, n = x.shape
    xc = x - x.mean(axis=1, keepdims=True)
    g0 = (xc**2).mean(axis=1)
    tot = g0.copy()
    for k in range(1, lag + 1):
        w = 1.0 - k / (lag + 1.0)
        tot += 2.0 * w * (xc[:, k:] * xc[:, :-k]).mean(axis=1)
    return tot


def t_stats(xb, s2, n):
    den = np.sqrt(s2[:, 1:] / n[1:] + s2[:, 0][:, None] / n[0])
    return (xb[:, 1:] - xb[:, 0][:, None]) / den


def run_ar1(
    phi_list=(0.0, 0.1, 0.2, 0.3, 0.4, 0.6),
    n_list=(36, 60),
    lam=(1.0, 4.0),
    alpha=0.05,
    B=2000,
    seed=777,
    workers=10,
    pool=None,
):
    rows = []
    for n_i in n_list:
        n = np.array([n_i] * 3, float)
        for phi in phi_list:
            rng = np.random.default_rng(seed + _int(1000 * phi) + n_i)
            xb = np.empty((B, 3))
            s2n = np.empty((B, 3))
            s2h = np.empty((B, 3))
            for i in range(3):
                x = ar1_batch(rng, phi, np.sqrt(lam[i - 1]) if i else 1.0, n_i, B)
                xb[:, i] = x.mean(axis=1)
                s2n[:, i] = x.var(axis=1, ddof=1)
                s2h[:, i] = long_run_var(x, nw_lag(n_i))
            lh_n = s2n[:, 1:] / s2n[:, 0][:, None]
            lh_h = s2h[:, 1:] / s2h[:, 0][:, None]
            t_n = t_stats(xb, s2n, n)
            t_h = (xb[:, 1:] - xb[:, 0][:, None]) / np.sqrt(
                s2h[:, 1:] / n_i + s2h[:, 0][:, None] / n_i
            )
            c_n = batch_c(n, lh_n, alpha=alpha, workers=workers, pool=pool)
            c_h = batch_c(n, lh_h, alpha=alpha, workers=workers, pool=pool)
            rej_n = f64(np.mean(t_n.max(axis=1) > c_n))
            rej_h = f64(np.mean(t_h.max(axis=1) > c_h))
            rows.append(
                {
                    "n": n_i,
                    "phi": phi,
                    "alpha": alpha,
                    "draws": B,
                    "fwer_naive": rej_n,
                    "mcse_naive": f64(np.sqrt(rej_n * (1 - rej_n) / B)),
                    "fwer_hac": rej_h,
                    "mcse_hac": f64(np.sqrt(rej_h * (1 - rej_h) / B)),
                    "mean_c_naive": f64(c_n.mean()),
                    "mean_c_hac": f64(c_h.mean()),
                }
            )
            print(
                f"n={n_i} phi={phi:.1f} naive {rej_n:.4f} "
                f"hac {rej_h:.4f} (nominal {alpha})",
                flush=True,
            )
    return rows


def run_factor(
    rho_list=(0.0, 0.25, 0.5, 0.75),
    K=2,
    n_i=36,
    alpha=0.025,
    B=2000,
    seed=555,
    workers=10,
    mode="same",
    pool=None,
):
    """One common factor with loadings set by `mode`.

    same       b = (+sqrt(rho), +sqrt(rho), ...)  benchmark moves with the
               candidates; differences less variable than the law assumes
    cand_only  b = (0, +sqrt(rho), ...)           candidates correlated with
               each other but not the benchmark
    sign_flip  b = (+sqrt(rho), -sqrt(rho), +sqrt(rho), ...)  opposite-sign
               exposures; differences MORE variable than the law assumes
    """
    rows = []
    n = np.array([n_i] * (K + 1), float)
    for rho in rho_list:
        rng = np.random.default_rng(seed + _int(100 * rho))
        b = np.full(K + 1, np.sqrt(rho))
        if mode == "cand_only":
            b[0] = 0.0
        elif mode == "sign_flip":
            b[0] = np.sqrt(rho)
            for i in range(1, K + 1):
                b[i] = (-1.0) ** i * np.sqrt(rho)
        se = np.sqrt(max(1e-12, 1 - rho))
        fac = rng.standard_normal((B, n_i))
        xb = np.empty((B, K + 1))
        s2 = np.empty((B, K + 1))
        for i in range(K + 1):
            x = b[i] * fac + se * rng.standard_normal((B, n_i))
            xb[:, i] = x.mean(axis=1)
            s2[:, i] = x.var(axis=1, ddof=1)
        lh = s2[:, 1:] / s2[:, 0][:, None]
        t = t_stats(xb, s2, n)
        c = batch_c(n, lh, alpha=alpha, workers=workers, pool=pool)
        rej = f64(np.mean(t.max(axis=1) > c))
        rows.append(
            {
                "K": K,
                "n": n_i,
                "rho": rho,
                "alpha": alpha,
                "draws": B,
                "mode": mode,
                "loadings": [f64(x) for x in b],
                "fwer": rej,
                "mcse": f64(np.sqrt(rej * (1 - rej) / B)),
                "mean_c": f64(c.mean()),
            }
        )
        print(f"mode={mode} rho={rho:.2f} fwer={rej:.4f} (nominal {alpha})", flush=True)
    return rows


def factor_conditional_check(K=2, n_i=36, rho=0.5, reps=400, seed=99):
    """Conditional on the factor path: is s_i^2 noncentral chi-square?

    X_i = b F + eps_i.  Given F, (n-1) s_i^2 / sig_eps^2 is noncentral chi-square
    with df n-1 and noncentrality b^2 sum_t (F_t - Fbar)^2 / sig_eps^2.
    Checked by the sample mean and variance of the ratio against the theoretical
    noncentral chi-square moments, and by comparing the empirical CDF at the
    theoretical quantiles.
    """
    from scipy.stats import ncx2

    rng = np.random.default_rng(seed)
    b = np.sqrt(rho)
    se = np.sqrt(max(1e-12, 1 - rho))
    ratios, nonc = [], []
    for _ in range(reps):
        f = rng.standard_normal(n_i)
        x = b * f + se * rng.standard_normal(n_i)
        s2 = x.var(ddof=1)
        lamnc = b**2 * ((f - f.mean()) ** 2).sum() / se**2
        ratios.append((n_i - 1) * s2 / se**2)
        nonc.append(lamnc)
    ratios = np.array(ratios)
    nonc = np.array(nonc)
    mean_theory = f64(np.mean((n_i - 1) + nonc))
    # Var(Y) = E[2(k+2lam)] + Var(lam) = 2k + 4E[lam] + Var(lam)
    var_theory = f64(2 * (n_i - 1) + 4 * np.mean(nonc) + np.var(nonc))
    # CDF comparison at the median quantile
    q = np.array([0.25, 0.5, 0.75])
    emp = np.quantile(ratios, q)
    theo = np.array([f64(ncx2.ppf(qq, n_i - 1, nonc.mean())) for qq in q])
    return {
        "reps": reps,
        "n": n_i,
        "rho": rho,
        "mean_ratio_empirical": f64(ratios.mean()),
        "mean_ratio_theory": f64(mean_theory),
        "var_ratio_empirical": f64(ratios.var(ddof=1)),
        "var_ratio_theory": f64(var_theory),
        "quantiles_empirical": [f64(x) for x in emp],
        "quantiles_theory_at_meannonc_jensen_biased": [f64(x) for x in theo],
        "noncentrality_mean": f64(nonc.mean()),
        "noncentrality_sd": f64(nonc.std(ddof=1)),
    }


def main():
    import multiprocessing as mp

    with mp.Pool(3) as pool:
        _main(pool)


def _main(pool):
    t0 = time.time()
    ar = run_ar1(pool=pool)
    dump_json(
        {
            "rows": ar,
            "note": "AR(1) errors; exact-law plug-in rule "
            "applied to ordinary and HAC denominators",
        },
        os.path.join(OUTDIR, "ar1.json"),
    )
    fac_same = run_factor(mode="same", pool=pool)
    fac_cand = run_factor(mode="cand_only", pool=pool)
    fac_flip = run_factor(mode="sign_flip", pool=pool)
    chk = factor_conditional_check()
    dump_json(
        {
            "same_sign": fac_same,
            "candidates_only": fac_cand,
            "sign_flip": fac_flip,
            "conditional_check": chk,
        },
        os.path.join(OUTDIR, "factor.json"),
    )
    print(f"done in {time.time() - t0:.0f}s")
    print("wrote", OUTDIR)


if __name__ == "__main__":
    main()
