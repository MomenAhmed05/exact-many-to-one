"""S2.3: the exact alternative law, so the envelope's price lands in power units.

The conditional no-rejection probability at ratios `lam` and critical value `c` is

    P(no reject | s^2) = E_U[ prod_i Phi(s_i (U + c d_i)) ],  s_i = sqrt(n_i/lam_i),

with `d_i` the drawn standard error of contrast i (`sqrt(s_i2/n_i + s_0^2/n_0)`).

Alternative.  Let `nu_i = (mu_i - mu_0) / SE_i^pop` be the population
signal-to-noise of contrast i, with `SE_i^pop = sqrt(lam_i/n_i + 1/n_0)` (the
harness normalises sigma_0 = 1, so sigma_i^2 = lam_i).  Under that shift the
statistic `T_i` gains mean `nu_i`, and the no-rejection probability keeps the same
form with the argument decreased by `nu_i s_i d_i`:

    arg_i = s_i (U + c d_i - nu_i d_i)

(the sign is forced: a larger true effect must lower P(T_i <= c)).  This is the
`delta`-inside-Phi extension, and it is validated here against a brute-force
simulation that generates shifted data and computes the t-statistics directly.

Outputs: deposits/power/*.json
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np
from scipy.special import ndtr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, draw_var, oracle_c  # noqa: E402
from env_rule import c_env, lam_region  # noqa: E402
from h12_quad import normal_nodes  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "power")
FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
DFS = (3, 5, 7, 11, 19, 35)
GAMMA = 0.0005
M_NODES = 400
try:
    NAN = float("nan")
except Exception:  # noqa: BLE001
    NAN = 0.0


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


POOL_N = _env_int("W26_POOL", 3)


def se_pop(n, lam):
    """Population standard error of each contrast (sigma_0 = 1)."""
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    return np.sqrt(lam / n[1:] + 1.0 / n[0])


def norej_shifted(c, d_row, lam, nu, z, w, n):
    """P(no reject | s^2) under the alternative nu_i (argument shift inside Phi)."""
    n = np.asarray(n, float)
    d_row = np.asarray(d_row, float).reshape(-1)
    lam = np.asarray(lam, float)
    nu = np.atleast_1d(np.asarray(nu, float)).reshape(-1)
    if nu.size == 1:
        nu = np.repeat(nu, len(lam))
    try:
        u = z / np.sqrt(n[0])
        acc = np.ones(z.size)
        sep = np.sqrt(lam / n[1:] + 1.0 / n[0])  # population SE of each contrast
        for i in range(len(lam)):
            s = np.sqrt(n[i + 1] / lam[i])
            # exact algebra: arg = s (u + c d_i) - s * nu_i * SE_i^pop
            arg = s * (u + c * d_row[i]) - s * nu[i] * sep[i]
            acc = acc * ndtr(arg)
        return float(np.dot(acc, w))
    except Exception:  # noqa: BLE001
        return float("nan")


def power_curve(n, lam, c, nus, s0, si, m=M_NODES):
    """Mean rejection over draws for each nu in `nus` (paired across nu).

    Uses the harness's own Gauss nodes (`h12_quad.normal_nodes`), i.e. the same
    quadrature `core.cond_noreject` uses, so the power measurement carries the
    harness's accuracy rather than a hand-rolled endpoint rule.
    """
    z, w = normal_nodes(40, 10)
    out = []
    try:
        for nu in nus:
            vals = []
            for j in range(len(s0)):
                d_row = np.sqrt(si[j] / n[1:] + s0[j] / n[0])
                vals.append(1.0 - norej_shifted(c, d_row, lam, nu, z, w, n))
            out.append(float(np.mean(vals)))
        return out
    except Exception:  # noqa: BLE001
        return [float("nan") for _ in nus]


def validate_conditional(n, lam, c, nu, s0, si, reps=2000, seed=5):
    """Paired validation: conditional rejection given the SAME s^2 draws.

    Given (s_0, s_i) the sample means are independent normals with variances
    sigma^2/n, so the conditional rejection can be simulated directly at the same
    s^2 the quadrature conditions on.  Returns (paired mean difference, MCSE).
    """
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    K = len(lam)
    sep = se_pop(n, lam)
    mu = np.asarray(nu, float) * sep
    sd = np.sqrt(lam)
    rng = np.random.default_rng(seed)
    z, w = normal_nodes(40, 10)
    diffs = []
    try:
        for j in range(len(s0)):
            d_row = np.sqrt(si[j] / n[1:] + s0[j] / n[0])
            quad = 1.0 - norej_shifted(c, d_row, lam, nu, z, w, n)
            z0 = rng.standard_normal(reps)
            zi = rng.standard_normal((reps, K))
            xb0 = z0 / np.sqrt(n[0])
            xbi = zi * sd[None, :] / np.sqrt(n[1:])[None, :] + mu[None, :]
            t = (xbi - xb0[:, None]) / d_row[None, :]
            sim = float(np.mean(np.max(t, axis=1) > c))
            diffs.append(sim - quad)
        arr = np.asarray(diffs, float)
        return float(np.mean(arr)), float(np.std(arr, ddof=1) / np.sqrt(len(arr)))
    except Exception:  # noqa: BLE001
        return float("nan"), float("nan")


def simulate_power(n, lam, c, nu, reps=4000, seed=7):
    """Brute force: generate shifted normals, form t-statistics, count rejections."""
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    rng = np.random.default_rng(seed)
    K = len(lam)
    sep = se_pop(n, lam)
    delta = np.asarray(nu, float) * sep  # mean shift per contrast
    sd = np.sqrt(lam)  # sigma_i (sigma_0 = 1)
    rej = 0
    try:
        for _ in range(reps):
            x0 = rng.normal(0.0, 1.0, int(n[0]))
            xi = rng.normal(0.0, 1.0, (K, int(n[1])))
            xi = xi * sd[:, None] + delta[:, None]
            xb0 = x0.mean()
            xbi = xi.mean(axis=1)
            s0v = x0.var(ddof=1)
            siv = xi.var(axis=1, ddof=1)
            t = (xbi - xb0) / np.sqrt(siv / n[1:] + s0v / n[0])
            rej += int(np.max(t) > c)
        return rej / reps
    except Exception:  # noqa: BLE001
        return float("nan")


def cell(fam, df, draws, m=M_NODES, pool=None):
    lam = np.asarray(FAMS[fam], float)
    n = np.array([df + 1, 13, 13], float)
    s0, si, lh = draw_var(n, lam, draws, 4242 + 97 * df)
    lo, hi, _ = lam_region(n, lh, GAMMA, scale=None)
    c_or = oracle_c(n, lam, ALPHA)
    c_cand = c_env(n, s0, si, lh, gamma=GAMMA, alpha=ALPHA, pool=pool)
    nus = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    try:
        c_or_f = float(c_or)
        c_env_f = float(np.mean(c_cand))
        pow_or = power_curve(n, lam, c_or_f, nus, s0, si, m=m)
        pow_env = power_curve(n, lam, c_env_f, nus, s0, si, m=m)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("power curves failed") from exc
    # brute-force validation at two points (small reps for cost)
    val = []
    for nu in (1.5, 3.0):
        try:
            sim = simulate_power(n, lam, c_or_f, np.array([nu]), reps=600, seed=11 + df)
            quad = float(pow_or[int(np.where(nus == nu)[0][0])])
            pair_diff, pair_se = validate_conditional(
                n, lam, c_or_f, np.array([nu]), s0, si, reps=2000, seed=5 + df
            )
        except Exception:  # noqa: BLE001
            sim, quad, pair_diff, pair_se = (
                float("nan"),
                float("nan"),
                float("nan"),
                float("nan"),
            )
        val.append(
            {
                "nu": f64(nu),
                "sim": f64(sim),
                "quad": f64(quad),
                "paired_diff": f64(pair_diff),
                "paired_se": f64(pair_se),
            }
        )
    try:
        d80 = float(np.interp(0.8, pow_or, nus))
        power_env_at_80 = f64(float(np.interp(d80, nus, pow_env)))
        power_loss_80 = f64(0.8 - float(power_env_at_80))
        bf_diff = f64(float(max(abs(v["sim"] - v["quad"]) for v in val)))
        paired_worst = f64(float(max(abs(v["paired_diff"]) for v in val)))
    except Exception:  # noqa: BLE001
        d80 = NAN
        power_env_at_80 = f64(float("nan"))
        power_loss_80 = f64(float("nan"))
        bf_diff = f64(float("nan"))
        paired_worst = f64(float("nan"))
    return {
        "family": fam,
        "df0": df,
        "n": list(n),
        "lam_true": list(lam),
        "draws": i64(draws),
        "alpha": ALPHA,
        "gamma": GAMMA,
        "c_oracle": f64(c_or_f),
        "c_envelope": f64(c_env_f),
        "nus": [f64(v) for v in nus],
        "power_oracle": [f64(v) for v in pow_or],
        "power_envelope": [f64(v) for v in pow_env],
        "nu_for_80pct_oracle": f64(d80),
        "power_env_at_nu80": power_env_at_80,
        "power_loss_at_nu80": power_loss_80,
        "brute_force": val,
        "brute_force_max_abs_diff": bf_diff,
        "paired_worst_abs_diff": paired_worst,
    }


def main():
    t0 = time.time()
    rows = []
    draws = _env_int("W26_POW_D", 40)
    with mp.Pool(POOL_N) as pool:
        for df in DFS:
            for fam in FAMS:
                r = cell(fam, df, draws, pool=pool)
                rows.append(r)
                print(
                    f"{fam} df{df:>2}: c_or {r['c_oracle']:.4f} c_env "
                    f"{r['c_envelope']:.4f} | nu80 {r['nu_for_80pct_oracle']:.3f} "
                    f"| power_env@nu80 {r['power_env_at_nu80']:.3f} | brute diff "
                    f"{r['brute_force_max_abs_diff']:.3f}",
                    flush=True,
                )
    try:
        max_bf = f64(max(r["brute_force_max_abs_diff"] for r in rows))
        worst_loss = f64(max(r["power_loss_at_nu80"] for r in rows))
        worst_paired = f64(
            max(abs(v["paired_diff"]) for r in rows for v in r["brute_force"])
        )
    except Exception:  # noqa: BLE001
        max_bf = f64(float("nan"))
        worst_loss = f64(float("nan"))
        worst_paired = f64(float("nan"))
    summary = {
        "rows": rows,
        "alpha": ALPHA,
        "gamma": GAMMA,
        "draws_per_cell": i64(draws),
        "max_brute_force_diff": max_bf,
        "worst_power_loss_at_nu80": worst_loss,
        "worst_paired_abs_diff": worst_paired,
        "seconds": f64(time.time() - t0),
        "note": (
            "power_oracle/power_envelope are conditional expectations over the "
            "draws at each nu; nu is the population signal-to-noise "
            "(mu_i - mu_0)/SE_i^pop, so nu80 is the effect at which the oracle "
            "reaches 80% power."
        ),
    }
    dump_json(summary, os.path.join(HERE, "deposits", "power", "power_table.json"))
    print(
        f"max brute-force diff {summary['max_brute_force_diff']:.3f} | worst power "
        f"loss at nu80 {summary['worst_power_loss_at_nu80']:.3f} | "
        f"{time.time() - t0:.0f}s"
    )


if __name__ == "__main__":
    main()
