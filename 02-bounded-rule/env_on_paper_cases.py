"""The envelope rule priced on the paper's own four cases.

Sections 7.4 and 7.5 of the paper compare the plug-in exact rule and Method 3 on
four designs (S1, S4, S5, S6). Those comparisons are at the plug-in, so they carry
the estimation error. This run adds the second arm the paper needs: the admission
region envelope, which is the rule whose familywise error rate is controlled by
construction, priced on exactly the same designs and the same nominal level so the
numbers sit beside the published ones without reinterpretation.

Per case it deposits, over the same draws and with the oracle as the paired
reference:

  oracle     c*_pop at the true ratios, and the rule that uses it
  plug-in    c*_pop at the estimated ratios
  envelope   max over the region's vertices of c*_pop, with the region built from
             the same estimates at the certified coverage 1 - gamma

plus the fraction of draws on which the true ratios fall outside the region, which
is the quantity the guarantee's additive term budgets for.

Outputs: deposits/env/paper_cases.json
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, batch_c, draw_var, oracle_c  # noqa: E402
from env_rule import c_env, fwer_of, lam_region  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "env")
GAMMA = 0.0005

# (label, n, sigma) exactly as the paper's Table 5 designs
CASES = (
    ("S1 balanced, equal sd", [10, 10, 10], [30.0, 30.0, 30.0]),
    ("S4 control smallest n and sd", [4, 13, 13], [10.0, 10.0, 50.0]),
    ("S5 trt smallest n, largest sd", [13, 13, 4], [10.0, 10.0, 50.0]),
    ("S6 control smallest n, largest sd", [4, 13, 13], [50.0, 10.0, 10.0]),
)


def _i(x, default=0):
    try:
        return int(x)
    except Exception:  # noqa: BLE001
        return int(default)


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


D = _env_int("W26_PC_D", 200)
SEED = _env_int("W26_PC_SEED", 20260918)
POOL_N = _env_int("W26_POOL", 4)


def cell(label, n_list, sig_list, pool):
    n = np.array(n_list, float)
    sig = np.array(sig_list, float)
    lam_true = (sig[1:] / sig[0]) ** 2  # ratios to the benchmark
    s0, si, lh = draw_var(n, lam_true, D, SEED + _i(7 * len(n_list)) + _i(sum(n_list)))
    lam_hat = np.atleast_2d(lh)
    c_or = oracle_c(n, lam_true, ALPHA)
    c_plug = batch_c(n, lam_hat, alpha=ALPHA, pool=pool)
    c_env_v = c_env(n, s0, si, lam_hat, alpha=ALPHA, gamma=GAMMA, pool=pool)
    r_or = fwer_of(n, lam_true, np.full(D, c_or), s0, si)
    r_pl = fwer_of(n, lam_true, c_plug, s0, si)
    r_ev = fwer_of(n, lam_true, c_env_v, s0, si)
    lam_lo, lam_hi, diag = lam_region(n, lam_hat, GAMMA)
    inside = np.all(
        (lam_lo <= lam_true[None, :]) & (lam_true[None, :] <= lam_hi), axis=1
    )
    return {
        "case": label,
        "n": [f64(x) for x in n],
        "sigma": [f64(x) for x in sig],
        "lam_true": [f64(x) for x in lam_true],
        "alpha": ALPHA,
        "gamma": GAMMA,
        "draws": i64(D),
        "c_oracle": f64(c_or),
        "c_plug_mean": f64(np.mean(c_plug)),
        "c_env_mean": f64(np.mean(c_env_v)),
        "fwer_oracle": r_or["fwer"],
        "fwer_plug": r_pl["fwer"],
        "fwer_env": r_ev["fwer"],
        "mcse_oracle": r_or["mcse"],
        "mcse_plug": r_pl["mcse"],
        "mcse_env": r_ev["mcse"],
        "env_excess_over_oracle": f64(r_ev["fwer"] - r_or["fwer"]),
        "plug_excess_over_oracle": f64(r_pl["fwer"] - r_or["fwer"]),
        "cost_env_pct": f64(100.0 * (np.mean(c_env_v) - c_or) / c_or),
        "cost_plug_pct": f64(100.0 * (np.mean(c_plug) - c_or) / c_or),
        "region_miss_frac": f64(1.0 - float(np.mean(inside))),
        "region_scale": f64(diag.get("scale")),
        "coverage_exact": f64(diag.get("coverage_exact")),
        "guarantee_holds": bool(r_ev["fwer"] <= ALPHA + GAMMA),
        "dominance_holds": bool(r_ev["fwer"] <= r_or["fwer"] + 3.0 * r_ev["mcse"]),
    }


def main():
    t0 = time.time()
    rows = []
    pool = mp.Pool(POOL_N) if POOL_N > 1 else None
    try:
        for label, n_list, sig_list in CASES:
            r = cell(label, n_list, sig_list, pool)
            rows.append(r)
            print(
                f"{label:<32s} env {r['fwer_env']:.4f} ({r['cost_env_pct']:+.2f}% in c) "
                f"| plug {r['fwer_plug']:.4f} ({r['cost_plug_pct']:+.2f}%) "
                f"| oracle {r['fwer_oracle']:.4f} | miss {r['region_miss_frac']:.4f} "
                f"| guarantee {r['guarantee_holds']}",
                flush=True,
            )
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    worst_env = max(r["fwer_env"] for r in rows)
    worst_plug = max(r["fwer_plug"] for r in rows)
    out = {
        "alpha": ALPHA,
        "gamma": GAMMA,
        "draws_per_case": i64(D),
        "seed": i64(SEED),
        "rows": rows,
        "n_cases": i64(len(rows)),
        "worst_fwer_env": f64(worst_env),
        "worst_fwer_plug": f64(worst_plug),
        "all_guarantees_hold": bool(all(r["guarantee_holds"] for r in rows)),
        "all_dominance_holds": bool(all(r["dominance_holds"] for r in rows)),
        "mean_cost_env_pct": f64(np.mean([r["cost_env_pct"] for r in rows])),
        "mean_cost_plug_pct": f64(np.mean([r["cost_plug_pct"] for r in rows])),
        "seconds": f64(time.time() - t0),
        "note": (
            "the envelope's guarantee is alpha + gamma by construction; the plug-in has no "
            "such guarantee and its excess over the oracle is the estimation error"
        ),
    }
    dump_json(out, os.path.join(OUT, "paper_cases.json"))
    print(
        f"\nworst env {worst_env:.4f} vs alpha+gamma {ALPHA + GAMMA:.4f} | worst plug "
        f"{worst_plug:.4f} | all guarantees {out['all_guarantees_hold']} | "
        f"{time.time() - t0:.0f}s"
    )


if __name__ == "__main__":
    main()
