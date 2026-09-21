"""W2 final acceptance table: multi-seed, both families, df 3..35.

Rules: plugin, cb0.9, cb0.7, cbdir0.9 (per-draw max of the lower/upper CB),
cap0.05, cap0.10 (conservative CB capped at +cap in the critical value).

3 seeds x 100 draws per cell; paired against the oracle; one batched solve per
rule set per cell.

Outputs: deposits/w2/final_df<df>_<fam>.json
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import (
    ALPHA,
    batch_c,
    cond_noreject,
    draw_var,
    f_quantile,  # noqa
    oracle_c,
)
from io_util import dump_json, f64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w2")
FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
RULES = ["plugin", "cb0.9", "cb0.7", "cbdir0.9", "cap0.05", "cap0.10"]
SEEDS = (11, 22, 33)
D = 100
CAPS = {"cap0.05": 0.05, "cap0.10": 0.10}


def _main(pool, dfs):
    for df in dfs:
        for fam, lam in FAMS.items():
            p = os.path.join(OUTDIR, f"final_df{df}_{fam}.json")
            if os.path.exists(p):
                continue
            n = np.asarray([df + 1, 13, 13], float)
            lam = np.asarray(lam, float)
            nu = n - 1.0
            c_or = oracle_c(n, lam, ALPHA)
            q09 = np.array(
                [f_quantile(nu[i + 1], nu[0], 0.9) for i in range(len(nu) - 1)]
            )
            q07 = np.array(
                [f_quantile(nu[i + 1], nu[0], 0.7) for i in range(len(nu) - 1)]
            )
            t0 = time.time()
            per_seed = {r: [] for r in RULES}
            for seed in SEEDS:
                s0, si, lh = draw_var(n, lam, D, seed)
                lams = [
                    lh,
                    lh / q09[None, :],
                    lh / q07[None, :],
                    lh / q09[None, :],
                    lh * q09[None, :],
                ]
                solved = batch_c(
                    n, np.vstack(lams), alpha=ALPHA, workers=3, pool=pool
                ).reshape(5, D)
                cplug, c09, c07, clo, chi = solved
                ccons = np.maximum(clo, chi)
                cvecs = {"plugin": cplug, "cb0.9": c09, "cb0.7": c07, "cbdir0.9": ccons}
                for name, cap in CAPS.items():
                    cvecs[name] = np.minimum(ccons, cplug + cap)
                nr_o = cond_noreject(np.full(D, c_or), s0, si, n, lam)
                for rule, cv in cvecs.items():
                    nr = cond_noreject(cv, s0, si, n, lam)
                    d = nr_o - nr
                    per_seed[rule].append(
                        (
                            ALPHA + f64(d.mean()),
                            f64(d.std(ddof=1) / np.sqrt(D)),
                            f64(cv.mean()),
                        )
                    )
            out = {
                "n": list(n),
                "lam": list(lam),
                "alpha": ALPHA,
                "D": D,
                "seeds": list(SEEDS),
                "c_oracle": c_or,
                "rules": {},
            }
            for rule in RULES:
                rows = np.array(per_seed[rule])
                f = rows[:, 0]
                out["rules"][rule] = {
                    "fwer": f64(f.mean()),
                    "seed_se_pct": f64(100 * f.std(ddof=1) / np.sqrt(len(f)) / ALPHA),
                    "rel_err_pct": 100 * (f64(f.mean()) - ALPHA) / ALPHA,
                    "mean_c": f64(rows[:, 2].mean()),
                    "power_cost_pct": 100 * (f64(rows[:, 2].mean()) - c_or) / c_or,
                    "seeds_fwer": [f64(x) for x in f],
                }
            out["seconds"] = time.time() - t0
            dump_json(out, p)
            line = " ".join(
                f"{k}:{v['rel_err_pct']:+.1f}%" for k, v in out["rules"].items()
            )
            print(f"df{df} {fam} {line}", flush=True)


def _cli_ints(argv):
    try:
        return tuple(int(x) for x in argv)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"bad df list: {argv}") from exc


if __name__ == "__main__":
    dfs = _cli_ints(sys.argv[1:]) if len(sys.argv) > 1 else (3, 5, 7, 11, 19, 35)
    with mp.Pool(3) as pool:
        _main(pool, dfs)
