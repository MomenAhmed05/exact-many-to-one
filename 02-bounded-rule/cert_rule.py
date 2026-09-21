"""S1.3: the certified rule's realised FWER, and the price of the proof.

The Stage-1 certificate (`cert_sup.c_rigorous`) is theorem-grade: for every draw
and every lambda in the region the conditional rejection at that c is at most
alpha, by construction, so the rule using it satisfies FWER <= alpha (no gamma
slack needed at all).  This script measures, per desk cell and over multi-seed
draws:

  * the realised FWER of the certified rule (`fwer_of` against the oracle);
  * the price of the proof in `c` (certified vs the candidate envelope);
  * and both together, so the report can state where the theorem is affordable.

Outputs: deposits/cert/certified_rule.json
"""

import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cert_sup import DFS, FAMS, GAMMA, certify_draws  # noqa: E402
from core import ALPHA, draw_var  # noqa: E402
from env_rule import c_env, fwer_of, lam_region  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "cert")


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


SEEDS = (101, 202, 303, 404)
D_PER_SEED = _env_int("W26_CR_D", 25)
POOL_N = _env_int("W26_POOL", 3)


def cell(fam, df, pool):
    lam = np.asarray(FAMS[fam], float)
    n = np.array([df + 1, 13, 13], float)
    s0s, sis, lhs = [], [], []
    for seed in SEEDS:
        s0, si, lh = draw_var(n, lam, D_PER_SEED, seed)
        s0s.append(s0)
        sis.append(si)
        lhs.append(lh)
    s0 = np.concatenate(s0s)
    si = np.concatenate(sis)
    lh = np.concatenate(lhs)
    lo, hi, _ = lam_region(n, lh, GAMMA, scale=None)
    c_rig = certify_draws(n, s0, si, lo, hi, alpha=ALPHA, pool=pool)
    c_cand = c_env(n, s0, si, lh, gamma=GAMMA, alpha=ALPHA, pool=pool)
    r_rig = fwer_of(n, lam, c_rig, s0, si)
    r_cand = fwer_of(n, lam, c_cand, s0, si)
    inside = np.all((lam[None, :] >= lo) & (lam[None, :] <= hi), axis=1)
    try:
        n_miss = i64(int(np.count_nonzero(~inside)))
        miss_rate = f64(float(np.mean(~inside)))
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("miss aggregation failed") from exc
    return {
        "family": fam,
        "df0": df,
        "n": list(n),
        "lam_true": list(lam),
        "seeds": [i64(s) for s in SEEDS],
        "draws_per_seed": i64(D_PER_SEED),
        "draws_total": i64(len(s0)),
        "certified": {**r_rig, "passes_alpha": bool(r_rig["fwer"] <= ALPHA + 1e-9)},
        "candidate": {
            **r_cand,
            "passes_guarantee": bool(r_cand["fwer"] <= ALPHA + GAMMA),
        },
        "price_of_proof_pct": f64(r_rig["power_cost_pct"] - r_cand["power_cost_pct"]),
        "n_region_miss": n_miss,
        "miss_rate": miss_rate,
        "gamma": GAMMA,
    }


def main():
    t0 = time.time()
    rows = []
    with mp.Pool(POOL_N) as pool:
        for df in DFS:
            for fam in FAMS:
                r = cell(fam, df, pool)
                rows.append(r)
                print(
                    f"{fam} df{df:>2}: certified FWER {r['certified']['fwer']:.4f} "
                    f"(<= {ALPHA}) {r['certified']['passes_alpha']} | candidate FWER "
                    f"{r['candidate']['fwer']:.4f} | price of proof "
                    f"{r['price_of_proof_pct']:+.1f}% | misses "
                    f"{r['n_region_miss']}/{r['draws_total']}",
                    flush=True,
                )
    try:
        prices = [r["price_of_proof_pct"] for r in rows]
        med_price = f64(float(np.median(prices)))
    except Exception:  # noqa: BLE001
        med_price = f64(float("nan"))
    summary = {
        "rows": rows,
        "alpha": ALPHA,
        "gamma": GAMMA,
        "seeds": [i64(s) for s in SEEDS],
        "draws_per_seed": i64(D_PER_SEED),
        "all_certified_pass_alpha": bool(
            all(r["certified"]["passes_alpha"] for r in rows)
        ),
        "all_candidate_pass_guarantee": bool(
            all(r["candidate"]["passes_guarantee"] for r in rows)
        ),
        "total_draws": i64(sum(r["draws_total"] for r in rows)),
        "total_region_misses": i64(sum(r["n_region_miss"] for r in rows)),
        "max_price_of_proof_pct": f64(max(r["price_of_proof_pct"] for r in rows)),
        "median_price_of_proof_pct": med_price,
        "seconds": f64(time.time() - t0),
        "note": (
            "The certified rule needs no gamma slack (its defining property is "
            "per-draw and per-lambda); the candidate rule is the desk rule and "
            "spends gamma on region misses."
        ),
    }
    dump_json(summary, os.path.join(OUT, "certified_rule.json"))
    print(
        f"certified passes alpha: {summary['all_certified_pass_alpha']} | candidate "
        f"passes guarantee: {summary['all_candidate_pass_guarantee']} | misses "
        f"{summary['total_region_misses']}/{summary['total_draws']} | median price "
        f"{summary['median_price_of_proof_pct']:+.1f}% | {time.time() - t0:.0f}s"
    )


if __name__ == "__main__":
    main()
