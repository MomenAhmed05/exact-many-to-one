"""W2: direction-aware, design-adaptive confidence-bound rules.

Designs (families A and B, benchmark df swept):
  A  n = (n0, 13, 13), sigma = (50, 10, 10)  ->  lam = (0.04, 0.04)
  B  n = (n0, 13, 13), sigma = (10, 10, 50)  ->  lam = (1.00, 25.0)

Rules evaluated (all as transformations of lam_hat, so they share the SAME
per-draw solves where possible):
  plug-in            c(lam_hat)
  CB-g               c(lam_hat / F_g)          lower confidence bound on lam
  CBdir-g            per draw max(c(lower), c(upper))
  CQ-g               g-quantile of c(lam_hat*), lam_hat* ~ F-law at lam_hat
  DEB               2 c(lam_hat) - mean_b c(lam_hat*)

Every rule is evaluated paired against the oracle constant with the true-ratio
conditional rejection probability (Rao-Blackwellised).

Usage:  python w2_cb.py <phase>      phase in {explore, final, mono}
"""

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import (
    ALPHA,
    batch_c,
    cond_noreject,
    draw_var,
    f_quantile,  # noqa
    lam_lower,
    lam_upper,
    oracle_c,
)
from io_util import dump_json, f64, i64  # noqa

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w2")


def bootstrap_c(n, lh, B, alpha, seed, workers=10, pool=None):
    """c(lam_hat*) for B F-law draws per outer draw.  Returns (B, D) array."""
    n = np.asarray(n, float)
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    D, K = lh.shape
    # F* ~ F(nu_i, nu_0); lam_hat* = lam_hat * F*
    fstar = np.empty((B, D, K))
    for i in range(K):
        fstar[:, :, i] = rng.f(nu[i + 1], nu[0], size=(B, D))
    lam_star = (lh[None, :, :] * fstar).reshape(B * D, K)
    cs = batch_c(n, lam_star, alpha=alpha, workers=workers, pool=pool).reshape(B, D)
    return cs


def evaluate(
    n, lam_true, alpha=ALPHA, D=200, seeds=(1,), B=25, workers=10, rules=None, pool=None
):
    n = np.asarray(n, float)
    lam_true = np.asarray(lam_true, float)
    out = {
        "n": list(n),
        "lam_true": list(lam_true),
        "alpha": alpha,
        "B_boot": B,
        "seeds": list(seeds),
        "draws_per_seed": D,
        "rules": {},
        "per_seed": {},
        "nan_solves": 0,
    }
    c_or = oracle_c(n, lam_true, alpha)
    out["c_oracle"] = c_or
    per_seed = {k: [] for k in rules}
    t0 = time.time()
    for seed in seeds:
        s0, si_, lh = draw_var(n, lam_true, D, seed)
        # --- one batched solve per rule family -------------------------
        lams, names = [lh], ["plugin"]
        for rule in rules:
            if rule.startswith("cbdir"):
                gam = f64(rule.split("-")[1])
                lams += [lam_lower(n, lh, gam), lam_upper(n, lh, gam)]
                names += [rule + ":lo", rule + ":hi"]
            elif rule.startswith("cb"):
                gam = f64(rule.split("-")[1])
                lams.append(lam_lower(n, lh, gam))
                names.append(rule)
        stack = np.vstack(lams)
        solved = batch_c(n, stack, alpha=alpha, workers=workers, pool=pool)
        out["nan_solves"] += i64(np.sum(~np.isfinite(solved)))
        solved = solved.reshape(len(lams), D)
        cmap = {nm: solved[j] for j, nm in enumerate(names)}
        cs_star = bootstrap_c(n, lh, B, alpha, seed + 500000, workers, pool=pool)
        cvecs = {}
        for rule in rules:
            if rule in cmap:
                cvecs[rule] = cmap[rule]
            elif rule.startswith("cbdir"):
                cvecs[rule] = np.maximum(cmap[rule + ":lo"], cmap[rule + ":hi"])
            elif rule.startswith("cq"):
                gam = f64(rule.split("-")[1])
                cvecs[rule] = np.quantile(cs_star, gam, axis=0)
            elif rule == "deb":
                cvecs[rule] = 2 * cmap["plugin"] - cs_star.mean(axis=0)
            elif rule == "oracle":
                cvecs[rule] = np.full(D, c_or)
            else:
                raise RuntimeError(f"unknown rule {rule}")
        nr_o = cond_noreject(np.full(D, c_or), s0, si_, n, lam_true)
        for rule, cv in cvecs.items():
            nr = cond_noreject(cv, s0, si_, n, lam_true)
            d = nr_o - nr
            per_seed[rule].append(
                (
                    alpha + f64(d.mean()),
                    f64(d.std(ddof=1) / np.sqrt(D)),
                    f64(np.mean(cv)),
                )
            )
    for rule in rules:
        rows = np.array(per_seed[rule])
        f = rows[:, 0]
        out["rules"][rule] = {
            "fwer_mean": f64(f.mean()),
            "seed_se": f64(f.std(ddof=1) / np.sqrt(len(f))) if len(f) > 1 else None,
            "fwer_seeds": [f64(x) for x in f],
            "rel_err_pct": 100 * (f64(f.mean()) - alpha) / alpha,
            "mcse_within_seed": f64(rows[:, 1].mean()),
            "mean_c": f64(rows[:, 2].mean()),
            "power_cost_pct": 100 * (f64(rows[:, 2].mean()) - c_or) / c_or,
        }
    out["seconds"] = time.time() - t0
    return out


def monotonicity_map(
    n0_list=(4, 6, 8, 12, 20, 36), lam_list=((0.04, 0.04), (1.0, 25.0), (0.25, 4.0))
):
    """Sign of dc/dlam at the true lam, as a function of benchmark df."""
    rows = []
    for lam in lam_list:
        for n0 in n0_list:
            n = [n0, 13, 13]
            lam = list(lam)
            q = [
                f_quantile(n[1] - 1, n[0] - 1, 0.9),
                f_quantile(n[2] - 1, n[0] - 1, 0.9),
            ]
            c0 = oracle_c(n, lam)
            c_lo = oracle_c(n, [lam[0] / q[0], lam[1] / q[1]])
            c_hi = oracle_c(n, [lam[0] * q[0], lam[1] * q[1]])
            rows.append(
                {
                    "n": n,
                    "lam": lam,
                    "c": c0,
                    "c_lower": c_lo,
                    "c_upper": c_hi,
                    "dc_dlam_sign": (
                        "negative (lower bound is conservative)"
                        if c_lo > c0 > c_hi
                        else "positive (upper bound is conservative)"
                        if c_hi > c0 > c_lo
                        else "mixed/non-monotone"
                    ),
                }
            )
    return rows


def main(phase):
    if phase == "mono":
        res = monotonicity_map()
        p = os.path.join(OUTDIR, "mono.json")
        dump_json(res, p)
        for r in res:
            print(
                f"n={r['n']} lam={r['lam']} c={r['c']:.4f} "
                f"lo={r['c_lower']:.4f} hi={r['c_upper']:.4f} "
                f"{r['dc_dlam_sign']}"
            )
        print("wrote", p)
        return
    fams = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
    if phase == "explore":
        rules = ["plugin", "cb-0.8", "cb-0.9", "cbdir-0.9", "cq-0.7", "cq-0.8", "deb"]
        cells = [(fam, n0) for fam in ("A", "B") for n0 in (4, 6, 8, 12, 20, 36)]
        D, seeds, B = 200, (11,), 25
    else:
        raise RuntimeError(
            'the final acceptance run lives in w2_final.py; w2_cb.py only '
            'runs the exploration phase'
        )
    import multiprocessing as mp

    with mp.Pool(6) as pool:
        for fam, n0 in cells:
            n = [n0, 13, 13]
            key = f"{fam}_n0{n0}"
            p = os.path.join(OUTDIR, f"explore_{key}.json")
            if os.path.exists(p):
                print("skip", p)
                continue
            res = evaluate(n, fams[fam], D=D, seeds=seeds, B=B, rules=rules, pool=pool)
            res["family"] = fam
            dump_json(res, p)
            print(
                key,
                json.dumps(
                    {
                        k: [round(v["rel_err_pct"], 2), round(v["power_cost_pct"], 2)]
                        for k, v in res["rules"].items()
                    }
                ),
            )
            print("wrote", p)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "explore")
