"""Close residual #1: the gamma-miss rate, certified by quadrature AND probed empirically.

Two independent quantities per desk cell:

  designed   miss_designed = 1 - coverage_exact, where coverage_exact is the region's
             exact 1-D quadrature coverage at the tight level scale.  This is NOT an
             estimate: the region is built so that P(lambda_true in Lambda) integrates
             to that value, so the designed miss rate is certified, not sampled.
  observed   miss_hat over 20,000 independent draws of s^2: form each draw's region
             (the same tight construction the rule uses) and test membership of the
             true lambda.  Purely arithmetic per draw -- no critical value is solved,
             so this is cheap and can be run at large D.

Agreement between the two closes the residual: the single observed event in the desk
sweep (~1200 draws) is then consistent with a designed rate that is itself certified,
and any discrepancy would indict the region machinery rather than the gamma budget.

Outputs: deposits/cert/miss_probe.json
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import draw_var  # noqa: E402
from env_rule import lam_region  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "cert")
GAMMA = 0.0005
DESK_FAMS = (("A", (0.04, 0.04)), ("B", (1.0, 25.0)))
GAMMA_DF = (3, 5, 7, 11, 19, 35)


def _f(x, default=float("nan")):
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return float(default)


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


D_PROBE = _env_int("W26_MISS_D", 20000)
SEED = _env_int("W26_MISS_SEED", 20260917)


def main():
    rows = []
    tot_draws = 0
    tot_miss = 0.0
    for fam, lam in DESK_FAMS:
        lam_true = np.asarray(lam, float)
        for df0 in GAMMA_DF:
            n = np.array([df0 + 1, 13, 13], float)
            _s0, _si, lh = draw_var(
                n, lam_true, D_PROBE, SEED + 7 * df0 + (131 if fam == "B" else 0)
            )
            lam_hat = np.atleast_2d(lh)
            lam_lo, lam_hi, diag = lam_region(n, lam_hat, GAMMA)
            inside = np.all(
                (lam_lo <= lam_true[None, :]) & (lam_true[None, :] <= lam_hi), axis=1
            )
            p = _f(1.0 - float(np.mean(inside)))
            cov_q = _f(diag.get("coverage_exact"))
            rows.append(
                {
                    "family": fam,
                    "df0": i64(df0),
                    "n": [f64(x) for x in n],
                    "lam_true": [f64(x) for x in lam_true],
                    "draws": i64(D_PROBE),
                    "miss_observed": f64(p),
                    "miss_observed_mcse": f64(
                        np.sqrt(max(1e-300, p * (1 - p) / D_PROBE))
                    ),
                    "miss_designed_quadrature": f64(1.0 - cov_q),
                    "coverage_exact": f64(cov_q),
                    "region_scale": f64(diag.get("scale")),
                    "abs_diff": f64(abs(p - (1.0 - cov_q))),
                }
            )
            tot_draws += D_PROBE
            tot_miss += p * D_PROBE
    worst = max((r["abs_diff"] for r in rows), default=float("nan"))
    zs = [
        (r["miss_observed"] - r["miss_designed_quadrature"]) / r["miss_observed_mcse"]
        for r in rows
        if r["miss_observed_mcse"] > 0
    ]
    out = {
        "gamma": GAMMA,
        "draws_per_cell": i64(D_PROBE),
        "seed": i64(SEED),
        "seconds_note": "membership tests only; no critical values solved",
        "rows": rows,
        "n_cells": i64(len(rows)),
        "total_draws": i64(tot_draws),
        "pooled_miss_observed": f64(tot_miss / max(1, tot_draws)),
        "worst_abs_diff": f64(worst),
        "max_abs_z": f64(max((abs(z) for z in zs), default=float("nan"))),
        "verdict": (
            "the observed miss rate matches the quadrature-certified designed rate "
            "in every cell within Monte-Carlo error"
            if worst == worst and all(abs(z) <= 4.0 for z in zs)
            else "DISCREPANCY: inspect the region machinery"
        ),
        "note": (
            "the designed miss rate is 1 - coverage_exact, an exact quadrature value "
            "(the region is constructed so that P(lambda_true in Lambda) equals it); "
            "the observed column is an independent end-to-end check of that machinery"
        ),
    }
    dump_json(out, os.path.join(OUT, "miss_probe.json"))
    print(
        f"cells {len(rows)} | draws {tot_draws} | worst |obs - designed| = {worst:.2e} "
        f"| max |z| = {out['max_abs_z']:.2f} | pooled miss {out['pooled_miss_observed']:.2e}"
    )
    print(f"designed gamma = {GAMMA:.2e} | verdict: {out['verdict']}")


if __name__ == "__main__":
    main()
