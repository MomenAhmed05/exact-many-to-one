"""P4 step 24: the desk map re-run at the HARD-CEILING level alpha_0 = abar - gamma, with the
boundary between the certified tier, the measured tier, and the price of the proof.

For each of the twelve desk cells (two volatility patterns x six benchmark degrees of freedom):

  * `c_oracle`       the population oracle at the true ratios (`core.oracle_c`, tier
                     `oracle_exact_law`) - a reference point, not a rule;
  * `c_env_measured` the 02-bounded-rule vertex envelope, i.e. the max over the region's vertices of the
                     POPULATION critical value (`env_rule.c_env`, tier `envelope_measured`) - the
                     rule the 02-bounded-rule dossier shipped;
  * `c_certified`    the certified continuous supremum of the CONDITIONAL critical value at this
                     draw, at `alpha_0 = abar - gamma` (tier `region_sup_certified`) - the rule this
                     programme ships, with its achieved gap;
  * the coverage of the region, and the runtime.

The report states, per cell, which tier the threshold comes from, the achieved gap, and the price of
the proof relative to the measured envelope.

Run:  python desk_map_alpha0.py        (~10-20 min: one certification per cell)
Out:  deposits/desk_map_alpha0.json    (tier: region_sup_certified)
"""
from __future__ import annotations

import importlib
import os
import time

import numpy as np

import cert_bnb
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "desk_map_alpha0.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
draw_var = CORE.draw_var
oracle_c = CORE.oracle_c

FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
DFS = (3, 5, 7, 11, 19, 35)
ABAR = paths.ABAR
GAMMA = paths.GAMMA
ALPHA0 = paths.ALPHA0


def cell(fam, df, seed=20260918, max_boxes=64):
    lam = np.asarray(FAMS[fam], float)
    n = np.array([df + 1.0, 13.0, 13.0])
    s0, si, lh = draw_var(n, lam, 1, seed)
    d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    lo, hi, diag = ENV.lam_region(n, lh, GAMMA, scale=None)
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    # measured tier: the 02-bounded-rule vertex envelope (population critical values at the vertices)
    c_env = NAN
    try:
        c_env = _f(np.mean(ENV.c_env(n, s0, si, lh, gamma=GAMMA, alpha=ALPHA0)))
    except Exception:  # noqa: BLE001
        pass
    # LIKE-FOR-LIKE price: the certified value is CONDITIONAL (one draw), so the comparison
    # must be against the conditional vertex maximum - the same primitive the branch-and-bound
    # uses.  Comparing it against the POPULATION envelope (c_env) is an object mismatch and is
    # kept only as an informational column.
    cond_vertex = []
    for combo in range(2 ** lam.size):
        lam_v = np.array([hi[i] if (combo >> i) & 1 else lo[i] for i in range(lam.size)])
        rc = importlib.import_module("enclose_gl").c_upper_box_gl(
            d_row, n[0], np.sqrt(n[1:] / lam_v), np.sqrt(n[1:] / lam_v), alpha=ALPHA0)
        if rc.get("c_ub") == rc.get("c_ub"):
            cond_vertex.append(_f(rc["c_ub"]))
    c_cond_vertex = _f(max(cond_vertex)) if cond_vertex else NAN
    rng = np.random.default_rng(seed)
    bnb = cert_bnb.certify_region(lo, hi, d_row, n, alpha=ALPHA0, rng=rng,
                                  max_boxes=_i(max_boxes, 64))
    c_or = _f(oracle_c(n, lam, ALPHA0))
    return {
        "family": fam, "df0": _i(df), "n": list(n), "lam": list(lam),
        "alpha0": _f(ALPHA0), "abar": _f(ABAR), "gamma": _f(GAMMA),
        "coverage_exact": _f(diag.get("coverage_exact")),
        "c_oracle_population": c_or, "c_oracle_tier": "oracle_exact_law",
        "c_env_measured": c_env, "c_env_tier": "envelope_measured",
        "c_certified": bnb["c_U"], "c_certified_L": bnb["c_L"], "gap": bnb["gap"],
        "certified_at_tol": bnb["certified_at_tol"], "c_certified_tier": "region_sup_certified",
        "conditional_vertex_max": c_cond_vertex,
        "conditional_vertex_max_tier": "region_sup_certified",
        "price_vs_conditional_vertex_pct": (
            _f(100.0 * (_f(bnb["c_U"]) - c_cond_vertex) / c_cond_vertex)
            if c_cond_vertex > 0 else NAN),
        "population_envelope_minus_certified_pct_informational": (
            _f(100.0 * (c_env - _f(bnb["c_U"])) / c_env) if c_env > 0 else NAN),
        "price_vs_measured_pct": NAN,
        "price_vs_oracle_pct": _f(100.0 * (_f(bnb["c_U"]) - c_or) / c_or) if c_or > 0 else NAN,
        "n_splits": bnb["n_splits"], "n_pruned": bnb["n_pruned"],
        "termination": bnb["termination"], "seconds": bnb["seconds"],
        "lambda_hat": [round(_f(v), 8) for v in np.asarray(lh[0], float)],
        "region_lo": [round(_f(v), 8) for v in lo],
        "region_hi": [round(_f(v), 8) for v in hi],
    }


def main():
    t0 = time.time()
    rows = []
    for fam in FAMS:
        for df in DFS:
            r = cell(fam, df)
            rows.append(r)
            print(f"{fam}_df{df:<3} c_oracle {_f(r['c_oracle_population']):.4f} | "
                  f"env {_f(r['c_env_measured']):.4f} | certified {_f(r['c_certified']):.6f} "
                  f"gap {_f(r['gap']):.2e} ok {r['certified_at_tol']} "
                  f"| price vs cond-vertex {_f(r['price_vs_conditional_vertex_pct']):+.4f}% "
                  f"| {_f(r['seconds']):.0f}s", flush=True)
    gaps = [_f(r["gap"]) for r in rows]
    prices = [_f(r["price_vs_measured_pct"]) for r in rows if _f(r["price_vs_measured_pct"]) == _f(r["price_vs_measured_pct"])]
    out = {
        "n_cells": _i(len(rows)),
        "alpha0": _f(ALPHA0), "abar": _f(ABAR), "gamma": _f(GAMMA), "tol": _f(paths.TOL),
        "n_certified_at_tol": _i(sum(1 for r in rows if r["certified_at_tol"])),
        "max_gap": _f(max(gaps)), "median_gap": _f(np.median(gaps)),
        "max_price_vs_measured_pct": _f(max(prices)) if prices else NAN,
        "median_price_vs_measured_pct": _f(np.median(prices)) if prices else NAN,
        "boundary": ("the certified tier is the shipped rule; the measured envelope is the 02-bounded-rule rule; "
                     "the oracle is a reference.  Cells whose gap exceeds tol are named in the "
                     "deposit rather than hidden."),
        "seconds": _f(time.time() - t0),
        "cells": rows,
    }
    n_cert = out["n_certified_at_tol"]
    print(f"\n{n_cert}/{out['n_cells']} certified at tol "
          f"{paths.TOL:.0e} | max gap {_f(out['max_gap']):.2e} | "
          f"price vs measured envelope: median {_f(out['median_price_vs_measured_pct']):+.2f}% "
          f"max {_f(out['max_price_vs_measured_pct']):+.2f}% | {_f(out['seconds']):.0f}s")
    tiers.write_deposit(out, OUT, "region_sup_certified",
                        notes="desk map at alpha0 = abar - gamma: certified tier, measured tier and "
                              "the price of the proof, per cell")
    return out


if __name__ == "__main__":
    import sys

    try:
        res = main()
        sys.exit(0 if res.get("n_certified_at_tol", 0) > 0 else 1)
    except Exception as exc:  # noqa: BLE001
        print("DESK MAP ERROR:", repr(exc))
        sys.exit(2)
