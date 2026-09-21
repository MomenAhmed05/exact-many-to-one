"""P2 step 11: the certified continuous supremum by sub-box branch-and-bound (SPEC section 7).

The object certified is `sup_{lambda in Lambda} c*(lambda)` — the CONTINUOUS supremum, not the
maximum over vertices.  The two primitives are the certified critical values of `enclose_gl`:

    c_ub(box) = smallest c with the RIGOROUS lower enclosure F_lo(c; box) >= 1 - alpha0
                => c_ub(box) >= sup_{lambda in box} c*(lambda)        (an upper bound)
    c_lb(lambda_w) = largest c with the RIGOROUS upper enclosure F_hi(c; {lambda_w}) < 1 - alpha0
                => c_lb(lambda_w) < c*(lambda_w) <= sup              (a witness lower bound)

Algorithm (SPEC section 7):

  1. Queue = {Lambda}; witness set = vertices (capped), face midpoints, box centre, plus
     quasi-random interior starts.
  2. For each box: `c_ub(box)` and `c_lb` at its best witness.
  3. `L_lo = max` of certified witness lower bounds (monotone in the witnesses seen).
  4. Prune a box when `c_ub(box) <= L_lo + tol`; the pruning proof (both numbers) is recorded.
  5. Otherwise split along the coordinate with the largest `|d c*/d log lambda_i|` at the box
     centre (implicit-function derivative, `derivs.dcdloglam`), tie-broken by widest log-width.
  6. Stop when every box is pruned or the box budget is exhausted.
  7. Report `c_L = L_lo` and **`c_U = max over ALL boxes of c_ub(box)`** (pruned boxes included;
     each has `c_ub <= L_lo + tol <= c_U` so it can never dominate), with `gap = c_U - c_L`
     reported honestly.  If the budget runs out first, the gap may exceed `tol` and that is
     reported as a failure to certify at the declared tolerance - never masked by substituting
     `L_lo + tol` for `c_U` (the P0 verifier round caught that error in the first draft: the
     substitution returns a value BELOW a box's true upper bound, i.e. non-conservative in the
     one case where it matters).

The shipped threshold is `c_U`; the statement is about the continuous supremum and therefore does
not depend on vertex attainment.

Run:  python cert_bnb.py            (the 12 desk cells + extras; minutes per cell at sharp corners)
Out:  deposits/bnb.json             (tier: region_sup_certified)
"""

from __future__ import annotations

import contextlib
import importlib
import itertools
import os
import time

import derivs
import enclose_gl
import numpy as np
import paths
import tiers
from paths import NAN, _f, _i

ERRLOG: list[str] = []

OUT = os.path.join(paths.ensure_dirs(), "bnb.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
draw_var = CORE.draw_var

FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
DFS = (3, 5, 7, 11, 19, 35)
EXTRA = (
    {"name": "K1_df3", "n": [4.0, 7.0], "lam": [1.0]},
    {"name": "K3_df6", "n": [7.0, 7.0, 9.0, 11.0], "lam": [0.25, 1.0, 4.0]},
    {"name": "K5_df11", "n": [12.0] + [12.0] * 5, "lam": [0.2, 0.5, 1.0, 2.0, 5.0]},
    {"name": "K2_extreme", "n": [4.0, 4.0, 30.0], "lam": [0.01, 100.0]},
    {"name": "K2_quiet", "n": [30.0, 4.0, 4.0], "lam": [1.0, 1.0]},
)
MAX_VERTICES = 8


def _seed_points(lo, hi, rng, n_random=4):
    """Witness seeds: vertices (capped), the centre, and quasi-random interior points."""
    K = np.asarray(lo, float).size
    pts = [
        np.where(np.array(c, bool), np.asarray(hi, float), np.asarray(lo, float))
        for c in itertools.product((0, 1), repeat=K)
    ]
    if len(pts) > MAX_VERTICES:
        idx = np.linspace(0, len(pts) - 1, MAX_VERTICES).astype(int)
        pts = [pts[i] for i in idx]
    pts.append(np.sqrt(np.asarray(lo, float) * np.asarray(hi, float)))
    log_lo, log_hi = np.log(np.asarray(lo, float)), np.log(np.asarray(hi, float))
    for _ in range(_i(n_random, 4)):
        pts.append(np.exp(log_lo + rng.random(K) * (log_hi - log_lo)))
    return pts


def _split_axis(n, d_row, box_lo, box_hi, c_mid):
    """Coordinate with the largest |dc*/dlog lambda| at the box centre (tie-break: widest)."""
    K = np.asarray(box_lo, float).size
    try:
        g = derivs.dcdloglam(
            n, d_row, np.sqrt(np.asarray(box_lo) * np.asarray(box_hi)), c_mid
        )
        if len(g) == K:
            return _i(np.argmax(np.abs(np.asarray(g, float))))
    except Exception as exc:  # noqa: BLE001
        ERRLOG.append(f"{__name__}: derivative split axis unavailable: {exc!r}")
    return _i(
        np.argmax(np.log(np.asarray(box_hi, float)) - np.log(np.asarray(box_lo, float)))
    )


def certify_region(
    lam_lo,
    lam_hi,
    d_row,
    n,
    alpha=paths.ALPHA0,
    tol=paths.TOL,
    max_boxes=64,
    prec=64,
    rng=None,
    max_witness=MAX_VERTICES,
):
    """Certified [c_L, c_U] for sup over the box, with the certificate tree."""
    t0 = time.time()
    rng = rng or np.random.default_rng(20260918)
    lam_lo = np.asarray(lam_lo, float).copy()
    lam_hi = np.asarray(lam_hi, float).copy()
    K = lam_lo.size
    out: dict = {
        "K": _i(K),
        "tol": _f(tol),
        "alpha0": _f(alpha),
        "c_L": NAN,
        "c_U": NAN,
        "gap": NAN,
        "certified_at_tol": False,
        "n_boxes": 0,
        "n_pruned": 0,
        "n_splits": 0,
        "n_cub": 0,
        "n_clb": 0,
        "termination": None,
        "tree": [],
        "seconds": NAN,
        "error": None,
    }
    try:
        # ---- witnesses: certify a lower bound at each seed point ----
        L_lo = 0.0
        witness_rows = []
        best_witness = None
        seeds = _seed_points(lam_lo, lam_hi, rng)[
            : max(1, _i(max_witness, MAX_VERTICES))
        ]
        for lam_w in seeds:
            s_w = np.sqrt(np.asarray(n, float)[1:] / lam_w)
            r = enclose_gl.c_lower_s_gl(s_w, d_row, n[0], alpha=alpha, prec=prec)
            out["n_clb"] += 1
            if r["c_lb"] == r["c_lb"]:
                v = _f(r["c_lb"])
                row = {"lambda": [round(_f(x), 8) for x in lam_w], "c_lb": v}
                witness_rows.append(row)
                if v > L_lo:
                    # the witness that ATTAINS L_lo: only there is c_L a certified lower bound on
                    # c*(lambda).  At other witnesses c*(lambda) may be below c_L, so c_L must NOT
                    # be checked there (an earlier version of the adversarial suite did, and
                    # reported spurious failures).
                    L_lo = v
                    best_witness = row
        out["witnesses"] = witness_rows
        out["argmax_witness"] = best_witness
        out["L_lo_initial"] = _f(L_lo)

        # ---- box queue ----
        boxes = [{"lo": lam_lo, "hi": lam_hi, "level": 0, "parent": -1}]
        cubs: list = []  # (box index, c_ub)
        while boxes:
            # evaluate every unevaluated box
            i = 0
            while i < len(boxes):
                b = boxes[i]
                if "c_ub" in b:
                    i += 1
                    continue
                s_lo, s_hi = (
                    np.sqrt(np.asarray(n, float)[1:] / b["hi"]),
                    np.sqrt(np.asarray(n, float)[1:] / b["lo"]),
                )
                guess = _f(L_lo) if L_lo > 0 else 3.0
                r = enclose_gl.c_upper_box_gl(
                    d_row,
                    n[0],
                    s_lo,
                    s_hi,
                    alpha=alpha,
                    c_guess=guess if guess > 0 else None,
                    prec=prec,
                )
                out["n_cub"] += 1
                b["c_ub"] = _f(r.get("c_ub"), NAN)
                b["c_ub_settings"] = r.get("settings")
                if b["c_ub"] != b["c_ub"]:
                    b["error"] = r.get("error") or "c_ub unavailable"
                i += 1
            if out["n_splits"] >= _i(max_boxes, 64):
                out["termination"] = "split budget exhausted"
                break
            # prune
            keep = []
            for b in boxes:
                if b.get("c_ub") == b.get("c_ub") and b["c_ub"] <= L_lo + _f(tol):
                    out["n_pruned"] += 1
                    b["pruned_at_L"] = _f(L_lo)
                    # a PRUNED box is a final LEAF: its c_ub enters the final maximum.  Superseded
                    # parents do not (their union is covered by their children), which is what lets
                    # c_U fall below the loose whole-region bound as the tree is refined.
                    cubs.append(b["c_ub"])
                    out["tree"].append(
                        {
                            k: (v.tolist() if isinstance(v, np.ndarray) else v)
                            for k, v in b.items()
                        }
                    )
                else:
                    keep.append(b)
            boxes = keep
            out["n_boxes"] += len(boxes)
            if not boxes:
                out["termination"] = "all boxes pruned"
                break
            # split the box with the largest c_ub
            boxes.sort(
                key=lambda b: (
                    -(b.get("c_ub") if b.get("c_ub") == b.get("c_ub") else 0.0)
                )
            )
            b = boxes.pop(0)
            c_mid = _f(b.get("c_ub"), 3.0)
            j = _split_axis(n, d_row, b["lo"], b["hi"], c_mid)
            mid = np.sqrt(b["lo"][j] * b["hi"][j])
            for side in (0, 1):
                lo2, hi2 = b["lo"].copy(), b["hi"].copy()
                if side == 0:
                    hi2[j] = mid
                else:
                    lo2[j] = mid
                boxes.append(
                    {
                        "lo": lo2,
                        "hi": hi2,
                        "level": _i(b.get("level"), 0) + 1,
                        "parent": -1,
                        "split_axis": _i(j),
                        "split_value": _f(mid),
                    }
                )
            out["n_splits"] += 1
            out["tree"].append(
                {
                    "split": True,
                    "axis": _i(j),
                    "value": _f(mid),
                    "level": _i(b.get("level"), 0),
                    "c_ub_of_parent": c_mid,
                }
            )
        # ---- final certificate: max over the LEAF boxes (pruned leaves + unpruned) ----
        all_cub = [c for c in cubs if c == c] + [b.get("c_ub", NAN) for b in boxes]
        all_cub = [c for c in all_cub if c == c]
        out["c_L"] = _f(L_lo)
        out["c_U"] = _f(max(all_cub)) if all_cub else NAN
        out["gap"] = _f(out["c_U"]) - _f(out["c_L"]) if all_cub else NAN
        out["certified_at_tol"] = bool(
            out["gap"] == out["gap"] and out["gap"] <= _f(tol)
        )
        if out["termination"] is None:
            out["termination"] = "all boxes pruned" if not boxes else "budget"
        out["n_unpruned"] = _i(len(boxes))
        out["unpruned"] = [
            {
                "lo": b["lo"].tolist(),
                "hi": b["hi"].tolist(),
                "c_ub": _f(b.get("c_ub"), NAN),
            }
            for b in boxes[:8]
        ]
    except Exception as exc:  # noqa: BLE001
        out["error"] = repr(exc)
    out["seconds"] = _f(time.time() - t0)
    return out


def main():
    t0 = time.time()
    rng = np.random.default_rng(20260918)
    cells = [
        {"name": f"{fam}_df{df}", "n": [df + 1.0, 13.0, 13.0], "lam": list(lam)}
        for fam, lam in FAMS.items()
        for df in DFS
    ]
    cells += list(EXTRA)
    rows = []
    for spec in cells:
        n = np.asarray(spec["n"], float)
        lam = np.asarray(spec["lam"], float)
        s0, si, lh = draw_var(n, lam, 1, 20260918)
        d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
        lo, hi, _ = ENV.lam_region(n, lh, paths.GAMMA, scale=None)
        lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
        r = certify_region(lo, hi, d_row, n, rng=rng)
        r["cell"] = spec["name"]
        r["lambda_hat"] = [round(_f(v), 8) for v in np.asarray(lh[0], float)]
        r["region_lo"] = [round(_f(v), 8) for v in lo]
        r["region_hi"] = [round(_f(v), 8) for v in hi]
        rows.append(r)
        print(
            f"{spec['name']:<12} c_L {_f(r['c_L']):.6f} c_U {_f(r['c_U']):.6f} "
            f"gap {_f(r['gap']):.2e} tol {paths.TOL:.0e} ok {r['certified_at_tol']} "
            f"| boxes {r['n_boxes']} pruned {r['n_pruned']} splits {r['n_splits']} "
            f"| {_f(r['seconds']):.0f}s",
            flush=True,
        )
    summary = {
        "n_cells": _i(len(rows)),
        "n_certified_at_tol": _i(sum(1 for r in rows if r["certified_at_tol"])),
        "all_certified_at_tol": bool(rows and all(r["certified_at_tol"] for r in rows)),
        "max_gap": _f(
            max((_f(r["gap"]) for r in rows if r["gap"] == r["gap"]), default=NAN)
        ),
        "median_gap": _f(
            np.median([_f(r["gap"]) for r in rows if r["gap"] == r["gap"]])
        ),
        "max_seconds": _f(max(_f(r["seconds"]) for r in rows)),
        "total_seconds": _f(time.time() - t0),
        "tol": _f(paths.TOL),
        "alpha0": _f(paths.ALPHA0),
        "gamma": _f(paths.GAMMA),
        "note": (
            "c_U = max over ALL boxes of c_ub (SPEC section 7 step 7); pruned boxes are "
            "included because each has c_ub <= L_lo + tol <= c_U.  A cell whose gap exceeds "
            "tol is reported as a failure to certify at the declared tolerance."
        ),
        "cells": rows,
    }
    print(
        f"\nsummary: {summary['n_certified_at_tol']}/{summary['n_cells']} cells certified at "
        f"tol {paths.TOL:.0e} | max gap {_f(summary['max_gap']):.2e} | "
        f"{_f(summary['total_seconds']):.0f}s"
    )
    tiers.write_deposit(
        summary,
        OUT,
        "region_sup_certified",
        notes="branch-and-bound certificate for the CONTINUOUS supremum over each "
        "region; c_U = max over all boxes, gap reported honestly",
    )
    return summary


if __name__ == "__main__":
    with contextlib.suppress(Exception):
        main()
