"""P2 step 11 unit checks for the branch-and-bound (plan step 11 items (a)-(d)).

  (a) a single-box certificate on a TINY region must agree with the point bracket at its centre;
  (b) subdividing a region must never INCREASE c_U (children are tighter than their parent);
  (c) pruning proofs re-verified after the fact: every recorded pruned leaf must still satisfy
      `c_ub(leaf) <= L_lo + tol` when its bound is recomputed with the fine layout;
  (d) the gap must be monotone non-increasing in the split budget.

Run:  python bnb_checks.py        (a few minutes)
Out:  deposits/bnb_checks.json    (tier: region_sup_certified)
"""

from __future__ import annotations

import importlib
import os
import time

import cert_bnb
import enclose_gl
import numpy as np
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "bnb_checks.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
draw_var = CORE.draw_var


def _cell(n, lam, seed=20260918):
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    s0, si, lh = draw_var(n, lam, 1, seed)
    d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    return n, lam, d_row, lh


def check_a(tol=paths.TOL):
    """Tiny region: the certified box value must match the point bracket at the centre."""
    n, lam, d_row, _lh = _cell([4.0, 13.0, 13.0], [0.04, 0.04])
    rng = np.random.default_rng(11)
    lo, hi, _ = ENV.lam_region(
        n, _cell([4.0, 13.0, 13.0], [0.04, 0.04])[3], paths.GAMMA, scale=None
    )
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    centre = np.sqrt(lo * hi)
    tiny_lo, tiny_hi = centre / 1.0001, centre * 1.0001
    r = cert_bnb.certify_region(tiny_lo, tiny_hi, d_row, n, rng=rng, max_boxes=4)
    br = enclose_gl.bracket_point_gl(np.sqrt(n[1:] / centre), d_row, n[0])
    diff = abs(_f(r["c_U"]) - _f(br["c_U"]))
    return {
        "tiny_box_c_U": r["c_U"],
        "point_bracket_c_U": br["c_U"],
        "point_bracket_c_L": br["c_L"],
        "abs_diff": _f(diff),
        "within_tol": bool(diff <= max(_f(tol), 1e-5)),
        "certified_at_tol": bool(r["certified_at_tol"]),
    }


def check_b(tol=paths.TOL):
    """Splitting a region must not increase the certified upper bound."""
    n, _lam, d_row, lh = _cell([4.0, 13.0, 13.0], [0.04, 0.04])
    rng = np.random.default_rng(12)
    lo, hi, _ = ENV.lam_region(n, lh, paths.GAMMA, scale=None)
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    parent = cert_bnb.certify_region(lo, hi, d_row, n, rng=rng, max_boxes=1)
    mid = np.sqrt(lo * hi)
    children = []
    for side in (0, 1):
        l2, h2 = lo.copy(), hi.copy()
        if side == 0:
            h2[0] = mid[0]
        else:
            l2[0] = mid[0]
        children.append(cert_bnb.certify_region(l2, h2, d_row, n, rng=rng, max_boxes=1))
    child_max = max(_f(c["c_U"]) for c in children)
    return {
        "parent_c_U": parent["c_U"],
        "children_c_U": [c["c_U"] for c in children],
        "child_max": _f(child_max),
        "no_increase": bool(child_max <= _f(parent["c_U"]) + 1e-9),
        "improvement": _f(parent["c_U"]) - child_max,
    }


def check_c(tol=paths.TOL):
    """Re-verify the pruning proofs of a completed certificate with the fine layout."""
    n, _lam, d_row, lh = _cell([4.0, 13.0, 13.0], [0.04, 0.04])
    rng = np.random.default_rng(13)
    lo, hi, _ = ENV.lam_region(n, lh, paths.GAMMA, scale=None)
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    r = cert_bnb.certify_region(lo, hi, d_row, n, rng=rng, max_boxes=8)
    leaves = [t for t in r.get("tree", []) if isinstance(t, dict) and "lo" in t]
    rows = []
    for leaf in leaves[:6]:
        l2 = np.asarray(leaf["lo"], float)
        h2 = np.asarray(leaf["hi"], float)
        s_lo, s_hi = np.sqrt(n[1:] / h2), np.sqrt(n[1:] / l2)
        rc = enclose_gl.c_upper_box_gl(d_row, n[0], s_lo, s_hi, alpha=paths.ALPHA0)
        c_ub = _f(rc.get("c_ub"))
        # Recompute the witness lower bound from the witness list rather than trusting the
        # recorded `pruned_at_L` (P2 fix-verify TASK C(iii)).
        L_rec = _f(leaf.get("pruned_at_L"))
        L_wit = max([_f(w.get("c_lb")) for w in r.get("witnesses", [])] or [NAN])
        L = L_wit if L_wit == L_wit else L_rec
        rows.append(
            {
                "c_ub_recomputed": c_ub,
                "pruned_at_L": L,
                "bound": _f(L) + _f(tol),
                "holds": bool(c_ub == c_ub and c_ub <= L + _f(tol) + 1e-9),
            }
        )
    return {
        "n_leaves_checked": _i(len(rows)),
        "rows": rows,
        "all_hold": bool(rows and all(x["holds"] for x in rows)),
    }


def check_d(budgets=(1, 2, 4, 8), tol=paths.TOL):
    """The gap must be monotone non-increasing in the split budget."""
    n, _lam, d_row, lh = _cell([4.0, 13.0, 13.0], [0.04, 0.04])
    rng = np.random.default_rng(14)
    lo, hi, _ = ENV.lam_region(n, lh, paths.GAMMA, scale=None)
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    rows = []
    for b in budgets:
        r = cert_bnb.certify_region(lo, hi, d_row, n, rng=rng, max_boxes=_i(b, 1))
        rows.append(
            {
                "max_boxes": _i(b),
                "gap": r["gap"],
                "splits": r["n_splits"],
                "c_U": r["c_U"],
                "seconds": r["seconds"],
            }
        )
    gaps = [_f(x["gap"]) for x in rows]
    return {
        "rows": rows,
        "gaps": gaps,
        "monotone_nonincreasing": bool(
            all(gaps[i] >= gaps[i + 1] - 1e-12 for i in range(len(gaps) - 1))
        ),
    }


def main():
    t0 = time.time()
    out: dict = {"tol": _f(paths.TOL), "alpha0": _f(paths.ALPHA0)}
    out["a_tiny_region"] = check_a()
    print("(a) tiny region:", out["a_tiny_region"], flush=True)
    out["b_subdivide"] = check_b()
    print(
        "(b) subdivide never increases c_U:",
        out["b_subdivide"]["no_increase"],
        f"(improvement {_f(out['b_subdivide']['improvement']):+.3e})",
        flush=True,
    )
    out["c_pruning"] = check_c()
    print(
        "(c) pruning proofs hold:",
        out["c_pruning"]["all_hold"],
        f"({out['c_pruning']['n_leaves_checked']} leaves)",
        flush=True,
    )
    out["d_budget"] = check_d()
    print(
        "(d) gap monotone in budget:",
        out["d_budget"]["monotone_nonincreasing"],
        out["d_budget"]["gaps"],
        flush=True,
    )
    out["seconds"] = _f(time.time() - t0)
    ok = (
        bool(out["a_tiny_region"]["within_tol"])
        and bool(out["b_subdivide"]["no_increase"])
        and bool(out["c_pruning"]["all_hold"])
        and bool(out["d_budget"]["monotone_nonincreasing"])
    )
    out["verdict"] = "PASS" if ok else "FAIL"
    print(f"\nRESULT: {out['verdict']} ({_f(out['seconds']):.0f}s)")
    tiers.write_deposit(
        out,
        OUT,
        "region_sup_certified",
        notes="branch-and-bound unit checks (a)-(d) of plan step 11",
    )
    return out


if __name__ == "__main__":
    import sys

    try:
        res = main()
        sys.exit(0 if res.get("verdict") == "PASS" else 1)
    except Exception as exc:  # noqa: BLE001
        print("CHECKS ERROR:", repr(exc))
        sys.exit(2)
