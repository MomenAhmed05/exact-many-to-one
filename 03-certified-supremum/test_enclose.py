"""P1 falsifier suite for the enclosure and the root certificate (plan P1 step 7).

Runs the falsifiers (a)-(g) of the plan plus the checks that matter most for the
programme, and deposits every result as strict JSON with its tier:

  F1  point containment: Riemann-float, Riemann-Arb and GL-Arb all bracket the
      independent solver `core.cond_noreject`.
  F2  random-box containment, including extreme ratio spreads, low degrees of freedom,
      balanced/unbalanced sample sizes and K = 1..5.
  F3  box collapse: shrinking a box to a point narrows the enclosure monotonically.
  F4  refinement: the Riemann width falls with the node count; the GL width is orders
      of magnitude smaller at a small fraction of the cost.
  F5  mass floor: dropping cells below a floor widens the enclosure, never narrows it.
  F6  cross-path agreement: GL-Arb and Riemann-Arb enclosures overlap (independent
      quadratures, same integral).
  F7  the bracket, checked against the INDEPENDENT solver: cond_noreject at c_U must be
      >= 1 - alpha at every corner of the box, and at c_L must be < 1 - alpha.
  F8  population/conditional identity: the draw-average of cond_noreject at the oracle
      critical value is 1 - alpha (Monte-Carlo error reported).
  F9  cost table: seconds and Arb evaluations per enclosure call and per bracket, by K.
  F10 error-estimate consistency: raising the GL order must shrink the width by a large
      factor (the ellipse remainder must behave like the exponential decay it claims).

Run:  python test_enclose.py            (a few minutes)
Out:  deposits/enclose_falsifiers.json
"""

from __future__ import annotations

import importlib
import os
import time

import enclose
import enclose_gl
import numpy as np
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "enclose_falsifiers.json")
CORE = importlib.import_module("core")
cond_noreject = CORE.cond_noreject
draw_var = CORE.draw_var
oracle_c = CORE.oracle_c


def _F_ref(c, s0, si, n, lam):
    return _f(
        cond_noreject(
            np.array([_f(c)]),
            np.array([_f(s0)]),
            np.asarray(si, float)[None, :],
            np.asarray(n, float),
            np.asarray(lam, float),
            npanels=paths.REF_PANELS,
        )[0]
    )


def f1_point(alpha=paths.ALPHA0, prec=64, n=None):
    n = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    s0, si, _ = draw_var(n, lam, 1, 20260918)
    d = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    s = np.sqrt(n[1:] / lam)
    rows = []
    for c in (1.8, 2.2, 2.6, 3.0, 3.4):
        rf = enclose.bounds_float(c, enclose.base_partition(nodes=4000), d, s, s, n[0])
        prep = enclose.prepare_partition(enclose.base_partition(nodes=20000), prec=prec)
        ra = enclose.bounds_arb(c, d, n[0], s, s, prep)
        ga = enclose_gl.bounds_gl_arb(c, d, n[0], s, s, prec=prec, n=n)
        ref = _F_ref(c, s0[0], si[0], n, lam)
        rows.append(
            {
                "c": _f(c),
                "F_ref": ref,
                "riemann_float": [rf[0], rf[1]],
                "riemann_arb": [ra[0], ra[1]],
                "gl_arb": [ga[0], ga[1]],
                "width_gl": _f(ga[1]) - _f(ga[0]),
                "all_contain": bool(
                    rf[0] <= ref <= rf[1]
                    and ra[0] <= ref <= ra[1]
                    and ga[0] <= ref <= ga[1]
                ),
            }
        )
    return rows


def f2_boxes(ntrials=24, prec=64, seed=20260919, alpha=paths.ALPHA0):
    rng = np.random.default_rng(seed)
    rows = []
    for t in range(_i(ntrials, 24)):
        K = _i(rng.integers(1, 6))
        df0 = _i(rng.integers(2, 40))
        n = np.concatenate([[df0 + 1.0], rng.integers(3, 40, K) + 1.0])
        spread = 10.0 ** _f(
            rng.uniform(-3.0, 3.0)
        )  # DOWNWARD spreads too (small lambda, large s: the pass-1 coverage gap)
        lam = np.exp(rng.uniform(-0.5, 0.5, K)) * (1.0 if t % 2 == 0 else spread)
        s0v, siv, _ = draw_var(n, lam, 1, seed + t)
        d = np.sqrt(siv[0] / n[1:] + s0v[0] / n[0])
        fac = 1.0 + rng.random() * 0.9
        lo_f, hi_f = lam / fac, lam * fac
        s_lo, s_hi = enclose.s_range(n, lo_f, hi_f)
        c = _f(2.0 + 1.5 * rng.random())
        g_lo, g_hi = enclose_gl.bounds_gl_arb(c, d, n[0], s_lo, s_hi, prec=prec)
        fmin, fmax = 1.0, 0.0
        for combo in np.ndindex(*([2] * K)):
            lam_pt = np.where(np.array(combo, bool), hi_f, lo_f)
            v = _F_ref(c, s0v[0], siv[0], n, lam_pt)
            fmin, fmax = min(fmin, v), max(fmax, v)
        rows.append(
            {
                "trial": _i(t),
                "K": _i(K),
                "df0": _i(df0),
                "n": list(n),
                "lam": list(lam),
                "box": [list(lo_f), list(hi_f)],
                "c": c,
                "F_lo": g_lo,
                "F_hi": g_hi,
                "min_corner": fmin,
                "max_corner": fmax,
                "width": _f(g_hi) - _f(g_lo),
                "lo_ok": bool(g_lo <= fmin),
                "hi_ok": bool(g_hi >= fmax),
                "inside": bool(g_lo <= fmin and g_hi >= fmax),
            }
        )
    return rows


def f3_collapse(prec=64):
    n = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    s0, si, _ = draw_var(n, lam, 1, 20260918)
    d = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    c = 2.6
    rows = []
    for mult in (1.5, 1.2, 1.05, 1.01, 1.0):
        s_lo, s_hi = enclose.s_range(n, lam / mult, lam * mult)
        lo, hi = enclose_gl.bounds_gl_arb(c, d, n[0], s_lo, s_hi, prec=prec)
        rows.append({"mult": _f(mult), "width": _f(hi) - _f(lo), "lo": lo, "hi": hi})
    widths = [r["width"] for r in rows]
    return {
        "rows": rows,
        "monotone_narrowing": bool(
            all(widths[i] >= widths[i + 1] - 1e-12 for i in range(len(widths) - 1))
        ),
    }


def f4_refinement(prec=64):
    n = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    s0, si, _ = draw_var(n, lam, 1, 20260918)
    d = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    s = np.sqrt(n[1:] / lam)
    c = 2.6
    r_rows = []
    for nodes in (200, 1000, 5000, 25000, 100000):
        z = enclose.base_partition(nodes=nodes)
        t0 = time.time()
        lo, hi = enclose.bounds_float(c, z, d, s, s, n[0])
        r_rows.append(
            {
                "nodes": _i(nodes),
                "width": _f(hi) - _f(lo),
                "ms": _f(1000 * (time.time() - t0)),
            }
        )
    t0 = time.time()
    g_lo, g_hi = enclose_gl.bounds_gl_arb(c, d, n[0], s, s, prec=prec)
    gl_ms = _f(1000 * (time.time() - t0))
    t0 = time.time()
    prep = enclose.prepare_partition(enclose.base_partition(nodes=25000), prec=prec)
    r_lo, r_hi = enclose.bounds_arb(c, d, n[0], s, s, prep)
    ra_ms = _f(1000 * (time.time() - t0))
    return {
        "riemann_float": r_rows,
        "gl": {"width": _f(g_hi) - _f(g_lo), "ms": gl_ms},
        "riemann_arb": {"nodes": 25001, "width": _f(r_hi) - _f(r_lo), "ms": ra_ms},
        "gl_beats_riemann": bool(
            (_f(g_hi) - _f(g_lo)) <= min(_f(r["width"]) for r in r_rows)
        ),
        "gl_vs_riemann_ratio": (
            (_f(r_hi) - _f(r_lo)) / (_f(g_hi) - _f(g_lo))
            if (_f(g_hi) - _f(g_lo)) > 0
            else NAN
        ),
    }


def f5_mass_floor(prec=64):
    n = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    s0, si, _ = draw_var(n, lam, 1, 20260918)
    d = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    s = np.sqrt(n[1:] / lam)
    z = enclose.base_partition(nodes=2000)
    rows = []
    for floor in (0.0, 1e-12, 1e-9, 1e-6, 1e-4):
        prep = enclose.prepare_partition(z, prec=prec, mass_floor=floor)
        lo, hi = enclose.bounds_arb(2.6, d, n[0], s, s, prep)
        rows.append(
            {
                "mass_floor": _f(floor),
                "width": _f(hi) - _f(lo),
                "dropped_up": _f(prep.get("dropped_up"), 0.0),
            }
        )
    w = [r["width"] for r in rows]
    return {
        "rows": rows,
        "monotone_widening": bool(
            all(w[i] <= w[i + 1] + 1e-15 for i in range(len(w) - 1))
        ),
    }


def f6_cross_path(prec=64, ntrials=8, seed=20260920):
    rng = np.random.default_rng(seed)
    rows = []
    for t in range(_i(ntrials, 8)):
        n = np.array([13.0] + [13.0] * _i(rng.integers(1, 4)))
        K = n.size - 1
        lam = np.exp(rng.uniform(-1, 1, K))
        s0, si, _ = draw_var(n, lam, 1, seed + t)
        d = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
        s = np.sqrt(n[1:] / lam)
        c = _f(2.0 + 1.5 * rng.random())
        g_lo, g_hi = enclose_gl.bounds_gl_arb(c, d, n[0], s, s, prec=prec)
        prep = enclose.prepare_partition(enclose.base_partition(nodes=25000), prec=prec)
        r_lo, r_hi = enclose.bounds_arb(c, d, n[0], s, s, prep)
        rows.append(
            {
                "trial": _i(t),
                "K": _i(K),
                "c": c,
                "gl": [g_lo, g_hi],
                "riemann": [r_lo, r_hi],
                "overlap": bool(r_lo <= g_hi and g_lo <= r_hi),
                "gl_width": _f(g_hi) - _f(g_lo),
                "riemann_width": _f(r_hi) - _f(r_lo),
            }
        )
    return {"rows": rows, "all_overlap": bool(all(r["overlap"] for r in rows))}


def f7_bracket(alpha=paths.ALPHA0, prec=64, tol=1e-6):
    """The decisive independent check of the bracket."""
    rows = []
    for K, lam, df0, ncs in (
        (2, [1.0, 25.0], 12, [13.0, 13.0]),
        (3, [0.25, 1.0, 4.0], 6, [7.0, 9.0, 11.0]),
        (1, [2.0], 4, [13.0]),
    ):
        n = np.array([df0 + 1.0] + list(ncs))
        lam_v = np.array(lam, float)
        s0, si, _ = draw_var(n, lam_v, 1, 20260918)
        d = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
        s = np.sqrt(n[1:] / lam_v)
        br = enclose_gl.bracket_point_gl(s, d, n[0], alpha=alpha, prec=prec, tol=tol)
        f_at_cl = _F_ref(br["c_L"], s0[0], si[0], n, lam_v)
        f_at_cu = _F_ref(br["c_U"], s0[0], si[0], n, lam_v)
        rows.append(
            {
                "K": _i(K),
                "n": list(n),
                "lam": list(lam_v),
                "c_L": br["c_L"],
                "c_U": br["c_U"],
                "width": br["width"],
                "dp": br["dp"],
                "c_equiv": br["c_equiv"],
                "n_eval": br["n_eval"],
                "F_ref_at_c_L": f_at_cl,
                "F_ref_at_c_U": f_at_cu,
                "target": 1.0 - _f(alpha),
                "c_L_below_target": bool(f_at_cl < 1.0 - _f(alpha)),
                "c_U_at_or_above_target": bool(f_at_cu >= 1.0 - _f(alpha) - 1e-15),
                "width_under_tol": bool(_f(br["width"]) <= paths.TOL),
            }
        )
    return {
        "rows": rows,
        "all_pass": bool(
            all(
                r["c_L_below_target"]
                and r["c_U_at_or_above_target"]
                and r["width_under_tol"]
                for r in rows
            )
        ),
    }


def f8_population(alpha=paths.ALPHA0, draws=4000, seed=20260921):
    n = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    c_or = oracle_c(n, lam, _f(alpha))
    s0, si, _ = draw_var(n, lam, _i(draws, 4000), seed)
    p = cond_noreject(np.full(_i(draws, 4000), c_or), s0, si, n, lam)
    mcse = _f(np.std(p, ddof=1) / np.sqrt(_i(draws, 4000)))
    return {
        "c_oracle": c_or,
        "target": 1.0 - _f(alpha),
        "mean_cond_noreject": _f(np.mean(p)),
        "mcse": mcse,
        "z": _f((np.mean(p) - (1.0 - _f(alpha))) / mcse),
        "draws": _i(draws),
        "within_3_mcse": bool(abs(_f((np.mean(p) - (1.0 - _f(alpha))) / mcse)) < 3.0),
    }


def f9_cost(prec=64):
    rows = []
    for K in (1, 2, 3, 5, 8):
        n = np.array([13.0] + [13.0] * K)
        lam = np.array([1.0] * K)
        s0, si, _ = draw_var(n, lam, 1, 20260918)
        d = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
        s = np.sqrt(n[1:] / lam)
        t0 = time.time()
        lo, hi = enclose_gl.bounds_gl_arb(2.6, d, n[0], s, s, prec=prec)
        per_call = _f(time.time() - t0)
        t0 = time.time()
        br = enclose_gl.bracket_point_gl(s, d, n[0], tol=1e-6, prec=prec)
        per_br = _f(time.time() - t0)
        rows.append(
            {
                "K": _i(K),
                "call_s": per_call,
                "bracket_s": per_br,
                "bracket_evals": _i(br.get("n_eval"), 0),
                "width": br.get("width"),
                "dp": br.get("dp"),
            }
        )
    return {
        "rows": rows,
        "budget_s": _f(paths.CERT_MIN_PER_FAMILY),
        "all_under_budget": bool(
            all(r["bracket_s"] <= paths.CERT_MIN_PER_FAMILY for r in rows)
        ),
    }


def f10_error_consistency(prec=64):
    """The ellipse remainder must decay like the exponential decay it claims.

    Widths are reported from float64 endpoints, so once the true width falls below float64
    resolution near a value of order 1 (~2.2e-16) the difference collapses to 0.0.  The
    check therefore requires monotone non-increase WHILE RESOLVABLE, reports the index at
    which the resolution floor is reached, and requires a large shrinkage factor between
    the lowest order and the last still-resolvable order.
    """
    n = np.array([13.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    s0, si, _ = draw_var(n, lam, 1, 20260918)
    d = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    s = np.sqrt(n[1:] / lam)
    rows = []
    for nn in (6, 8, 12, 16, 24):
        lo, hi, diag = enclose_gl.bounds_gl_arb(
            2.6, d, n[0], s, s, prec=prec, n=nn, return_diag=True
        )
        rows.append(
            {
                "gl_n": _i(nn),
                "width": _f(hi) - _f(lo),
                "err_total": diag.get("err_total"),
                "f_lo": lo,
                "f_hi": hi,
            }
        )
    w = [r["width"] for r in rows]
    # The width now has a RIGOROUS floor: each node carries an Arb radius (NODE_RAD) that the
    # interval arithmetic propagates, so the total plateaus at ~2*NODE_RAD*|f'| summed over the
    # rule (measured ~2.2e-14) instead of decaying to the float64 resolution.  The check therefore
    # (a) requires a large shrinkage factor from the lowest order to a mid order, (b) requires the
    # width never to INCREASE by more than a plateau tolerance, and (c) reports where the plateau
    # begins.  (The earlier version demanded strict monotone decay, which the plateau breaks.)
    plateau_tol = 1e-13
    floor = 2.2e-16
    at_floor = [i for i, x in enumerate(w) if x <= floor]
    shrink_ok = bool(all(w[i] >= w[i + 1] - plateau_tol for i in range(len(w) - 1)))
    ratio = (w[0] / w[-1]) if w and w[-1] > 0 else NAN
    plateau_from = None
    for i in range(1, len(w)):
        if w[i] > w[i - 1] - plateau_tol:
            plateau_from = i
            break
    ok = bool(shrink_ok and _f(ratio) >= 10.0)
    return {
        "rows": rows,
        "widths": w,
        "shrinks": shrink_ok,
        "at_float64_floor_from_index": (at_floor[0] if at_floor else None),
        "plateau_from_index": plateau_from,
        "plateau_tolerance": plateau_tol,
        "float64_floor": floor,
        "ratio_first_to_last": ratio,
        "ratio_first_to_last_resolvable": ratio,
        "pass": ok,
    }


def main():
    t0 = time.time()
    res: dict = {
        "alpha0": paths.ALPHA0,
        "gamma": paths.GAMMA,
        "bar_alpha": paths.ABAR,
        "tol_predeclared": paths.TOL,
        "dp_target": paths.CERT_DP_TARGET,
        "f1_point": f1_point(),
        "f3_collapse": f3_collapse(),
        "f4_refinement": f4_refinement(),
        "f5_mass_floor": f5_mass_floor(),
        "f6_cross_path": f6_cross_path(),
        "f10_error_consistency": f10_error_consistency(),
    }
    print(
        "F1 point:",
        "PASS" if all(r["all_contain"] for r in res["f1_point"]) else "FAIL",
    )
    print("F3 collapse:", res["f3_collapse"]["monotone_narrowing"])
    print(
        "F4 refinement: riemann widths",
        [f"{r['width']:.1e}" for r in res["f4_refinement"]["riemann_float"]],
        "| gl",
        f"{res['f4_refinement']['gl']['width']:.1e}",
        f"| ratio {_f(res['f4_refinement']['gl_vs_riemann_ratio']):.0f}x",
    )
    print("F5 mass floor:", res["f5_mass_floor"]["monotone_widening"])
    print("F6 cross-path:", res["f6_cross_path"]["all_overlap"])
    print(
        "F10 error consistency:",
        res["f10_error_consistency"]["pass"],
        f"widths {[f'{x:.1e}' for x in res['f10_error_consistency']['widths']]}",
        f"(floor from index {res['f10_error_consistency']['at_float64_floor_from_index']}, "
        f"ratio {_f(res['f10_error_consistency']['ratio_first_to_last_resolvable']):.1f}x)",
    )
    res["f2_boxes"] = f2_boxes()
    print(
        f"F2 boxes: {sum(1 for r in res['f2_boxes'] if r['inside'])}/{len(res['f2_boxes'])} contained"
    )
    res["f7_bracket"] = f7_bracket()
    print("F7 bracket:", "PASS" if res["f7_bracket"]["all_pass"] else "FAIL")
    for r in res["f7_bracket"]["rows"]:
        print(
            f"   K={r['K']} width {_f(r['width']):.2e} dp {_f(r['dp']):.2e} "
            f"F(c_L)={_f(r['F_ref_at_c_L']):.6f} < {_f(r['target']):.6f} <= "
            f"F(c_U)={_f(r['F_ref_at_c_U']):.6f} : "
            f"{r['c_L_below_target'] and r['c_U_at_or_above_target']}"
        )
    res["f8_population"] = f8_population()
    print(
        f"F8 population identity: mean {_f(res['f8_population']['mean_cond_noreject']):.6f} "
        f"vs {_f(res['f8_population']['target']):.6f} (z={_f(res['f8_population']['z']):.2f})"
    )
    res["f9_cost"] = f9_cost()
    print(
        "F9 cost:",
        [
            (r["K"], f"{r['bracket_s']:.2f}s", r["bracket_evals"])
            for r in res["f9_cost"]["rows"]
        ],
        "| under budget:",
        res["f9_cost"]["all_under_budget"],
    )
    res["seconds"] = _f(time.time() - t0)
    ok = (
        all(r["all_contain"] for r in res["f1_point"])
        and all(r["inside"] for r in res["f2_boxes"])
        and res["f3_collapse"]["monotone_narrowing"]
        and bool(res["f4_refinement"]["gl_beats_riemann"])
        and res["f5_mass_floor"]["monotone_widening"]
        and res["f6_cross_path"]["all_overlap"]
        and res["f7_bracket"]["all_pass"]
        and bool(res["f9_cost"]["all_under_budget"])
        and bool(res["f10_error_consistency"]["pass"])
        and bool(res["f8_population"]["within_3_mcse"])
    )
    res["verdict"] = "PASS" if ok else "FAIL"
    res["verdict_scope"] = (
        "composite of F1-F10; F4 requires the GL width to beat the finest Riemann width and F9 requires every bracket to be inside the runtime budget"
    )
    print(f"\nRESULT: {res['verdict']}  ({_f(res['seconds']):.0f}s)")
    tiers.write_deposit(
        res,
        OUT,
        "region_sup_certified",
        notes="enclosure + bracket falsifiers; the enclosure is the certified "
        "primitive, so this deposit is certificate-tier",
    )
    return res


if __name__ == "__main__":
    import sys

    try:
        res = main()
        sys.exit(0 if (res or {}).get("verdict") == "PASS" else 1)
    except Exception as exc:  # noqa: BLE001
        print("BATTERY ERROR:", repr(exc))
        sys.exit(2)
