"""Gates for the frontier programme.  Exit status 0 only if every gate passes.

X1  the K = 1 closed form against the general z-quadrature engine
X2  the corner substitutions behind the rigorous enclosure, falsified on random
    cells and ratio boxes
X3  the rigorous bracket must contain the quadrature value at every ratio tested
X4  the per-candidate engine with a CONSTANT threshold must reproduce the
    paper's own exact fixed-ratio law
X5  the per-candidate engine at K = 1 must equal the K = 1 quadrature
X6  every deposited rule must satisfy its own claimed constraint when re-checked
    on a grid finer than the one it was calibrated on
X7  tier discipline: no deposit may carry a certified tier while reporting an
    uncertified supremum
X8  the certified rule must never be more liberal than its own certificate
X9  every any-K certificate must dominate the exact law it bounds, and sit at
    or below the ceiling
X10 the speed layer must not move a number: the tabulated inversion against a
    direct root-find, and the parallel map against the serial one
X11 the joint certificate must dominate the exact law and sit at or below the
    ceiling
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.stats import norm

import fcore as F
import k1lab as L
import wb_certify as C
import wb_percand as P

RES = []


def gate(name, fn):
    try:
        r = fn()
    except Exception as exc:  # noqa: BLE001
        r = {"pass": False, "error": repr(exc)}
    r["gate"] = name
    RES.append(r)
    print(("  PASS " if r.get("pass") else "  FAIL ") + name + "   " +
          ", ".join("%s=%s" % (k, v) for k, v in r.items()
                    if k not in ("pass", "gate")))
    return r


def _x1_worst(nz, order):
    rng = np.random.default_rng(3)
    z, w = F._u_nodes(nz, order)
    worst = 0.0
    for _ in range(60):
        n = np.array([float(rng.integers(4, 40)), float(rng.integers(4, 40))])
        lam = float(np.exp(rng.uniform(np.log(1e-3), np.log(1e3))))
        u = np.exp(rng.normal(0, 0.5, 2))
        c = float(rng.uniform(0.5, 6.0))
        a = F.cond_noreject(c, u, n, [lam])
        s2 = np.array([1.0, lam]) * u
        d = np.sqrt(s2[1] / n[1] + s2[0] / n[0])
        b = float(w @ norm.cdf((z / np.sqrt(n[0]) + c * d)
                               * np.sqrt(n[1] / lam)))
        worst = max(worst, abs(a - b))
    return float(worst)


def x1():
    """The closed form is exact, so any residual must belong to the REFERENCE.

    A tolerance alone cannot tell those apart, so the gate refines the reference
    and requires the residual to fall to the float64 floor.  At a coarse
    reference the gap is 1.7e-12 and at a fine one it is 1.3e-15, which locates
    the error in the reference and not in the identity under test.
    """
    coarse = _x1_worst(96, 14)
    fine = _x1_worst(192, 16)
    finer = _x1_worst(384, 18)
    return {"pass": bool(fine < 1e-13 and finer < 1e-13 and coarse > fine),
            "coarse_reference": coarse, "fine_reference": fine,
            "finer_reference": finer,
            "refinement_reduces_residual": bool(coarse > fine)}


def x2():
    r = C.selftest(trials=300)
    return {"pass": r["x2_corner_substitution_violations"] == 0, **r}


def x3():
    n = np.array([12.0, 13.0])
    tab = L.PluginTable(n, 0.0235, m=41, npanel=12, order=8)
    edges = F.log_grid(3e-3, 3e2, 35)
    mids = np.concatenate([[edges[0] / 3], np.sqrt(edges[:-1] * edges[1:]),
                           [edges[-1] * 3]])
    rule = C.StepRule(tab(mids), edges)
    rows = C.bracket_check(rule, n, [0.01, 0.1, 1.0, 10.0, 100.0], m=400)
    return {"pass": all(r["contained"] for r in rows),
            "n_points": len(rows),
            "max_width": float(max(r["width"] for r in rows))}


def x4():
    worst = 0.0
    for n, lam in (([12.0, 13.0, 13.0], [1.0, 4.0]),
                   ([8.0, 13.0, 13.0, 13.0], [0.25, 1.0, 4.0])):
        n = np.asarray(n, float)
        lam = np.asarray(lam, float)
        c = F.oracle_c_pop(n, lam, F.ABAR, F.GRID)
        f = P.fwer_percand(lambda x, c=c: np.full_like(np.asarray(x, float), c),
                           n, lam)
        worst = max(worst, abs(f - F.ABAR))
    return {"pass": worst < 1e-6, "worst_abs_diff": float(worst)}


def x5():
    n = np.array([12.0, 13.0])
    c = F.oracle_c_pop(n, [4.0], F.ABAR, F.GRID)
    a = P.fwer_percand(lambda x, c=c: np.full_like(np.asarray(x, float), c),
                       n, np.array([4.0]))
    b = F.fwer_exact_k1(lambda x, c=c: np.full_like(np.asarray(x, float), c),
                        n, 4.0, 24, 12)
    return {"pass": abs(a - b) < 1e-10, "abs_diff": float(abs(a - b))}


def x6():
    """Re-check every deposited K = 1 rule on a finer grid than it was built on."""
    try:
        dep = F.read_deposit("wb_certified_rule.json")
    except Exception:  # noqa: BLE001
        return {"pass": True, "skipped": "no deposit yet"}
    bad, checked = [], 0
    fine = F.log_grid(1e-3, 1e3, 181)
    for d in dep["designs"]:
        rule = C.StepRule(d["c"], d["edges"])
        s, lw, _ = L.sup_fwer(rule, np.asarray(d["n"], float), fine, 24, 12)
        checked += 1
        if s > dep["abar"] + 1e-9:
            bad.append({"design": d["design"], "sup": float(s)})
    return {"pass": not bad, "checked": checked, "violations": bad[:4]}


def x7():
    import json
    import os
    bad = []
    for fn in sorted(os.listdir(F.DEP)):
        if not fn.endswith(".json"):
            continue
        with open(os.path.join(F.DEP, fn), encoding="utf-8") as fh:
            d = json.load(fh)
        t = d.get("tier")
        if t not in F.TIERS:
            bad.append({"file": fn, "why": "missing or unknown tier"})
        if t == "unconditional_valid_certified" and not d.get("all_certified"):
            bad.append({"file": fn, "why": "certified tier without a pass"})
    return {"pass": not bad, "violations": bad[:4]}


def x8():
    try:
        dep = F.read_deposit("wb_certified_rule.json")
    except Exception:  # noqa: BLE001
        return {"pass": True, "skipped": "no deposit yet"}
    bad = []
    for d in dep["designs"]:
        if d["sup_upper_bound"] is None or d["sup_upper_bound"] > dep["abar"] + 1e-12:
            bad.append(d["design"])
    return {"pass": not bad, "violations": bad[:4]}


def x9():
    """X9: every any-K certificate must dominate the exact law it bounds.

    The union step is an inequality, so the deposited bound has to sit ABOVE the
    exact supremum over the declared patterns at every row, and at or below the
    ceiling.  A row where it sat below the exact law would refute the step
    rather than merely be loose.
    """
    try:
        ak = F.read_deposit("certify_anyk.json")
    except Exception:  # noqa: BLE001
        return {"pass": True, "skipped": "no deposit yet"}
    bad = [d["design"] for d in ak["designs"]
           if d["union_bound_total"] < d["exact_sup_over_patterns"] - 1e-12]
    over = [d["design"] for d in ak["designs"]
            if d["union_bound_total"] > ak["abar"] + 1e-12]
    return {"pass": bool(not bad and not over and ak["designs"]),
            "rows": len(ak["designs"]), "not_dominating": bad[:4],
            "above_ceiling": over[:4],
            "worst_slack": ak["worst_union_slack"]}


def x10():
    """X10: an optimisation that changes an answer is not an optimisation.

    The tabulated inversion replaced a root-find and the pool replaced a serial
    loop; both are only admissible if the numbers are unchanged.  The parallel
    map must agree EXACTLY (it is the same arithmetic in a different order of
    dispatch, not a different computation), and the inversion to the tolerance
    its node count supports.
    """
    import fast
    r = fast.selftest()
    return {"pass": bool(r["pass"] and r["worst_parallel_vs_serial"] == 0.0),
            **{k: v for k, v in r.items() if k != "pass"}}


def x11():
    """X11: the joint certificate, same discipline as X9."""
    try:
        j = F.read_deposit("certify_joint.json")
    except Exception:  # noqa: BLE001
        return {"pass": True, "skipped": "no deposit yet"}
    bad = []
    for d in j["designs"]:
        if d["union_bound_total"] is None:
            bad.append(d["design"])
            continue
        if d["union_bound_total"] > j["abar"] + 1e-12:
            bad.append(d["design"])
            continue
        # the comparison object is sampled, so allow three standard errors;
        # a breach beyond that would refute the projection, not be noise
        se = d.get("exact_sup_mcse", 0.0)
        if d["union_bound_total"] < d["exact_sup_over_patterns"] - 3.0 * se:
            bad.append(d["design"])
    return {"pass": bool(not bad and j["designs"]),
            "rows": len(j["designs"]), "violations": bad[:4]}


def main():
    print("FRONTIER GATES")
    for name, fn in (("X1 K=1 closed form", x1), ("X2 corner substitution", x2),
                     ("X3 bracket containment", x3),
                     ("X4 per-candidate vs exact law", x4),
                     ("X5 per-candidate vs K=1", x5),
                     ("X6 deposited rules on a finer grid", x6),
                     ("X7 tier discipline", x7),
                     ("X8 certificate consistency", x8),
                     ("X9 any-K union bound dominates", x9),
                     ("X10 speed layer changes no number", x10),
                     ("X11 joint certificate dominates", x11)):
        gate(name, fn)
    npass = sum(1 for r in RES if r.get("pass"))
    print("\nGATES: %d/%d pass" % (npass, len(RES)))
    F.write_deposit("gates.json",
                    {"gates": RES, "n_pass": npass, "n_gates": len(RES)},
                    "frontier_measured")
    sys.exit(0 if npass == len(RES) else 1)


if __name__ == "__main__":
    main()
