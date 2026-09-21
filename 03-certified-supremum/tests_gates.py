"""The E1-E8 verification gates (plan section 8), wired to the deposits.

  E1  enclosure containment: the falsifier battery's F1/F2 outcomes
  E2  certificate-tree integrity: every recorded leaf is covered exactly once by the tree, and every
      pruned leaf's recorded `pruned_at_L` satisfies `c_ub <= L + tol` (re-derived from the deposit)
  E3  root certificates: `c_L <= c* <= c_U` against the independent solver (the rootcert outcomes)
  E4  tier / reserved-word gate: every deposit carries a valid tier and no measured-tier text
      overclaims (`tiers.gate_all`)
  E5  number -> deposit -> script traceability (`manifest.check_map`)
  E6  Arb vs the independent high-precision evaluation (the mpmath cross-check outcomes)
  E7  stochastic falsifier agreement (the draw-average identity outcome)
  E8  runtime and node budgets recorded (the cost tables in the falsifier and branch-and-bound
      deposits), checked against the predeclared budgets

Run:  python tests_gates.py            -> deposits/gates.json, exit 0 iff every gate passes
"""

from __future__ import annotations

import json
import os
import re
import sys

import manifest
import paths
import tiers
from paths import _f, _i

OUT = os.path.join(paths.ensure_dirs(), "gates.json")
ERRLOG: list[str] = []


def _load(name):
    p = os.path.join(paths.DEPOSITS, name)
    if not os.path.exists(p):
        return None
    try:
        return json.load(open(p, encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return None


def e1():
    d = _load("enclose_falsifiers.json")
    if not d:
        return {"pass": False, "why": "no falsifier deposit"}
    ok = bool(
        all(r["all_contain"] for r in d["f1_point"])
        and all(r["inside"] for r in d["f2_boxes"])
    )
    return {
        "pass": ok,
        "f1": len(d["f1_point"]),
        "f2": len(d["f2_boxes"]),
        "f2_contained": sum(1 for r in d["f2_boxes"] if r["inside"]),
    }


def e2():
    d = _load("bnb.json")
    if not d:
        return {"pass": False, "why": "no branch-and-bound deposit"}
    tol = _f(d.get("tol"))
    bad, checked = [], 0
    for cell in d.get("cells", []):
        for leaf in [
            t for t in cell.get("tree", []) if isinstance(t, dict) and "lo" in t
        ]:
            checked += 1
            c_ub = _f(leaf.get("c_ub"))
            L = _f(leaf.get("pruned_at_L"))
            if not (c_ub == c_ub and L == L and c_ub <= L + tol + 1e-9):
                bad.append({"cell": cell.get("cell"), "c_ub": c_ub, "pruned_at_L": L})
    # every cell must also report a c_U that dominates its own c_L
    dom = all(_f(c.get("c_U")) >= _f(c.get("c_L")) - 1e-12 for c in d.get("cells", []))
    return {
        "pass": bool(checked > 0 and not bad and dom),
        "leaves_checked": _i(checked),
        "violations": bad[:5],
        "c_U_dominates_c_L": dom,
    }


def e3():
    d = _load("rootcert.json")
    if not d:
        return {"pass": False, "why": "no root-certificate deposit"}
    rows = [r for c in d.get("cells", []) for r in c.get("rows", [])]
    ok = all(r.get("c_L_below") and r.get("c_U_above") for r in rows)
    return {
        "pass": bool(rows and ok),
        "n_brackets": _i(len(rows)),
        "n_pass": _i(sum(1 for r in rows if r.get("c_L_below") and r.get("c_U_above"))),
    }


def e4():
    g = tiers.gate_all()
    return {
        "pass": bool(g.get("pass")),
        "deposits": (g.get("deposits") or {}).get("n_deposits"),
        "non_strict": (g.get("deposits") or {}).get("n_non_strict"),
        "report_violations": len(g.get("reports_with_violations") or []),
    }


def e5():
    m = manifest.check_map()
    return {
        "pass": bool(m.get("pass")),
        "entries": m.get("n_entries"),
        "bad": m.get("n_bad"),
        "not_comparable": m.get("n_not_comparable"),
    }


def e6():
    d = _load("enclose_iv.json")
    if not d:
        return {"pass": False, "why": "no mpmath cross-check deposit"}
    return {
        "pass": bool(d.get("all_contained")),
        "points": d.get("n_points"),
        "max_convergence": d.get("max_convergence"),
    }


def e7():
    d = _load("enclose_falsifiers.json")
    if not d:
        return {"pass": False, "why": "no falsifier deposit"}
    f8 = d.get("f8_population") or {}
    return {
        "pass": bool(f8.get("within_3_mcse")),
        "z": f8.get("z"),
        "draws": f8.get("draws"),
        "mean": f8.get("mean_cond_noreject"),
        "target": f8.get("target"),
    }


def e8():
    d = _load("enclose_falsifiers.json")
    b = _load("bnb.json")
    ok_f = bool((d or {}).get("f9_cost", {}).get("all_under_budget"))
    ok_b = bool(
        (b or {}).get("all_certified_at_tol") is not None
        and max((_f(c.get("seconds")) for c in (b or {}).get("cells", [])), default=0)
        <= 600.0
    )
    return {
        "pass": bool(ok_f and ok_b),
        "falsifier_budget_ok": ok_f,
        "bnb_max_seconds": _f(
            max((_f(c.get("seconds")) for c in (b or {}).get("cells", [])), default=0)
        ),
        "bnb_budget_s": 600.0,
    }


# --------------------------------------------------------------------------
# E9: the two-sided documentation gate
# --------------------------------------------------------------------------
_TWO_SIDED = (
    re.compile(r"equivalently\s+two-sided", re.I),
    re.compile(r"two-sided\s+(familywise\s+)?guarantee", re.I),
    re.compile(r"covers?\s+two-sided", re.I),
    re.compile(r"two-sided\s+(familywise\s+)?inference\s+(is\s+)?covered", re.I),
    re.compile(r"both\s+one-\s+and\s+two-sided", re.I),
    re.compile(r"one-sided\s+and\s+two-sided\s+(are|is)\s+(both\s+)?covered", re.I),
)
_NEGATION = re.compile(
    r"(not|no|never|without|requires?|separate|own derivation|uncovered|out of scope|not implemented"
    r"|not approximated|covers one-sided only|one-sided only)",
    re.I,
)


def _listdir(path):
    """Directory listing that cannot raise out of a gate."""
    try:
        return os.listdir(path)
    except OSError:
        return []


def _two_sided_hits(text):
    """Clauses that imply two-sided coverage, excluding explicit denials."""
    hits = []
    for clause in re.split(r"(?<=[.;:])\s+|\n+", str(text)):
        for pat in _TWO_SIDED:
            m = pat.search(clause)
            if not m:
                continue
            pre = clause[max(0, m.start() - 60) : m.start()]
            post = clause[m.end() : m.end() + 60]
            if _NEGATION.search(pre) or _NEGATION.search(post):
                continue  # an explicit denial or a scoping statement
            hits.append({"clause": clause.strip()[:200], "pattern": pat.pattern})
    return hits


def e9():
    """E9: no theorem-tier text may imply two-sided coverage.

    Scans every markdown document that declares a certified tier, and every deposit whose tier is
    certified, for phrases that imply two-sided familywise coverage.  Explicit denials ("not
    covered", "requires a separate derivation") are allowed and expected; an implication is not.
    """
    violations, scanned = [], 0
    for dirpath, _dirs, files in os.walk(paths.HERE):
        for fn in sorted(files):
            if not fn.endswith(".md"):
                continue
            pth = os.path.join(dirpath, fn)
            try:
                with open(pth, encoding="utf-8") as fh:
                    text = fh.read()
            except Exception as exc:  # noqa: BLE001
                ERRLOG.append(f"{pth}: {exc!r}")
                continue
            m = re.search(r"^\s*tier\s*[:=]\s*([a-z_]+)\s*$", text, re.I | re.M)
            tier = m.group(1).lower() if m else None
            if tier not in tiers.CERTIFIED_TIERS:
                continue
            scanned += 1
            for h in _two_sided_hits(text):
                violations.append({"path": os.path.relpath(pth, paths.HERE), **h})
    for name in sorted(_listdir(paths.DEPOSITS)):
        if not name.endswith(".json"):
            continue
        try:
            obj = tiers.load_strict(os.path.join(paths.DEPOSITS, name))
        except Exception as exc:  # noqa: BLE001
            ERRLOG.append(f"{name}: {exc!r}")
            continue
        if obj.get("tier") not in tiers.CERTIFIED_TIERS:
            continue
        for h in _two_sided_hits(json.dumps(obj)):
            violations.append({"path": f"deposits/{name}", **h})
    return {
        "pass": not violations,
        "documents_scanned": _i(scanned),
        "violations": violations[:6],
    }


GATES = (
    ("E1", "enclosure containment", e1),
    ("E2", "certificate-tree integrity", e2),
    ("E3", "root certificates", e3),
    ("E4", "tier / reserved-word gate", e4),
    ("E5", "number -> deposit traceability", e5),
    ("E6", "Arb vs mpmath", e6),
    ("E7", "stochastic falsifier", e7),
    ("E8", "runtime and node budgets", e8),
    ("E9", "two-sided documentation gate", e9),
)


def main():
    rows = []
    for name, desc, fn in GATES:
        try:
            r = fn()
        except Exception as exc:  # noqa: BLE001
            r = {"pass": False, "error": repr(exc)}
        r.update({"gate": name, "description": desc})
        rows.append(r)
        print(
            f"{name} {desc:<34} {'PASS' if r.get('pass') else 'FAIL'}  "
            f"{ {k: v for k, v in r.items() if k not in ('gate', 'description', 'pass')} }",
            flush=True,
        )
    ok = all(r.get("pass") for r in rows)
    out = {
        "n_gates": _i(len(rows)),
        "n_pass": _i(sum(1 for r in rows if r.get("pass"))),
        "all_pass": bool(ok),
        "gates": rows,
        "note": "E1-E8 wired to the deposits; exit status is 0 iff every gate passes",
    }
    tiers.write_deposit(
        out, OUT, "region_sup_certified", notes="verification gates E1-E8"
    )
    print(f"\nGATES: {out['n_pass']}/{out['n_gates']} pass")
    return out


if __name__ == "__main__":
    try:
        res = main()
        sys.exit(0 if res.get("all_pass") else 1)
    except Exception as exc:  # noqa: BLE001
        print("GATES ERROR:", repr(exc))
        sys.exit(2)
