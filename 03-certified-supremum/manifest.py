"""Hashes, staleness verification, and the number -> deposit -> script map (gate E5).

Rewritten after the P0 verifier round, which found that the first version could silently
drop a file whose hash failed, never verified a stored manifest against the tree, and never
compared a claimed number to the value actually deposited.  Now:

* `scan()` walks with a visited-realpath guard (junction cycles), records unreadable files
  explicitly in `failures` instead of dropping them, covers `.py/.md/.json/.txt/.csv` plus
  the `release/` artefacts, and can include the read-only `02-bounded-rule/` input tree.
* `verify()` recomputes the scan and reports changed / added / removed against the stored
  manifest, so a stale manifest cannot pass unnoticed.
* `check_map()` requires a concrete `field`, and (for numeric claims) compares the claimed
  number with the deposited value; non-numeric claims such as `12/12` are recorded as
  `not_comparable` rather than silently passing.
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys

import paths
import tiers
from paths import _i

MANIFEST = os.path.join(paths.RELEASE, "manifest.json")
NUMBER_MAP = os.path.join(paths.RELEASE, "number_map.json")
EXTS = (".py", ".md", ".json", ".txt", ".csv", ".yml", ".yaml")

# Seeded from the 02-bounded-rule dossier (REPORT-PROVEN-BOUNDARY.md, STAGE1-CERTIFICATE.md).
P0_SEED = [
    {
        "number": "12/12",
        "claim": "desk cells inside alpha + gamma (envelope); the desk sweep summary "
        "records all_pass=True over twelve rows",
        "deposit": "../02-bounded-rule/deposits/env/sweep_summary.json",
        "field": "all_pass",
        "script": "02-bounded-rule/w5_env.py",
        "report": "02-bounded-rule/REPORT-PROVEN-BOUNDARY.md",
        "tier": "envelope_measured",
    },
    {
        "number": "4.05e-07",
        "claim": "the K=8 pointwise artifact (quadrature level), carried forward "
        "rather than suppressed",
        "deposit": "../02-bounded-rule/deposits/env/w5env_tight_summary.json",
        "field": "conditional_violations.0.max_excess_vs_oracle_inside",
        "script": "02-bounded-rule/w5_env.py",
        "report": "02-bounded-rule/REPORT-PROVEN-BOUNDARY.md",
        "tier": "envelope_measured",
    },
    {
        "number": "0.0000",
        "claim": "worst corner deficit vs dense-grid sup (12/12 cells)",
        "deposit": "../02-bounded-rule/deposits/cert/sup_falsifier.json",
        "field": "worst_corner_deficit",
        "script": "02-bounded-rule/cert_falsify.py",
        "report": "02-bounded-rule/REPORT-PROVEN-BOUNDARY.md",
        "tier": "envelope_measured",
    },
    {
        "number": "+0.01180",
        "claim": "the one dominance violation (a gamma event)",
        "deposit": "../02-bounded-rule/deposits/cert/leak_probe.json",
        "field": "worst_excess_rej",
        "script": "02-bounded-rule/cert_leak.py",
        "report": "02-bounded-rule/REPORT-PROVEN-BOUNDARY.md",
        "tier": "envelope_measured",
    },
    {
        "number": "4.96e-04",
        "claim": "observed region miss rate (pooled 240k draws)",
        "deposit": "../02-bounded-rule/deposits/cert/miss_probe.json",
        "field": "pooled_miss_observed",
        "script": "02-bounded-rule/miss_probe.py",
        "report": "02-bounded-rule/REPORT-PROVEN-BOUNDARY.md",
        "tier": "envelope_measured",
    },
    {
        "number": "+422%",
        "claim": "max price of the proven (whole-region) tier",
        "deposit": "../02-bounded-rule/deposits/cert/certified_rule.json",
        "field": "max_price_of_proof_pct",
        "script": "02-bounded-rule/cert_rule.py",
        "report": "02-bounded-rule/REPORT-PROVEN-BOUNDARY.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "0.024964",
        "claim": "worst conditional rejection at the certified constant",
        "deposit": "../02-bounded-rule/deposits/cert/rigorous.json",
        "field": "worst_falsify_rej",
        "script": "02-bounded-rule/cert_sup.py",
        "report": "02-bounded-rule/STAGE1-CERTIFICATE.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "12/17",
        "claim": "desk cells certified at tol = 1e-4 by the branch-and-bound",
        "deposit": "deposits/bnb.json",
        "field": "n_certified_at_tol",
        "script": "cert_bnb.py",
        "report": "PHASE2-CERTIFIED-SUP.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "3.17e-03",
        "claim": "largest achieved branch-and-bound gap (family A df11)",
        "deposit": "deposits/bnb.json",
        "field": "max_gap",
        "script": "cert_bnb.py",
        "report": "PHASE2-CERTIFIED-SUP.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "91/91",
        "claim": "point root certificates passing the independent check",
        "deposit": "deposits/rootcert.json",
        "field": "n_brackets",
        "script": "rootcert.py",
        "report": "PHASE2-CERTIFIED-SUP.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "16/20",
        "claim": "adversarial designs certified at tol",
        "deposit": "deposits/adversarial.json",
        "field": "n_certified",
        "script": "adversarial.py",
        "report": "PHASE2-CERTIFIED-SUP.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "2.88e-02",
        "claim": "structural box slack at a factor-3 box (c-units)",
        "deposit": "deposits/box_slack.json",
        "field": "slack_max",
        "script": "box_slack.py",
        "report": "PHASE2-CERTIFIED-SUP.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "8.38e-14",
        "claim": "largest measured Delta p across the error-budget points",
        "deposit": "deposits/error_budget.json",
        "field": "max_Delta_p",
        "script": "error_budget.py",
        "report": "PHASE3-NUMERICAL-CORE.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "5/5",
        "claim": "mpmath cross-check points contained in the Arb enclosure",
        "deposit": "deposits/enclose_iv.json",
        "field": "n_points",
        "script": "enclose_iv.py",
        "report": "PHASE3-NUMERICAL-CORE.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "9",
        "claim": "verification gates defined (the PASS COUNT is read from the gate's own run; mapping it here would be self-referential)",
        "deposit": "deposits/gates.json",
        "field": "n_gates",
        "script": "tests_gates.py",
        "report": "REPORT.md",
        "tier": "region_sup_certified",
    },
    {
        "number": "0.0245",
        "claim": "the hard-ceiling internal level alpha0 = abar - gamma",
        "deposit": "deposits/rule_certified.json",
        "field": "alpha0",
        "script": "rule_certified.py",
        "report": "PHASE4-RULE.md",
        "tier": "region_sup_certified",
    },
]

_NUM_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def sha256(path):
    """SHA-256 of a file, or None if it cannot be read."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:  # noqa: BLE001
        return None


def _git_rev():
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=paths.HERE,
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


def scan(root, exts=EXTS, respect_skip=True):
    """{relpath: sha256} plus a `failures` list; cycle-guarded walk."""
    out: dict = {}
    failures: list = []
    seen = set()
    for dirpath, dirs, files in os.walk(root, followlinks=False):
        try:
            real = os.path.realpath(dirpath)
        except OSError:
            real = dirpath
        if real in seen:
            dirs[:] = []
            continue
        seen.add(real)
        if respect_skip:
            dirs[:] = [
                d for d in dirs if d not in {"__pycache__", ".git", ".pytest_cache"}
            ]
        for fn in sorted(files):
            if not fn.endswith(exts):
                continue
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, root).replace("\\", "/")
            h = sha256(p)
            if h is None:
                failures.append(rel)
            else:
                out[rel] = h
    return out, failures


def build(save=True, include_w26=False):
    pins: dict[str, object] = {"python": sys.version.split()[0]}
    for name, mod in (
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("mpmath", "mpmath"),
        ("python-flint", "flint"),
    ):
        try:
            pins[name] = __import__(mod).__version__
        except Exception:  # noqa: BLE001
            pins[name] = None
    scripts, s_fail = scan(paths.HERE, (".py",))
    reports, r_fail = scan(paths.HERE, (".md",))
    deposits, d_fail = scan(paths.DEPOSITS, (".json",))
    obj = {
        "git_rev": _git_rev(),
        "interpreter": sys.executable,
        "pins": pins,
        "scripts": scripts,
        "reports": reports,
        "deposits": deposits,
        "failures": {"scripts": s_fail, "reports": r_fail, "deposits": d_fail},
        "counts": {
            "scripts": len(scripts),
            "reports": len(reports),
            "deposits": len(deposits),
            "unreadable": len(s_fail) + len(r_fail) + len(d_fail),
        },
    }
    if include_w26:
        w26, w_fail = scan(paths.W26, (".py", ".md"))
        obj["w26_inputs"] = w26
        obj["failures"]["w26"] = w_fail
        obj["counts"]["w26_inputs"] = len(w26)
    try:
        import datetime

        obj["generated_at"] = datetime.datetime.now().isoformat(timespec="seconds")
    except Exception as exc:  # noqa: BLE001
        obj["generated_at"] = None
        obj.setdefault("notes", []).append(f"timestamp unavailable: {exc!r}")
    if save:
        tiers.dump_strict_atomic(obj, MANIFEST)
    return obj


def verify(path=MANIFEST):
    """Compare the stored manifest with the current tree: changed/added/removed."""
    try:
        old = tiers.load_strict(path)
    except Exception as exc:  # noqa: BLE001
        return {"path": path, "error": repr(exc), "pass": False}
    now = build(save=False)
    res = {"pass": True, "changed": [], "added": [], "removed": [], "unreadable": []}
    for key in ("scripts", "reports", "deposits"):
        o, n = old.get(key, {}), now.get(key, {})
        res["changed"] += [
            f"{key}/{k}" for k in sorted(set(o) & set(n)) if o[k] != n[k]
        ]
        res["added"] += [f"{key}/{k}" for k in sorted(set(n) - set(o))]
        res["removed"] += [f"{key}/{k}" for k in sorted(set(o) - set(n))]
    for group, fails in (now.get("failures") or {}).items():
        res["unreadable"] += [f"{group}/{f}" for f in fails]
    res["pass"] = not (
        res["changed"] or res["added"] or res["removed"] or res["unreadable"]
    )
    res["counts"] = now.get("counts")
    res["generated_at"] = old.get("generated_at")
    return res


def write_map(entries=None):
    ents = P0_SEED if entries is None else entries
    tiers.dump_strict_atomic(
        {
            "entries": ents,
            "n": len(ents),
            "note": "number -> deposit -> script traceability (gate E5)",
        },
        NUMBER_MAP,
    )
    return ents


def _resolve(deposit):
    p = str(deposit or "")
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(paths.HERE, p))


def _dig(obj, field):
    cur = obj
    for part in str(field).split("."):
        if isinstance(cur, list):
            try:
                cur = cur[_i(part)]
            except Exception:  # noqa: BLE001
                return False, None
        elif isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            return False, None
    return True, cur


def _numbers_match(claim, value):
    """(comparable, ok) for a claimed number string against a deposited value."""
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        return False, True
    m = _NUM_RE.search(str(claim))
    if not m:
        return False, True
    try:
        cv = float(m.group(0))
    except Exception:  # noqa: BLE001
        return False, True
    try:
        vv = float(value)
    except Exception:  # noqa: BLE001
        return False, True
    if vv != vv:
        return False, True
    scale = 1.0
    if "%" in str(claim):
        scale = 1.0 if abs(vv) > 5 else 100.0
    return True, bool(abs(cv - vv * scale) <= max(5e-4, 0.05 * abs(cv)))


def check_map(entries=None, path=NUMBER_MAP):
    """E5 gate: every mapped number resolves to a concrete deposit field and value."""
    if entries is None:
        try:
            entries = tiers.load_strict(path).get("entries", [])
        except Exception:  # noqa: BLE001
            entries = P0_SEED
    bad, not_comparable = [], []
    for e in entries:
        deposit, field = e.get("deposit", ""), e.get("field")
        if field in (None, "", "-"):
            bad.append(
                {
                    "number": e.get("number"),
                    "problem": "no field named",
                    "deposit": deposit,
                }
            )
            continue
        p = _resolve(deposit)
        if not os.path.exists(p):
            bad.append(
                {"number": e.get("number"), "problem": "missing deposit", "path": p}
            )
            continue
        try:
            obj = tiers.load_strict(p)
        except Exception as exc:  # noqa: BLE001
            bad.append({"number": e.get("number"), "problem": repr(exc), "path": p})
            continue
        ok, val = _dig(obj, field)
        if not ok:
            bad.append(
                {
                    "number": e.get("number"),
                    "problem": "missing field",
                    "field": field,
                    "path": p,
                }
            )
            continue
        comparable, match = _numbers_match(e.get("number"), val)
        if not comparable:
            not_comparable.append(
                {"number": e.get("number"), "field": field, "value": val}
            )
        elif not match:
            bad.append(
                {
                    "number": e.get("number"),
                    "problem": "value mismatch",
                    "field": field,
                    "deposited": val,
                }
            )
    return {
        "n_entries": len(entries),
        "n_bad": len(bad),
        "n_not_comparable": len(not_comparable),
        "bad": bad,
        "not_comparable": not_comparable[:20],
        "pass": not bad,
    }


if __name__ == "__main__":
    import json

    m = build(include_w26=False)
    write_map()
    v = verify()
    print(
        json.dumps(
            {
                "counts": m["counts"],
                "git_rev": m["git_rev"],
                "map": check_map(),
                "verify_next_run_will_report_changes": len(v.get("changed", []))
                + len(v.get("added", [])),
            },
            indent=2,
        )[:1500]
    )
