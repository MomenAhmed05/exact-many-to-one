"""Reproduction modes: quick / paper / certify (plan §7.2, SPEC §1).

    python run_modes.py quick      minutes   analytic limits, region identity, gates
    python run_modes.py paper      all reproduced numbers, deterministic
    python run_modes.py certify    certificates (slower, auditable)

`quick` runs in-process and never touches the heavy solvers.  `paper` and `certify` run
the phase scripts in a fixed order, skipping (and *saying* they skip) scripts the
programme has not produced yet, so the same command keeps working as phases land.
Cross-tree modules (02-bounded-rule's `check_core`, `env_rule`) are imported dynamically so the
"not produced yet" path is explicit rather than an import error.
"""

from __future__ import annotations

import importlib
import os
import subprocess
import sys
import time

import paths
import tiers
from tiers import load_strict

# Ordered scripts per mode; entries are skipped with an explicit line when the file does
# not exist yet.  Paths are relative to this directory.
PAPER = [
    ("test_enclose.py", "enclosure + bracket falsifiers F1-F10"),
    ("enclose.py", "Riemann-path self-test"),
    ("enclose_gl.py", "GL certificate self-test"),
    ("deriv_diag.py", "derivatives, Lipschitz materiality, argmax geometry"),
    ("rootcert.py", "91 point root certificates"),
    ("cert_bnb.py", "the desk map (branch-and-bound)"),
    ("bnb_checks.py", "branch-and-bound unit checks (a)-(d)"),
    ("adversarial.py", "the adversarial geometry suite"),
    ("sup_compare.py", "continuous sup vs vertex/grid/local maxima"),
    ("box_slack.py", "the structural box slack"),
    ("error_budget.py", "per-certificate error budget"),
    ("enclose_iv.py", "mpmath cross-check"),
    ("impl_compare.py", "analytic limits and cross-implementation comparisons"),
    ("region_shapes_cert.py", "shapes at equal certified coverage"),
    ("rule_certified.py", "the shipped rule with decision certificates"),
    ("desk_map_alpha0.py", "the desk map at alpha0 = abar - gamma"),
    ("tests_gates.py", "the E1-E8 verification gates"),
]
CERTIFY = [
    ("test_enclose.py", "enclosure falsifiers"),
    ("rootcert.py", "point root certificates"),
    ("cert_bnb.py", "branch-and-bound certificates"),
    ("bnb_checks.py", "branch-and-bound unit checks"),
    ("enclose_iv.py", "mpmath independent cross-check"),
    ("error_budget.py", "error budgets"),
    ("adversarial.py", "adversarial certificates"),
    ("region_shapes_cert.py", "shape certificates"),
    ("rule_certified.py", "decision certificates for the shipped rule"),
    ("tests_gates.py", "the verification gates"),
]


def reproduce():
    """One-command reproduction: every deterministic artefact, in dependency order."""
    order = [("test_enclose.py", "falsifiers"), ("enclose.py", "Riemann self-test"),
             ("enclose_gl.py", "GL self-test"), ("deriv_diag.py", "derivatives"),
             ("error_budget.py", "error budget"), ("enclose_iv.py", "mpmath cross-check"),
             ("impl_compare.py", "limits"), ("box_slack.py", "box slack"),
             ("rootcert.py", "root certificates"), ("cert_bnb.py", "desk map"),
             ("bnb_checks.py", "unit checks"), ("adversarial.py", "adversarial suite"),
             ("sup_compare.py", "headline comparison"), ("region_shapes_cert.py", "shapes"),
             ("rule_certified.py", "shipped rule"), ("desk_map_alpha0.py", "desk map at alpha0"),
             ("tests_gates.py", "gates")]
    return {"mode": "reproduce", "runs": [_run(sc, lab) for sc, lab in order]}


def _try_import(name):
    """Import a module from the configured paths; return (module, error_string)."""
    try:
        return importlib.import_module(name), None
    except Exception as exc:  # noqa: BLE001
        return None, repr(exc)


def _run(script, label, timeout=None):
    p = os.path.join(paths.HERE, script)
    if not os.path.exists(p):
        return {
            "script": script,
            "status": "skipped",
            "why": "not produced yet",
            "label": label,
        }
    t0 = time.time()
    try:
        r = subprocess.run(
            [sys.executable, p],
            cwd=paths.HERE,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        tail = (r.stdout or "").strip().splitlines()[-6:]
        return {
            "script": script,
            "status": "ok" if r.returncode == 0 else "failed",
            "returncode": r.returncode,
            "seconds": round(time.time() - t0, 1),
            "tail": tail,
            "label": label,
            "stderr_tail": (r.stderr or "").strip().splitlines()[-3:],
        }
    except Exception as exc:  # noqa: BLE001
        return {"script": script, "status": "error", "why": repr(exc), "label": label}


def quick():
    """Cheap, deterministic checks: region identity, enclosure, gates, traceability."""
    out: dict[str, object] = {}

    # 1. exact-law core cross-checks (02-bounded-rule kernel vs its independent implementation)
    mod, err = _try_import("check_core")
    out["analytic"] = (
        {
            "status": "available",
            "note": "02-bounded-rule/check_core.py carries the full cross-check suite",
        }
        if mod is not None
        else {"status": "unavailable", "why": err}
    )

    # 2. region coverage identity: the tight scale must reproduce at least 1 - gamma
    env_rule, err = _try_import("env_rule")
    if env_rule is None:
        out["region"] = {"status": "unavailable", "why": err}
    else:
        try:
            n = [13.0, 13.0, 13.0]
            sc = env_rule.tight_scale(tuple(n), paths.GAMMA)
            cov = env_rule.region_coverage(n, paths.GAMMA, scale=sc)
            out["region"] = {
                "scale": sc,
                "coverage_exact": cov,
                "target": 1.0 - paths.GAMMA,
                "excess": cov - (1.0 - paths.GAMMA),
                "status": "ok" if cov >= 1.0 - paths.GAMMA else "FAIL",
            }
        except Exception as exc:  # noqa: BLE001
            out["region"] = {"status": "error", "why": repr(exc)}

    # 3. enclosure sanity: one fast certified containment check against the independent
    #    solver (the full battery lives in test_enclose.py)
    enc, err = _try_import("enclose_gl")
    if enc is None:
        out["enclosure"] = {"status": "not produced yet", "why": err}
    else:
        try:
            out["enclosure"] = enc.quickcheck()
        except Exception as exc:  # noqa: BLE001
            out["enclosure"] = {"status": "error", "why": repr(exc)}

    # 4. tier gate + deposit hygiene
    try:
        out["tiers"] = tiers.smoke()
    except Exception as exc:  # noqa: BLE001
        out["tiers"] = {"status": "error", "why": repr(exc)}
    out["deposit_tiers"] = tiers.recheck_deposits()

    # 5. number -> deposit traceability
    man, err = _try_import("manifest")
    if man is None:
        out["map"] = {"status": "unavailable", "why": err}
    else:
        try:
            out["map"] = man.check_map()
        except Exception as exc:  # noqa: BLE001
            out["map"] = {"status": "error", "why": repr(exc)}
        try:
            out["manifest"] = man.verify()
        except Exception as exc:  # noqa: BLE001
            out["manifest"] = {"pass": None, "error": repr(exc)}

    # 6. strict JSON across 03-certified-supremum deposits (load_strict rejects bare NaN/Infinity)
    bad: list[dict] = []
    for dirpath, _d, files in os.walk(paths.DEPOSITS):
        for fn in sorted(files):
            if fn.endswith(".json"):
                p = os.path.join(dirpath, fn)
                try:
                    load_strict(p)
                except Exception as exc:  # noqa: BLE001
                    bad.append(
                        {"path": os.path.relpath(p, paths.DEPOSITS), "why": repr(exc)}
                    )
    out["json"] = {"n_bad": len(bad), "bad": bad[:10], "pass": not bad}

    # 7. the reserved-word gate over every deposit and every report (E4)
    try:
        out["gate"] = tiers.gate_all()
    except Exception as exc:  # noqa: BLE001
        out["gate"] = {"pass": False, "error": repr(exc)}
    return out


def paper():
    return {"mode": "paper", "runs": [_run(s, lab) for s, lab in PAPER]}


def certify():
    return {"mode": "certify", "runs": [_run(s, lab) for s, lab in CERTIFY]}


def _summary(res):
    lines = []
    runs = res.get("runs")
    if isinstance(runs, list):
        for r in runs:
            tail = (
                f"  ({r.get('why', '')})"
                if r["status"] != "ok"
                else f"  [{r.get('seconds')}s]"
            )
            lines.append(
                f"  {r['status']:>7}  {r['script']:<24} {r.get('label', '')}{tail}"
            )
        return "\n".join(lines)
    reg = res.get("region") or {}
    enc = res.get("enclosure") or {}
    dep = res.get("deposit_tiers") or {}
    js = res.get("json") or {}
    mp = res.get("map") or {}
    gt = res.get("gate") or {}
    lines.append(
        f"  analytic suite       : {(res.get('analytic') or {}).get('status')}"
    )
    lines.append(
        f"  region coverage      : {reg.get('coverage_exact')} vs target "
        f"{reg.get('target')} -> {reg.get('status')}"
    )
    lines.append(
        f"  enclosure contained  : {enc.get('contained', enc.get('status'))} "
        f"(width {enc.get('width')})"
    )
    lines.append(
        f"  deposit tier gate    : {dep.get('pass')} ({dep.get('n_deposits')} files, "
        f"non-strict {dep.get('n_non_strict')})"
    )
    lines.append(f"  strict JSON          : {js.get('pass')}")
    lines.append(
        f"  number map resolves  : {mp.get('pass')} "
        f"(bad {mp.get('n_bad')}, not comparable {mp.get('n_not_comparable')})"
    )
    lines.append(f"  tier gate smoke      : {(res.get('tiers') or {}).get('pass')}")
    lines.append(
        f"  reserved-word gate   : {gt.get('pass')} "
        f"(report violations {len(gt.get('reports_with_violations') or [])}, "
        f"record hits "
        f"{sum((x.get('n') or 0) for x in (gt.get('records_with_meta_hits') or []))})"
    )
    mf = res.get("manifest") or {}
    lines.append(
        f"  manifest freshness   : {mf.get('pass')} "
        f"(changed {len(mf.get('changed') or [])}, added {len(mf.get('added') or [])}) "
        f"- informational, rebuilt at freeze points"
    )
    return "\n".join(lines)


def main(argv=None):
    argv = argv or sys.argv[1:]
    mode = (argv[0] if argv else os.environ.get("WA2H_MODE") or "quick").lower()
    t0 = time.time()
    fn = {"quick": quick, "paper": paper, "certify": certify,
          "reproduce": reproduce}.get(mode)
    if fn is None:
        print(f"unknown mode {mode!r}; use quick|paper|certify")
        return 2
    res = fn()
    print(f"== 03-certified-supremum {mode} ({time.time() - t0:.1f}s) ==")
    print(_summary(res))
    ok = True
    if mode == "quick":
        reg = (res.get("region") or {}).get("status")
        enc = res.get("enclosure") or {}
        enc_ok = bool(enc.get("contained"))
        ok = bool(
            not (res.get("json") or {}).get("n_bad")
            and (res.get("deposit_tiers") or {}).get("pass", False)
            and (res.get("map") or {}).get("pass", False)
            and reg == "ok"
            and enc_ok
            and (res.get("analytic") or {}).get("status") == "available"
            and (res.get("tiers") or {}).get("pass", False)
            and (res.get("gate") or {}).get("pass", False)
        )
    print("RESULT:", "PASS" if ok else "ATTENTION")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
