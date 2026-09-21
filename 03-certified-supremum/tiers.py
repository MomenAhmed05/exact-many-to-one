"""Decision tiers, atomic deposit stamping, strict reading, and the reserved-word gate.

This module is the enforcement layer for the plan's claim discipline (SPEC section 8) and
for deposit hygiene.  It was rewritten after the P0 verifier round found (and this file now
fixes) these defects:

* `write_deposit` could leave a TRUNCATED deposit when the serialisation failed mid-write
  (the verifier reproduced 25- and 10-byte corpses).  Writes are now atomic: build the
  payload in memory, write to a temporary file in the same directory, fsync, then
  `os.replace`.  A failed write leaves the previous file untouched.
* the reader-side "strict JSON" check used `json.load`, which ACCEPTS bare `NaN`/`Infinity`,
  so the check could not detect what it existed to detect.  `load_strict` now rejects
  non-finite constants via `parse_constant`.
* the reserved-word gate had false negatives ("we certify the envelope", because `certify`
  was also listed as a reference marker; "controls the FWER"; anything but the literal
  "controls fwer") and false positives (honest denials such as "not certified").  The gate
  is now clause-wise with three precise exemptions: the reserved token is *inside* a
  reference phrase, a negation precedes the token, or the token is part of a word that only
  contains it.  Overclaim words outside the three reserved terms (`proven`, `theorem`, ...)
  are reported separately as WARNINGS, not violations.
* `gate_file` had no call sites.  `gate_all` now scans every deposit and every markdown
  report in the programme and is called by `run_modes.quick`.
"""

from __future__ import annotations

import contextlib
import json
import math
import os
import re
import tempfile

import paths  # noqa: F401  (wires 02-bounded-rule/ and 03-certified-supremum/ onto sys.path)

TIERS = (
    "oracle_exact_law",
    "plugin_measured",
    "calibrated_measured",
    "envelope_measured",
    "region_sup_certified",
    "global_sup_certified",
)
CERTIFIED_TIERS = ("region_sup_certified", "global_sup_certified")

# --- reserved words (violations when used of a measured-tier object) -------------
RESERVED_PATTERNS = (
    ("guarantee", re.compile(r"\bguarantee(?:s|d)?\b", re.IGNORECASE)),
    ("certify", re.compile(r"\bcertif(?:y|ies|ied|ication)\b", re.IGNORECASE)),
    (
        "control_fwer",
        re.compile(
            r"\bcontrol(?:s|led)?\s+(?:the\s+)?(?:fwer|family-?wise(?:\s+error)?"
            r"|familywise(?:\s+error)?)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "control_fwer_of",
        re.compile(
            r"\b(?:fwer|family-?wise(?:\s+error)?|familywise(?:\s+error)?)\s+control\b",
            re.IGNORECASE,
        ),
    ),
)
# phrases that legitimately NAME the certified tier: a reserved token inside one is exempt
REFERENCE_PATTERN = re.compile(
    r"(?:region_sup_certified|global_sup_certified"
    r"|certified\s+(?:tier|constant|threshold|column|rule|envelope|supremum|value|critical"
    r"|object|programme|program|path|region|upper\s+bound|interval|bound)"
    r"|price\s+of\s+the\s+proof"
    r"|certification\s+tier"
    r"|certificate[s]?"
    r"|certify\s+mode"
    r"|certifying)",
    re.IGNORECASE,
)
NEGATION_PATTERN = re.compile(
    r"(?:\bnot\b|\bno\b|\bnever\b|\bwithout\b|\bcannot\b|\bcan't\b|\bisn't\b|\baren't\b"
    r"|\blacks?\b|\buncertified\b|\bun-?proven\b|\bfails?\s+to\b|\brather\s+than\b)",
    re.IGNORECASE,
)
# words that are strong claims but not on the reserved list: reported as warnings
WARN_PATTERNS = (
    ("proven", re.compile(r"\bprove[sdn]?\b|\bproof\b", re.IGNORECASE)),
    ("theorem", re.compile(r"\btheorem\b", re.IGNORECASE)),
)

# sentences that are ABOUT the gate / claims / defects (or quote an overclaim as a defect) are not
# themselves overclaims; and documentation records have a different register from result reports
META_PATTERN = re.compile(
    r"(?:gate|verifier|audit|defect|finding|overclaim|false\s+negative|false\s+positive"
    r"|reserved\s+word|counterexample|quoted|flag(?:s|ged)?|must\s+not\s+be\s+called"
    r"|tier\s+label|claim\s+discipline)",
    re.IGNORECASE,
)
RECORD_NAMES = (
    "SPEC.md",
    "PROGRESS.md",
    "ENVIRONMENT.md",
    "README.md",
    "PRIOR-ART.md",
    "ESTIMAND-DECISION.md",
    "GOVERNANCE.md",
    "METHOD.md",
    "FALSIFIERS.md",
)
NEG_WINDOW = 25


def _negated_before(clause, pos, window=NEG_WINDOW):
    """Is the token at `pos` directly negated?

    Pass-4 finding N4: a 45-character lookback exempted evasions such as
    'not that it matters, we guarantee FWER control'.  The window is now 25 characters and a
    comma, semicolon or full stop inside the window breaks the link (a new sub-clause).
    """
    pre = clause[max(0, pos - window) : pos]
    for sep in (",", ";", ".", "!"):
        if sep in pre:
            pre = pre.rsplit(sep, 1)[-1]
    return bool(NEGATION_PATTERN.search(pre))


_CLAUSE_SPLIT = re.compile(r"(?<=[.;:!?])\s+|\n+|,\s+(?=(?:unlike|whereas|while)\b)")


def _f(x, default=None):
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return float("nan") if default is None else default


def is_certified(tier):
    """True for tiers that carry the finite-sample theorem."""
    return str(tier) in CERTIFIED_TIERS


# --------------------------------------------------------------------------
# strict, atomic IO
# --------------------------------------------------------------------------
def _reject_constant(name):
    raise ValueError(f"non-finite JSON constant {name!r} (strict JSON forbids it)")


def _sanitize(obj, path="", nulled=None):
    """Recursively replace non-finite floats by None and numpy scalars by Python ones.

    Returns (clean_object, [paths_that_were_nulled]).  The P0 verifier showed that a
    stray NaN reaching `json.dump(..., allow_nan=False)` raises *after* partial output,
    which is how a truncated deposit appeared; sanitising first makes that impossible.
    """
    if nulled is None:
        nulled = []
    if isinstance(obj, dict):
        clean = {}
        for k, v in obj.items():
            clean[k], nulled = _sanitize(v, f"{path}.{k}" if path else str(k), nulled)
        return clean, nulled
    if isinstance(obj, (list, tuple)):
        items = []
        for i, v in enumerate(obj):
            cv, nulled = _sanitize(v, f"{path}[{i}]", nulled)
            items.append(cv)
        return items, nulled
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        try:
            obj = obj.item()
        except Exception:  # noqa: BLE001
            return obj, nulled
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            nulled.append(path)
            return None, nulled
        return obj, nulled
    return obj, nulled


def load_strict(path):
    """Load JSON, refusing bare `NaN`/`Infinity`/`-Infinity`."""
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh, parse_constant=_reject_constant)
    except ValueError:
        raise
    except OSError as exc:
        raise RuntimeError(f"could not read {path}") from exc


def dump_strict_atomic(obj, path):
    """Write strict JSON atomically: temp file in the same directory, then replace."""
    d = os.path.dirname(os.path.abspath(path))
    try:
        os.makedirs(d, exist_ok=True)
    except OSError as exc:
        raise RuntimeError(f"could not create deposit directory {d}") from exc
    fd, tmp = tempfile.mkstemp(prefix=".tmp-", suffix=".json", dir=d)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2, allow_nan=False)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except Exception:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise
    return path


def write_deposit(obj, path, tier, notes=None):
    """Stamp a deposit with its tier and write it atomically as strict JSON.

    `notes` is appended as a list under `notes` (the place for honest qualifications such
    as null-valued fields, which strict JSON forces us to state rather than emit as NaN).
    """
    if tier not in TIERS:
        raise ValueError(f"unknown tier {tier!r}; known: {TIERS}")
    out: dict[str, object] = dict(obj) if isinstance(obj, dict) else {"value": obj}
    out["tier"] = tier
    out["tier_claim"] = (
        "carries the alpha_0 + gamma finite-sample theorem"
        if is_certified(tier)
        else "measured; no finite-sample guarantee"
    )
    if notes is not None:
        prev = out.get("notes")
        add = list(notes) if isinstance(notes, list) else [str(notes)]
        if isinstance(prev, list):
            out["notes"] = prev + add
        elif prev:
            out["notes"] = [str(prev)] + add
        else:
            out["notes"] = add
    _clean, nulled = _sanitize(out)
    out = _clean if isinstance(_clean, dict) else {"value": _clean}
    if nulled:
        note = f"non-finite values replaced by null at: {', '.join(nulled[:20])}"
        prev = out.get("notes")
        if isinstance(prev, list):
            out["notes"] = [*prev, note]
        elif isinstance(prev, str):
            out["notes"] = [prev, note]
        else:
            out["notes"] = [note]
    return dump_strict_atomic(out, path)


# --------------------------------------------------------------------------
# the gate
# --------------------------------------------------------------------------
def _spans(pattern, text):
    return [(m.start(), m.end(), m.group(0)) for m in pattern.finditer(text)]


def scan_text(text):
    """All reserved occurrences and warnings, clause by clause.

    Each record: {kind: 'violation'|'warning'|'error', word, clause}.  A reserved occurrence is
    exempt when it lies INSIDE a reference phrase (which names the certified tier) or when a
    negation directly precedes it.  There is deliberately NO clause-level "meta" exemption: the
    pass-2 re-verifier showed that one produced six false negatives such as 'we certify the
    envelope, a finding of this audit'.  Document-level register (record vs report) is handled in
    `gate_file`.
    """
    out = []
    try:
        for clause in _CLAUSE_SPLIT.split(str(text)):
            cl = clause.strip()
            if not cl:
                continue
            # pass-4 finding N1: reference spans must be computed PER CLAUSE, because the
            # reserved-match offsets used below are clause-local.  Comparing clause-local
            # offsets against whole-text spans produced 3 false negatives and 4 false positives.
            refs = _spans(REFERENCE_PATTERN, cl)
            for word, pat in RESERVED_PATTERNS:
                for m in pat.finditer(cl):
                    if any(s <= m.start() and m.end() <= e for s, e, _t in refs):
                        continue  # inside a reference phrase
                    if _negated_before(cl, m.start()):
                        continue  # an honest denial
                    out.append({"kind": "violation", "word": word, "clause": cl[:400]})
            for word, pat in WARN_PATTERNS:
                if pat.search(cl) and not any(
                    s <= m.start() and m.end() <= e
                    for s, e, _t in refs
                    for m in pat.finditer(cl)
                ):
                    out.append({"kind": "warning", "word": word, "clause": cl[:400]})
    except Exception as exc:  # noqa: BLE001
        out.append({"kind": "error", "word": "scan", "clause": repr(exc)})
    return out


def is_record(path):
    """Documentation records (spec, progress, environment, phase/verifier reports) vs reports.

    Path-form independent (pass-4 finding N2: a relative 'verify/...' path was not recognised
    while an absolute one was), and PHASE*/REPORT* documents are records.
    """
    try:
        p = os.path.normpath(str(path)).replace("\\", "/")
    except Exception:  # noqa: BLE001
        p = str(path).replace("\\", "/")
    parts = [x for x in p.split("/") if x and x not in ("..", ".")]
    base = parts[-1] if parts else ""
    return bool(
        base in RECORD_NAMES
        or "verify" in parts
        or base.startswith(("PHASE", "REPORT"))
    )


def gate_text(text, tier):
    """Violations for a text belonging to `tier` (certified tiers are never gated)."""
    if is_certified(tier):
        return []
    return [r for r in scan_text(text) if r.get("kind") == "violation"]


_TIER_LINE = re.compile(
    r"^\s*tier\s*[:=]\s*([a-z_]+)\s*$", re.IGNORECASE | re.MULTILINE
)


def gate_file(path, tier=None):
    """Gate one markdown/text/JSON file.

    Tier comes from the argument or a standalone `tier: x` line.  An UNDECLARED documentation
    record (`SPEC.md`, `PROGRESS.md`, `ENVIRONMENT.md`, `README.md`, `PRIOR-ART.md`, `PHASE*`,
    `REPORT*`, anything under `verify/`) reports reserved words as `record_hits` rather than
    violations, because those documents legitimately name the certified tier and quote the claims
    the gate polices.  An explicit `tier:` line DEFEATS that demotion, so a report declaring a
    measured tier is still gated (pass-4 finding N3).
    """
    try:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    except Exception as exc:  # noqa: BLE001
        return {
            "path": path,
            "tier": tier,
            "error": repr(exc),
            "violations": [],
            "warnings": [],
            "record_hits": [],
        }
    m = _TIER_LINE.search(text)
    t = tier or (m.group(1).lower() if m else "envelope_measured")
    explicit = tier is not None or bool(m)
    recs = scan_text(text)
    record = is_record(path)
    violations = (
        [] if is_certified(t) else [r for r in recs if r.get("kind") == "violation"]
    )
    meta: list = []
    if record and not explicit:
        meta = violations
        violations = []
    return {
        "path": path,
        "tier": t,
        "record": record,
        "violations": violations,
        "record_hits": meta,
        "warnings": [r for r in recs if r.get("kind") == "warning"],
        "errors": [r for r in recs if r.get("kind") == "error"],
    }


def recheck_deposits(root=None, deep=False):
    """Deposit hygiene: strict JSON, valid tier, and (deep) the gate on string fields."""
    root = root or paths.DEPOSITS
    missing, invalid, nonstrit, violations, n = [], [], [], [], 0
    for dirpath, _dirs, files in os.walk(root):
        for fn in sorted(files):
            if not fn.endswith(".json"):
                continue
            n += 1
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, root)
            try:
                obj = load_strict(p)
            except ValueError as exc:
                nonstrit.append({"path": rel, "error": repr(exc)})
                continue
            except Exception as exc:  # noqa: BLE001
                invalid.append({"path": rel, "error": repr(exc)})
                continue
            tier = obj.get("tier")
            if tier is None:
                missing.append(rel)
            elif tier not in TIERS:
                invalid.append({"path": rel, "tier": str(tier)})
            if deep and not is_certified(tier):
                for v in gate_text(json.dumps(obj), tier):
                    violations.append({"path": rel, "tier": tier, **v})
    return {
        "n_deposits": n,
        "n_missing_tier": len(missing),
        "n_invalid": len(invalid),
        "n_non_strict": len(nonstrit),
        "n_violations": len(violations),
        "missing": missing[:50],
        "invalid": invalid[:50],
        "non_strict": nonstrit[:20],
        "violations": violations[:20],
        "pass": not (missing or invalid or nonstrit or violations),
    }


def gate_all(root=None, include_reports=True):
    """E4 gate over the whole programme: deposits (deep) and markdown reports."""
    root = root or paths.HERE
    rep_viol, rep_meta, rep_err = [], [], []
    if include_reports:
        for dirpath, _dirs, files in os.walk(root):
            for fn in sorted(files):
                if not fn.endswith(".md"):
                    continue
                p = os.path.join(dirpath, fn)
                g = gate_file(p)
                rel = os.path.relpath(p, root)
                if g["violations"]:
                    rep_viol.append(
                        {
                            "path": rel,
                            "tier": g["tier"],
                            "n": len(g["violations"]),
                            "first": g["violations"][0],
                        }
                    )
                if g.get("record_hits"):
                    rep_meta.append({"path": rel, "n": len(g["record_hits"])})
                if g.get("errors"):
                    rep_err.append({"path": rel, "errors": g["errors"]})
    dep = recheck_deposits(deep=True)
    return {
        "deposits": dep,
        "reports_with_violations": rep_viol,
        "records_with_meta_hits": rep_meta,
        "report_errors": rep_err,
        "pass": bool(dep["pass"] and not rep_viol and not rep_err),
    }


def smoke():
    """Self-check of the gate on adversarial cases (used by tests_gates.py)."""
    cases = [
        ("the envelope guarantees familywise error control", "envelope_measured", True),
        ("we certify the envelope", "envelope_measured", True),
        ("the rule controls the FWER", "envelope_measured", True),
        ("it controls family-wise error", "envelope_measured", True),
        (
            "this guarantees control of the familywise error rate",
            "envelope_measured",
            True,
        ),
        ("the certified tier removes the additive term", "envelope_measured", False),
        ("the certified constant is what to quote", "envelope_measured", False),
        ("this is not a guarantee of FWER control", "envelope_measured", False),
        ("no finite-sample guarantee is claimed here", "envelope_measured", False),
        ("the envelope is uncertified", "envelope_measured", False),
        ("the envelope is measured, not certified", "envelope_measured", False),
        (
            "the constant is guaranteed by the certified rule",
            "region_sup_certified",
            False,
        ),
    ]
    rows, bad = [], []
    for text, tier, should_flag in cases:
        got = bool(gate_text(text, tier))
        rows.append(
            {
                "text": text,
                "tier": tier,
                "flagged": got,
                "expected": should_flag,
                "ok": got == should_flag,
            }
        )
        if got != should_flag:
            bad.append(rows[-1])
    return {"rows": rows, "n_wrong": len(bad), "wrong": bad, "pass": not bad}


if __name__ == "__main__":
    print(
        json.dumps(
            {"smoke": smoke(), "deposits": recheck_deposits(deep=True)}, indent=2
        )[:2000]
    )
