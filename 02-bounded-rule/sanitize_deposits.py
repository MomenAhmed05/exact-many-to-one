"""Make every deposit strict-JSON: non-finite floats become null, with notes.

Bare NaN/Infinity is invalid strict JSON and was how a polluted deposit arose
earlier.  Each non-finite slot in this project falls into one of three benign
classes, and the sanitiser records which:

  * calib_vr/*.json  -- control_r2 for arms that use no control variates
                        (not applicable by construction);
  * w3/route1_*.json -- the oracle_cal rule row, computed only when the
                        benchmark group has n[0] <= 6 (calibrator cost), so it
                        is null at larger benchmark df rather than a failure;
  * ar2/ar1.json     -- bootstrap_boundary_phi_max, null when no phi falls inside
                        the 5% band.

Run:  python sanitize_deposits.py            (report only)
      python sanitize_deposits.py --apply    (rewrite the affected files)
"""

import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

NOTES = {
    "calib_vr/variance.json": {
        "field": "control_r2",
        "note": (
            "null where the arm uses no control variates (plain/antithetic); "
            "not applicable, not a failed measurement"
        ),
    },
    "calib_vr/accept.json": {
        "field": "control_r2",
        "note": (
            "null where the arm uses no control variates (plain); "
            "not applicable, not a failed measurement"
        ),
    },
    "w3/route1_df7_A.json": {"field": "oracle_cal", "note": None},
    "w3/route1_df7_B.json": {"field": "oracle_cal", "note": None},
    "w3/route1_df35_A.json": {"field": "oracle_cal", "note": None},
    "w3/route1_df35_B.json": {"field": "oracle_cal", "note": None},
    "ar2/ar1.json": {"field": "bootstrap_boundary_phi_max", "note": None},
}

W3_NOTE = (
    "oracle_cal is computed only when the benchmark group has n[0] <= 6 "
    "(calibrator cost); it is null at larger benchmark df by design, and "
    "the other rules in this file are unaffected"
)
AR2_NOTE = "null when no phi falls inside the 5% band (empty set), not a failure"


def _clean(obj):
    """(cleaned object, count of non-finite slots replaced)."""
    n = 0
    if isinstance(obj, float):
        if math.isfinite(obj):
            return obj, 0
        return None, 1
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            cv, c = _clean(v)
            out[k] = cv
            n += c
        return out, n
    if isinstance(obj, list):
        out = []
        for v in obj:
            cv, c = _clean(v)
            out.append(cv)
            n += c
        return out, n
    return obj, 0


def _add_notes(rel, obj):
    added = []
    spec = NOTES.get(rel.replace(os.sep, "/"))
    if not spec:
        return added
    field = spec["field"]
    if rel.endswith("w3/route1_df7_A.json") or "route1_df" in rel:
        obj.setdefault("notes", {})["oracle_cal"] = W3_NOTE
        added.append("notes.oracle_cal")
    elif rel.endswith("ar2/ar1.json"):
        obj.setdefault("notes", {})["bootstrap_boundary_phi_max"] = AR2_NOTE
        added.append("notes.bootstrap_boundary_phi_max")
    elif spec["note"]:
        obj.setdefault("notes", {})[field] = spec["note"]
        added.append(f"notes.{field}")
    return added


def main():
    apply = "--apply" in sys.argv
    total_files = total_vals = 0
    for p in sorted(
        glob.glob(os.path.join(HERE, "deposits", "**", "*.json"), recursive=True)
    ):
        rel = os.path.relpath(p, os.path.join(HERE, "deposits"))
        try:
            with open(p, encoding="utf-8") as fh:
                obj = json.load(fh)
        except Exception as exc:  # noqa: BLE001
            print(f"UNREADABLE {rel}: {exc!r}")
            continue
        cleaned, n = _clean(obj)
        if not n:
            continue
        total_files += 1
        total_vals += n
        notes = _add_notes(rel, cleaned)
        print(
            f"{rel}: {n} non-finite -> null"
            + (f" (+{','.join(notes)})" if notes else "")
        )
        if apply:
            try:
                with open(p, "w", encoding="utf-8") as fh:
                    json.dump(cleaned, fh, indent=2, allow_nan=False)
            except Exception as exc:  # noqa: BLE001
                raise SystemExit(f"cannot write {p}: {exc!r}") from exc
    print(
        f"{'APPLIED' if apply else 'DRY RUN'}: {total_files} files, {total_vals} values"
    )


if __name__ == "__main__":
    main()
