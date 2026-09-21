"""Write release/expected.json: the canonical values and the deposit hashes.

The plan asks for "expected hashes or canonical formatted outputs" so a replicator can compare
without being told the answers in advance.  This script writes the canonical fields of the key
deposits (not their full contents) together with the SHA-256 of each deposit file, so a fresh
implementation can be checked field by field.

    python make_expected.py        -> release/expected.json
"""
from __future__ import annotations

import json
import os

import manifest
import paths
import tiers

FIELDS = {
    "bnb.json": ("n_cells", "n_certified_at_tol", "max_gap", "tol", "alpha0"),
    "bnb_checks.json": ("verdict",),
    "rootcert.json": ("n_brackets", "all_independent_check_pass", "max_width", "max_dp"),
    "adversarial.json": ("n_designs", "n_certified", "all_coverage_ok",
                         "all_witness_checks_pass", "max_gap"),
    "box_slack.json": ("slack_max", "slack_at_factor_1p001", "tol"),
    "error_budget.json": ("n_points", "max_Delta_p", "all_dp_under_target"),
    "enclose_iv.json": ("n_points", "all_contained", "max_convergence"),
    "impl_compare.json": ("verdict",),
    "rule_certified.json": ("n_families", "all_coverage_certified",
                            "all_under_runtime_target", "alpha0"),
    "region_shapes_cert.json": ("n_cells", "all_coverage_equal", "max_c_U_spread"),
    "desk_map_alpha0.json": ("n_cells", "n_certified_at_tol", "max_gap", "alpha0"),
    "gates.json": ("n_gates", "n_pass", "all_pass"),
    "enclose_falsifiers.json": ("verdict",),
}


def main():
    out = {"fields": {}, "hashes": {}, "note": "canonical fields + deposit hashes for a blind replica"}
    for name, keys in FIELDS.items():
        p = os.path.join(paths.DEPOSITS, name)
        if not os.path.exists(p):
            out["fields"][name] = None
            continue
        try:
            d = tiers.load_strict(p)
        except Exception as exc:  # noqa: BLE001
            out["fields"][name] = {"error": repr(exc)}
            continue
        out["fields"][name] = {k: d.get(k) for k in keys}
        out["hashes"][name] = manifest.sha256(p)
    tiers.write_deposit(out, os.path.join(paths.RELEASE, "expected.json"),
                        "region_sup_certified",
                        notes="expected canonical values and hashes; used by the blind replication kit")
    print(json.dumps({k: v for k, v in out["fields"].items() if v}, indent=2)[:1800])
    return out


if __name__ == "__main__":
    main()
