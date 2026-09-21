"""W4: the promotion-gate verifier.

The familywise guarantee is a statement about a procedure that was fixed before
anyone looked at the outcome.  Nothing in the mathematics can check that; it is
a fact about a workflow.  This module makes it a fact about a FILE.

`check(manifest)` returns a verdict with one of three labels:

    GUARANTEED   every condition the certificate needs is recorded and
                 internally consistent;
    CONDITIONAL  the procedure is valid given something the manifest declares
                 but cannot itself evidence, and the verdict names it;
    SCREEN       at least one condition is violated.  The test may still be
                 worth running, but it is a screen and must be reported as one.

A run refuses on: a missing required field, a family locked after the evaluation
window closed, a ratio domain inferred from the returns under test, a
calibration window overlapping the evaluation window, a dependence tier of C
presented as a certificate, or a threshold table whose hash does not match.

The refusals are deliberately unforgiving, because the failure mode this guards
against is not dishonesty, it is a desk remembering the order of events
generously six months later.

Usage:
    python w4_manifest.py example        # write and check a valid example
    python w4_manifest.py check FILE     # check a manifest
Deposit: deposits/w4_manifest.json
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import io
import json
import os
import sys
import time

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA = os.path.join(HERE, "manifest.schema.json")

REQUIRED_TOP = ("manifest_version", "locked_at", "locked_by", "benchmark",
                "candidates", "returns", "window", "hypothesis", "dependence",
                "nuisance_domain", "threshold_table", "specifications")

CIRCULAR_SOURCE = "same_returns_being_tested"


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _date(s):
    return _dt.date.fromisoformat(s)


def _dtime(s):
    return _dt.datetime.fromisoformat(s.replace("Z", "+00:00"))


def check(man, table_path=None):
    """Return {verdict, refusals, conditions, notes}."""
    refuse, cond, notes = [], [], []

    # ---- structural -----------------------------------------------------
    for k in REQUIRED_TOP:
        if k not in man:
            refuse.append("missing required field: %s" % k)
    if refuse:
        return {"verdict": "SCREEN", "refusals": refuse, "conditions": cond,
                "notes": notes}

    if man.get("manifest_version") != 1:
        refuse.append("unsupported manifest_version %r"
                      % man.get("manifest_version"))

    w = man["window"]
    try:
        locked = _dtime(man["locked_at"]).date()
        ev0, ev1 = _date(w["evaluation_start"]), _date(w["evaluation_end"])
    except Exception as exc:                      # noqa: BLE001
        return {"verdict": "SCREEN",
                "refusals": ["unparseable dates: %s" % exc],
                "conditions": cond, "notes": notes}

    # ---- the ordering the guarantee depends on --------------------------
    if locked > ev1:
        refuse.append(
            "family locked %s, after the evaluation window closed %s: the "
            "candidate set could have been chosen on the outcome"
            % (locked, ev1))
    elif locked > ev0:
        cond.append(
            "family locked %s, inside the evaluation window (%s to %s): valid "
            "only for the part of the window after the lock" % (locked, ev0, ev1))

    if "research_end" in w and _date(w["research_end"]) > ev0:
        refuse.append(
            "research window ends %s, after the evaluation window opens %s: "
            "hypotheses were generated on evaluation data"
            % (w["research_end"], ev0))

    # ---- the ratio domain must not be circular --------------------------
    nd = man["nuisance_domain"]
    if nd.get("source") == CIRCULAR_SOURCE:
        refuse.append(
            "variance-ratio domain inferred from the same returns being "
            "tested: the domain is then a function of the data and the "
            "uniform bound over it is not a bound")
    lo, hi = nd["variance_ratio"]
    if not (0 < lo < hi):
        refuse.append("variance_ratio domain %r is not a positive interval"
                      % (nd["variance_ratio"],))

    # ---- dependence tiers ------------------------------------------------
    dep = man["dependence"]
    for key in ("serial", "cross_sectional"):
        tier = dep[key]["tier"]
        if tier == "C_same_period":
            cond.append(
                "%s nuisance estimated on the evaluation period: the result is "
                "a plug-in, not a certificate, and must be reported as such"
                % key)
        if tier == "B_calibration_period":
            if "calibration_end" not in w or "calibration_start" not in w:
                refuse.append(
                    "%s tier is B but no calibration window is declared" % key)
            else:
                c0, c1 = _date(w["calibration_start"]), _date(w["calibration_end"])
                if c1 > ev0:
                    refuse.append(
                        "%s calibration window ends %s, overlapping the "
                        "evaluation window that opens %s" % (key, c1, ev0))
                else:
                    cond.append(
                        "%s nuisance frozen from %s to %s: valid conditionally "
                        "on that estimate, and the transform's error is not in "
                        "the certificate" % (key, c0, c1))

    if dep["serial"]["tier"] in ("A_declared", "B_calibration_period"):
        if dep["serial"].get("transform") == "prais_winsten":
            refuse.append(
                "serial transform is prais_winsten: it rescales the first "
                "observation, so the whitened mean depends on the group's "
                "length and the rule drifts off its ceiling at unequal "
                "lengths.  Use cochrane_orcutt")

    if dep["cross_sectional"]["model"] == "independent_groups":
        cond.append(
            "cross-sectional model is independent_groups: this asserts the "
            "benchmark and candidates do not share dates or shocks, which for "
            "a panel observed on common dates is usually false")
    if (dep["cross_sectional"]["model"] in ("observed_factor", "latent_factor")
            and man["hypothesis"].get("estimand") == "raw_mean_difference"
            and dep["cross_sectional"]["model"] == "observed_factor"):
        cond.append(
            "residualising against an observed factor changes the estimand to "
            "a factor alpha; estimand is declared as raw_mean_difference")

    # ---- the threshold table --------------------------------------------
    tt = man["threshold_table"]
    if tt["alpha0"] > man["hypothesis"]["fwer_target"]:
        refuse.append(
            "internal level %.5f exceeds the familywise target %.5f"
            % (tt["alpha0"], man["hypothesis"]["fwer_target"]))
    if table_path is not None:
        got = sha256_of(table_path)
        if got != tt["sha256"]:
            refuse.append("threshold table hash mismatch: manifest says %s, "
                          "file is %s" % (tt["sha256"][:12], got[:12]))
        else:
            notes.append("threshold table hash verified")

    # ---- family integrity ------------------------------------------------
    ids = [c["id"] for c in man["candidates"]]
    if len(set(ids)) != len(ids):
        refuse.append("duplicate candidate ids in the declared family")
    if man["benchmark"]["id"] in ids:
        refuse.append("the benchmark is also listed as a candidate")
    for a in man.get("amendments", []):
        if _dtime(a["at"]).date() > ev0:
            cond.append("manifest amended %s, inside the evaluation window: %s"
                        % (a["at"][:10], a["change"]))

    verdict = "SCREEN" if refuse else ("CONDITIONAL" if cond else "GUARANTEED")
    return {"verdict": verdict, "refusals": refuse, "conditions": cond,
            "notes": notes, "n_candidates": len(ids)}


def example_manifest(alpha0=0.02377):
    return {
        "manifest_version": 1,
        "locked_at": "2026-01-05T09:00:00+00:00",
        "locked_by": "desk-research",
        "note": "Example manifest; the shape a real one takes.",
        "benchmark": {"id": "live_book_v7",
                      "code_sha256": "a" * 64,
                      "description": "the deployed multi-sleeve book"},
        "candidates": [
            {"id": "carry_v3", "code_sha256": "b" * 64},
            {"id": "trend_v2", "code_sha256": "c" * 64},
            {"id": "fx_fix_v1", "code_sha256": "d" * 64},
        ],
        "returns": {
            "frequency": "monthly",
            "cost_model": {"commission_bp": 2.0, "slippage_bp": 4.0,
                           "financing": "SOFR + 90bp on gross",
                           "source": "broker schedule 2025-11"},
            "rebalance": "monthly, last business day, close",
            "missing_data_policy": "reject_run",
        },
        "window": {
            "research_end": "2025-12-31",
            "calibration_start": "2021-01-01",
            "calibration_end": "2025-12-31",
            "evaluation_start": "2026-01-01",
            "evaluation_end": "2028-12-31",
        },
        "hypothesis": {"direction": "one_sided_greater", "fwer_target": 0.025,
                       "statistic": "welch_separate_variance",
                       "estimand": "raw_mean_difference"},
        "dependence": {
            "serial": {"tier": "A_declared", "transform": "cochrane_orcutt",
                       "phi_max": 0.2,
                       "source": "declared bound, whitened at the top"},
            "cross_sectional": {"tier": "none",
                                "model": "independent_groups"},
        },
        "nuisance_domain": {
            "variance_ratio": [0.001, 1000.0],
            "source": "published_volatility_target",
            "evidence": "book targets 8% annualised; sleeve mandates cap 25%",
        },
        "threshold_table": {"version": "wb_certified_rule/df35",
                            "sha256": "e" * 64, "alpha0": alpha0,
                            "certificate_deposit": "wb_certified_rule.json"},
        "specifications": {
            "primary": "monthly net returns, Welch, one-sided, FWER 0.025",
            "sensitivity": ["gross of costs", "weekly frequency"],
        },
    }


BAD_CASES = [
    ("locked after the window closed",
     lambda m: m.update({"locked_at": "2029-06-01T00:00:00+00:00"})),
    ("ratio domain from the tested returns",
     lambda m: m["nuisance_domain"].update({"source": CIRCULAR_SOURCE})),
    ("research window overlaps evaluation",
     lambda m: m["window"].update({"research_end": "2026-06-30"})),
    ("calibration overlaps evaluation",
     lambda m: (m["dependence"]["serial"].update(
         {"tier": "B_calibration_period"}),
         m["window"].update({"calibration_end": "2026-06-30"}))),
    ("prais-winsten whitening",
     lambda m: m["dependence"]["serial"].update(
         {"transform": "prais_winsten"})),
    ("internal level above the target",
     lambda m: m["threshold_table"].update({"alpha0": 0.05})),
    ("benchmark listed as a candidate",
     lambda m: m["candidates"].append(
         {"id": "live_book_v7", "code_sha256": "f" * 64})),
    ("missing nuisance domain",
     lambda m: m.pop("nuisance_domain")),
]


def main():
    t0 = time.time()
    arg = sys.argv[1] if len(sys.argv) > 1 else "example"
    if arg == "check":
        man = json.load(io.open(sys.argv[2], encoding="utf-8"))
        v = check(man)
        print(json.dumps(v, indent=1))
        return

    good = example_manifest()
    path = os.path.join(HERE, "manifest.example.json")
    io.open(path, "w", encoding="utf-8").write(
        json.dumps(good, indent=1) + "\n")
    v = check(good)
    print("example manifest -> %s" % v["verdict"])
    for c in v["conditions"]:
        print("   condition: %s" % c)

    print("\nREFUSAL PROBES (each must be refused on its own ground)")
    probes = []
    for tag, mut in BAD_CASES:
        m = json.loads(json.dumps(good))
        mut(m)
        r = check(m)
        ok = r["verdict"] == "SCREEN" and len(r["refusals"]) >= 1
        probes.append({"probe": tag, "verdict": r["verdict"],
                       "refusals": r["refusals"], "refused": bool(ok)})
        print("  %-38s %-11s %s" % (tag, r["verdict"],
                                    r["refusals"][0][:70] if r["refusals"]
                                    else ""))

    # the tier-C case is a CONDITIONAL, not a refusal, and that is the point
    m = json.loads(json.dumps(good))
    m["dependence"]["cross_sectional"] = {"tier": "C_same_period",
                                          "model": "observed_factor"}
    tc = check(m)

    out = {
        "schema": os.path.basename(SCHEMA),
        "example_verdict": v["verdict"],
        "example_conditions": v["conditions"],
        "probes": probes,
        "n_probes": len(probes),
        "n_refused": int(sum(1 for p in probes if p["refused"])),
        "all_probes_refused": bool(all(p["refused"] for p in probes)),
        "tier_c_verdict": tc["verdict"],
        "tier_c_conditions": tc["conditions"],
        "verdict_labels": ["GUARANTEED", "CONDITIONAL", "SCREEN"],
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w4_manifest.json", out, "diagnostic")
    print("\n%d of %d refusal probes refused; tier-C case -> %s"
          % (out["n_refused"], out["n_probes"], out["tier_c_verdict"]))


if __name__ == "__main__":
    main()
