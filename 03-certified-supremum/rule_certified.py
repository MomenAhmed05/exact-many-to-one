"""P4 step 20: the shipped theorem-grade rule (plan section 4.3), with a decision certificate.

The eight steps of the plan:

  1. inputs: `n_0..n_K`, sample means, sample variances, the hard ceiling `abar`, the miss budget
     `gamma`;
  2. `alpha0 = abar - gamma` (the HARD-CEILING parameterisation: the inner level is reduced so the
     guarantee reads `FWER <= alpha0 + gamma = abar`, never "alpha + gamma described as a hard abar");
  3. `lambda_hat_i = S_i^2 / S_0^2`;
  4. the ratio region `Lambda(lambda_hat)` with certified coverage at least `1 - gamma`;
  5. a certified `c_hat >= sup_{lambda in Lambda} c*(lambda; alpha0)` from the branch-and-bound,
     with its achieved gap reported;
  6. the Welch statistics `T_i`;
  7. promote candidate i iff `T_i > c_hat`;
  8. a decision certificate recording inputs, region, coverage lower bound, the `c` bracket, the
     branch-and-bound certificate identity, the error budget, the software version and the tier -
     plus the required theorem statement and the model exclusions.

The runtime per comparison family is measured against the predeclared `<= 60 s` target (the
certification tier; the fast vertex rule remains the research-triage path).

Run:  python rule_certified.py        (a few minutes: one certification per family)
Out:  deposits/rule_certified.json    (tier: region_sup_certified)
"""

from __future__ import annotations

import hashlib
import importlib
import os
import time

import cert_bnb
import numpy as np
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "rule_certified.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
draw_var = CORE.draw_var

ABAR = paths.ABAR
GAMMA = paths.GAMMA
ALPHA0 = paths.ALPHA0

FAMILIES = (
    ("desk_B_df3", [4.0, 13.0, 13.0], [1.0, 25.0]),
    ("desk_A_df3", [4.0, 13.0, 13.0], [0.04, 0.04]),
    ("desk_B_df35", [36.0, 13.0, 13.0], [1.0, 25.0]),
    ("K5_df11", [12.0] + [12.0] * 5, [0.2, 0.5, 1.0, 2.0, 5.0]),
)


def decide(n, means, variances, abar=ABAR, gamma=GAMMA, seed=20260918, max_boxes=64):
    """The full eight-step rule.  Returns the decision certificate."""
    t0 = time.time()
    n = np.asarray(n, float)
    means = np.asarray(means, float)
    variances = np.asarray(variances, float)
    alpha0 = _f(abar) - _f(gamma)
    cert: dict = {
        "tier": "region_sup_certified",
        "abar": _f(abar),
        "gamma": _f(gamma),
        "alpha0": alpha0,
        "n": list(n),
        "means": list(means),
        "variances": list(variances),
        "step": {},
        "decisions": None,
    }
    # step 3: ratio estimates
    lam_hat = variances[1:] / variances[0]
    cert["lambda_hat"] = [round(_f(v), 10) for v in lam_hat]
    # step 4: the region (coverage >= 1 - gamma)
    lo, hi, diag = ENV.lam_region(n, lam_hat, _f(gamma), scale=None)
    lo, hi = np.asarray(lo[0], float), np.asarray(hi[0], float)
    cert["region"] = {
        "lo": [round(_f(v), 10) for v in lo],
        "hi": [round(_f(v), 10) for v in hi],
        "scale": _f(diag.get("scale")),
        "coverage_exact": _f(diag.get("coverage_exact")),
        "coverage_ge_1_minus_gamma": bool(
            _f(diag.get("coverage_exact")) >= 1.0 - _f(gamma)
        ),
    }
    # step 5: the certified threshold (conditional on this draw's statistics)
    s0_var = variances[0]
    si_var = variances[1:]
    d_row = np.sqrt(si_var / n[1:] + s0_var / n[0])
    rng = np.random.default_rng(seed)
    bnb = cert_bnb.certify_region(
        lo, hi, d_row, n, alpha=alpha0, rng=rng, max_boxes=_i(max_boxes, 64)
    )
    c_hat = _f(bnb["c_U"])
    cert["threshold"] = {
        "c_hat": c_hat,
        "c_L": bnb["c_L"],
        "gap": bnb["gap"],
        "certified_at_tol": bnb["certified_at_tol"],
        "tol": _f(paths.TOL),
        "n_splits": bnb["n_splits"],
        "n_pruned": bnb["n_pruned"],
        "termination": bnb["termination"],
        "bnb_settings": bnb.get("witnesses") is not None,
    }
    cert["bnb_certificate_id"] = hashlib.sha256(
        repr([bnb["c_L"], bnb["c_U"], bnb["n_splits"], bnb["n_pruned"]]).encode()
    ).hexdigest()[:16]
    # step 6: the statistics
    se = np.sqrt(si_var / n[1:] + s0_var / n[0])
    T = (means[1:] - means[0]) / se
    cert["statistics"] = [round(_f(v), 6) for v in T]
    # step 7: the decisions
    cert["decisions"] = [bool(_f(v) > c_hat) for v in T]
    cert["n_promoted"] = _i(np.sum(cert["decisions"]))
    # step 8: the error budget and the statement
    cert["error_budget"] = {
        "Delta_p": _f(bnb.get("gap")) if bnb.get("gap") == bnb.get("gap") else NAN,
        "note": "the gap is the certified supremum enclosure; the numerical "
        "enclosure width is ~1e-13 (see deposits/error_budget.json)",
    }
    cert["software_version"] = "03-certified-supremum P4 rule_certified.py"
    cert["seconds"] = _f(time.time() - t0)
    cert["runtime_under_target"] = bool(cert["seconds"] <= paths.CERT_MIN_PER_FAMILY)
    return cert


def theorem_statement():
    return (
        "Under independent Gaussian samples across groups, with the stated sample sizes and the "
        "specified separate-variance Welch statistics, this decision uses internal level alpha0 = "
        "abar - gamma and a ratio-region miss budget gamma.  The certified region coverage is at "
        "least 1 - gamma, and the reported threshold is a certified upper bound on the oracle "
        "critical value throughout that region.  Therefore the familywise error rate is at most "
        "alpha0 + gamma = abar."
    )


def exclusions():
    return [
        "serial dependence",
        "cross-group dependence beyond the model",
        "factor exposure",
        "non-Gaussianity",
        "rolling re-estimation",
        "research-process multiplicity",
    ]


def main():
    t0 = time.time()
    out: dict = {
        "abar": _f(ABAR),
        "gamma": _f(GAMMA),
        "alpha0": _f(ALPHA0),
        "tol": _f(paths.TOL),
        "runtime_target_s": _f(paths.CERT_MIN_PER_FAMILY),
        "theorem_statement": theorem_statement(),
        "model_exclusions": exclusions(),
        "families": [],
    }
    for name, n, lam in FAMILIES:
        n_a = np.asarray(n, float)
        lam_a = np.asarray(lam, float)
        # a synthetic desk draw with a known truth (sigma_0 = 1), so the decision can be exhibited
        s0, si, _lh = draw_var(n_a, lam_a, 1, 20260918)
        rng = np.random.default_rng(4242)
        means = np.concatenate([[0.0], rng.normal(0.0, np.sqrt(lam_a / n_a[1:]))])
        variances = np.concatenate([[s0[0]], si[0]])
        cert = decide(n_a, means, variances)
        cert["family"] = name
        cert["lambda_true"] = list(lam_a)
        out["families"].append(cert)
        print(
            f"{name:<12} alpha0 {_f(ALPHA0):.4f} | c_hat {_f(cert['threshold']['c_hat']):.6f} "
            f"gap {_f(cert['threshold']['gap']):.2e} ok {cert['threshold']['certified_at_tol']} "
            f"| T {cert['statistics']} -> promoted {cert['n_promoted']} "
            f"| {_f(cert['seconds']):.0f}s (target {paths.CERT_MIN_PER_FAMILY:.0f}s)",
            flush=True,
        )
    secs = [_f(c["seconds"]) for c in out["families"]]
    out["n_families"] = _i(len(out["families"]))
    out["max_seconds"] = _f(max(secs))
    out["all_under_runtime_target"] = bool(
        all(s <= paths.CERT_MIN_PER_FAMILY for s in secs)
    )
    out["all_coverage_certified"] = bool(
        all(c["region"]["coverage_ge_1_minus_gamma"] for c in out["families"])
    )
    out["n_certified_at_tol"] = _i(
        sum(1 for c in out["families"] if c["threshold"]["certified_at_tol"])
    )
    out["seconds"] = _f(time.time() - t0)
    print(
        f"\nfamilies {out['n_families']} | certified at tol {out['n_certified_at_tol']} | "
        f"coverage ok {out['all_coverage_certified']} | runtime ok "
        f"{out['all_under_runtime_target']} (max {_f(out['max_seconds']):.0f}s)"
    )
    tiers.write_deposit(
        out,
        OUT,
        "region_sup_certified",
        notes="shipped rule at alpha0 = abar - gamma with a decision certificate, "
        "the required theorem statement and the model exclusions",
    )
    return out


if __name__ == "__main__":
    import sys

    try:
        res = main()
        sys.exit(0 if res.get("all_coverage_certified") else 1)
    except Exception as exc:  # noqa: BLE001
        print("RULE ERROR:", repr(exc))
        sys.exit(2)
