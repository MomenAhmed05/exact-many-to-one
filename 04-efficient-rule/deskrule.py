"""THE DESK PROCEDURE.

What a desk runs, end to end, with no region, no simulation, no seed, and no
branch-and-bound at decision time:

    from deskrule import promote
    out = promote(bench_returns, [cand1, cand2, ...], abar=0.025)
    out["decisions"]      # one boolean per candidate
    out["threshold"]      # the constant it compared against
    out["certificate"]    # what was proven, and over what domain

THE RULE.

    1. lambda_hat_i = S_i^2 / S_0^2 .
    2. c = c*( lambda_hat ; alpha0(n, abar) ), the paper's exact fixed-ratio
       constant evaluated at the estimated ratios, at an INTERNAL level alpha0.
    3. promote candidate i iff T_i > c.

Only step 2's internal level is new, and it is the whole of the construction:
alpha0 is the largest level whose UNCONDITIONAL familywise error rate stays at
or below the ceiling for every ratio vector in the declared domain,

    sup_{lambda in [1e-3, 1e3]^K}  E_S[ P_lambda( max_i T_i > c(S) | S ) ]  <=  abar .

WHY NOT THE REGION RULE.  The certified region rule of `03-certified-supremum` protects, on every
draw, against every ratio in an admission region.  The guarantee it has to meet
is an average over draws, so that is strictly more than required, and the price
is severe: at a 36-month benchmark its threshold is about three times the oracle
and its measured power is a couple of per cent.  Calibrating one internal level
against the exact unconditional constraint costs a few per cent of the level
instead, and keeps the power.

LEVELS OF ASSURANCE, and they are not the same:

  K = 1, joint          CERTIFIED over the CONTINUUM of ratios by the enclosure
                        in `wb_certify.py`.
  any K, per-candidate  CERTIFIED over the continuum by `certify_anyk.py`: the
                        supremum over the ratio cube is at most K times the
                        one-candidate supremum, which the K = 1 enclosure bounds
                        tightly, so certifying at abar/K suffices.  Two-sided
                        follows at abar/(2K).
  K > 1, joint          MEASURED.  The joint threshold couples all the
                        coordinates, so neither the union step nor the enclosure
                        applies; it is the more powerful rule and the one that
                        cannot be certified here.

`promote(..., require_certificate=True)` returns the per-candidate rule wherever
a certificate exists, and raises rather than quietly serving a measured rule.
"""

from __future__ import annotations

import os

import numpy as np

import fcore as F
import k1lab as L
import wb_certify as C

_CACHE = {}

MODEL_EXCLUSIONS = [
    "serial dependence, UNLESS the series is whitened with a declared "
    "autocorrelation first (see serial.py); an estimated one is a plug-in",
    "factor exposure, UNLESS the returns are residualised with declared "
    "loadings first (see factor.py); estimated loadings are a plug-in",
    "non-Gaussian returns",
    "rolling re-estimation and research-process multiplicity",
]


def _key(n, abar):
    return (tuple(float(x) for x in n), float(abar))


def _certified_table(name, n, abar):
    """Load the deposited certificate for this design, if one exists."""
    try:
        dep = F.read_deposit("wb_certified_rule.json")
    except Exception:  # noqa: BLE001
        return None
    if abs(dep.get("abar", -1) - abar) > 1e-12:
        return None
    for d in dep["designs"]:
        if (len(d["n"]) == n.size and np.allclose(d["n"], n)
                and d.get("certified")):
            return d
    return None


def threshold_rule(n, abar=F.ABAR, force_recompute=False):
    """Return (callable c(lambda_hat), certificate dict) for this design."""
    n = np.asarray(n, float)
    k = _key(n, abar)
    if k in _CACHE and not force_recompute:
        return _CACHE[k]
    K = n.size - 1

    if K == 1:
        rec = None if force_recompute else _certified_table(None, n, abar)
        if rec is not None:
            rule = C.StepRule(rec["c"], rec["edges"])
            cert = {
                "tier": "unconditional_valid_certified",
                "statement": (
                    "sup over lambda in [%g, %g] of the unconditional FWER is "
                    "at most %.6f, proven by rigorous enclosure over the "
                    "continuum." % (1e-3, 1e3, rec["sup_upper_bound"])),
                "alpha0": rec["alpha0"],
                "sup_upper_bound": rec["sup_upper_bound"],
                "boxes": rec["certificate"]["n_boxes_accepted"],
                "source": "deposits/wb_certified_rule.json",
            }
            _CACHE[k] = (rule, cert)
            return _CACHE[k]
        # compute one on the spot, measured tier only
        cal = L.calibrate_level(n, F.log_grid(1e-3, 1e3, 41), abar,
                                npanel=12, order=8, iters=11)
        cert = {"tier": "unconditional_candidate_measured",
                "statement": ("sup over a 41-point ratio grid of the "
                              "unconditional FWER is %.6f <= %.6f; NOT "
                              "certified over the continuum."
                              % (cal["sup_fwer"], abar)),
                "alpha0": cal["alpha0"], "source": "computed on demand"}
        _CACHE[k] = (cal["table"], cert)
        return _CACHE[k]

    # K > 1: the joint plug-in at its calibrated level, measured tier
    dep = None
    try:
        dep = F.read_deposit("wb_general.json")
    except Exception:  # noqa: BLE001
        pass
    a0 = None
    if dep:
        for d in dep["designs"]:
            if len(d["n"]) == n.size and np.allclose(d["n"], n):
                a0 = d["alpha0"]
                break
    if a0 is None:
        raise RuntimeError(
            "no calibrated level deposited for n=%s at abar=%g; run "
            "wb_general.py for this design first" % (n.tolist(), abar))
    cert = {"tier": "unconditional_candidate_measured",
            "statement": ("sup over the declared ratio pattern set of the "
                          "unconditional FWER is at most %.6f; the supremum "
                          "is over patterns, not over the continuum, so this "
                          "is measured and not certified." % abar),
            "alpha0": a0, "source": "deposits/wb_general.json"}
    _CACHE[k] = (a0, cert)
    return _CACHE[k]


def certified_percand(n, abar=F.ABAR, two_sided=False):
    """The any-K certified per-candidate rule for this design, or None."""
    n = np.asarray(n, float)
    try:
        dep = F.read_deposit("certify_anyk.json")
    except Exception:  # noqa: BLE001
        return None
    if abs(dep.get("abar", -1) - abar) > 1e-12:
        return None
    for d in dep["designs"]:
        if (len(d["n"]) == n.size and np.allclose(d["n"], n)
                and bool(d["two_sided"]) == bool(two_sided)
                and d.get("certified")):
            return d
    return None


def promote(bench, candidates, abar=F.ABAR, require_certificate=False,
            two_sided=False):
    """Run the rule on real return series.

    `bench` and each entry of `candidates` are 1-D arrays of returns.  Sample
    sizes may differ across groups.
    """
    b = np.asarray(bench, float).ravel()
    cs = [np.asarray(c, float).ravel() for c in candidates]
    n = np.array([b.size] + [c.size for c in cs], float)
    K = len(cs)
    s2 = np.array([b.var(ddof=1)] + [c.var(ddof=1) for c in cs])
    m = np.array([b.mean()] + [c.mean() for c in cs])
    lam_hat = s2[1:] / s2[0]
    d = np.sqrt(s2[1:] / n[1:] + s2[0] / n[0])
    t = (m[1:] - m[0]) / d

    if require_certificate or two_sided:
        rec = certified_percand(n, abar, two_sided)
        if rec is None:
            raise RuntimeError(
                "no certificate deposited for n=%s at abar=%g, %s-sided; run "
                "certify_anyk.py for this design, or call with "
                "require_certificate=False to use the measured joint rule"
                % (n.tolist(), abar, "two" if two_sided else "one"))
        rule = C.StepRule(rec["c"], rec["edges"])
        c = float(np.asarray(rule(lam_hat)).ravel().max())
        cv = np.asarray(rule(lam_hat)).ravel()
        stat = np.abs(t) if two_sided else t
        return {
            "n": n.tolist(), "K": int(K), "abar": float(abar),
            "two_sided": bool(two_sided),
            "lambda_hat": lam_hat.tolist(), "t": t.tolist(),
            "threshold": [float(x) for x in cv],
            "decisions": [bool(x) for x in (stat > cv)],
            "n_promoted": int((stat > cv).sum()),
            "certificate": {
                "tier": "unconditional_valid_certified",
                "rule": "per-candidate",
                "statement": (
                    "sup over the ratio cube [1e-3, 1e3]^K of the "
                    "unconditional %s-sided FWER is at most %.6f, by the union "
                    "step over a one-candidate supremum certified by rigorous "
                    "enclosure over the continuum."
                    % ("two" if two_sided else "one",
                       rec["union_bound_total"])),
                "alpha0": rec["alpha0"],
                "union_bound_total": rec["union_bound_total"],
                "source": "deposits/certify_anyk.json"},
            "model_exclusions": MODEL_EXCLUSIONS,
        }

    rule, cert = threshold_rule(n, abar)
    if K == 1:
        c = float(np.asarray(rule(np.array([lam_hat[0]]))).ravel()[0])
    else:
        c = F.oracle_c_pop(n, lam_hat, cert["alpha0"], F.FAST)
    return {
        "n": n.tolist(), "K": int(K), "abar": float(abar),
        "lambda_hat": lam_hat.tolist(),
        "t": t.tolist(), "threshold": float(c),
        "decisions": [bool(x) for x in (t > c)],
        "n_promoted": int((t > c).sum()),
        "certificate": cert,
        "model_exclusions": MODEL_EXCLUSIONS,
    }


def _demo():
    """Three worked decisions, one per level of assurance.

    The one-candidate case takes the CERTIFIED path; the two-candidate case
    takes the measured one, and says so.
    """
    rng = np.random.default_rng(1)
    print("--- K = 1, a design with a deposited certificate (n = 12, 13) ---")
    b = rng.normal(0.0, 1.0, 12)
    c1 = rng.normal(1.4, 2.0, 13)
    out = promote(b, [c1])
    for k in ("n", "lambda_hat", "t", "threshold", "decisions"):
        print("  %-12s %s" % (k, out[k]))
    print("  tier         %s" % out["certificate"]["tier"])
    print("  statement    %s" % out["certificate"]["statement"])

    print("")
    print("--- K = 2, measured tier (n = 12, 13, 13) ---")
    b = rng.normal(0.0, 1.0, 12)
    c1 = rng.normal(0.0, 2.0, 13)
    c2 = rng.normal(1.6, 1.5, 13)
    out = promote(b, [c1, c2])
    for k in ("n", "lambda_hat", "t", "threshold", "decisions"):
        print("  %-12s %s" % (k, out[k]))
    print("  tier         %s" % out["certificate"]["tier"])
    print("  statement    %s" % out["certificate"]["statement"])

    print("")
    print("--- K = 3, CERTIFIED at any K, two-sided (n = 36, 13, 13, 13) ---")
    b = rng.normal(0.0, 1.0, 36)
    cs = [rng.normal(0.0, 2.0, 13), rng.normal(-2.2, 1.0, 13),
          rng.normal(2.4, 1.0, 13)]
    out = promote(b, cs, require_certificate=True, two_sided=True)
    for k in ("n", "t", "decisions"):
        print("  %-12s %s" % (k, out[k]))
    print("  threshold    %s" % [round(x, 3) for x in out["threshold"]])
    print("  tier         %s" % out["certificate"]["tier"])
    print("  statement    %s" % out["certificate"]["statement"])


if __name__ == "__main__":
    _demo()
