"""What the Proposition actually rests on, checked directly.

The linear-program falsifier was the wrong instrument and it is recorded as
unbuilt.  This is the right one, and it is cheaper.

A theorem is only as good as its hypotheses.  The Proposition has three that are
specific to this problem, as opposed to textbook, and all three are numerically
checkable:

  H1  the maximal invariant is (T, lambda_hat)      Lemma 1
  H2  the null and the alternative are G-invariant  so the averaging argument
                                                    stays inside the class
  H3  the LOSS is invariant, which means the oracle power beta* depends on the
      parameter only through the maximal invariant

H3 is the load-bearing one.  It is what lets the invariant mean pass through
the regret in the proof, and it is the step that fails for weighted average
power.  If it were false the Proposition would be false, and nothing else in
the argument would show it.

PART 1 proves Lemma 1 SYMBOLICALLY rather than at 1e-11.  The recovery formulas
are algebra, so sympy can close them exactly and the numerical check in
`w1_maximal_k.py` becomes a corroboration rather than the evidence.

PART 2 checks H2 and H3 by applying random group elements to random parameter
points and comparing.

PART 3 tests the CONCLUSION without a linear program.  Rather than optimise over
a mesh, it constructs explicit non-invariant rules, reading `s_0^2` as well as
`(T, lambda_hat)`, calibrates each one against the uniform constraint over
sigma_0^2 by exact quadrature, and measures what it achieves.  The Proposition
says none of them can beat the invariant floor; a mesh was never needed to ask
that.

Deposit: deposits/w1_hypotheses.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import norm

import fcore as F
import k1lab as L

ABAR = F.ABAR
ORACLE = 0.80
LAM_DOMAIN = (1e-3, 1e3)
ALT_LAMS = (0.04, 0.25, 1.0, 4.0, 25.0)


# ------------------------------------------------------------ part 1, symbolic
def lemma1_symbolic():
    """Close Lemma 1's four identities in exact algebra, at K = 1 and K = 2."""
    import sympy as sp

    out = []
    for K in (1, 2):
        n = sp.symbols("n0:%d" % (K + 1), positive=True)
        v = sp.symbols("v0:%d" % (K + 1), positive=True)   # s_i^2 of point 1
        m = sp.symbols("m0:%d" % (K + 1), real=True)       # Xbar_i of point 1
        w0, p0 = sp.symbols("w0 p0", positive=True), sp.symbols("p0_", real=True)
        w0 = sp.Symbol("w0", positive=True)                # s_0^2 of point 2
        p0 = sp.Symbol("p0", real=True)                    # Xbar_0 of point 2

        # point 2 is built to share the invariant, so its coordinates are
        # DETERMINED by point 1's invariant and its own (w0, p0)
        lam = [v[i] / v[0] for i in range(1, K + 1)]
        w = [w0] + [lam[i - 1] * w0 for i in range(1, K + 1)]
        d1 = [sp.sqrt(v[i] / n[i] + v[0] / n[0]) for i in range(1, K + 1)]
        d2 = [sp.sqrt(w[i] / n[i] + w[0] / n[0]) for i in range(1, K + 1)]
        T = [(m[i] - m[0]) / d1[i - 1] for i in range(1, K + 1)]
        p = [p0] + [p0 + T[i - 1] * d2[i - 1] for i in range(1, K + 1)]

        # the element, fitted on the benchmark pair alone
        b = sp.sqrt(w0 / v[0])
        a = p0 - b * m[0]

        ids = []
        for i in range(K + 1):
            ids.append(sp.simplify(sp.expand(a + b * m[i] - p[i])))
            ids.append(sp.simplify(sp.expand(b ** 2 * v[i] - w[i])))
        closed = all(x == 0 for x in ids)
        out.append({"K": K, "identities": len(ids), "all_zero": bool(closed),
                    "residuals": [str(x) for x in ids]})
        print("  K=%d  %d identities, all identically zero: %s"
              % (K, len(ids), closed))
    return out


# --------------------------------------------------- part 2, the hypotheses
def invariants_of_parameter(mu, sig2, n):
    """(eta_i, lambda_i), the maximal invariant on the parameter."""
    mu = np.asarray(mu, float)
    sig2 = np.asarray(sig2, float)
    n = np.asarray(n, float)
    eta = (mu[1:] - mu[0]) / np.sqrt(sig2[1:] / n[1:] + sig2[0] / n[0])
    return eta, sig2[1:] / sig2[0]


def beta_star_from_theta(mu, sig2, n, alpha=ABAR):
    """beta*, the oracle's power, computed from the RAW parameter.

    Nothing here is written in terms of the invariant, so the check that this
    is unmoved by the group is a check and not a restatement.  The oracle is
    the test that knows the true ratio, so its constant is `c_pop_k1` at that
    ratio, and its power is the exact quadrature the paper uses.
    """
    mu = np.asarray(mu, float)
    sig2 = np.asarray(sig2, float)
    n = np.asarray(n, float)
    lam = sig2[1] / sig2[0]
    delta = (mu[1] - mu[0]) / np.sqrt(sig2[0])   # the quadrature is in sigma_0
    c = L.c_pop_k1(n, float(lam), alpha)
    return F.power_exact_k1(lambda x: np.full_like(np.asarray(x, float), c),
                            n, float(lam), float(delta), 20, 10)


def hypothesis_checks(n, D=4000, seed=404):
    """Apply random group elements to random parameters; the invariant, the
    null, the alternative and beta* must all be unmoved."""
    rng = np.random.default_rng(int(seed))
    n = np.asarray(n, float)
    mu = rng.normal(0.0, 4.0, size=(D, 2))
    s2 = np.exp(rng.normal(0.0, 1.0, size=(D, 2)))
    a = rng.normal(0.0, 8.0, D)
    b = np.exp(rng.normal(0.0, 1.2, D))

    e0, l0 = invariants_of_parameter(mu.T, s2.T, n)
    mu2 = a[:, None] + b[:, None] * mu
    s22 = (b ** 2)[:, None] * s2
    e1, l1 = invariants_of_parameter(mu2.T, s22.T, n)
    drift_eta = float(np.max(np.abs(e1 - e0)))
    drift_lam = float(np.max(np.abs(l1 / l0 - 1.0)))

    # the null: equal means stay equal, and the ratio is untouched
    mun = np.repeat(rng.normal(0.0, 4.0, D)[:, None], 2, axis=1)
    mun2 = a[:, None] + b[:, None] * mun
    null_preserved = bool(np.max(np.abs(mun2[:, 1] - mun2[:, 0])) < 1e-9)

    # the alternative: sign of the shift is preserved because b > 0
    alt_preserved = bool(np.all(np.sign(e1) == np.sign(e0)))

    # H3: beta*, computed from the raw parameter, must be unmoved by the
    # group.  This is the load-bearing hypothesis: it is what lets the
    # invariant mean pass through the regret, and it is the step that fails for
    # weighted average power.
    #
    # The shift is drawn to put beta* in the INTERIOR.  Drawing the mean freely
    # saturates it at 0 or 1 on almost every draw, and a quantity pinned at an
    # endpoint is unmoved by anything, so that version of the check passes
    # whether or not the hypothesis holds.
    rows = []
    for i in range(48):
        sg = s2[i]
        sig = np.sqrt((sg[1] / sg[0]) / n[1] + 1.0 / n[0])
        target = 0.05 + 0.90 * ((i + 0.5) / 48.0)
        eta = norm.ppf(target) + norm.ppf(1.0 - ABAR)
        m0 = np.array([mu[i, 0], mu[i, 0] + eta * sig * np.sqrt(sg[0])])
        base = beta_star_from_theta(m0, sg, n)
        moved = beta_star_from_theta(a[i] + b[i] * m0, (b[i] ** 2) * sg, n)
        rows.append({"a": float(a[i]), "b": float(b[i]),
                     "beta_star": float(base),
                     "beta_star_after_g": float(moved),
                     "abs_diff": float(abs(base - moved))})
    interior = [r for r in rows if 0.02 < r["beta_star"] < 0.98]
    return {"drift_eta": drift_eta, "drift_lambda": drift_lam,
            "null_preserved": null_preserved,
            "alternative_preserved": alt_preserved,
            "beta_star": rows,
            "beta_star_range": [float(min(r["beta_star"] for r in rows)),
                                float(max(r["beta_star"] for r in rows))],
            "n_beta_star_draws": len(rows),
            "n_beta_star_interior": len(interior),
            "worst_beta_star_drift": float(max(r["abs_diff"] for r in rows)),
            "worst_beta_star_drift_interior": float(
                max(r["abs_diff"] for r in interior)) if interior else None}


# ----------------------------------------- part 3, explicit non-invariant rules
def rule_fwer(cfun, n, lam, v, npanel=20, order=10):
    """Exact FWER of a rule c(lambda_hat, s_0^2) at (lambda, sigma_0^2 = v)."""
    n = np.asarray(n, float)
    q0, w0 = F.state_nodes(n[0] - 1.0, npanel, order)
    q1, w1 = F.state_nodes(n[1] - 1.0, npanel, order)
    U0, U1 = np.meshgrid(q0, q1, indexing="ij")
    W = np.outer(w0, w1)
    lh = lam * U1 / U0
    d = np.sqrt(lam * U1 / n[1] + U0 / n[0])
    sig = np.sqrt(lam / n[1] + 1.0 / n[0])
    c = cfun(lh, v * U0)
    return float((W * norm.cdf(-c * d / sig)).sum() / W.sum())


def rule_power(cfun, n, lam, v, delta, npanel=20, order=10):
    n = np.asarray(n, float)
    q0, w0 = F.state_nodes(n[0] - 1.0, npanel, order)
    q1, w1 = F.state_nodes(n[1] - 1.0, npanel, order)
    U0, U1 = np.meshgrid(q0, q1, indexing="ij")
    W = np.outer(w0, w1)
    lh = lam * U1 / U0
    d = np.sqrt(lam * U1 / n[1] + U0 / n[0])
    sig = np.sqrt(lam / n[1] + 1.0 / n[0])
    c = cfun(lh, v * U0)
    return float((W * norm.cdf((delta - c * d) / sig)).sum() / W.sum())


def tilt(tab, gamma, kind):
    """A rule that reads s_0^2 as well as the estimated ratio.

    The `deproj` family is the DIRECTED probe, and it is the rule a sceptic
    would actually propose.  Since s_0^2 = sigma_0^2 u_0 and the estimated
    ratio is lambda u_1 / u_0, a rule that guesses sigma_0^2 = v0 recovers an
    estimate of u_0 and multiplies the ratio back by it, which de-biases the
    ratio when the guess is good.  At gamma = 0 it IS the deployed rule, so
    that arm is a sanity check that has to return the baseline; at gamma = 1
    it is the full correction.  The generic shapes elsewhere in this file ask
    whether ANY tilt helps; this one asks whether the tilt that should help
    does.
    """
    if kind.startswith("deproj"):
        v0 = float(kind.split("@")[1])

        def f(lh, s0):
            return tab(lh * (np.maximum(s0, 1e-300) / v0) ** gamma)
        return f

    def f(lh, s0):
        base = tab(lh)
        x = np.log(np.maximum(s0, 1e-300))
        if kind == "bump":
            pert = np.exp(-0.5 * x ** 2)
        elif kind == "step":
            pert = np.where(x > 0.0, 1.0, -1.0)
        elif kind.startswith("band"):
            # conservative baseline, LIBERAL in a band: the only shape that can
            # buy anything, since being liberal somewhere has to be paid for
            x0 = float(kind.split("@")[1])
            pert = 1.0 - 2.0 * np.exp(-0.5 * ((x - x0) / 0.8) ** 2)
        else:                                    # linear, clipped
            pert = np.clip(x / 3.0, -1.0, 1.0)
        return base * (1.0 + gamma * pert)
    return f


def main():
    t0 = time.time()
    n = np.array([36.0, 13.0])

    print("PART 1  Lemma 1 in exact algebra")
    sym = lemma1_symbolic()

    print("\nPART 2  the hypotheses the Proposition rests on")
    hyp = hypothesis_checks(n)
    print("  invariant unmoved by the group: eta %.2e, lambda %.2e"
          % (hyp["drift_eta"], hyp["drift_lambda"]))
    print("  null preserved: %s ; alternative preserved: %s"
          % (hyp["null_preserved"], hyp["alternative_preserved"]))
    print("  beta* spans %.4f-%.4f, %d of %d draws interior, and the group "
          "moves it by at most %.2e"
          % (hyp["beta_star_range"][0], hyp["beta_star_range"][1],
             hyp["n_beta_star_interior"], hyp["n_beta_star_draws"],
             hyp["worst_beta_star_drift"]))

    print("\nPART 3  explicit non-invariant rules, no mesh and no program")
    cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                            npanel=12, order=8, iters=10)
    tab = cal["table"]
    lam_grid = F.log_grid(*LAM_DOMAIN, 41)
    v_grid = np.exp(np.linspace(-6.0, 6.0, 25))
    alts = [(l, F.delta_for_oracle_power_k1(n, [l], ORACLE))
            for l in ALT_LAMS]

    def worst_fwer(cfun):
        return max(rule_fwer(cfun, n, l, v)
                   for l in lam_grid for v in v_grid)

    def worst_regret(cfun):
        return max(ORACLE - rule_power(cfun, n, l, v, d)
                   for l, d in alts for v in v_grid)

    flat = lambda lh, s0: tab(lh)                              # noqa: E731
    base_f, base_r = worst_fwer(flat), worst_regret(flat)
    print("  the invariant rule: worst FWER %.5f, worst regret %.2f pp"
          % (base_f, 100 * base_r))

    rows = []
    for kind in ("deproj@1.0", "deproj@0.3", "deproj@3.0"):
        for g in (0.0, 0.25, 0.5, 1.0):
            r = tilt(tab, g, kind)
            f = worst_fwer(r)
            # calibrate: scale the whole rule until it meets the ceiling
            lo, hi = -2.0, 40.0
            for _ in range(26):
                mid = 0.5 * (lo + hi)
                if worst_fwer(lambda lh, s0, t=mid, rr=r: rr(lh, s0) + t) \
                        > ABAR:
                    lo = mid
                else:
                    hi = mid
            k = hi
            cal_rule = (lambda lh, s0, t=k, rr=r: rr(lh, s0) + t)
            ach = worst_fwer(cal_rule)
            reg = worst_regret(cal_rule)
            gain = 100.0 * (base_r - reg)
            rows.append({"kind": kind, "gamma": float(g),
                     "liberal_somewhere": True,
                     "calibrated": bool(ach <= ABAR + 1e-9),
                     "achieved_worst_fwer": float(ach),
                         "uncalibrated_worst_fwer": float(f),
                         "calibration_scale": float(k),
                         "worst_regret_pp": float(100 * reg),
                         "gain_over_invariant_pp": float(gain),
                         "beats_the_floor": bool(gain > 0.01)})
            print("  %-10s gamma %.2f  shift %+.4f  worst regret %.2f pp  "
                  "gain %+6.3f pp  %s"
                  % (kind, g, k, 100 * reg, gain,
                     "BEATS THE FLOOR" if gain > 0.01 else ""))

    sanity = [r for r in rows if r["gamma"] == 0.0]
    sanity_ok = bool(sanity and all(
        abs(r["worst_regret_pp"] - 100 * base_r) < 0.01 for r in sanity))
    print("  sanity arm (gamma = 0) reproduces the deployed rule: %s"
          % sanity_ok)
    infeasible = [r for r in rows if not r["calibrated"]]
    print("  %d of %d rules could not be calibrated and are not losses"
          % (len(infeasible), len(rows)))
    beat = [r for r in rows if r["beats_the_floor"] and r["calibrated"]]
    out = {
        "abar": ABAR, "oracle_power": ORACLE, "n": n.tolist(),
        "lemma1_symbolic": sym,
        "lemma1_closed_symbolically": bool(all(s["all_zero"] for s in sym)),
        "hypotheses": hyp,
        "H2_null_and_alternative_invariant": bool(
            hyp["null_preserved"] and hyp["alternative_preserved"]),
        "H3_loss_invariant": bool(hyp["worst_beta_star_drift"] < 1e-10
                              and hyp["n_beta_star_interior"] >= 40),
        "invariant_rule_worst_fwer": float(base_f),
        "invariant_rule_worst_regret_pp": float(100 * base_r),
        "rules": rows,
        "sanity_arm_reproduces_baseline": sanity_ok,
        "n_infeasible": len(infeasible),
        "calibration": "additive shift on the threshold",
        "n_rules": len(rows),
        "n_rules_liberal_somewhere": sum(1 for r in rows
                                         if r["liberal_somewhere"]),
        "n_that_beat_the_floor": len(beat),
        "conclusion_survives": bool(not beat),
        "worst_gain_pp": float(max(r["gain_over_invariant_pp"]
                                   for r in rows)),
        "note": ("Part 3 replaces the linear program that failed.  Each rule "
                 "is written down rather than optimised, calibrated against "
                 "the uniform constraint over sigma_0^2 by the same exact "
                 "quadrature the paper uses, and compared with the invariant "
                 "rule on the same instrument.  No mesh sits between the two, "
                 "so no mesh mismatch can be mistaken for an effect."),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w1_hypotheses_directed.json", out, "exact_in_model")
    print("\nLemma 1 closes symbolically: %s" % out["lemma1_closed_symbolically"])
    print("H2 holds: %s ; H3 holds: %s"
          % (out["H2_null_and_alternative_invariant"], out["H3_loss_invariant"]))
    print("%d of %d explicit non-invariant rules beat the floor; worst gain "
          "%+.3f points" % (len(beat), len(rows), out["worst_gain_pp"]))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
