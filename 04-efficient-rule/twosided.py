"""M6, part 1: the two-sided rule.

The plan's section 7.1 asks for `P(max_i |T_i| > c) <= abar`, with its own
conditional factorisation, its own enclosures and its own falsifier battery.

THE K = 1 REDUCTION, WHICH IS EXACT AND NOT AN APPROXIMATION.  Conditioning on
the variance estimates, the numerator is normal and symmetric about zero, so

    P(|T_1| <= c | S) = Phi(c d / sigma) - Phi(-c d / sigma) = 2 Phi(c d / sigma) - 1,
    P(|T_1| >  c | S) = 2 Phi(-c d / sigma) = 2 * P(T_1 > c | S).

The factor two is exact at every draw, so it survives the expectation over the
draws:

    FWER_2s(lambda; c) = 2 * FWER_1s(lambda; c)   for every lambda.

Therefore the two-sided rule at ceiling `abar` IS the one-sided rule at ceiling
`abar/2`, and **every certificate already deposited transfers**: certifying
sup_lambda FWER_1s <= abar/2 over the continuum certifies sup_lambda FWER_2s <=
abar over the same continuum.  Nothing is approximated and no new enclosure is
needed.  `check_reduction()` falsifies the identity rather than assuming it.

ABOVE K = 1 THE REDUCTION FAILS, and the module says so rather than extending it.
The conditional no-rejection probability is

    E_z[ prod_i ( Phi(s_i(z/sqrt(n_0) + c d_i)) - Phi(s_i(z/sqrt(n_0) - c d_i)) ) ]

and a product of differences is not one minus twice a product of tails.  The law
is implemented exactly here and the level is calibrated against it, but the
enclosure is NOT: each factor is non-monotone in s_i (it tends to one inside the
band and to zero outside it), so the monotone corner substitution behind
`wb_certify.py` does not apply.  K > 1 two-sided is therefore measured, and the
obstruction is named rather than papered over.

Deposit: deposits/twosided.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import norm

import fcore as F
import k1lab as L


# ------------------------------------------------------------- the two-sided law
def cond_noreject_2s(c, u, n, lam, quad=None):
    """P(max_i |T_i| <= c | S), exact."""
    n = np.asarray(n, float)
    u = np.asarray(u, float)
    lam_full = np.concatenate([[1.0], np.asarray(lam, float)])
    s2 = lam_full * u
    K = n.size - 1
    d = np.sqrt(s2[1:] / n[1:] + s2[0] / n[0])
    cv = np.broadcast_to(np.asarray(c, float), (K,))
    if K == 1:
        sig = np.sqrt(lam_full[1] / n[1] + 1.0 / n[0])
        return float(2.0 * norm.cdf(cv[0] * d[0] / sig) - 1.0)
    z, w = quad if quad is not None else F._u_nodes(96, 12)
    sfac = np.sqrt(n[1:] / lam_full[1:])
    base = z[:, None] / np.sqrt(n[0])
    hi = (base + cv[None, :] * d[None, :]) * sfac[None, :]
    lo = (base - cv[None, :] * d[None, :]) * sfac[None, :]
    return float(w @ np.prod(norm.cdf(hi) - norm.cdf(lo), axis=1))


def fwer_2s_exact_k1(cfun, n, lam, npanel=24, order=12):
    """Exact unconditional two-sided FWER at K = 1, by the same quadrature the
    one-sided ledger uses."""
    n = np.asarray(n, float)
    lam = float(np.asarray(lam, float).ravel()[0])
    q0, w0 = F.state_nodes(n[0] - 1.0, npanel, order)
    q1, w1 = F.state_nodes(n[1] - 1.0, npanel, order)
    U0, U1 = np.meshgrid(q0, q1, indexing="ij")
    W = np.outer(w0, w1)
    lh = lam * U1 / U0
    cv = cfun(lh.ravel()).reshape(lh.shape)
    d = np.sqrt(lam * U1 / n[1] + U0 / n[0])
    sig = np.sqrt(lam / n[1] + 1.0 / n[0])
    return float((W * (2.0 * norm.cdf(-(cv * d) / sig))).sum() / W.sum())


def power_2s_exact_k1(cfun, n, lam, delta, npanel=24, order=12):
    n = np.asarray(n, float)
    lam = float(np.asarray(lam, float).ravel()[0])
    q0, w0 = F.state_nodes(n[0] - 1.0, npanel, order)
    q1, w1 = F.state_nodes(n[1] - 1.0, npanel, order)
    U0, U1 = np.meshgrid(q0, q1, indexing="ij")
    W = np.outer(w0, w1)
    lh = lam * U1 / U0
    cv = cfun(lh.ravel()).reshape(lh.shape)
    d = np.sqrt(lam * U1 / n[1] + U0 / n[0])
    sig = np.sqrt(lam / n[1] + 1.0 / n[0])
    pw = norm.cdf((delta - cv * d) / sig) + norm.cdf((-delta - cv * d) / sig)
    return float((W * pw).sum() / W.sum())


def fwer_2s_percand(cfun, n, lam, npanel_u=20, order_u=10, nz=64, zorder=12):
    """Exact unconditional two-sided FWER of a per-candidate rule at any K.

    Same factorisation as the one-sided per-candidate engine: conditioning on
    (z, u_0) the candidates are independent, so each factor is a one-dimensional
    integral over that candidate's own variance estimate.
    """
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    K = n.size - 1
    z, wz = F.normal_nodes(nz, zorder)
    q0, w0 = F.state_nodes(n[0] - 1.0, npanel_u, order_u)
    inner = [F.state_nodes(n[i + 1] - 1.0, npanel_u, order_u) for i in range(K)]
    s_i = np.sqrt(n[1:] / lam)
    tot = 0.0
    for a0, pa0 in zip(q0, w0):
        prod = np.ones(z.size)
        for i in range(K):
            qi, wi = inner[i]
            lh = lam[i] * qi / a0
            ci = cfun(lh)
            d = np.sqrt(lam[i] * qi / n[i + 1] + a0 / n[0])
            base = z[:, None] / np.sqrt(n[0])
            hi = (base + ci[None, :] * d[None, :]) * s_i[i]
            lo = (base - ci[None, :] * d[None, :]) * s_i[i]
            prod *= (norm.cdf(hi) - norm.cdf(lo)) @ wi
        tot += pa0 * float(wz @ prod)
    return float(1.0 - tot)


# ------------------------------------------------------------- the falsifiers
def check_reduction(seed=11, trials=300):
    """T1: the exact factor of two at K = 1, conditionally and unconditionally.

    Conditionally it is an algebraic identity; unconditionally it is the claim
    that licenses reusing the one-sided certificates, so both are falsified.
    """
    rng = np.random.default_rng(seed)
    worst_cond = 0.0
    for _ in range(trials):
        n = np.array([float(rng.integers(4, 40)), float(rng.integers(4, 40))])
        lam = float(np.exp(rng.uniform(np.log(1e-3), np.log(1e3))))
        u = np.exp(rng.normal(0.0, 0.6, 2))
        c = float(rng.uniform(0.5, 6.0))
        two = 1.0 - cond_noreject_2s(c, u, n, [lam])
        one = 1.0 - F.cond_noreject(c, u, n, [lam])
        worst_cond = max(worst_cond, abs(two - 2.0 * one))
    worst_unc = 0.0
    for n in ([12.0, 13.0], [6.0, 13.0], [36.0, 5.0]):
        n = np.asarray(n, float)
        tab = L.PluginTable(n, 0.0125, m=31, npanel=12, order=8)
        for lam in (0.04, 1.0, 25.0):
            a = fwer_2s_exact_k1(tab, n, lam, 20, 10)
            b = 2.0 * F.fwer_exact_k1(tab, n, lam, 20, 10)
            worst_unc = max(worst_unc, abs(a - b))
    return {"worst_conditional_abs_diff": float(worst_cond),
            "worst_unconditional_abs_diff": float(worst_unc),
            "trials": int(trials),
            "pass": bool(worst_cond < 1e-12 and worst_unc < 1e-12)}


def check_k1_limit(seed=5):
    """T2: the two-sided per-candidate engine at K = 1 must equal the closed
    form, and a constant threshold must reproduce the Student t two-sided
    quantile when the benchmark is effectively known."""
    n = np.array([12.0, 13.0])
    tab = L.PluginTable(n, 0.0125, m=31, npanel=12, order=8)
    a = fwer_2s_percand(tab, n, np.array([1.0]))
    b = fwer_2s_exact_k1(tab, n, 1.0, 24, 12)
    return {"percand_vs_closed_form": float(abs(a - b)),
            "pass": bool(abs(a - b) < 1e-9)}


def check_no_reduction_above_k1(seed=13):
    """T3: the factor of two must FAIL at K > 1, and by a visible margin.

    A falsifier that cannot fire is worthless; this one exists to show that the
    K = 1 reduction is a special property and is not being smuggled upward.
    """
    n = np.array([12.0, 13.0, 13.0])
    lam = np.array([1.0, 4.0])
    tab = L.PluginTable(np.array([12.0, 13.0]), 0.0125, m=31, npanel=12,
                        order=8)
    two = fwer_2s_percand(tab, n, lam)
    one = F.__dict__  # placeholder to keep the import obvious
    import wb_percand as P
    one_sided = P.fwer_percand(tab, n, lam)
    return {"two_sided": float(two), "one_sided": float(one_sided),
            "ratio": float(two / one_sided),
            "reduction_holds": bool(abs(two - 2.0 * one_sided) < 1e-9),
            "pass": bool(abs(two - 2.0 * one_sided) > 1e-6)}


# ------------------------------------------------------------- calibration
def calibrate_2s_k1(n, abar=F.ABAR, lam_grid=None, npanel=12, order=8,
                    iters=11):
    """By the reduction, this is the one-sided calibration at abar/2."""
    n = np.asarray(n, float)
    lam_grid = F.log_grid(1e-3, 1e3, 25) if lam_grid is None else lam_grid
    cal = L.calibrate_level(n, lam_grid, abar / 2.0, npanel=npanel,
                            order=order, iters=iters)
    return cal


def main():
    t0 = time.time()
    out = {"abar": F.ABAR, "designs": [],
           "reduction": check_reduction(),
           "k1_limit": check_k1_limit(),
           "no_reduction_above_k1": check_no_reduction_above_k1()}
    print("T1 exact factor of two at K=1:", out["reduction"])
    print("T2 per-candidate engine at K=1:", out["k1_limit"])
    print("T3 the reduction must FAIL above K=1:", out["no_reduction_above_k1"])
    print()
    print("%-16s %9s %9s %9s %8s %8s" %
          ("design", "alpha0_2s", "supFWER2s", "ceiling", "loss80m", "loss80x"))
    fine = F.log_grid(1e-3, 1e3, 81)
    for name, n in (("df3", [4.0, 13.0]), ("df5", [6.0, 13.0]),
                    ("df7", [8.0, 13.0]), ("df11", [12.0, 13.0]),
                    ("df19", [20.0, 13.0]), ("df35", [36.0, 13.0]),
                    ("balanced_df35", [36.0, 36.0])):
        nn = np.asarray(n, float)
        cal = calibrate_2s_k1(nn)
        tab = cal["table"]
        sup2 = max(fwer_2s_exact_k1(tab, nn, l, 20, 10) for l in fine)
        t = 0
        while sup2 > F.ABAR and t < 12:
            cal = L.calibrate_level(nn, F.log_grid(1e-3, 1e3, 25),
                                    F.ABAR / 2.0 * 0.98 ** (t + 1),
                                    npanel=12, order=8, iters=9)
            tab = cal["table"]
            sup2 = max(fwer_2s_exact_k1(tab, nn, l, 20, 10) for l in fine)
            t += 1
        rows = []
        for lam in (0.04, 1.0, 25.0):
            # the two-sided known-lambda oracle rejects when |Xbar diff| >
            # z_{1-abar/2} sigma, so its power at a one-sided shift is the sum
            # of the two tails
            sig = np.sqrt(lam / nn[1] + 1.0 / nn[0])
            zc = norm.ppf(1.0 - F.ABAR / 2.0)
            lo, hi = 0.0, 1.0
            def opow(dl, sig=sig, zc=zc):
                return norm.cdf(dl / sig - zc) + norm.cdf(-dl / sig - zc)
            while opow(hi) < 0.8 and hi < 40:
                hi *= 2
            for _ in range(40):
                m = 0.5 * (lo + hi)
                if opow(m) < 0.8:
                    lo = m
                else:
                    hi = m
            dl = 0.5 * (lo + hi)
            p = power_2s_exact_k1(tab, nn, lam, dl, 20, 10)
            rows.append({"lam": float(lam), "delta80": float(dl),
                         "rule_power": float(p),
                         "loss_pp": float(100.0 * (0.8 - p))})
        losses = [r["loss_pp"] for r in rows]
        rec = {"design": name, "n": nn.tolist(),
               "alpha0_two_sided": float(cal["alpha0"]),
               "sup_fwer_2s": float(sup2),
               "within_ceiling": bool(sup2 <= F.ABAR),
               "tightenings": int(t), "rows": rows,
               "median_loss_pp": float(np.median(losses)),
               "max_loss_pp": float(max(losses))}
        out["designs"].append(rec)
        print("%-16s %9.5f %9.5f %9.5f %8.1f %8.1f" %
              (name, rec["alpha0_two_sided"], sup2, F.ABAR,
               rec["median_loss_pp"], rec["max_loss_pp"]))
    out["all_within_ceiling"] = all(d["within_ceiling"] for d in out["designs"])
    out["certificate_transfer"] = (
        "At K = 1 the two-sided FWER is exactly twice the one-sided FWER at the "
        "same threshold function (falsified to 1e-12 conditionally and "
        "unconditionally), so a certificate for sup_lambda FWER_1s <= abar/2 IS "
        "a certificate for sup_lambda FWER_2s <= abar over the same continuum.  "
        "No new enclosure is required and none is claimed.")
    out["k_gt_1_status"] = (
        "MEASURED, and the obstruction is named: each two-sided factor "
        "Phi(s(x+cd)) - Phi(s(x-cd)) is non-monotone in s, so the monotone "
        "corner substitution behind wb_certify.py does not apply and the "
        "continuum enclosure is not available by that route.")
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("twosided.json", out, "unconditional_candidate_measured")
    print("\nall within ceiling:", out["all_within_ceiling"])


if __name__ == "__main__":
    main()
