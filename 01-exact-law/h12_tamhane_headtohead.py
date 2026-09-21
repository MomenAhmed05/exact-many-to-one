"""Head-to-head against Tamhane and Xi (2023) on their own fourteen scenarios.

WHY THIS EXISTS. Tamhane and Xi, Biometrical Journal 65, e2200300, report the
familywise error rate of four methods across seven designs with two treatments
and seven with three, each estimated from 100,000 simulated data sets. Their
Method 3 is the state of the art and is the paper's contribution. Every number in
their tables is a Monte Carlo estimate of a quantity that, with the conditional
reduction in h12_hetero_exact.py, can be computed instead.

So this recomputes the same quantity exactly, on the same designs, and adds the
exact procedure as a fifth arm. It is the comparison that decides whether the
reduction is worth anything against the best published alternative.

WHAT THEIR METHODS ARE, from the paper, so the reimplementation is checkable.
All of Methods 1, 2 and 3 use the same Welch statistic
T_i = (Xbar_i - Xbar_0)/sqrt(S_i^2/n_i + S_0^2/n_0), written by them as
T_i = Z_i / sqrt(U_i/nuhat_i).

  Method 1  Tamhane's (1977, 1979) multiplicative inequality:
            ptilde_i <= 1 - prod_j F_{nuhat_j}(t_i), so one common critical value
            c solves prod_j F_{nuhat_j}(c) = 1 - alpha. Conservative by
            construction, since the inequality is an upper bound on the p-value.

  Method 2  Hasler and Hothorn (2008), the "PI" method: a k-variate t with the
            estimated product correlation rhohat_ij = lambdahat_i lambdahat_j and
            degrees of freedom nuhat_i, giving a DIFFERENT critical value per
            comparison. Their own remark on it: "this method does not take into
            account the correlations among the U_i's".

  Method 3  Theirs. Simulate the joint null by simulating Z (correlated normals)
            and U (correlated chi-squares). The U's are built as U_i = V_i + W_i
            with a nested ladder of chi-squares sharing V-components, and
            c_i = eta_i^2 nuhat_i, d_i = (1 - eta_i^2) nuhat_i.

  Method 0  assumes homoskedasticity and uses a POOLED variance, so it is a
            different statistic. It is excluded here rather than misrepresented;
            both Hasler and Hothorn and Tamhane and Xi already report it as
            non-robust, and their tables show it reaching 0.197 against a nominal
            0.025.

WHERE THE APPROXIMATIONS ARE IN METHOD 3, which is the whole point. Four layers,
all stated in their own text. The paper enumerates the same four in Section 2.3;
layer 2 below was previously folded into layer 3 here, which undercounted them:
  1. Welch-Satterthwaite. U_i is a weighted sum of two chi-squares; they replace
     it by a single chi-square with nuhat_i d.f. by matching two moments.
  2. True degrees of freedom for the denominator correlation. Their Appendix
     derives corr(U_i,U_j) after replacing the ESTIMATED d.f. by the true ones,
     justified as tending to one as the sample sizes grow.
  3. The correlation cannot be reproduced. Their nested construction gives
     corr(U_i,U_j) = c_i/sqrt(nuhat_i nuhat_j), and they write: "This formula
     does not have a product structure like (8), and so, it is not possible to
     use it to derive expressions for c_i and d_i. Instead, we propose to use
     c_i = eta_i^2 nuhat_i". In their worked example: "This correlation matrix
     cannot be exactly reproduced via chi-squared distributed U_i's; however, we
     approximate it".
  4. Monte Carlo error on top of 1, 2 and 3.

The conditional reduction has none of the four. It never forms the joint law of
the denominators, so it never approximates it, and it uses no simulation.

HOW FWER IS MEASURED HERE, which is itself a by-product of the reduction. Not by
simulating data sets and counting rejections, which is what produces their
+/- 0.0005 at 100,000 replications. Two devices:

  Rao-Blackwellisation. Only the variance estimates are drawn. Every method's
  critical value is computed from those estimates exactly as it would be on real
  data, and the rejection probability CONDITIONAL on the drawn variances is then
  evaluated in closed form by the same conditioning argument that gives the exact
  law. The normal draws are integrated out rather than sampled.

  A control variate whose true value is known EXACTLY. The oracle arm uses the
  exact critical value at the TRUE variances, so its familywise error rate is
  alpha by construction, verified analytically per scenario as 1 - cdf(c) = alpha
  to six decimals. Every arm is run on the identical draws, so

      FWER(method) = alpha + mean( reject_method - reject_oracle )

  and the standard error is that of the PAIRED difference. Measured, the paired
  standard error at 1000 draws is 0.00006 to 0.00038, against 0.00049 for their
  100,000 data sets. The reduction therefore buys precision even where simulation
  is still used, which is a second and separate consequence of it.

FALSIFIERS. Two, both of which can fail.
  (a) The DIRECT, unpaired oracle mean must agree with alpha within its own Monte
      Carlo error. This tests the Rao-Blackwellised estimator itself.
  (b) The control-variate estimate and the direct estimate of each method must
      agree within the direct estimate's error. If the control variate were
      mis-paired, they would not.
"""

import json
import os
from math import sqrt

import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import t as student_t

from h12_hetero_compare import c_onefactor_t
from h12_hetero_exact import HeteroMaxT
from h12_quad import normal_nodes

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "series", "h12_tamhane_headtohead.json")

ALPHA = 0.025                 # one-sided, as in their tables
FAST = {"u": (8, 10), "v": (8, 10), "s": (3, 10)}
FAST_1F = dict(n_panel=12, order=10)
M3_SIMS = 40000               # inner simulation size for their Method 3

# --------------------------------------------------------------------------
# their scenarios, verbatim from Sections 5.1 and 5.2
# --------------------------------------------------------------------------
SCEN_K2 = [
    ("S1 balanced, equal sd", [10, 10, 10], [30., 30., 30.]),
    ("S2 balanced, control largest sd", [10, 10, 10], [50., 10., 10.]),
    ("S3 balanced, one trt largest sd", [10, 10, 10], [10., 10., 50.]),
    ("S4 control smallest n and sd", [4, 13, 13], [10., 10., 50.]),
    ("S5 trt smallest n, largest sd", [13, 13, 4], [10., 10., 50.]),
    ("S6 control smallest n, largest sd", [4, 13, 13], [50., 10., 10.]),
    ("S7 trt smallest n and sd", [13, 13, 4], [50., 10., 10.]),
]
SCEN_K3 = [
    ("S1 balanced, equal sd", [10, 10, 10, 10], [30., 30., 30., 30.]),
    ("S2 balanced, control largest sd", [10, 10, 10, 10], [50., 10., 10., 10.]),
    ("S3 balanced, one trt largest sd", [10, 10, 10, 10], [10., 10., 10., 50.]),
    ("S4 control smallest n and sd", [4, 12, 12, 12], [10., 10., 10., 50.]),
    ("S5 trt smallest n, largest sd", [12, 12, 12, 4], [10., 10., 10., 50.]),
    ("S6 control smallest n, largest sd", [4, 12, 12, 12], [50., 10., 10., 10.]),
    ("S7 trt smallest n and sd", [12, 12, 12, 4], [50., 10., 10., 10.]),
]
# their published Monte Carlo values, global null, one-sided alpha = 0.025
# Table 1 (k=2) and Table 3 (k=3); Method 0 kept for the record, not recomputed
PUBLISHED = {
    2: {0: [0.02507, 0.04760, 0.03967, 0.00450, 0.17313, 0.19703, 0.02195],
        1: [0.02307, 0.01552, 0.02397, 0.02613, 0.02651, 0.01713, 0.01753],
        2: [0.02516, 0.02343, 0.02507, 0.02914, 0.02715, 0.02436, 0.02509],
        3: [0.02455, 0.02315, 0.02419, 0.02646, 0.02664, 0.02501, 0.02455]},
    3: {0: [0.02501, 0.06547, 0.04882, 0.01063, 0.18052, 0.21530, 0.04250],
        1: [0.02104, 0.01234, 0.02325, 0.02390, 0.02636, 0.01301, 0.01305],
        2: [0.02441, 0.02410, 0.02573, 0.02880, 0.02800, 0.02323, 0.02426],
        3: [0.02347, 0.02371, 0.02427, 0.02596, 0.02660, 0.02341, 0.02369]},
}
PUBLISHED_SIMS = 100000


# --------------------------------------------------------------------------
# plug-in quantities, all from the sample variances as on real data
# --------------------------------------------------------------------------
def welch_df(n, s2):
    """nuhat_i, their equation (2)."""
    w0 = s2[0] / n[0]
    w = s2[1:] / n[1:]
    return (w + w0) ** 2 / (w ** 2 / (n[1:] - 1) + w0 ** 2 / (n[0] - 1))


def lambdas(n, s2):
    """lambdahat_i, their equation (7): corr(Z_i,Z_j) = lambda_i lambda_j."""
    w0 = s2[0] / n[0]
    return np.sqrt(w0 / (s2[1:] / n[1:] + w0))


def etas(n, s2):
    """etahat_i, their equation (9): corr(U_i,U_j) = eta_i eta_j."""
    num = (s2[0] / n[0]) / sqrt(n[0] - 1.0)
    den = np.sqrt(s2[0] ** 2 / (n[0] ** 2 * (n[0] - 1.0))
                  + s2[1:] ** 2 / (n[1:] ** 2 * (n[1:] - 1.0)))
    return num / den


# --------------------------------------------------------------------------
# the four critical-value rules
# --------------------------------------------------------------------------
def c_method1(dfs, alpha=ALPHA):
    """prod_j F_{nuhat_j}(c) = 1 - alpha."""
    def f(c):
        return float(np.prod(student_t.cdf(c, dfs))) - (1.0 - alpha)
    return brentq(f, 0.2, 25.0, xtol=1e-10, rtol=1e-13, maxiter=200)


def c_method2(lam, dfs, alpha=ALPHA):
    """Per-comparison: k-variate t with loadings lam and this comparison's df."""
    return np.array([c_onefactor_t(lam, float(d), alpha, **FAST_1F) for d in dfs])


def c_method3(n, s2, rng, alpha=ALPHA, sims=M3_SIMS):
    """Their algorithm, steps 1-7, reimplemented from the paper."""
    K = n.size - 1
    nuh = welch_df(n, s2)
    lam = lambdas(n, s2)
    eta = etas(n, s2)
    ci = eta ** 2 * nuh
    di = (1.0 - eta ** 2) * nuh
    order = np.argsort(ci)
    c_ord, d_ord, nu_ord = ci[order], di[order], nuh[order]

    # step 1: Z with corr = lam_i lam_j, i.e. one-factor
    y = rng.standard_normal(sims)
    z = (lam[None, :] * y[:, None]
         + np.sqrt(np.clip(1 - lam ** 2, 0, None))[None, :]
         * rng.standard_normal((sims, K)))
    # step 2: nested ladder V_[i] = V_[i-1] + chi2_{c[i]-c[i-1]}
    V = np.empty((sims, K))
    prev = np.zeros(sims)
    prev_c = 0.0
    for j in range(K):
        inc = c_ord[j] - prev_c
        if inc > 1e-12:
            prev = prev + rng.chisquare(inc, sims)
        V[:, j] = prev
        prev_c = c_ord[j]
    # step 3-4: W and U, then restore the original comparison order
    W = np.stack([rng.chisquare(max(d, 1e-12), sims) for d in d_ord], axis=1)
    U_ord = V + W
    U = np.empty_like(U_ord)
    U[:, order] = U_ord
    nu_back = np.empty(K)
    nu_back[order] = nu_ord
    # step 5-7
    T = z / np.sqrt(U / nu_back[None, :])
    return float(np.quantile(T.max(axis=1), 1.0 - alpha))


def c_exact(n, s2, alpha=ALPHA):
    """The conditional reduction at the plug-in variance ratios."""
    sig_hat = np.sqrt(s2 / s2[0])
    sig_hat[0] = 1.0
    return HeteroMaxT(n, sig_hat, grid=FAST).critical_value(alpha)


# --------------------------------------------------------------------------
# exact conditional rejection probability, given the drawn variances
# --------------------------------------------------------------------------
def cond_reject(n, sigma_true, s2, c, u_nodes):
    u, wu = u_nodes
    d = np.sqrt(s2[1:] / n[1:] + s2[0] / n[0])
    scale = np.sqrt(n[1:]) / sigma_true[1:]
    c = np.atleast_1d(np.asarray(c, float))
    if c.size == 1:
        c = np.repeat(c, n.size - 1)
    return 1.0 - float(wu @ ndtr((u[:, None] + (c * d)[None, :]) * scale[None, :])
                       .prod(axis=1))


def run_scenario(label, n_l, sig_l, n_draws, seed):
    n = np.asarray(n_l, float)
    sig = np.asarray(sig_l, float)
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    un = normal_nodes(20, 12)
    u_nodes = (un[0] * sig[0] / sqrt(n[0]), un[1])

    c_oracle = HeteroMaxT(n, sig, grid=FAST).critical_value(ALPHA)
    arms = {k: [] for k in ("oracle", "exact", "method1", "method2", "method3")}
    for _ in range(n_draws):
        s2 = (rng.chisquare(nu, size=n.size) / nu) * sig ** 2
        dfs = welch_df(n, s2)
        lam = lambdas(n, s2)
        cs = {"oracle": c_oracle,
              "exact": c_exact(n, s2),
              "method1": c_method1(dfs),
              "method2": c_method2(lam, dfs),
              "method3": c_method3(n, s2, rng)}
        for k, c in cs.items():
            arms[k].append(cond_reject(n, sig, s2, c, u_nodes))

    out = {"scenario": label, "n": n_l, "sigma": sig_l, "draws": n_draws,
           "c_oracle": c_oracle,
           "oracle_analytic_fwer": 1.0 - HeteroMaxT(n, sig, grid=FAST).cdf(c_oracle)}
    orc = np.asarray(arms["oracle"])
    for k, v in arms.items():
        v = np.asarray(v)
        direct = float(v.mean())
        dse = float(v.std(ddof=1) / sqrt(len(v)))
        if k == "oracle":
            out[k] = {"fwer": ALPHA, "mcse": 0.0, "direct": direct,
                      "direct_mcse": dse,
                      "direct_z_vs_alpha": (direct - ALPHA) / dse}
        else:
            d = v - orc
            cv = ALPHA + float(d.mean())
            out[k] = {"fwer": cv,
                      "mcse": float(d.std(ddof=1) / sqrt(len(d))),
                      "direct": direct, "direct_mcse": dse,
                      "cv_vs_direct_z": (cv - direct) / dse}
    return out


def main(n_draws=1000):
    res = {"alpha": ALPHA, "draws": n_draws, "m3_sims": M3_SIMS,
           "published_sims": PUBLISHED_SIMS, "blocks": {}}
    pub_se = sqrt(ALPHA * (1 - ALPHA) / PUBLISHED_SIMS)
    print("Head-to-head with Tamhane and Xi (2023), global null, one-sided "
          "alpha = %.3f" % ALPHA)
    print("Their tables: 100,000 simulated data sets, standard error %.5f."
          % pub_se)
    print("Here: %d Rao-Blackwellised draws, all arms on identical draws.\n"
          % n_draws)

    for K, scen in ((2, SCEN_K2), (3, SCEN_K3)):
        print("=" * 104)
        print("k = %d treatments versus control   (their Table %d)"
              % (K, 1 if K == 2 else 3))
        print("%-34s | %8s %8s %8s %8s | %8s %8s %8s"
              % ("scenario", "EXACT", "M1", "M2", "M3",
                 "M1 pub", "M2 pub", "M3 pub"))
        rows = []
        for i, (label, n_l, sig_l) in enumerate(scen):
            r = run_scenario(label, n_l, sig_l, n_draws, seed=4200 + 31 * K + i)
            r["published"] = {str(m): PUBLISHED[K][m][i] for m in (0, 1, 2, 3)}
            m3_pub = PUBLISHED[K][3][i]
            r["m3_published_minus_recomputed"] = m3_pub - r["method3"]["fwer"]
            rows.append(r)
            print("%-34s | %8.5f %8.5f %8.5f %8.5f | %8.5f %8.5f %8.5f"
                  % (label, r["exact"]["fwer"], r["method1"]["fwer"],
                     r["method2"]["fwer"], r["method3"]["fwer"],
                     PUBLISHED[K][1][i], PUBLISHED[K][2][i], m3_pub))
        res["blocks"][str(K)] = rows

        za = max(abs(r["oracle"]["direct_z_vs_alpha"]) for r in rows)
        zb = max(abs(r[k]["cv_vs_direct_z"]) for r in rows
                 for k in ("exact", "method1", "method2", "method3"))
        aw = max(abs(r["oracle_analytic_fwer"] - ALPHA) for r in rows)
        pse = max(r[k]["mcse"] for r in rows
                  for k in ("exact", "method1", "method2", "method3"))
        print("  FALSIFIER (a) direct oracle mean vs alpha: worst |z| = %.2f" % za)
        print("  FALSIFIER (b) control-variate vs direct:   worst |z| = %.2f" % zb)
        print("  oracle exact by construction: worst |1 - cdf(c) - alpha| = %.2e"
              % aw)
        print("  worst paired standard error %.5f, against %.5f for their 100,000"
              % (pse, pub_se))
        for k, name in (("exact", "EXACT"), ("method1", "Method 1"),
                        ("method2", "Method 2"), ("method3", "Method 3")):
            f = [r[k]["fwer"] for r in rows]
            worst = max(f, key=lambda x: abs(x - ALPHA))
            print("  %-9s FWER %.5f to %.5f   worst relative error %+.1f%%"
                  % (name, min(f), max(f), 100 * (worst - ALPHA) / ALPHA))

    summ = {}
    for K in (2, 3):
        rows = res["blocks"][str(K)]
        for k in ("oracle", "exact", "method1", "method2", "method3"):
            f = [r[k]["fwer"] for r in rows]
            summ.setdefault(k, {})[str(K)] = {
                "min": min(f), "max": max(f),
                "worst_rel_err_pct": 100 * (max(f, key=lambda x: abs(x - ALPHA))
                                            - ALPHA) / ALPHA}
    res["summary"] = summ
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    main()
