"""What the shipped approximations cost, measured against the exact law.

h12_hetero_exact.py computes the exact null distribution of the heteroscedastic
many-to-one max-t by conditioning on the control. That makes a question
answerable that previously could only be simulated: how far wrong is each
approximation the field actually ships? Every approximate rule below is fed back
through the EXACT law to get the familywise error rate it really delivers. Those
numbers are the finding; everything else is machinery.

THE RULES COMPARED, all on the same Welch statistic
t_i = (Xbar_i - Xbar_0)/sqrt(S_i^2/n_i + S_0^2/n_0):

  EXACT        one constant c from the conditional reduction.
  T3           Dunnett's T3 / Games-Howell: compare each t_i against the
               studentized maximum modulus at that comparison's own
               Welch-Satterthwaite df. The SMM treats the K statistics as
               independent normals over ONE common scale.
  SIDAK-t      the fully independent reading: c_i = F_{nu_i*}^{-1}((1-a)^{1/K}).
               This is the rule H12's own deposit computed, and the contrast
               with T3 is worth stating: a shared denominator leaves the
               statistics dependent even when their numerators are uncorrelated,
               so SMM and Sidak-t are different objects and neither is the other.
  ONE-FACTOR t Dunnett's own computation carried across unchanged: keep the
               correlation, force a single df. The numerator correlation here is
               EXACTLY one-factor, rho_ij = h_i h_j with h_i = sqrt(w_0/(w_i+w_0))
               and w_i = sigma_i^2/n_i, so this is the honest generalization of
               Dunnett's rho = 1/2 and not a straw version of it. Two df
               conventions are priced: the pooled df sum(nu_i) the homoscedastic
               procedure uses, and the mean Welch df.

CHECKS ON THE COMPARISON ARMS, so the contrast is not an artefact of a bad rival:
  C1 one-factor t at UNEQUAL loadings against direct Monte Carlo of the same
     multivariate t. This is the case C3 and C4 do not reach.
  C2 one-factor t at K = 1 must reproduce t_nu exactly.
  C3 one-factor t at equal h = sqrt(1/2), two-sided, must reproduce the
     published Dunnett critical value 2.901 for k = 6 groups including the
     control at df = 12 and alpha = .05 (Dunnett 1964 table, reproduced and
     confirmed in text at real-statistics.com/statistics-tables/dunnetts-table).
  C4 one-factor t must reproduce the p-values of scipy.stats.dunnett, an
     independent implementation, to that routine's own quasi-Monte-Carlo error.

A NOTE ON WHAT IS NOT CLAIMED. The exact law depends on the variance ratios, so
it is exact GIVEN them. Every rule here, the exact one included, needs those
ratios plugged in to be used on data. The plug-in residual is a separate
question and is priced separately; it is not swept into these numbers.
"""

import json
import os
from math import sqrt

import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import t as student_t

from h12_hetero_exact import HeteroMaxT
from h12_quad import chi2_over_df_nodes, normal_nodes

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "series", "h12_hetero_compare.json")

# the one published value the source states in text rather than inside an image:
# k = 6 groups INCLUDING the control (so K = 5 treatments), df = 12,
# alpha = .05, TWO-SIDED.
DUNNETT_PUBLISHED_TWO_SIDED = {"K": 5, "nu": 12, "alpha": 0.05, "c": 2.901}
ALPHA = 0.05


def welch_df(n, sigma, i):
    """Welch-Satterthwaite df for the comparison of treatment i with control."""
    w0 = sigma[0] ** 2 / n[0]
    wi = sigma[i] ** 2 / n[i]
    return (wi + w0) ** 2 / (wi ** 2 / (n[i] - 1) + w0 ** 2 / (n[0] - 1))


def factor_loadings(n, sigma):
    """h_i with rho_ij = h_i h_j: the exact numerator correlation."""
    w0 = sigma[0] ** 2 / n[0]
    w = sigma[1:] ** 2 / n[1:]
    return np.sqrt(w0 / (w + w0))


# --------------------------------------------------------------------------
# the rival critical values
# --------------------------------------------------------------------------
def smm_cdf(c, K, df, n_panel=40, order=14):
    """One-sided studentized maximum modulus: E_s[Phi(c s)^K], s^2 = chi2_df/df."""
    s2, ws = chi2_over_df_nodes(df, n_panel, order)
    return float(ws @ ndtr(c * np.sqrt(s2)) ** K)


def c_smm(K, df, alpha=ALPHA, n_panel=40, order=14):
    return brentq(lambda c: smm_cdf(c, K, df, n_panel, order) - (1 - alpha),
                  0.2, 25.0, xtol=1e-12, rtol=1e-14, maxiter=300)


def c_sidak_t(df, K, alpha=ALPHA):
    return float(student_t.ppf((1 - alpha) ** (1.0 / K), df))


def onefactor_t_cdf(c, h, df, n_panel=40, order=14, two_sided=False):
    """P(max_i t_i <= c) for a one-factor multivariate t with loadings h."""
    y, wy = normal_nodes(n_panel, order)
    s2, ws = chi2_over_df_nodes(df, n_panel, order)
    s = np.sqrt(s2)
    h = np.asarray(h, float)
    root = np.sqrt(np.clip(1.0 - h ** 2, 1e-300, None))
    base = -h[None, None, :] * y[None, :, None]
    hi = ndtr((c * s[:, None, None] + base) / root)
    core = hi - ndtr((-c * s[:, None, None] + base) / root) if two_sided else hi
    return float(ws @ core.prod(axis=2) @ wy)


def c_onefactor_t(h, df, alpha=ALPHA, two_sided=False, n_panel=40, order=14):
    return brentq(
        lambda c: onefactor_t_cdf(c, h, df, n_panel, order, two_sided) - (1 - alpha),
        0.2, 25.0, xtol=1e-12, rtol=1e-14, maxiter=300)


# --------------------------------------------------------------------------
# checks on the comparison arms
# --------------------------------------------------------------------------
def c1_unequal_loadings(n_draws=4_000_000, seed=451):
    """Unequal h against direct Monte Carlo of the same multivariate t."""
    rng = np.random.default_rng(seed)
    rows = []
    for h_l, df in (([0.2, 0.5, 0.8, 0.95], 17.0),
                    ([0.9, 0.9, 0.3], 9.0),
                    ([0.1, 0.4, 0.6, 0.7, 0.99], 40.0)):
        h = np.asarray(h_l)
        c = c_onefactor_t(h, df)
        y = rng.standard_normal(n_draws)
        z = rng.standard_normal((n_draws, h.size))
        s = np.sqrt(rng.chisquare(df, n_draws) / df)
        num = h[None, :] * y[:, None] + np.sqrt(1 - h ** 2)[None, :] * z
        p = float((( num / s[:, None]).max(axis=1) > c).mean())
        se = sqrt(p * (1 - p) / n_draws)
        rows.append({"h": h_l, "df": df, "c": c, "mc_tail": p, "mc_se": se,
                     "nominal": ALPHA, "z": (p - ALPHA) / se})
    return rows


def c2_single():
    rows = []
    for nu in (8, 20, 60):
        got = c_onefactor_t(np.array([0.7]), nu)
        want = float(student_t.ppf(1 - ALPHA, nu))
        rows.append({"nu": nu, "onefactor": got, "student_t": want,
                     "abs_err": abs(got - want)})
    return rows


def c3_dunnett_published():
    d = DUNNETT_PUBLISHED_TWO_SIDED
    got = c_onefactor_t(np.full(d["K"], sqrt(0.5)), float(d["nu"]), d["alpha"],
                        two_sided=True)
    return {"K": d["K"], "nu": d["nu"], "alpha": d["alpha"], "onefactor": got,
            "published": d["c"], "abs_err": abs(got - d["c"]),
            "rounds_to_published": round(got, 3) == d["c"]}


def c4_scipy_dunnett(seed=7):
    from scipy.stats import dunnett
    rng = np.random.default_rng(seed)
    rows = []
    for K, n in ((4, 6), (3, 12), (6, 9)):
        ctrl = rng.standard_normal(n)
        grps = [rng.standard_normal(n) + 0.3 * j for j in range(K)]
        r = dunnett(*grps, control=ctrl, alternative="greater")
        nu = float(sum(len(g) - 1 for g in [ctrl] + grps))
        sp2 = sum((len(g) - 1) * np.var(g, ddof=1) for g in [ctrl] + grps) / nu
        tst = [(g.mean() - ctrl.mean()) / sqrt(sp2 * (1 / len(g) + 1 / n))
               for g in grps]
        mine = [1 - onefactor_t_cdf(t, np.full(K, sqrt(0.5)), nu) for t in tst]
        rows.append({"K": K, "n": n, "scipy_p": [float(x) for x in r.pvalue],
                     "onefactor_p": [float(x) for x in mine],
                     "max_abs_p_diff": max(abs(a - b)
                                           for a, b in zip(r.pvalue, mine))})
    return rows


# --------------------------------------------------------------------------
# the comparison
# --------------------------------------------------------------------------
DESIGNS = [
    ("balanced n=10, homoscedastic", [10, 10, 10, 10], [1.0, 1.0, 1.0, 1.0]),
    ("balanced n=10, sd 1:1:2:4", [10, 10, 10, 10], [1.0, 1.0, 2.0, 4.0]),
    ("balanced n=20, sd 1:3:3:3", [20, 20, 20, 20], [1.0, 3.0, 3.0, 3.0]),
    ("balanced n=20, sd 3:1:1:1 noisy control", [20, 20, 20, 20],
     [3.0, 1.0, 1.0, 1.0]),
    ("small control n0=6, sd 1:.25:4:1", [6, 21, 6, 6], [1.0, 0.25, 4.0, 1.0]),
    ("unbalanced K=4, sd spread 6x", [31, 11, 11, 51, 11],
     [1.0, 3.0, 0.5, 1.0, 2.0]),
    ("K=8 balanced n=15, graded sd", [15] * 9,
     [1.0, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0]),
    ("large df n0=400, sd 1:2:.5:3", [400, 60, 200, 60], [1.0, 2.0, 0.5, 3.0]),
]


def compare(alpha=ALPHA):
    rows = []
    for label, n_l, sig_l in DESIGNS:
        n = np.asarray(n_l, float)
        sig = np.asarray(sig_l, float)
        K = n.size - 1
        mod = HeteroMaxT(n, sig)
        c_ex = mod.critical_value(alpha)
        dfs = np.array([welch_df(n, sig, i) for i in range(1, K + 1)])
        h = factor_loadings(n, sig)

        rules = {
            "T3_smm_welch_df": np.array([c_smm(K, d, alpha) for d in dfs]),
            "sidak_t_welch_df": np.array([c_sidak_t(d, K, alpha) for d in dfs]),
            "onefactor_pooled_df": c_onefactor_t(h, float((n - 1).sum()), alpha),
            "onefactor_mean_welch_df": c_onefactor_t(h, float(dfs.mean()), alpha),
        }
        out = {"design": label, "n": n_l, "sigma": sig_l, "K": K,
               "welch_df": [round(float(d), 3) for d in dfs],
               "h": [round(float(x), 4) for x in h], "c_exact": c_ex}
        for name, c in rules.items():
            fwer = 1.0 - mod.cdf(c)
            cc = np.atleast_1d(c)
            out[name] = {
                "c": [round(float(x), 6) for x in cc] if cc.size > 1
                     else round(float(cc[0]), 6),
                "true_fwer": fwer,
                "fwer_rel_err_pct": 100 * (fwer - alpha) / alpha}
        rows.append(out)
    return rows


def main():
    res = {"alpha": ALPHA}

    print("C1  one-factor t at UNEQUAL loadings vs direct Monte Carlo (4e6)")
    res["c1_unequal_loadings"] = c1_unequal_loadings()
    for r in res["c1_unequal_loadings"]:
        print("    h=%-28s df=%-5.1f c=%.7f  MC tail %.5f (se %.5f)  z=%+.2f"
              % (r["h"], r["df"], r["c"], r["mc_tail"], r["mc_se"], r["z"]))

    print("\nC2  one-factor t at K=1 must be t_nu")
    res["c2_single"] = c2_single()
    for r in res["c2_single"]:
        print("    nu=%-3d  %.13f vs %.13f   err %.2e"
              % (r["nu"], r["onefactor"], r["student_t"], r["abs_err"]))

    print("\nC3  one-factor t vs the published Dunnett value (two-sided)")
    r3 = c3_dunnett_published()
    res["c3_dunnett_published"] = r3
    print("    K=%d nu=%d a=%.2f  computed %.6f  published %.3f  diff %.5f  "
          "rounds to published: %s"
          % (r3["K"], r3["nu"], r3["alpha"], r3["onefactor"], r3["published"],
             r3["abs_err"], r3["rounds_to_published"]))

    print("\nC4  one-factor t vs scipy.stats.dunnett (independent implementation)")
    res["c4_scipy_dunnett"] = c4_scipy_dunnett()
    for r in res["c4_scipy_dunnett"]:
        print("    K=%d n=%-3d  max |p difference| = %.2e"
              % (r["K"], r["n"], r["max_abs_p_diff"]))

    print("\nTRUE FWER OF EACH RULE AT NOMINAL 0.05 (computed exactly)")
    res["comparison"] = compare()
    names = ["T3_smm_welch_df", "sidak_t_welch_df", "onefactor_pooled_df",
             "onefactor_mean_welch_df"]
    print("%-42s %8s %9s %9s %9s %9s"
          % ("design", "c_exact", "T3", "Sidak-t", "1f pooled", "1f Welch"))
    for r in res["comparison"]:
        print("%-42s %8.5f %9.5f %9.5f %9.5f %9.5f"
              % (r["design"], r["c_exact"], *[r[k]["true_fwer"] for k in names]))

    res["summary"] = {
        "c1_worst_abs_z": max(abs(r["z"]) for r in res["c1_unequal_loadings"]),
        "c2_worst": max(r["abs_err"] for r in res["c2_single"]),
        "c3_abs_err_vs_published": r3["abs_err"],
        "c4_worst_p_diff": max(r["max_abs_p_diff"] for r in res["c4_scipy_dunnett"]),
    }
    for k in names:
        f = [r[k]["true_fwer"] for r in res["comparison"]]
        rel = [100 * (x - ALPHA) / ALPHA for x in f]
        worst = max(rel, key=abs)
        res["summary"][k] = {"min_fwer": min(f), "max_fwer": max(f),
                             "min_rel_err_pct": min(rel), "max_rel_err_pct": max(rel),
                             "worst_signed_rel_err_pct": worst}
        print("  %-26s FWER %.5f to %.5f   relative error %+.1f%% to %+.1f%%"
              % (k, min(f), max(f), min(rel), max(rel)))

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    main()
