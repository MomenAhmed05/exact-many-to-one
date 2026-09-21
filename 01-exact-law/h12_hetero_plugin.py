"""Does the exact law survive the plug-in step? The honest scope test.

The exact null distribution from h12_hetero_exact.py is exact GIVEN the variance
ratios lambda_i = sigma_i^2/sigma_0^2. On data those are not given. Every rule in
this family, the exact one included, must estimate them, and the estimates are
built from the SAME sample variances that appear in the statistics. So the
question that decides whether the exact law is useful rather than merely true is:

    once lambda is estimated, does the exact rule still beat the approximations,
    or does plug-in noise swallow the difference?

If it swallows it, the exact law is a theoretical result with no practical gain,
and this file says so.

METHOD. Rao-Blackwellised, not raw simulation. For each draw of the scale vector
(S_0^2, ..., S_K^2), the critical values are computed from the estimates, and the
rejection probability CONDITIONAL on those scales is then evaluated in closed
form by the same conditioning argument that drives the whole route:

    P(reject | S_0..S_K) = 1 - E_U [ prod_i Phi( sqrt(n_i)/sigma_i * (U + c_i d_i) ) ],
    d_i = sqrt(S_i^2/n_i + S_0^2/n_0),   U ~ N(0, sigma_0^2/n_0),

with the TRUE sigma_i in the evaluation and the ESTIMATED ones inside c_i. Only
the scales are simulated; the normal draws are integrated out exactly. That cuts
the Monte Carlo error by roughly an order of magnitude for the same draw count,
which matters because the differences being measured are a few thousandths.

ARMS
  oracle_exact   exact critical value at the TRUE ratios. Must return alpha; this
                 is the falsifier for the machinery, not a result.
  plugin_exact   exact critical value at the ESTIMATED ratios.
  plugin_T3      studentized maximum modulus at the ESTIMATED Welch df, which is
                 what T3 actually does on data.
  plugin_1f_pool one-factor t at the estimated loadings with pooled df.
"""

import json
import os
from math import sqrt

import numpy as np
from scipy.special import ndtr

from h12_hetero_compare import ALPHA, c_onefactor_t, c_smm
from h12_hetero_exact import HeteroMaxT
from h12_quad import normal_nodes

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "series", "h12_hetero_plugin.json")

# coarse grid for the inner critical-value solves: measured against the
# reference grid at 5.2e-6 in c, which moves a size by about 5e-7.
FAST = {"u": (8, 10), "v": (8, 10), "s": (3, 10)}
FAST_1F = dict(n_panel=12, order=10)

DESIGNS = [
    ("balanced n=10, sd 1:1:2:4", [10, 10, 10, 10], [1.0, 1.0, 2.0, 4.0]),
    ("balanced n=20, sd 3:1:1:1 noisy control", [20, 20, 20, 20],
     [3.0, 1.0, 1.0, 1.0]),
    ("small control n0=6, sd 1:.25:4:1", [6, 21, 6, 6], [1.0, 0.25, 4.0, 1.0]),
]


def cond_reject(n, sigma_true, s2, c, u_nodes):
    """P(reject | scales), exactly, for a per-comparison critical vector c."""
    u, wu = u_nodes
    d = np.sqrt(s2[1:] / n[1:] + s2[0] / n[0])
    scale = np.sqrt(n[1:]) / sigma_true[1:]
    c = np.atleast_1d(np.asarray(c, float))
    if c.size == 1:
        c = np.repeat(c, n.size - 1)
    arg = (u[:, None] + (c * d)[None, :]) * scale[None, :]
    return 1.0 - float(wu @ ndtr(arg).prod(axis=1))


def welch_df_from(n, s2):
    w0 = s2[0] / n[0]
    w = s2[1:] / n[1:]
    return (w + w0) ** 2 / (w ** 2 / (n[1:] - 1) + w0 ** 2 / (n[0] - 1))


def run_design(label, n_l, sig_l, n_draws=1200, seed=90913, alpha=ALPHA):
    n = np.asarray(n_l, float)
    sig = np.asarray(sig_l, float)
    K = n.size - 1
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    u_nodes = normal_nodes(20, 12)
    u_nodes = (u_nodes[0] * sig[0] / sqrt(n[0]), u_nodes[1])

    c_oracle = HeteroMaxT(n, sig, grid=FAST).critical_value(alpha)
    pooled_df = float(nu.sum())

    acc = {k: [] for k in ("oracle_exact", "plugin_exact", "plugin_T3",
                           "plugin_1f_pool")}
    c_spread = []
    for _ in range(n_draws):
        s2 = (rng.chisquare(nu, size=n.size) / nu) * sig ** 2
        sig_hat = np.sqrt(s2 / s2[0])           # scale-invariant: sigma_0_hat = 1
        sig_hat[0] = 1.0
        c_pl = HeteroMaxT(n, sig_hat, grid=FAST).critical_value(alpha)
        dfh = welch_df_from(n, s2)
        c_t3 = np.array([c_smm(K, float(d), alpha, **FAST_1F) for d in dfh])
        w0h = s2[0] / n[0]
        h_hat = np.sqrt(w0h / (s2[1:] / n[1:] + w0h))
        c_1f = c_onefactor_t(h_hat, pooled_df, alpha, **FAST_1F)
        c_spread.append(c_pl)
        for k, c in (("oracle_exact", c_oracle), ("plugin_exact", c_pl),
                     ("plugin_T3", c_t3), ("plugin_1f_pool", c_1f)):
            acc[k].append(cond_reject(n, sig, s2, c, u_nodes))

    out = {"design": label, "n": n_l, "sigma": sig_l, "K": K, "draws": n_draws,
           "alpha": alpha, "c_oracle": c_oracle,
           "c_plugin_mean": float(np.mean(c_spread)),
           "c_plugin_sd": float(np.std(c_spread, ddof=1))}
    for k, v in acc.items():
        v = np.asarray(v)
        out[k] = {"fwer": float(v.mean()),
                  "mcse": float(v.std(ddof=1) / sqrt(len(v))),
                  "rel_err_pct": float(100 * (v.mean() - alpha) / alpha)}
    return out


def main():
    res = {"alpha": ALPHA, "grid": FAST, "designs": []}
    hdr = ("%-40s %9s %9s %9s %9s" %
           ("design", "oracle", "plug exact", "plug T3", "plug 1f"))
    print(hdr)
    for label, n_l, sig_l in DESIGNS:
        r = run_design(label, n_l, sig_l)
        res["designs"].append(r)
        print("%-40s %9.5f %9.5f %9.5f %9.5f"
              % (label, r["oracle_exact"]["fwer"], r["plugin_exact"]["fwer"],
                 r["plugin_T3"]["fwer"], r["plugin_1f_pool"]["fwer"]))
        print("%-40s %9.5f %9.5f %9.5f %9.5f   (MCSE)"
              % ("", r["oracle_exact"]["mcse"], r["plugin_exact"]["mcse"],
                 r["plugin_T3"]["mcse"], r["plugin_1f_pool"]["mcse"]))
        print("     plug-in critical value: mean %.5f  sd %.5f  (oracle %.5f)"
              % (r["c_plugin_mean"], r["c_plugin_sd"], r["c_oracle"]))

    worst_oracle = max(abs(d["oracle_exact"]["fwer"] - ALPHA) /
                       d["oracle_exact"]["mcse"] for d in res["designs"])
    res["summary"] = {
        "oracle_worst_abs_z": worst_oracle,
        "plugin_exact_fwer": [d["plugin_exact"]["fwer"] for d in res["designs"]],
        "plugin_T3_fwer": [d["plugin_T3"]["fwer"] for d in res["designs"]],
        "plugin_1f_pool_fwer": [d["plugin_1f_pool"]["fwer"] for d in res["designs"]],
    }
    print("\nFALSIFIER: the oracle arm must return alpha. Worst |z| = %.2f"
          % worst_oracle)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
