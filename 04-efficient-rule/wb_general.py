"""Workstream B at general K: the uniformly calibrated plug-in rule.

Same construction as `wb_k1.py`, carried to the candidate counts a desk uses.
The rule is one scalar:

    reject candidate i  iff  T_i > c*( lambda_hat ; alpha0 ) ,

where c* is the paper's exact fixed-ratio constant and alpha0 is the largest
internal level whose UNCONDITIONAL familywise error rate stays at or below the
ceiling for every ratio vector in the declared domain.

Above K = 1 the outer expectation over the variance estimates is sampled rather
than quadratured, so the FWER here is Rao-Blackwellised (exact inner, sampled
outer) and the tier stays measured.  Two devices keep that honest:

  * common random numbers: one state sample is shared by every lambda and every
    candidate level, so the bisection sees a smooth deterministic surface;
  * the constraint is applied to FWER + 2 * mcse, not to the point estimate.

COST.  The expensive object is c*(lambda_hat) per draw per lambda.  Rather than
root-find it for every candidate alpha0, the law is built ONCE per (lambda,
draw) and its distribution function is tabulated on a grid of c; the critical
value at any level is then an interpolation, so the bisection over alpha0 is
free.  That is the only reason this sweep is affordable.

Deposit: deposits/wb_general.json
"""

from __future__ import annotations

import time

import numpy as np

import fcore as F
from h12_hetero_exact import HeteroMaxT

BOOT = {"u": (4, 8), "v": (4, 8), "s": (3, 8)}
CGRID = np.linspace(1.4, 9.0, 18)
D_DRAWS = 500
SEED = 20260919

# the declared nuisance domain: patterns spanning three decades of ratio in
# every direction a desk can meet, declared before the run
def lam_patterns(K):
    out = []
    for lev in (0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0):
        out.append(np.full(K, lev))
    if K >= 2:
        for lo, hi in ((0.04, 25.0), (0.1, 10.0), (0.01, 1.0), (1.0, 100.0)):
            out.append(np.exp(np.linspace(np.log(lo), np.log(hi), K)))
        v = np.ones(K)
        v[0] = 25.0
        out.append(v.copy())
        v = np.ones(K)
        v[-1] = 0.04
        out.append(v.copy())
    return [np.asarray(x, float) for x in out]


DESIGNS = [
    ("K2_df3", [4.0, 13.0, 13.0]),
    ("K2_df5", [6.0, 13.0, 13.0]),
    ("K2_df7", [8.0, 13.0, 13.0]),
    ("K2_df11", [12.0, 13.0, 13.0]),
    ("K2_df19", [20.0, 13.0, 13.0]),
    ("K2_df35", [36.0, 13.0, 13.0]),
    ("K3_df11", [12.0, 13.0, 13.0, 13.0]),
    ("K3_df35", [36.0, 13.0, 13.0, 13.0]),
    ("K5_df11", [12.0, 12.0, 12.0, 12.0, 12.0, 12.0]),
    ("K5_df35", [36.0, 13.0, 13.0, 13.0, 13.0, 13.0]),
]


def _c_of_alpha(cdf_row, alpha):
    """Invert a tabulated distribution function for the level-alpha constant."""
    y = 1.0 - np.asarray(cdf_row, float)      # rejection probability, decreasing
    if y[0] <= alpha:
        return float(CGRID[0])
    if y[-1] >= alpha:
        return float(CGRID[-1])
    j = int(np.searchsorted(-y, -alpha))
    j = max(1, min(j, y.size - 1))
    y0, y1 = y[j - 1], y[j]
    t = 0.0 if y0 == y1 else (y0 - alpha) / (y0 - y1)
    return float(CGRID[j - 1] + t * (CGRID[j] - CGRID[j - 1]))


def build_tables(n, lams, u, quad):
    """cdf[l][d, :] = P(max T <= CGRID | law at that draw's estimated ratios)."""
    n = np.asarray(n, float)
    out = []
    for lam in lams:
        lh = lam[None, :] * u[:, 1:] / u[:, :1]
        rows = np.empty((u.shape[0], CGRID.size))
        for d in range(u.shape[0]):
            m = HeteroMaxT(n, np.sqrt(np.concatenate([[1.0], lh[d]])),
                           grid=BOOT)
            rows[d] = [m.cdf(float(c)) for c in CGRID]
        out.append(rows)
    return out


def fwer_at(alpha0, n, lam, rows, u, quad):
    cs = np.array([_c_of_alpha(rows[d], alpha0) for d in range(u.shape[0])])
    r = np.array([F.cond_reject(cs[d], u[d], n, lam, quad)
                  for d in range(u.shape[0])])
    return float(r.mean()), float(r.std(ddof=1) / np.sqrt(r.size)), cs


def run_design(name, n, abar=F.ABAR, D=D_DRAWS, seed=SEED):
    t0 = time.time()
    n = np.asarray(n, float)
    K = n.size - 1
    quad = F._u_nodes(48, 10)
    u = F.draw_states(n, D, seed)
    lams = lam_patterns(K)
    tabs = build_tables(n, lams, u, quad)

    def sup_fwer(a0):
        best, rec = -1.0, None
        for lam, rows in zip(lams, tabs):
            f, se, _ = fwer_at(a0, n, lam, rows, u, quad)
            if f + 2.0 * se > best:
                best, rec = f + 2.0 * se, (lam, f, se)
        return best, rec

    s0, rec0 = sup_fwer(abar)
    lo, hi = 0.004, abar
    if s0 <= abar:
        a0 = abar
    else:
        for _ in range(16):
            m = 0.5 * (lo + hi)
            if sup_fwer(m)[0] <= abar:
                lo = m
            else:
                hi = m
        a0 = lo
    s1, rec1 = sup_fwer(a0)

    # price it in power at the effect where the population oracle reaches 0.8
    rows_out = []
    for idx, lam in enumerate(lams[:7]):
        c_or = F.oracle_c_pop(n, lam, abar, F.FAST)

        def opow(dl):
            return float(np.mean([F.cond_power(c_or, u[d], n, lam, dl, quad)
                                  for d in range(D)]))

        a, b = 0.0, 1.0
        while opow(b) < 0.8 and b < 40:
            b *= 2
        for _ in range(28):
            m = 0.5 * (a + b)
            if opow(m) < 0.8:
                a = m
            else:
                b = m
        dl = 0.5 * (a + b)
        cs = np.array([_c_of_alpha(tabs[idx][d], a0) for d in range(D)])
        rp = float(np.mean([F.cond_power(cs[d], u[d], n, lam, dl, quad)
                            for d in range(D)]))
        rows_out.append({"lam": lam.tolist(), "delta80": float(dl),
                         "oracle_power": 0.8, "rule_power": rp,
                         "loss_pp": float(100.0 * (0.8 - rp)),
                         "mean_c": float(cs.mean()),
                         "oracle_c_pop": float(c_or)})
    losses = [r["loss_pp"] for r in rows_out]
    return {
        "design": name, "n": n.tolist(), "K": int(K), "draws": int(D),
        "seed": int(seed), "abar": abar,
        "plugin_at_abar_sup_fwer": float(s0),
        "plugin_at_abar_rel_pct": float(100.0 * (s0 - abar) / abar),
        "alpha0": float(a0),
        "level_price_pct": float(100.0 * (abar - a0) / abar),
        "sup_fwer_at_alpha0": float(s1),
        "worst_lam": rec1[0].tolist(), "worst_fwer": float(rec1[1]),
        "worst_mcse": float(rec1[2]),
        "rows": rows_out,
        "median_loss_pp": float(np.median(losses)),
        "max_loss_pp": float(max(losses)),
        "seconds": float(time.time() - t0),
    }


def main():
    out = {"abar": F.ABAR, "draws": D_DRAWS, "seed": SEED,
           "c_grid": CGRID.tolist(), "inner_grid": BOOT, "designs": []}
    print("%-10s %9s %9s %9s %8s %8s %6s" %
          ("design", "supF@abar", "alpha0", "supF@a0", "lossMed", "lossMax",
           "sec"))
    for name, n in DESIGNS:
        r = run_design(name, n)
        out["designs"].append(r)
        print("%-10s %9.5f %9.5f %9.5f %8.1f %8.1f %6.0f" %
              (name, r["plugin_at_abar_sup_fwer"], r["alpha0"],
               r["sup_fwer_at_alpha0"], r["median_loss_pp"], r["max_loss_pp"],
               r["seconds"]))
    out["all_within_ceiling"] = all(d["sup_fwer_at_alpha0"] <= F.ABAR + 1e-9
                                    for d in out["designs"])
    out["worst_loss_pp"] = max(d["max_loss_pp"] for d in out["designs"])
    out["note"] = (
        "Rao-Blackwellised outer expectation on common random numbers; the "
        "constraint is applied to FWER + 2*mcse.  Measured tier: the supremum "
        "is over a declared pattern set, not over the continuum."
    )
    F.write_deposit("wb_general.json", out, "unconditional_candidate_measured")
    print("\nall within ceiling:", out["all_within_ceiling"],
          "| worst loss %.1f pp" % out["worst_loss_pp"])


if __name__ == "__main__":
    main()
