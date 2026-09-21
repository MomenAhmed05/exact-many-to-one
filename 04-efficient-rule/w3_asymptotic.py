"""W3: what the asymptotic and resampling alternatives actually deliver.

The roadmap asks for a comparison against HAC standard errors and block
bootstraps.  The value of that comparison is not to pick a winner; it is the
same instrument the preceding paper used on the 2023 method.  These procedures
are the ones a desk would otherwise reach for, they are justified by limit
theorems, and at the sample sizes a monthly strategy panel actually has, the
limit has not arrived.  Measuring by how much is a result.

Data: AR(1) with a lag-one coefficient the procedure does not know, under the
null of equal means, one candidate against one benchmark.  Target size 0.025,
one-sided.

  exact, declared phi   the paper's rule after Cochrane-Orcutt whitening at a
                        DECLARED coefficient.  This is the reference: it is
                        entitled to be exact and the measurement checks that it
                        is.  It is NOT a fair contest against the others, which
                        are not given the coefficient.
  exact, phi estimated  the same rule with the coefficient estimated on the
                        evaluation data, separately in each group and pooled
                        across them.  THIS is the fair comparison, because it
                        is what a desk that does not know the coefficient
                        actually has to run.
  Newey-West, normal    HAC long-run variance in the denominator, standard
                        normal critical value.
  Newey-West, t         the same with a t_{n-1} critical value, which is what
                        most software does.
  moving block          nonparametric block bootstrap of the statistic on
                        centred data, block length chosen by the usual n^(1/3).
  stationary bootstrap  geometric block lengths, same mean block length.

Deposit: deposits/w3_asymptotic.json
"""

from __future__ import annotations

import time

import numpy as np
from scipy.stats import norm, t as tdist

import fcore as F
import k1lab as L

ABAR = F.ABAR
LAM_DOMAIN = (1e-3, 1e3)
MU = 2.0     # the common mean; the estimated-coefficient variants only differ
             # from each other when it is non-zero


def ar1(D, n, phi, sig, rng):
    e = rng.normal(size=(D, n))
    x = np.empty((D, n))
    x[:, 0] = e[:, 0]
    for t in range(1, n):
        x[:, t] = phi * x[:, t - 1] + np.sqrt(1.0 - phi ** 2) * e[:, t]
    return sig * x


def nw_lrv(x, L_lag):
    """Newey-West long-run variance, per row."""
    z = x - x.mean(1, keepdims=True)
    n = x.shape[1]
    g = (z * z).sum(1) / n
    for j in range(1, L_lag + 1):
        gj = (z[:, j:] * z[:, :-j]).sum(1) / n
        g = g + 2.0 * (1.0 - j / (L_lag + 1.0)) * gj
    return np.maximum(g, 1e-12)


def nw_lag(n):
    return max(1, int(4.0 * (n / 100.0) ** (2.0 / 9.0)))


def co_whiten(x, phi):
    return (x[:, 1:] - phi * x[:, :-1]) / np.sqrt(1.0 - phi ** 2)


def acf1(x):
    """Lag-one autocorrelation per row, clipped to a stationary range."""
    z = x - x.mean(1, keepdims=True)
    return np.clip((z[:, 1:] * z[:, :-1]).sum(1)
                   / np.maximum((z * z).sum(1), 1e-12), -0.95, 0.95)


def co_whiten_rows(x, phi):
    phi = np.asarray(phi, float)
    if phi.ndim == 1:
        phi = phi[:, None]
    return (x[:, 1:] - phi * x[:, :-1]) / np.sqrt(1.0 - phi ** 2)


def welch_t(b, c):
    n0, n1 = b.shape[1], c.shape[1]
    s0, s1 = b.var(1, ddof=1), c.var(1, ddof=1)
    return (c.mean(1) - b.mean(1)) / np.sqrt(s1 / n1 + s0 / n0), s1 / s0


def hac_t(b, c):
    n0, n1 = b.shape[1], c.shape[1]
    v = nw_lrv(b, nw_lag(n0)) / n0 + nw_lrv(c, nw_lag(n1)) / n1
    return (c.mean(1) - b.mean(1)) / np.sqrt(v)


def block_indices(D, n, bl, rng, stationary=False):
    """Index array (D, n) of moving-block or stationary-bootstrap draws."""
    if stationary:
        # vectorised: a geometric restart mask, then each position takes its
        # own block's start plus its offset within that block, wrapped
        newblk = rng.random((D, n)) < (1.0 / bl)
        newblk[:, 0] = True
        pos = np.arange(n)[None, :]
        start_pos = np.maximum.accumulate(np.where(newblk, pos, 0), axis=1)
        fresh = rng.integers(0, n, size=(D, n))
        s_at = np.take_along_axis(fresh, start_pos, axis=1)
        return (s_at + (pos - start_pos)) % n
    nb = int(np.ceil(n / bl))
    starts = rng.integers(0, n - bl + 1, size=(D, nb))
    off = np.arange(bl)
    idx = (starts[:, :, None] + off[None, None, :]).reshape(D, nb * bl)
    return idx[:, :n]


def boot_size(phi, n0, n1, lam, D, B, rng, stationary=False, alpha=ABAR,
              chunk=600):
    bl0 = max(2, int(round(n0 ** (1.0 / 3.0))))
    bl1 = max(2, int(round(n1 ** (1.0 / 3.0))))
    hits = done = 0
    while done < D:
        m = min(chunk, D - done)
        b = ar1(m, n0, phi, 1.0, rng)
        c = ar1(m, n1, phi, np.sqrt(lam), rng)
        t_obs, _ = welch_t(b, c)
        bz = b - b.mean(1, keepdims=True)
        cz = c - c.mean(1, keepdims=True)
        tb = np.empty((m, B))
        for r in range(B):
            i0 = block_indices(m, n0, bl0, rng, stationary)
            i1 = block_indices(m, n1, bl1, rng, stationary)
            bb = np.take_along_axis(bz, i0, axis=1)
            cc = np.take_along_axis(cz, i1, axis=1)
            tb[:, r], _ = welch_t(bb, cc)
        crit = np.quantile(tb, 1.0 - alpha, axis=1)
        hits += int((t_obs > crit).sum())
        done += m
    p = hits / D
    return p, float(np.sqrt(p * (1 - p) / D))


def main():
    t0 = time.time()
    rng = np.random.default_rng(31337)
    designs = [("n36", 36, 36), ("n120", 120, 120), ("n36_13", 36, 13)]
    phis = (0.0, 0.2, 0.4, 0.6)
    lam = 1.0
    D_fast, D_boot, B = 400000, 12000, 299

    rows = []
    print("%-8s %5s %-22s %9s %9s %8s" %
          ("design", "phi", "method", "size", "mcse", "z"))
    for dtag, n0, n1 in designs:
        n = np.array([float(n0), float(n1)])
        cal = L.calibrate_level(n, F.log_grid(*LAM_DOMAIN, 25), ABAR,
                                npanel=12, order=8, iters=10)
        tab = cal["table"]
        # the declared-phi reference needs one extra observation, since
        # Cochrane-Orcutt discards the first
        for phi in phis:
            def add(method, p, se):
                z = (p - ABAR) / se
                rows.append({"design": dtag, "n": [n0, n1], "phi": float(phi),
                             "method": method, "size": p, "mcse": se,
                             "z": float(z), "breach": bool(z > 3.0)})
                print("%-8s %5.1f %-22s %9.5f %9.5f %+8.1f %s"
                      % (dtag, phi, method, p, se, z,
                         "BREACH" if z > 3.0 else ""))

            hits = hitsn = hitst = hsep = hpool = done = 0
            while done < D_fast:
                m = min(100000, D_fast - done)
                b = ar1(m, n0 + 1, phi, 1.0, rng)
                c = ar1(m, n1 + 1, phi, np.sqrt(lam), rng)
                tw, lh = welch_t(co_whiten(b, phi), co_whiten(c, phi))
                hits += int((tw > tab(lh)).sum())
                th = hac_t(b[:, 1:], c[:, 1:])
                hitsn += int((th > norm.isf(ABAR)).sum())
                hitst += int((th > tdist.isf(ABAR, min(n0, n1) - 1)).sum())
                # the estimated-coefficient variants, on the SAME draws, with
                # a non-zero common mean since that is what separates them
                bb, cc = b + MU, c + MU
                p0, p1 = acf1(bb), acf1(cc)
                ts, ls = welch_t(co_whiten_rows(bb, p0),
                                 co_whiten_rows(cc, p1))
                hsep += int((ts > tab(ls)).sum())
                pp = np.clip(0.5 * (p0 + p1), -0.95, 0.95)
                tp, lp = welch_t(co_whiten_rows(bb, pp),
                                 co_whiten_rows(cc, pp))
                hpool += int((tp > tab(lp)).sum())
                done += m
            for method, h in (("the exact rule, coefficient declared", hits),
                              ("the exact rule, coefficient per group", hsep),
                              ("the exact rule, coefficient pooled", hpool),
                              ("Newey-West, normal critical value", hitsn),
                              ("Newey-West, t critical value", hitst)):
                p = h / D_fast
                add(method, p, float(np.sqrt(p * (1 - p) / D_fast)))

            p, se = boot_size(phi, n0, n1, lam, D_boot, B, rng, False)
            add("moving-block bootstrap", p, se)
            p, se = boot_size(phi, n0, n1, lam, D_boot, B, rng, True)
            add("stationary bootstrap", p, se)

    def worst(method):
        v = [r for r in rows if r["method"] == method]
        w = max(v, key=lambda r: r["size"])
        return w

    methods = ["the exact rule, coefficient declared", "the exact rule, coefficient per group",
               "the exact rule, coefficient pooled", "Newey-West, normal critical value",
               "Newey-West, t critical value", "moving-block bootstrap", "stationary bootstrap"]
    out = {
        "abar": ABAR, "draws_fast": D_fast, "draws_boot": D_boot,
        "boot_reps": B, "phis": list(phis), "rows": rows,
        "methods": methods,
        "worst_size": dict((mth, float(worst(mth)["size"])) for mth in methods),
        "worst_ratio_to_abar": dict(
            (mth, float(worst(mth)["size"] / ABAR)) for mth in methods),
        "breaches": dict(
            (mth, int(sum(1 for r in rows
                          if r["method"] == mth and r["breach"])))
            for mth in methods),
        "cells_per_method": int(len(rows) / len(methods)),
        "mu": MU,
        "note": ("The declared-coefficient row is a REFERENCE, not a "
                 "competitor: it is given something the others are not, and "
                 "its job is to confirm the instrument reads correctly.  The "
                 "contest is between the rows that all have to estimate the "
                 "coefficient from the same data, and the common mean is set "
                 "to %.1f because the separate and pooled variants coincide "
                 "when it is zero." % MU),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("w3_asymptotic.json", out, "measured")
    print("\nworst size by method:")
    for mth in methods:
        print("  %-22s %.5f  (%.2fx abar, %d of %d cells breach)"
              % (mth, out["worst_size"][mth], out["worst_ratio_to_abar"][mth],
                 out["breaches"][mth], out["cells_per_method"]))
    print("%.0f s" % out["seconds"])


if __name__ == "__main__":
    main()
