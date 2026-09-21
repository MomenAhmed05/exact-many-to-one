"""Stage 1: the sup-at-candidate step, replaced by a rigorous enclosure.

The conditional no-rejection probability at a draw with denominators `d` and true
ratios `lam` is

    F(c, lam) = E_U[ prod_i Phi(s_i(lam) (U + c d_i)) ],   s_i = sqrt(n_i / lam_i),

with `U ~ N(0, 1/n_0)` (`u = z / sqrt(n_0)` in `env_rule._rej_at_lam`) and
`d_i` the drawn standard error of contrast i.

Over a region `Lambda = prod_i [lam_lo_i, lam_hi_i]` the *s-range* is
`s_i in [s_lo_i, s_hi_i] = [sqrt(n_i/lam_hi_i), sqrt(n_i/lam_lo_i)]`.

Monotone pointwise minimum.  For a fixed `U` and `c`, `Phi(s x)` with `x = U + c d_i`
is increasing in `s` when `x >= 0` and decreasing in `s` when `x < 0`.  Hence

    Psi_i(U; c) = Phi(s_lo_i x) if x >= 0 else Phi(s_hi_i x)

is the pointwise minimum of `Phi(s x)` over the region's `s`-range, and
`G(U; c) = prod_i Psi_i(U; c)` is continuous and INCREASING in U (each `Psi_i` is a
continuity-matched pair of increasing branches: `Phi(s_hi x) < 1/2` for `x < 0`,
`Phi(s_lo x) > 1/2` for `x > 0`, and both branches equal `1/2` at `x = 0`).
Therefore `G(U; c) <= F(c, lam)` for EVERY `lam` in the region, and the monotone
step bound

    F_lo(c) = sum_j G(U_j; c) [Phi(z_{j+1}) - Phi(z_j)] + G(U_M; c) [1 - Phi(z_M)]

satisfies `F_lo(c) <= E[G(.; c)] <= F(c, lam)` (a Riemann lower sum for an increasing
integrand, with the upper tail bounded by its value at the last node; the lower tail
contributes at least 0).  Consequently the smallest `c` with `F_lo(c) >= 1 - alpha`
satisfies `c >= c_star(lam)` for EVERY `lam` in the region: a certified envelope over
the continuum, with no candidate set and no Lipschitz constant.

This is the "interval arithmetic over the quadrature itself" route, and it is strictly
stronger than a covering-radius bound because it never needs the supremum to be
attained anywhere.

Also provided: `dcdloglam`, the implicit-function derivative of `c_star` with respect
to `log lam_i` (the quantity the Lipschitz route would have bounded), used only as a
diagnostic for how far the envelope has to reach.
"""
import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ALPHA, cond_noreject, draw_var, oracle_c  # noqa: E402
from env_rule import lam_region  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "cert")
FAMS = {"A": (0.04, 0.04), "B": (1.0, 25.0)}
DFS = (3, 5, 7, 11, 19, 35)
GAMMA = 0.0005
try:
    NAN = float("nan")
except Exception:  # noqa: BLE001
    NAN = 0.0


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


M_NODES = _env_int("W26_CERT_M", 800)
POOL_N = _env_int("W26_POOL", 4)


def s_bounds(n, lam_lo, lam_hi):
    """Region s-range: s_i = sqrt(n_i / lam_i), so lam_hi gives s_lo."""
    n = np.asarray(n, float)
    lam_lo = np.asarray(lam_lo, float)
    lam_hi = np.asarray(lam_hi, float)
    return np.sqrt(n[1:] / lam_hi), np.sqrt(n[1:] / lam_lo)


def _z_grid(m):
    """Ascending z-nodes (uniform; the weight of each cell is a normal increment)."""
    return np.linspace(-9.0, 9.0, m + 1)


def _G_at_c(z, c, d_row, s_lo, s_hi, n0):
    """G(U_j; c) = prod_i min_{s in [s_lo,s_hi]} Phi(s (U_j + c d_i)), vectorised."""
    from scipy.special import ndtr

    u = z / np.sqrt(n0)
    acc = np.ones(z.size)
    for i in range(len(s_lo)):
        x = u + c * d_row[i]
        s = np.where(x >= 0.0, s_lo[i], s_hi[i])
        acc = acc * ndtr(s * x)
    return acc


def F_lower(c, z, d_row, s_lo, s_hi, n0, weight_pairs):
    """Monotone step lower bound on min over the region of F(c, lam)."""
    from scipy.special import ndtr

    try:
        g = _G_at_c(z, c, d_row, s_lo, s_hi, n0)
        mid = float(np.dot(g[:-1], weight_pairs))
        tail = float(g[-1] * (1.0 - ndtr(z[-1])))
        return mid + tail
    except Exception:  # noqa: BLE001
        return float("nan")


def _weights(z):
    from scipy.special import ndtr

    return ndtr(z[1:]) - ndtr(z[:-1])


def c_rigorous(n, d_row, lam_lo, lam_hi, alpha=ALPHA, m=M_NODES, c_lo=0.05,
               c_hi=200.0, tol=1e-8, grow=30, z=None, w=None):
    """Smallest c with F_lo(c) >= 1 - alpha: a certified envelope over the region."""
    n = np.asarray(n, float)
    d_row = np.asarray(d_row, float).reshape(-1)
    if z is None:
        z = _z_grid(m)
    if w is None:
        w = _weights(z)
    s_lo, s_hi = s_bounds(n, lam_lo, lam_hi)
    try:
        f_lo = F_lower(c_lo, z, d_row, s_lo, s_hi, n[0], w)
        if f_lo >= 1.0 - alpha:
            return float(c_lo)
        hi = c_hi
        for _ in range(grow):
            if F_lower(hi, z, d_row, s_lo, s_hi, n[0], w) >= 1.0 - alpha:
                break
            hi *= 2.0
        else:
            return float("nan")
        lo = c_lo
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if F_lower(mid, z, d_row, s_lo, s_hi, n[0], w) >= 1.0 - alpha:
                hi = mid
            else:
                lo = mid
            if hi - lo <= tol:
                break
        return float(hi)
    except Exception:  # noqa: BLE001
        return float("nan")


def dcdloglam(n, d_row, lam, c, m=200):
    """Implicit-function diagnostic: dc*/d log lam_i at a point (node quadrature).

    dF/dc        = E_U[ P * sum_i (phi/Phi)(a_i) * s_i d_i ]                ( > 0 )
    dF/dloglam_i = E_U[ P * (phi/Phi)(a_i) * (U + c d_i) * (-s_i/2) ]
    dc*/dloglam_i = -(dF/dloglam_i) / (dF/dc)
    """
    from scipy.special import ndtr

    z = _z_grid(m)
    w = _weights(z)
    n = np.asarray(n, float)
    d_row = np.asarray(d_row, float).reshape(-1)
    lam = np.asarray(lam, float)
    u = z / np.sqrt(n[0])
    P = np.ones(z.size)
    terms_dc = np.zeros(z.size)
    terms_dlam = np.zeros((len(lam), z.size))
    for i in range(len(lam)):
        s = np.sqrt(n[i + 1] / lam[i])
        a = s * (u + c * d_row[i])
        phi = np.exp(-0.5 * a * a) / np.sqrt(2.0 * np.pi)
        ph = ndtr(a)
        ratio = np.where(ph > 0.0, phi / np.maximum(ph, 1e-300), 0.0)
        P = P * ph
        terms_dc += ratio * s * d_row[i]
        terms_dlam[i] = ratio * (u + c * d_row[i]) * (-s / 2.0)
    try:
        num_dc = float(np.dot(P * terms_dc, w))
        den = num_dc if abs(num_dc) > 1e-300 else float("nan")
        return [float(-np.dot(P * terms_dlam[i], w) / den)
                for i in range(len(lam))]
    except Exception:  # noqa: BLE001
        return []


def certify_draws(n, s0, si, lam_lo, lam_hi, alpha=ALPHA, pool=None, m=M_NODES):
    """Certified c for each draw (per-draw denominators)."""
    z = _z_grid(m)
    w = _weights(z)
    out = np.empty(len(s0))
    for j in range(len(s0)):
        d_row = np.sqrt(si[j] / n[1:] + s0[j] / n[0])
        out[j] = c_rigorous(n, d_row, lam_lo[j], lam_hi[j], alpha=alpha, m=m,
                            z=z, w=w)
    return out


def cell(fam, df, draws, alpha=ALPHA, m=M_NODES, pool=None):
    """One desk cell: candidate envelope vs the certified envelope."""
    lam = np.asarray(FAMS[fam], float)
    n = np.array([df + 1, 13, 13], float)
    s0, si, lh = draw_var(n, lam, draws, 4242 + 97 * df)
    lo, hi, diag = lam_region(n, lh, GAMMA, scale=None)
    c_cand = None
    try:
        from env_rule import c_env

        c_cand = c_env(n, s0, si, lh, gamma=GAMMA, alpha=alpha, pool=pool)
    except Exception as exc:  # noqa: BLE001
        print(f"  c_env failed: {exc!r}", flush=True)
    c_rg = certify_draws(n, s0, si, lo, hi, alpha=alpha, pool=pool, m=m)
    c_or = oracle_c(n, lam, alpha)
    try:
        gap = float(np.nanmean(c_rg) - np.nanmean(c_cand))
        worst_gap = float(np.nanmax(c_rg - c_cand))
        mean_cand = float(np.nanmean(c_cand))
        mean_rg = float(np.nanmean(c_rg))
        min_rg = float(np.nanmin(c_rg))
        max_rg = float(np.nanmax(c_rg))
        n_nan = int(np.sum(~np.isfinite(c_rg)))
        deriv = dcdloglam(n, np.sqrt(si[0] / n[1:] + s0[0] / n[0]), lam, c_or)
    except Exception:  # noqa: BLE001
        gap = worst_gap = NAN
        mean_cand = mean_rg = min_rg = max_rg = NAN
        n_nan = 0
        deriv = []
    return {
        "family": fam, "df0": df, "lam": list(lam), "n": list(n),
        "draws": i64(draws), "m_nodes": i64(m), "alpha": alpha, "gamma": GAMMA,
        "c_oracle": f64(c_or),
        "c_candidate_mean": f64(mean_cand),
        "c_rigorous_mean": f64(mean_rg),
        "c_rigorous_min": f64(min_rg),
        "c_rigorous_max": f64(max_rg),
        "n_nan": i64(n_nan),
        "mean_gap": f64(gap), "worst_gap": f64(worst_gap),
        "mean_gap_pct_of_oracle": f64(100 * gap / c_or) if np.isfinite(gap) else NAN,
        "dcdloglam_at_oracle": [f64(v) for v in deriv],
        "region_scale": f64(diag.get("scale", 1.0)),
        "region_coverage": f64(diag.get("coverage_exact", NAN)),
    }


def falsify_conditional(n, s0, si, lam_lo, lam_hi, alpha=ALPHA, m=M_NODES,
                        n_lam=5):
    """Direct check of the certified rule's defining property.

    `c_rigorous` is built so that for every draw and every `lam` in the region the
    CONDITIONAL rejection probability at `c_rigorous` is at most `alpha` (up to the
    monotone step bound's slack).  This measures that directly: for each draw,
    compute `c_rigorous`, then evaluate the conditional rejection on a dense
    log-`lam` grid inside the region and report the worst excess over `alpha`.
    """
    n = np.asarray(n, float)
    out = []
    for j in range(len(s0)):
        try:
            d_row = np.sqrt(si[j] / n[1:] + s0[j] / n[0])
            c_rg = c_rigorous(n, d_row, lam_lo[j], lam_hi[j], alpha=alpha, m=m)
            grids = [np.geomspace(lam_lo[j][i], lam_hi[j][i], n_lam)
                     for i in range(len(lam_lo[j]))]
            pts = []
            for combo in np.ndindex(*([n_lam] * len(grids))):
                pts.append([grids[i][combo[i]] for i in range(len(grids))])
            pts = np.asarray(pts, float)
            worst = -np.inf
            for lam_pt in pts:
                rej = 1.0 - float(cond_noreject(
                    np.array([c_rg]), np.array([s0[j]]), si[j:j + 1], n, lam_pt)[0])
                worst = max(worst, rej)
            out.append((f64(c_rg), f64(worst), i64(len(pts))))
        except Exception:  # noqa: BLE001
            out.append((NAN, NAN, i64(0)))
    return out


def main():
    t0 = time.time()
    rows = []
    draws = _env_int("W26_CERT_D", 60)
    fals_draws = _env_int("W26_CERT_FALS", 4)
    with mp.Pool(POOL_N) as pool:
        for df in DFS:
            for fam in FAMS:
                r = cell(fam, df, draws, pool=pool)
                rows.append(r)
                lam = np.asarray(FAMS[fam], float)
                n = np.array([df + 1, 13, 13], float)
                s0, si, lh = draw_var(n, lam, fals_draws, 991 + 31 * df)
                lo, hi, _ = lam_region(n, lh, GAMMA, scale=None)
                fals = falsify_conditional(n, s0, si, lo, hi)
                try:
                    worst_rej = float(np.nanmax([f[1] for f in fals]))
                except Exception:  # noqa: BLE001
                    worst_rej = NAN
                r["falsify_max_rej"] = f64(worst_rej)
                try:
                    r["falsify_n_lam_points"] = i64(sum(int(f[2]) for f in fals))
                except Exception:  # noqa: BLE001
                    r["falsify_n_lam_points"] = i64(0)
                r["falsify_draws"] = i64(fals_draws)
                r["falsify_passes"] = bool(worst_rej <= ALPHA + 1e-9)
                print(f"{fam} df{df:>2}: oracle {r['c_oracle']:.4f} | cand "
                      f"{r['c_candidate_mean']:.4f} | certified {r['c_rigorous_mean']:.4f} "
                      f"| price {r['mean_gap_pct_of_oracle']:+.1f}% | maxcond "
                      f"{worst_rej:.4f} <= {ALPHA} {r['falsify_passes']}",
                      flush=True)
    try:
        max_gap_pct = float(max(r["mean_gap_pct_of_oracle"] for r in rows))
        dominates = bool(all(r["worst_gap"] >= -1e-9 for r in rows))
        worst_falsify = float(max(r.get("falsify_max_rej", NAN) for r in rows))
    except Exception:  # noqa: BLE001
        max_gap_pct = NAN
        dominates = False
        worst_falsify = NAN
    summary = {
        "rows": rows, "alpha": ALPHA, "gamma": GAMMA, "m_nodes": i64(M_NODES),
        "draws_per_cell": i64(draws),
        "max_mean_gap_pct": f64(max_gap_pct),
        "all_rigorous_dominate_candidate": dominates,
        "all_falsify_pass": bool(all(r.get("falsify_passes", False) for r in rows)),
        "worst_falsify_rej": f64(worst_falsify),
        "falsify_draws_per_cell": i64(fals_draws),
        "seconds": f64(time.time() - t0),
        "note": ("c_rigorous is the smallest c whose monotone step lower bound on "
                 "min over the region of F(c, lam) reaches 1 - alpha; it dominates "
                 "c_star(lam) for EVERY lam in the region, so it needs no candidate "
                 "set and no Lipschitz constant."),
    }
    dump_json(summary, os.path.join(OUT, "rigorous.json"))
    print(f"max mean gap {summary['max_mean_gap_pct']:+.2f}% of oracle | "
          f"dominates candidate: {summary['all_rigorous_dominate_candidate']} | "
          f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
