"""Independent replication of the desk-guarantee method (METHOD.md sections 1-5).

Written from METHOD.md and replicate/inputs/designs.json only.

The one substantive ambiguity resolved here (see REPLICATION.md): the formula
displayed in section 1 is the max-t law *conditional on the variance ratios*.
The exact law -- the one for which F1 (Student t identity) and F2 (known-control
limit -> product of t's) are stated to hold -- additionally integrates the
chi-square variance ratios W0 = s0^2 and Y_i = s_i^2/lambda_i out.  The literal
conditional formula reproduces the standard-normal quantile, not t_{2n-2}; the
variance-integrated law reproduces t to ~1e-11.  We implement the latter.
"""

from __future__ import annotations

import concurrent.futures as cf
import functools
import itertools
import json
import math
import os
import time

import numpy as np
import numpy.polynomial.hermite as npherm
import numpy.polynomial.legendre as npleg
from scipy.special import gammaln, ndtr
from scipy.stats import chi2, t


# ----------------------------------------------------------------------------
# quadrature
# ----------------------------------------------------------------------------


@functools.cache
def gh_nodes(q):
    """E_{Z~N(0,1)}[h(Z)] = sum_j w_j h(z_j)."""
    z, w = npherm.hermgauss(q)
    return math.sqrt(2.0) * z, w / math.sqrt(math.pi)


@functools.cache
def chi2_nodes(panels, order, nu, gamma=2.5):
    """E_{X~chi2_nu}[h(X)] = sum_j w_j h(x_j).

    Composite Gauss-Legendre on v in [0,1] with the chi-square quantile map
    y = sqrt(chi2.ppf(...)) = ymax * v**gamma.  The power map concentrates
    panels near y=0, where the integrand of interest has its transition for
    small variance ratios.  Integrands are analytic in y, so this is spectrally
    accurate (plain genlaguerre / uniform panels are not).
    """
    ymax = math.sqrt(chi2.ppf(1.0 - 1e-13, nu))
    gx, gw = npleg.leggauss(order)
    edges = np.linspace(0.0, 1.0, panels + 1)
    const = 2.0 / (2.0 ** (nu / 2.0) * math.exp(gammaln(nu / 2.0)))
    xs, ws = [], []
    for a, b in zip(edges[:-1], edges[1:], strict=True):
        v = 0.5 * (b - a) * gx + 0.5 * (a + b)
        wj = 0.5 * (b - a) * gw
        y = ymax * np.power(v, gamma)
        dy = ymax * gamma * np.power(v, gamma - 1.0)
        xs.append(y**2)
        ws.append(const * np.power(y, nu - 1.0) * np.exp(-(y**2) / 2.0) * dy * wj)
    return np.concatenate(xs), np.concatenate(ws)


class Quad:
    """qu Gauss-Hermite nodes for U; panels x order chi-square nodes per
    variance axis (W0 and each Y_i)."""

    def __init__(self, qu=24, panels=8, order=10, gamma=2.5):
        self.qu = qu
        self.panels = panels
        self.order = order
        self.gamma = gamma

    def __repr__(self):
        return (
            f"Quad(qu={self.qu},panels={self.panels},order={self.order},"
            f"gamma={self.gamma})"
        )


# ----------------------------------------------------------------------------
# section 1: the exact law of max_i T_i
# ----------------------------------------------------------------------------


def exact_F_batch(c, LAM, n, Q):
    """Exact null CDF F(c, lambda) of max_i T_i, vectorised over M lambda rows."""
    LAM = np.atleast_2d(np.asarray(LAM, float))
    M, K = LAM.shape
    n = np.asarray(n, float)
    n0 = n[0]
    nu0 = n0 - 1.0
    c = np.broadcast_to(np.asarray(c, float), (M,))

    zu, wu = gh_nodes(Q.qu)
    u = zu / math.sqrt(n0)  # U0 ~ N(0,1/n0)
    W0v, ww = chi2_nodes(Q.panels, Q.order, nu0, Q.gamma)
    W0 = W0v / nu0  # W ~ chi2_nu0/nu0
    s = np.sqrt(n[1:][None, :] / LAM)  # (M,K)

    G = np.ones((M, Q.qu, W0.size))
    for i in range(K):
        nu = n[i + 1] - 1.0
        Yv, wy = chi2_nodes(Q.panels, Q.order, nu, Q.gamma)
        Y = Yv / nu  # Y_i ~ chi2_nu_i/nu_i
        D = np.sqrt(
            LAM[:, i][:, None, None] * Y[None, None, :] / n[i + 1]
            + W0[None, :, None] / n0
        )  # (M,qw,qy)
        term = s[:, i][:, None, None, None] * (
            u[None, :, None, None] + c[:, None, None, None] * D[:, None, :, :]
        )
        G *= np.tensordot(ndtr(term), wy, axes=([3], [0]))
    return np.einsum("mij,i,j->m", G, wu, ww)


def bisect_scalar(Ffun, target, lo=0.0, hi=40.0, iters=80):
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if Ffun(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cond_F_1d(c, lam, n, qu=64):
    """The literal 1-D U integral of METHOD.md section 1 (variance ratios held
    fixed at lambda).  Kept only to document the ambiguity; it returns the
    standard-normal, not the Student-t, quantile for K=1."""
    lam = np.asarray(lam, float)
    n = np.asarray(n, float)
    zu, wu = gh_nodes(qu)
    u = zu / math.sqrt(n[0])
    s = np.sqrt(n[1:] / lam)
    d = np.sqrt(lam / n[1:] + 1.0 / n[0])
    term = s * (u[:, None] + c * d[None, :])
    return float(np.sum(wu * np.prod(ndtr(term), axis=1)))


def solve_roots(LAM, n, target, Q, lo=0.0, hi=40.0, iters=32, chunk=48):
    """Vectorised bisection for roots of exact_F_batch(c, LAM, n) = target."""
    LAM = np.atleast_2d(np.asarray(LAM, float))
    M = LAM.shape[0]
    c_lo = np.full(M, float(lo))
    c_hi = np.full(M, float(hi))
    for _ in range(iters):
        mid = 0.5 * (c_lo + c_hi)
        for k in range(0, M, chunk):
            sl = slice(k, min(k + chunk, M))
            Fm = exact_F_batch(mid[sl], LAM[sl], n, Q) - target
            go_up = Fm < 0
            c_lo[sl] = np.where(go_up, mid[sl], c_lo[sl])
            c_hi[sl] = np.where(go_up, c_hi[sl], mid[sl])
    return 0.5 * (c_lo + c_hi)


# ----------------------------------------------------------------------------
# section 2: region and exact coverage
# ----------------------------------------------------------------------------


def region_bounds(lam_hat, n, scale, gamma, K):
    lam_hat = np.asarray(lam_hat, float)
    nu0 = n[0] - 1.0
    w_lo = chi2.ppf((scale * gamma / 2.0) / 2.0, nu0) / nu0
    w_hi = chi2.ppf(1.0 - (scale * gamma / 2.0) / 2.0, nu0) / nu0
    lo = np.empty(K)
    hi = np.empty(K)
    for i in range(K):
        nu = n[i + 1] - 1.0
        y_lo = chi2.ppf((scale * gamma / (2.0 * K)) / 2.0, nu) / nu
        y_hi = chi2.ppf(1.0 - (scale * gamma / (2.0 * K)) / 2.0, nu) / nu
        lo[i] = lam_hat[i] * w_lo / y_hi
        hi[i] = lam_hat[i] * w_hi / y_lo
    return lo, hi, w_lo, w_hi


def coverage(scale, n, gamma, K, panels=12, order=14):
    """Exact cov(s) (independent of lambda_hat)."""
    nu0 = n[0] - 1.0
    nu = [n[i + 1] - 1.0 for i in range(K)]
    Wv, ww = chi2_nodes(panels, order, nu0)
    W = Wv / nu0

    w_lo = chi2.ppf((scale * gamma / 2.0) / 2.0, nu0) / nu0
    w_hi = chi2.ppf(1.0 - (scale * gamma / 2.0) / 2.0, nu0) / nu0
    integrand = np.ones_like(W)
    for i in range(K):
        y_lo = chi2.ppf((scale * gamma / (2.0 * K)) / 2.0, nu[i]) / nu[i]
        y_hi = chi2.ppf(1.0 - (scale * gamma / (2.0 * K)) / 2.0, nu[i]) / nu[i]
        integrand *= chi2.cdf(nu[i] * W * (y_hi / w_lo), nu[i]) - chi2.cdf(
            nu[i] * W * (y_lo / w_hi), nu[i]
        )
    return float(np.sum(ww * integrand))


def find_scale(n, gamma, K, margin=1e-6):
    """Largest scale s with cov(s) >= 1 - gamma + margin."""
    target = 1.0 - gamma + margin
    s_hi = 1.0
    while coverage(s_hi, n, gamma, K) > target and s_hi < 1e6:
        s_hi *= 2.0
    lo, hi = 0.0, s_hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if coverage(mid, n, gamma, K) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ----------------------------------------------------------------------------
# section 3/4
# ----------------------------------------------------------------------------


def envelope_vertices(lam_hat, lo, hi):
    grid = [np.array([lo[i], lam_hat[i], hi[i]]) for i in range(len(lam_hat))]
    return np.array(list(itertools.product(*grid)))


def draw_cell(seeds, draws_per_seed, n, lam_true):
    D = len(seeds) * draws_per_seed
    K = len(lam_true)
    nu0 = n[0] - 1
    nu = np.array([n[i + 1] - 1 for i in range(K)])
    s0sq = np.empty(D)
    sisq = np.empty((D, K))
    xbd = np.empty((D, K))
    idx = 0
    for seed in seeds:
        rng = np.random.default_rng(seed)
        for _ in range(draws_per_seed):
            s0sq[idx] = rng.chisquare(nu0) / nu0
            for i in range(K):
                sisq[idx, i] = lam_true[i] * rng.chisquare(nu[i]) / nu[i]
            U0 = rng.normal(0.0, math.sqrt(1.0 / n[0]))
            for i in range(K):
                X = rng.normal(0.0, math.sqrt(lam_true[i] / n[i + 1]))
                xbd[idx, i] = X - U0
            idx += 1
    return s0sq, sisq, xbd


def run_cell(cell, consts, qparams):
    Q = Quad(**qparams)
    gamma = consts["gamma"]
    alpha = consts["alpha"]
    per_seed = consts["draws_per_seed"]
    n = np.array(cell["n"], float)
    K = len(cell["lam_true"])
    lam_true = np.array(cell["lam_true"], float)

    scale = find_scale(n, gamma, K)
    cov = coverage(scale, n, gamma, K)
    covs = [coverage(f * scale, n, gamma, K) for f in (0.25, 0.5, 1.0, 2.0, 4.0)]
    monotone = all(covs[i] >= covs[i + 1] - 1e-12 for i in range(len(covs) - 1))

    s0sq, sisq, xbd = draw_cell(cell["seeds"], per_seed, n, lam_true)
    D = len(s0sq)
    lam_hat = sisq / s0sq[:, None]
    T = xbd / np.sqrt(sisq / n[1:][None, :] + s0sq[:, None] / n[0])
    maxT = T.max(axis=1)

    c_oracle = solve_roots(lam_true[None, :], n, 1 - alpha, Q)[0]
    c_plug = solve_roots(lam_hat, n, 1 - alpha, Q)

    verts = np.empty((D, 3**K, K))
    outside = np.zeros(D, dtype=bool)
    for d in range(D):
        lo, hi, _, _ = region_bounds(lam_hat[d], n, scale, gamma, K)
        verts[d] = envelope_vertices(lam_hat[d], lo, hi)
        outside[d] = np.any((lam_true < lo) | (lam_true > hi))
    c_env = (
        solve_roots(verts.reshape(D * 3**K, K), n, 1 - alpha, Q)
        .reshape(D, 3**K)
        .max(axis=1)
    )

    rej_env = (maxT > c_env).astype(float)
    rej_plug = (maxT > c_plug).astype(float)
    rej_or = (maxT > c_oracle).astype(float)
    diff = rej_env - rej_or
    pd_mean = float(diff.mean())
    pd_se = float(diff.std(ddof=1) / math.sqrt(D))
    return {
        "family": cell["family"],
        "df0": cell["df0"],
        "n": cell["n"],
        "lam_true": cell["lam_true"],
        "scale": float(scale),
        "coverage_exact": float(cov),
        "coverage_monotone_decreasing": bool(monotone),
        "coverage_samples": [float(x) for x in covs],
        "c_env_mean": float(c_env.mean()),
        "c_plug_mean": float(c_plug.mean()),
        "c_oracle_population": float(c_oracle),
        "cost_pct": float(c_env.mean() / c_oracle - 1.0) * 100.0,
        "fwer_envelope": float(rej_env.mean()),
        "fwer_plugin": float(rej_plug.mean()),
        "fwer_oracle": float(rej_or.mean()),
        "paired_dominance_mean": pd_mean,
        "paired_dominance_se": pd_se,
        "paired_dominance_ok": bool(pd_mean + 2 * pd_se <= 0),
        "frac_lambda_true_outside_region": float(outside.mean()),
    }


# ----------------------------------------------------------------------------
# section 5
# ----------------------------------------------------------------------------


def falsifier_F1(Q):
    worst, rows = 0.0, []
    for n in (10, 31, 61):
        for alpha in (0.05, 0.01):
            tq = t.ppf(1 - alpha, 2 * n - 2)
            r = solve_roots(np.array([[1.0]]), np.array([n, n]), 1 - alpha, Q)[0]
            dev = abs(r - tq)
            worst = max(worst, dev)
            rows.append(
                {
                    "n": n,
                    "alpha": alpha,
                    "t_quantile": float(tq),
                    "root": float(r),
                    "abs_dev": float(dev),
                }
            )
    return worst, rows


def falsifier_F2(Q, cs=(2.0, 2.5, 3.0)):
    sigma0 = 1e-5
    worst, rows = 0.0, []
    n = np.array([7.0] + [4.0 + i for i in range(8)])
    for K in (3, 5, 8):
        nn = np.concatenate([n[:1], n[1 : 1 + K]])
        lam = np.full(K, 1.0 / sigma0**2)
        for c in cs:
            Fv = exact_F_batch(np.array([c]), lam[None, :], nn, Q)[0]
            ref = float(np.prod([t.cdf(c, nn[i + 1] - 1.0) for i in range(K)]))
            worst = max(worst, abs(Fv - ref))
            rows.append(
                {
                    "K": K,
                    "c": c,
                    "F": float(Fv),
                    "ref_t": ref,
                    "abs_dev": float(abs(Fv - ref)),
                }
            )
    return worst, rows


def falsifier_F3(n, lam, Q, target=0.975):
    roots_u = {}
    for qu in (10, 32, 48):
        roots_u[qu] = float(
            solve_roots(
                np.array([lam]),
                np.array(n),
                target,
                Quad(qu=qu, panels=Q.panels, order=Q.order, gamma=Q.gamma),
            )[0]
        )
    roots_p = {}
    for p in (10, 32, 48):
        roots_p[p] = float(
            solve_roots(
                np.array([lam]),
                np.array(n),
                target,
                Quad(qu=Q.qu, panels=p, order=Q.order, gamma=Q.gamma),
            )[0]
        )
    return (
        abs(roots_u[32] - roots_u[48]),
        abs(roots_p[32] - roots_p[48]),
        roots_u,
        roots_p,
    )


def falsifier_F4(n, lam_true, N, alpha, Q, seed=20240607):
    c_star = solve_roots(np.array([lam_true]), np.array(n), 1 - alpha, Q)[0]
    K = len(lam_true)
    nu0 = n[0] - 1
    nu = np.array([n[i + 1] - 1 for i in range(K)])
    rng = np.random.default_rng(seed)
    s0sq = rng.chisquare(nu0, N) / nu0
    sisq = np.empty((N, K))
    for i in range(K):
        sisq[:, i] = lam_true[i] * rng.chisquare(nu[i], N) / nu[i]
    U0 = rng.normal(0.0, math.sqrt(1.0 / n[0]), N)
    T = np.empty((N, K))
    for i in range(K):
        X = rng.normal(0.0, math.sqrt(lam_true[i] / n[i + 1]), N)
        T[:, i] = (X - U0) / np.sqrt(sisq[:, i] / n[i + 1] + s0sq / n[0])
    maxT = T.max(axis=1)
    phat = float((maxT > c_star).mean())
    return {
        "c_star": float(c_star),
        "N": int(N),
        "phat": phat,
        "alpha": alpha,
        "z_score": float((phat - alpha) / math.sqrt(alpha * (1 - alpha) / N)),
    }


# ----------------------------------------------------------------------------


def _worker(args):
    cell, consts, qparams = args
    t0 = time.time()
    rec = run_cell(cell, consts, qparams)
    rec["_seconds"] = time.time() - t0
    return rec


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    with open(os.path.join(root, "inputs", "designs.json")) as fh:
        designs = json.load(fh)
    consts = designs["constants"]
    qparams = {"qu": 24, "panels": 8, "order": 10, "gamma": 2.5}
    Q = Quad(**qparams)
    QF1 = Quad(qu=48, panels=10, order=12, gamma=2.5)
    t_start = time.time()
    print(f"desk quadrature {Q}; falsifier quadrature {QF1}", flush=True)

    print("== section 4: desk (parallel over cells) ==", flush=True)
    workers = min(12, len(designs["desk_cells"]))
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        cells = list(
            ex.map(_worker, [(c, consts, qparams) for c in designs["desk_cells"]])
        )
    for c in cells:
        print(
            f"  {c['family']} df0={c['df0']} scale={c['scale']:.2f} "
            f"cov={c['coverage_exact']:.7f} c_env={c['c_env_mean']:.3f} "
            f"c_plug={c['c_plug_mean']:.3f} c_or={c['c_oracle_population']:.3f} "
            f"fwer_env={c['fwer_envelope']:.3f} fwer_pl={c['fwer_plugin']:.3f} "
            f"pd={c['paired_dominance_mean']:+.4f}+-{c['paired_dominance_se']:.4f} "
            f"out={c['frac_lambda_true_outside_region']:.3f} "
            f"({c['_seconds']:.0f}s)",
            flush=True,
        )

    print("== section 5: falsifiers ==", flush=True)
    f1_worst, f1_rows = falsifier_F1(QF1)
    f1_deskq, _ = falsifier_F1(Q)
    print(
        f"F1 worst abs dev = {f1_worst:.3e} (desk quadrature {f1_deskq:.3e})",
        flush=True,
    )

    lit = []
    for n in (10, 31, 61):
        for alpha in (0.05, 0.01):
            tq = t.ppf(1 - alpha, 2 * n - 2)
            r = bisect_scalar(lambda c, n=n: cond_F_1d(c, [1.0], [n, n]), 1 - alpha)
            lit.append(
                {
                    "n": n,
                    "alpha": alpha,
                    "literal_root": float(r),
                    "t_quantile": float(tq),
                    "abs_dev": float(abs(r - tq)),
                }
            )

    f2_worst, f2_rows = falsifier_F2(Q)
    print(f"F2 worst abs dev = {f2_worst:.3e}", flush=True)

    f3_u, f3_p, f3_roots_u, f3_roots_p = falsifier_F3([4, 13, 13], [1.0, 25.0], Q)
    print(
        f"F3 U-axis movement 32->48 = {f3_u:.3e}; chi-panel movement = {f3_p:.3e}",
        flush=True,
    )

    f4 = falsifier_F4([6, 13, 13], [1.0, 25.0], 2_000_000, consts["alpha"], Q)
    print(
        f"F4 phat={f4['phat']:.5f} z={f4['z_score']:.3f} c*={f4['c_star']:.4f}",
        flush=True,
    )

    lit_f2 = []
    nn = np.array([7, 4, 5, 6], float)
    for c in (2.0, 2.5, 3.0):
        Fv = cond_F_1d(c, np.full(3, 1e10), nn)
        ref = float(np.prod([t.cdf(c, nn[i + 1] - 1) for i in range(3)]))
        lit_f2.append({"c": c, "F_literal": Fv, "ref_t": ref})

    results = {
        "meta": {
            "law": (
                "exact max-t null law: N(0,1/n0) benchmark-mean integral x "
                "chi-square variance-ratio (W0,Y_i) integrals"
            ),
            "quadrature_desk": repr(Q),
            "quadrature_falsifiers": repr(QF1),
            "candidate_set": "3^K product of {lo,mid,hi}; K=2 in every desk cell",
            "draws_per_cell": consts["draws_per_cell"],
        },
        "per_cell": cells,
        "falsifiers": {
            "F1_worst_abs_dev_vs_student_t": f1_worst,
            "F1_worst_abs_dev_desk_quadrature": f1_deskq,
            "F1_rows": f1_rows,
            "F1_literal_1d_formula": lit,
            "F2_worst_abs_dev_vs_independent_t": f2_worst,
            "F2_rows": f2_rows,
            "F2_literal_1d_formula_L": lit_f2,
            "F3_movement_32_to_48_Uaxis": f3_u,
            "F3_movement_32_to_48_chi_panels": f3_p,
            "F3_roots_Uaxis": f3_roots_u,
            "F3_roots_chi_panels": f3_roots_p,
            "F4_monte_carlo_size_z": f4["z_score"],
            "F4_detail": f4,
        },
    }
    out = os.path.join(here, "replication_results.json")
    with open(out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"wrote {out}  ({time.time() - t_start:.1f}s total)", flush=True)


if __name__ == "__main__":
    main()
