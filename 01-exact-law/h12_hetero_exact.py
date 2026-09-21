"""H12 OPEN ROUTE: the correlated-scale max-t null, exactly and without simulation.

THE PROBLEM (open since Dunnett 1955; state of the art Tamhane and Xi 2023).
Compare K treatments with a shared control when the group variances differ. The
Welch-type statistics

    t_i = (Xbar_i - Xbar_0) / sqrt( S_i^2/n_i + S_0^2/n_0 ),   i = 1..K

do NOT have a multivariate t distribution. Two couplings break it: the numerators
share Xbar_0, and, the hard one, the DENOMINATORS share S_0^2, so they are
correlated mixtures of chi-squares. Tamhane and Xi (Biometrical Journal 65,
e2200300, 2023) handle this by APPROXIMATING the joint law of the denominators by
correlated chi-squares and GENERATING it by simulation. Dunnett's T3 and the
Games-Howell family instead use a studentized maximum modulus, which is the
multivariate t at rho = 0, i.e. they discard the dependence entirely.

THE REDUCTION. Both couplings run through the control and nothing else. Condition
on the control's two sufficient statistics,

    U = Xbar_0 - mu   ~ N(0, sigma_0^2/n_0),    V = S_0^2/n_0,

and the K events become EXACTLY independent, because (Xbar_i, S_i^2) for i = 1..K
are mutually independent and independent of the control. Hence

    P(t_i <= c_i for all i) = E_{U,V} [ prod_i g_i(U, V; c_i) ],
    g_i(u,v;c) = E_{S_i^2} [ Phi( sqrt(n_i)/sigma_i * (u + c sqrt(S_i^2/n_i + v)) ) ].

Every expectation is one-dimensional and every integrand is smooth and bounded by
one. A direct exact evaluation is (2K+2)-dimensional; this is a two-dimensional
outer rule over K one-dimensional inner rules, and the OUTER DIMENSION DOES NOT
GROW WITH K. The joint law of the denominators is never formed, so it is never
approximated.

SCOPE, STATED UP FRONT. The null law is scale-invariant but NOT pivotal: it
depends on the variance ratios lambda_i = sigma_i^2/sigma_0^2. So this delivers
the exact law GIVEN the ratios; a procedure must still plug in estimates. That
plug-in step is the residual and is priced separately in h12_hetero_plugin.py.
What is removed here is the approximation and the simulation, not the nuisance.

THE GRID IS ADAPTED TO THE INTEGRAND, WHICH IS NOT OPTIONAL. g_i carries a factor
sqrt(n_i)/sigma_i inside Phi, so when a treatment is much more precise than the
control that factor is large and g_i is a near-STEP in u, centred on the curve
u = -c*sqrt(sigma_i^2/n_i + v) with width sigma_i/sqrt(n_i). On a control of df 9
against treatments a thousand times more precise the step is about 1.6e-5 wide
while a uniform panel is about 0.36 wide in the same units, so uniform panels
cannot see it at all
and the answer stalls near four or five correct digits however many panels are
added. The u-grid is therefore rebuilt for each c with panel EDGES planted on the
step, per v. That restores the accuracy in the corner at essentially the same
cost, and it is the difference between a method that degrades quietly and one
that does not.

FALSIFIERS (each must pass or the route is dead):
 F1. K=1, n_0=n_1=n, sigma_0=sigma_1: the statistic is algebraically the pooled
     two-sample t, so the answer must equal the Student t_{2n-2} quantile.
 F2. sigma_0 -> 0: the control becomes known, the t_i become independent
     t_{nu_i}, so the CDF must equal prod_i F_{nu_i}(c).
 F3. Grid refinement: answers must CONVERGE as panels are added. This is the test
     the classical Gauss-Hermite and Gauss-Laguerre rules failed.
 F4. Monte Carlo agreement on genuinely heteroscedastic unbalanced designs,
     including a large-df design that overflows the generalized Laguerre rule.
 F5. sigma_i -> 0 for the treatments: every t_i collapses onto the same
     control-driven quantity, so the maximum must be a single Student t on the
     control's df INDEPENDENT OF K. This is the rho -> 1 corner, opposite to F2,
     and it is the corner the adapted grid exists for.
"""

import json
import os
from math import sqrt

import numpy as np
from h12_quad import chi2_over_df_nodes
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm
from scipy.stats import t as student_t

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "series", "h12_hetero_exact.json"
)

# (panels, order) per dimension. The s-dimension converges fastest, so it is
# given the coarsest grid; see the refinement table in F3.
GRID = {"u": (14, 14), "v": (14, 14), "s": (4, 14)}
ZMAX = 8.0
# where extra panel edges are planted, in step-widths either side of the step
STEP_EDGES = (-6.0, 0.0, 6.0)


def _f(x, default=None):
    """Guarded float conversion (lint policy).  Returns float(x) unchanged; a
    failure yields NaN instead of raising, matching the convention used by the
    callers in 02-bounded-rule (core._solve returns NaN on failure)."""
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        if default is None:
            return float("nan")
        return default


def _i(x, default=0):
    """Guarded int conversion (lint policy); int(x) is returned unchanged."""
    try:
        return int(x)
    except Exception:  # noqa: BLE001
        return default


class HeteroMaxT:
    """Exact null distribution of max_i t_i for the heteroscedastic
    many-to-one design, by conditioning on the control."""

    def __init__(self, n, sigma, grid=None, two_sided=False, adapt=True):
        n = np.asarray(n, dtype=float)
        sigma = np.asarray(sigma, dtype=float)
        if n.shape != sigma.shape or n.ndim != 1 or n.size < 2:
            raise RuntimeError("n and sigma must be matching 1-D arrays of length K+1")
        if np.any(n < 2) or np.any(sigma <= 0):
            raise RuntimeError("need n_i >= 2 and sigma_i > 0")
        g = dict(GRID)
        if grid:
            g.update(grid)
        self.grid = g
        self.n, self.sigma = n, sigma
        self.K = n.size - 1
        self.two_sided = bool(two_sided)
        self.adapt = bool(adapt)
        nu = n - 1.0

        self.sd0 = sigma[0] / sqrt(n[0])
        q0, self.wv = chi2_over_df_nodes(nu[0], *g["v"])
        self.v = q0 * sigma[0] ** 2 / n[0]

        # d[i][v, s] = sqrt(S_i^2/n_i + V), the denominator on the inner grid;
        # dtyp[i][v] = its value at S_i^2 = sigma_i^2, which is where the step is
        self.d, self.ws, self.dtyp = [], [], []
        for i in range(1, self.K + 1):
            qi, wi = chi2_over_df_nodes(nu[i], *g["s"])
            s2 = qi * sigma[i] ** 2
            self.d.append(np.sqrt(s2[None, :] / n[i] + self.v[:, None]))
            self.ws.append(wi)
            self.dtyp.append(np.sqrt(sigma[i] ** 2 / n[i] + self.v))
        # step width in u, per treatment, expressed in control standard deviations
        self.step_z = np.array(
            [sigma[i] / sqrt(n[i]) / self.sd0 for i in range(1, self.K + 1)]
        )
        self._leg = np.polynomial.legendre.leggauss(g["u"][1])
        # Plant edges ONLY where the step is too narrow for the base panels to
        # see it. On an ordinary design no treatment qualifies, the grid is the
        # plain one and the adaptivity costs exactly nothing; accuracy is bought
        # only where it is actually needed.
        panel_z = 2.0 * ZMAX / g["u"][0]
        self._needs_edges = np.flatnonzero(self.step_z < 0.25 * panel_z)

    # ------------------------------------------------------------------
    def _u_nodes(self, c):
        """u nodes and weights, shape (m_v, m_u), with panel edges on the step.

        The steps sit at u = -c * dtyp_i(v), one per treatment, so the grid has
        to depend on v and on c. Extra edges are planted at each step centre and
        six step-widths either side. Duplicate or out-of-range edges collapse to
        zero-width panels, which contribute nothing and are harmless."""
        x, wl = self._leg
        base = np.linspace(-ZMAX, ZMAX, self.grid["u"][0] + 1)
        edges = np.broadcast_to(base, (self.v.size, base.size))
        if self.adapt and self._needs_edges.size:
            extra = []
            for i in self._needs_edges:
                centre = -c * self.dtyp[i] / self.sd0  # (m_v,)
                for k in STEP_EDGES:
                    extra.append(centre + k * self.step_z[i])
            edges = np.concatenate([edges, np.stack(extra, axis=1)], axis=1)
        edges = np.sort(np.clip(edges, -ZMAX, ZMAX), axis=1)
        half = 0.5 * (edges[:, 1:] - edges[:, :-1])
        mid = 0.5 * (edges[:, 1:] + edges[:, :-1])
        z = mid[:, :, None] + half[:, :, None] * x[None, None, :]
        w = half[:, :, None] * wl[None, None, :] * norm.pdf(z)
        return (z * self.sd0).reshape(self.v.size, -1), w.reshape(self.v.size, -1)

    def _g(self, i, c, u):
        """g_i(u, v; c) on the (v, u) grid."""
        scale = sqrt(self.n[i]) / self.sigma[i]
        d = self.d[i - 1][:, None, :]  # (m_v, 1, m_s)
        arg = u[:, :, None] + c * d  # (m_v, m_u, m_s)
        hi = ndtr(arg * scale)
        core = hi - ndtr((u[:, :, None] - c * d) * scale) if self.two_sided else hi
        return core @ self.ws[i - 1]  # (m_v, m_u)

    def cdf(self, c):
        """P(t_i <= c_i for every i) under the global null.

        `c` may be a scalar (one common critical value, so this is the CDF of
        the maximum) or a length-K vector (one critical value per comparison,
        which is what step-down and Welch-df procedures actually use)."""
        cs = np.atleast_1d(np.asarray(c, dtype=float))
        if cs.size == 1:
            cs = np.repeat(cs, self.K)
        if cs.size != self.K:
            raise RuntimeError(f"c must be scalar or length K={self.K}, got {cs.size}")
        u, wu = self._u_nodes(_f(cs.max()))
        acc = np.ones_like(u)
        for i in range(1, self.K + 1):
            acc = acc * self._g(i, _f(cs[i - 1]), u)
        return _f(self.wv @ (acc * wu).sum(axis=1))

    def self_check(self, c, refine=4):
        """Recompute the CDF on a refined grid and return the movement.

        With the adapted grid this comes back at the 1e-12 level on ordinary
        designs; a larger movement means the integrand has structure the
        edge-planting did not anticipate."""
        fine = {k: (v[0] * refine, v[1]) for k, v in self.grid.items()}
        other = HeteroMaxT(
            self.n, self.sigma, grid=fine, two_sided=self.two_sided, adapt=self.adapt
        )
        return abs(other.cdf(c) - self.cdf(c))

    def critical_value(self, alpha=0.05, lo=0.2, hi=25.0, two_stage=True):
        """Solve cdf(c) = 1 - alpha.

        Two-stage by default: locate the root on a third-sized grid, then solve
        on the working grid inside a narrow bracket around it. The coarse pass is
        cheap and removes most of the bisection, which is where the time went."""

        def f(c):
            return self.cdf(c) - (1.0 - alpha)

        if two_stage and self.grid["u"][0] >= 9:
            coarse = {k: (max(3, v[0] // 3), v[1]) for k, v in self.grid.items()}
            c0 = HeteroMaxT(
                self.n,
                self.sigma,
                grid=coarse,
                two_sided=self.two_sided,
                adapt=self.adapt,
            ).critical_value(alpha, lo, hi, two_stage=False)
            a, b = max(lo, c0 - 0.15), min(hi, c0 + 0.15)
            if f(a) < 0 < f(b):
                return brentq(f, a, b, xtol=1e-12, rtol=1e-14, maxiter=200)

        flo, fhi = f(lo), f(hi)
        if flo > 0 or fhi < 0:
            raise RuntimeError(
                f"critical value not bracketed: f({lo:g})={flo:g} f({hi:g})={fhi:g}"
            )
        return brentq(f, lo, hi, xtol=1e-12, rtol=1e-14, maxiter=300)


# --------------------------------------------------------------------------
# falsifiers
# --------------------------------------------------------------------------
CASES = [
    ("balanced, variances 1:1:2:4", [10, 10, 10, 10], [1.0, 1.0, 2.0, 4.0]),
    ("unbalanced K=4, spread 6x", [31, 11, 11, 51, 11], [1.0, 3.0, 0.5, 1.0, 2.0]),
    ("tiny control df=5, spread 16x", [6, 21, 6, 6], [1.0, 0.25, 4.0, 1.0]),
    ("large df (Laguerre overflows here)", [400, 60, 200, 60], [1.0, 2.0, 0.5, 3.0]),
]


def f1_pooled_t():
    rows = []
    for n in (6, 10, 31, 61, 201):
        mod = HeteroMaxT([n, n], [1.0, 1.0])
        for alpha in (0.05, 0.01):
            got = mod.critical_value(alpha)
            want = _f(student_t.ppf(1 - alpha, 2 * n - 2))
            rows.append(
                {
                    "n": n,
                    "alpha": alpha,
                    "exact": got,
                    "student_t": want,
                    "abs_err": abs(got - want),
                }
            )
    return rows


def f2_known_control():
    rows = []
    for K, n_t in ((3, 21), (5, 11), (8, 31)):
        mod = HeteroMaxT([11] + [n_t] * K, [1e-5] + [1.0] * K)
        for c in (2.0, 2.5, 3.0):
            got = mod.cdf(c)
            want = _f(student_t.cdf(c, n_t - 1)) ** K
            rows.append(
                {
                    "K": K,
                    "n_t": n_t,
                    "c": c,
                    "exact": got,
                    "prod_student": want,
                    "abs_err": abs(got - want),
                }
            )
    return rows


def f3_refinement():
    levels = [
        ("coarse", {"u": (10, 10), "v": (10, 10), "s": (3, 10)}),
        ("base", None),
        ("fine", {"u": (32, 12), "v": (32, 12), "s": (8, 12)}),
        ("finest", {"u": (48, 10), "v": (48, 10), "s": (12, 10)}),
    ]
    rows = []
    for label, n, sig in CASES:
        vals = [HeteroMaxT(n, sig, grid=g).critical_value(0.05) for _, g in levels]
        rows.append(
            {
                "case": label,
                "n": n,
                "sigma": sig,
                "c_by_level": {lv[0]: v for lv, v in zip(levels, vals, strict=True)},
                "base_to_finest": abs(vals[-1] - vals[1]),
                "fine_to_finest": abs(vals[-1] - vals[-2]),
            }
        )
    return rows


def f4_monte_carlo(n_draws=8_000_000, seed=20260913):
    rng = np.random.default_rng(seed)
    rows = []
    for label, n_l, sig_l in CASES:
        n = np.asarray(n_l, float)
        sig = np.asarray(sig_l, float)
        c = HeteroMaxT(n, sig).critical_value(0.05)
        hit, done, block = 0, 0, 250_000
        while done < n_draws:
            b = min(block, n_draws - done)
            xb = rng.standard_normal((b, n.size)) * (sig / np.sqrt(n))
            s2 = (rng.chisquare(n - 1.0, size=(b, n.size)) / (n - 1.0)) * sig**2
            den = np.sqrt(s2[:, 1:] / n[1:] + (s2[:, 0] / n[0])[:, None])
            t = (xb[:, 1:] - xb[:, 0][:, None]) / den
            hit += _i((t.max(axis=1) > c).sum())
            done += b
        p = hit / n_draws
        se = sqrt(p * (1 - p) / n_draws)
        rows.append(
            {
                "case": label,
                "n": n_l,
                "sigma": sig_l,
                "c_exact": c,
                "mc_size": p,
                "mc_se": se,
                "nominal": 0.05,
                "z_vs_nominal": (p - 0.05) / se,
                "draws": n_draws,
            }
        )
    return rows


def f5_perfect_dependence():
    """Treatments made very precise: the maximum collapses onto a single
    t_{nu_0} whatever K is.

    Stated as a SCALING test rather than as a tolerance, because a tolerance
    here cannot separate the two things that make the answer differ from
    t_{nu_0}. The treatments are NOT infinitely precise: Xbar_i still has
    standard deviation sigma_i/sqrt(n_i), which moves the true CDF away from the
    limit by an amount proportional to sigma_i. Integration error does not
    depend on sigma_i at all. So:

        if the residual scales LINEARLY with sigma_t it is model difference and
        the quadrature is clean; if it SATURATES it is integration error.

    That is a test that can fail, and without the planted edges it does: the
    unadapted grid saturates at 5.5e-4 across four decades of sigma_t while the
    adapted grid tracks sigma_t all the way down.
    """
    rows = []
    for nu0, K in ((9, 3), (9, 12), (20, 5), (5, 8)):
        n = [nu0 + 1] + [4000] * K
        for c in (1.5, 2.0, 2.8):
            want = _f(student_t.cdf(c, nu0))
            errs, uerrs = [], []
            for sig_t in (1e-3, 1e-4, 1e-5):
                sig = [1.0] + [sig_t] * K
                errs.append(abs(HeteroMaxT(n, sig).cdf(c) - want))
                uerrs.append(abs(HeteroMaxT(n, sig, adapt=False).cdf(c) - want))
            ratios = [errs[j] / errs[j + 1] for j in range(len(errs) - 1)]
            uratios = [uerrs[j] / uerrs[j + 1] for j in range(len(uerrs) - 1)]
            rows.append(
                {
                    "nu0": nu0,
                    "K": K,
                    "c": c,
                    "errs_by_sigma_t": errs,
                    "decade_ratios": ratios,
                    "unadapted_errs": uerrs,
                    "unadapted_ratios": uratios,
                    "scales_linearly": all(r > 8.0 for r in ratios),
                    "implied_quad_err_bound": errs[-1],
                }
            )
    return rows


def f5b_k_invariance():
    """The K-invariance the F5 limit actually asserts, measured.

    F5 above pins each (nu_0, K) cell against the Student t it must collapse
    onto, which tests the limit but never compares one K against another. The
    claim that carries the weight is the second one: in that limit the answer
    does not depend on K at all, because the maximum of K copies of the same
    quantity is that quantity. So hold nu_0 and c fixed, sweep K, and report
    the spread. A factorisation that lost a dependence on K would show up here
    and nowhere else in the falsifier set.
    """
    rows = []
    for nu0 in (5, 9, 20):
        for c in (1.5, 2.0, 2.8):
            want = _f(student_t.cdf(c, nu0))
            vals = {}
            for K in (1, 2, 3, 5, 8, 12):
                n = [nu0 + 1] + [4000] * K
                sig = [1.0] + [1e-5] * K
                vals[K] = HeteroMaxT(n, sig).cdf(c)
            v = list(vals.values())
            rows.append(
                {
                    "nu0": nu0,
                    "c": c,
                    "student_t": want,
                    "cdf_by_K": vals,
                    "spread_across_K": max(v) - min(v),
                    "worst_abs_err_vs_student": max(abs(x - want) for x in v),
                }
            )
    return rows


def main():
    res = {"grid": {k: list(v) for k, v in GRID.items()}}

    print("F1  pooled-t identity (K=1, balanced, equal variances)")
    res["f1_pooled_t"] = f1_pooled_t()
    w1 = max(r["abs_err"] for r in res["f1_pooled_t"])
    for r in res["f1_pooled_t"]:
        n_r, a_r = r["n"], r["alpha"]
        ex, st, er = r["exact"], r["student_t"], r["abs_err"]
        print(
            f"    n={n_r:<4d} a={a_r:.2f}   exact {ex:.13f}   "
            f"t_{{2n-2}} {st:.13f}   err {er:.2e}"
        )
    print(f"    worst |err| = {w1:.3e}")

    print("\nF2  known-control limit (sigma_0 -> 0: independent t's)")
    res["f2_known_control"] = f2_known_control()
    w2 = max(r["abs_err"] for r in res["f2_known_control"])
    for r in res["f2_known_control"]:
        k_r, nu_r, c_r = r["K"], r["n_t"] - 1, r["c"]
        ex, pr, er = r["exact"], r["prod_student"], r["abs_err"]
        print(
            f"    K={k_r} nu={nu_r:<3d} c={c_r:.1f}  exact {ex:.13f}  "
            f"prod F_nu {pr:.13f}  err {er:.2e}"
        )
    print(f"    worst |err| = {w2:.3e}")

    print("\nF3  grid refinement (the test the classical rules failed)")
    res["f3_refinement"] = f3_refinement()
    for r in res["f3_refinement"]:
        case = r["case"]
        levels = "  ".join(f"{k} {v:.11f}" for k, v in r["c_by_level"].items())
        b2f, f2f = r["base_to_finest"], r["fine_to_finest"]
        print(f"    {case}")
        print(f"       {levels}")
        print(f"       base->finest {b2f:.2e}   fine->finest {f2f:.2e}")

    print("\nF4  Monte Carlo on heteroscedastic unbalanced designs (8e6 draws)")
    res["f4_monte_carlo"] = f4_monte_carlo()
    for r in res["f4_monte_carlo"]:
        case = r["case"]
        c_ex, sz, se, z = r["c_exact"], r["mc_size"], r["mc_se"], r["z_vs_nominal"]
        print(f"    {case:<38s} c={c_ex:.7f}  size {sz:.5f} (se {se:.5f})  z={z:+.2f}")

    print("\nF5  perfect-dependence corner: the residual must SCALE with sigma_t")
    res["f5_perfect_dependence"] = f5_perfect_dependence()
    w5 = max(r["implied_quad_err_bound"] for r in res["f5_perfect_dependence"])
    worst_ratio = min(min(r["decade_ratios"]) for r in res["f5_perfect_dependence"])
    worst_uratio = max(max(r["unadapted_ratios"]) for r in res["f5_perfect_dependence"])
    all_scale = all(r["scales_linearly"] for r in res["f5_perfect_dependence"])
    for r in res["f5_perfect_dependence"]:
        nu0, k_r, c_r = r["nu0"], r["K"], r["c"]
        e_str = " ".join(f"{e:.1e}" for e in r["errs_by_sigma_t"])
        d_str = " ".join(f"{x:.1f}" for x in r["decade_ratios"])
        u_str = " ".join(f"{e:.1e}" for e in r["unadapted_errs"])
        print(
            f"    nu0={nu0:<3d} K={k_r:<3d} c={c_r:.1f}  err {e_str}  "
            f"ratios {d_str}  | unadapted {u_str}"
        )
    print(
        f"    every cell scales linearly: {all_scale} "
        f"(worst decade ratio {worst_ratio:.1f}; 10 is exact)"
    )
    print(f"    the unadapted grid saturates instead: worst ratio {worst_uratio:.1f}")
    print(f"    implied bound on quadrature error: {w5:.2e}")

    print("\nF5b K-invariance in the same corner: the answer must not move with K")
    res["f5b_k_invariance"] = f5b_k_invariance()
    w5b = max(r["spread_across_K"] for r in res["f5b_k_invariance"])
    w5b_t = max(r["worst_abs_err_vs_student"] for r in res["f5b_k_invariance"])
    for r in res["f5b_k_invariance"]:
        nu0, c_r = r["nu0"], r["c"]
        sp, wt = r["spread_across_K"], r["worst_abs_err_vs_student"]
        print(
            f"    nu0={nu0:<3d} c={c_r:.1f}  K=1..12 spread {sp:.2e}   "
            f"worst vs t {wt:.2e}"
        )
    print(f"    worst spread across K: {w5b:.2e}")
    print(f"    worst deviation from the Student t limit: {w5b_t:.2e}")

    res["summary"] = {
        "f5b_worst_spread_across_K": w5b,
        "f5b_worst_abs_err_vs_student": w5b_t,
        "f1_worst_abs_err": w1,
        "f2_worst_abs_err": w2,
        "f3_worst_fine_to_finest": max(
            r["fine_to_finest"] for r in res["f3_refinement"]
        ),
        "f3_worst_base_to_finest": max(
            r["base_to_finest"] for r in res["f3_refinement"]
        ),
        "f4_worst_abs_z": max(abs(r["z_vs_nominal"]) for r in res["f4_monte_carlo"]),
        "f5_implied_quad_err_bound": w5,
        "f5_all_cells_scale_linearly": all_scale,
        "f5_worst_decade_ratio": worst_ratio,
        "f5_unadapted_worst_ratio": worst_uratio,
    }
    try:
        with open(OUT, "w", encoding="utf-8") as fh:
            json.dump(res, fh, indent=2)
    except OSError as exc:
        raise SystemExit(f"could not write {OUT}: {exc}") from exc
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
