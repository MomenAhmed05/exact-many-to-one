"""Figures for the H12 paper. Every series is read from a deposited JSON or
recomputed from the deposited code, so no figure carries a number that is not
reproducible from the repository."""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtr

from h12_hetero_compare import factor_loadings
from h12_hetero_exact import HeteroMaxT
from h12_quad import chi2_over_df_nodes

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(HERE, "series")
FIG = os.path.join(HERE, "figs_h12")
os.makedirs(FIG, exist_ok=True)

GREY = "#9a9a9a"
RED = "#c24f34"
DARK = "#4d4d4d"
BLUE = "#4c6f9c"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.edgecolor": "#cccccc", "axes.labelcolor": DARK,
    "axes.labelsize": 12, "xtick.color": DARK, "ytick.color": DARK,
    "axes.grid": True, "grid.color": "#e8e8e8", "grid.linewidth": 0.9,
    "axes.axisbelow": True, "figure.dpi": 110,
})


def frame(ax, title, ylab=None, xlab=None, grid_axis="y"):
    ax.set_title(title, color="#666666", fontsize=12, pad=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.grid(grid_axis in ("x", "both"))
    ax.yaxis.grid(grid_axis in ("y", "both"))
    if ylab:
        ax.set_ylabel(ylab)
    if xlab:
        ax.set_xlabel(xlab)


def save(fig, name):
    fig.tight_layout()
    p = os.path.join(FIG, name)
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  wrote %s" % p)


def load(n):
    return json.load(open(os.path.join(S, n)))


# --------------------------------------------------------------- 1 loadings
def fig_loadings():
    n = np.array([36.] * 6)
    sig = np.array([12., 8., 10., 15., 20., 30.])
    h = factor_loadings(n, sig)
    lab = ["8%", "10%", "15%", "20%", "30%"]
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    x = np.arange(len(h))
    ax.bar(x, h, color=GREY, width=0.6, label="true loading")
    ax.axhline(np.sqrt(0.5), color=RED, lw=2.2, ls="--",
               label="value the homoscedastic textbook assumes (0.707)")
    for xi, hi in zip(x, h):
        ax.text(xi, hi + 0.018, "%.3f" % hi, ha="center", fontsize=10.5,
                color=DARK, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(lab)
    ax.set_ylim(0, 1.0)
    frame(ax, "one benchmark, five candidates: the comparisons are not equally "
               "correlated, and the textbook assumes they are",
          ylab="loading $h_i$", xlab="candidate annualised volatility")
    ax.legend(frameon=False, loc="upper right", fontsize=10.5)
    save(fig, "fig1_loadings.png")


# ------------------------------------------------------------------- 2 step
def fig_step():
    n = np.array([10., 4000., 4000., 4000.])
    sig = np.array([1., 1e-3, 1e-3, 1e-3])
    c, i = 2.0, 1
    v = sig[0] ** 2 / n[0]
    q, w = chi2_over_df_nodes(n[i] - 1.0, 8, 12)
    s2 = q * sig[i] ** 2
    scale = np.sqrt(n[i]) / sig[i]

    def g(u):
        u = np.atleast_1d(u)
        d = np.sqrt(s2[None, :] / n[i] + v)
        return ndtr((u[:, None] + c * d) * scale) @ w

    centre = -c * np.sqrt(sig[i] ** 2 / n[i] + v)
    width = sig[i] / np.sqrt(n[i])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))

    u = np.linspace(-2.0, 1.0, 4000)
    axes[0].plot(u, g(u), color=DARK, lw=2)
    for e in np.linspace(-8, 8, 15) * (sig[0] / np.sqrt(n[0])):
        if -2.0 <= e <= 1.0:
            axes[0].axvline(e, color=GREY, lw=0.8, alpha=0.65)
    axes[0].axvline(centre, color=RED, lw=1.6)
    panel_u = (2 * 8.0 / 14) * (sig[0] / np.sqrt(n[0]))
    frame(axes[0], "as the quadrature sees it: the transition (red) is about\n"
                   "%s times narrower than a panel (grey edges)"
                   % format(int(round(panel_u / width, -2)), ","),
          ylab="$g_i(u)$", xlab="$u$")

    uz = np.linspace(centre - 8 * width, centre + 8 * width, 1200)
    axes[1].plot(uz, g(uz), color=DARK, lw=2)
    for k in (-6, 0, 6):
        axes[1].axvline(centre + k * width, color=RED, lw=1.6,
                        ls="--" if k else "-")
    axes[1].ticklabel_format(useOffset=False, style="plain", axis="x")
    axes[1].set_xticks([centre - 6 * width, centre, centre + 6 * width])
    axes[1].set_xticklabels(["$-6\\delta_i$", "step", "$+6\\delta_i$"])
    frame(axes[1], "the same function over a window %s times narrower:\n"
                   "the planted edges resolve it"
                   % format(int(round(3.0 / (16 * width), -2)), ","),
          xlab="$u$ near the step, width $\\delta_i=\\sigma_i/\\sqrt{n_i}$")
    save(fig, "fig2_step.png")


# ---------------------------------------------------------------- 3 scaling
def fig_scaling():
    d = load("h12_hetero_exact.json")["f5_perfect_dependence"]
    r = [x for x in d if x["nu0"] == 9 and x["K"] == 3 and x["c"] == 2.0][0]
    st = np.array([1e-3, 1e-4, 1e-5])
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    ax.loglog(st, r["errs_by_sigma_t"], "o-", color=RED, lw=2.2, ms=8,
              label="edges planted on the step")
    ax.loglog(st, r["unadapted_errs"], "s-", color=GREY, lw=2.2, ms=8,
              label="uniform panels")
    ax.loglog(st, r["errs_by_sigma_t"][0] * (st / st[0]), ":", color=DARK,
              lw=1.4, label="exact linear scaling")
    frame(ax, "the falsifier that separates model difference from integration "
              "error:\nonly one of these scales with the candidate volatility",
          ylab="|CDF $-$ Student $t_{\\nu_0}$|",
          xlab="candidate volatility $\\sigma_t$")
    ax.legend(frameon=False, fontsize=10.5)
    save(fig, "fig3_scaling.png")


# ------------------------------------------------------------------ 4 rules
def fig_rules():
    d = load("h12_hetero_compare.json")
    rows = d["comparison"]
    names = [("T3_smm_welch_df", "T3 / studentised max modulus", GREY),
             ("sidak_t_welch_df", "independent scales", "#7a7a7a"),
             ("onefactor_pooled_df", "Dunnett constant carried across", RED),
             ("onefactor_mean_welch_df", "one-factor t, mean Welch df", BLUE)]
    labs = [r["design"].replace(", ", ",\n") for r in rows]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(11.6, 5.6))
    for j, (k, nm, col) in enumerate(names):
        ax.plot(x + (j - 1.5) * 0.13, [r[k]["true_fwer"] for r in rows], "o",
                ms=9, color=col, label=nm)
    ax.axhline(0.05, color=DARK, lw=1.6, ls="--", label="nominal 0.05")
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=8.6)
    frame(ax, "the true familywise error rate each shipped rule delivers, "
              "computed exactly rather than simulated",
          ylab="true FWER")
    ax.legend(frameon=False, fontsize=10, loc="upper left")
    save(fig, "fig4_rules.png")


# -------------------------------------------------------------- 5 crossover
def fig_crossover():
    d = load("h12_crossover.json")["d1_crossover"]
    df = [r["df0"] for r in d]
    e = [r["uncal_rel_err_pct"] for r in d]
    m = [r["m3_rel_err_pct"] for r in d]
    fig, ax = plt.subplots(figsize=(9.8, 5.2))
    ax.plot(df, e, "o-", color=RED, lw=2.4, ms=8, label="exact law, plug-in ratios")
    ax.plot(df, m, "s-", color=GREY, lw=2.4, ms=8, label="Tamhane and Xi Method 3")
    ax.axhline(0, color=DARK, lw=1.4)
    ax.axhspan(-2, 2, color="#eeeeee", zorder=0)
    ax.set_xscale("log")
    ax.set_xticks(df)
    ax.set_xticklabels([str(v) for v in df])
    ax.annotate("clinical designs", xy=(3, 13.8), xytext=(3.3, 11.2),
                color=DARK, fontsize=10.5)
    ax.annotate("a desk with 36 monthly returns", xy=(35, 1.9),
                xytext=(9, 6.4), color=DARK, fontsize=10.5,
                arrowprops=dict(arrowstyle="->", color=DARK, lw=1.2))
    frame(ax, "the whole disadvantage lives below five degrees of freedom, "
              "and is gone by seven",
          ylab="relative error in FWER, per cent",
          xlab="degrees of freedom on the benchmark")
    ax.legend(frameon=False, fontsize=10.5)
    save(fig, "fig5_crossover.png")


# ------------------------------------------------------------------ 6 decomp
def fig_decomp():
    d = load("h12_crossover.json")["d2_decompose"]
    labs = [r["case"].split(" ", 1)[0] + "\n" + r["case"].split(" ", 1)[1]
            for r in d]
    ap = [r["m3_approximation_bias"] for r in d]
    es = [r["m3_estimation_bias"] for r in d]
    ex = [r["exact_plugin_bias"] for r in d]
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(10.6, 5.4))
    ax.bar(x - 0.22, ap, 0.2, color=RED, label="Method 3 approximation bias")
    ax.bar(x, es, 0.2, color=GREY, label="Method 3 estimation bias")
    ax.bar(x + 0.22, ex, 0.2, color=BLUE, label="exact law, estimation bias only")
    ax.axhline(0, color=DARK, lw=1.4)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=8.8)
    frame(ax, "Method 3 carries a large positive approximation bias and a "
              "negative estimation bias; the exact law carries only the latter",
          ylab="bias in the critical value")
    ax.legend(frameon=False, fontsize=10)
    save(fig, "fig6_decomposition.png")


# ------------------------------------------------------------- 7 calibration
def fig_calibration():
    d = load("h12_calibrated.json")["summary"]
    arms = [("uncalibrated", "plug-in, uncalibrated", GREY),
            ("calibrated", "plug-in, calibrated", RED),
            ("method3", "Method 3", BLUE)]
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    x = np.arange(2)
    for j, (k, nm, col) in enumerate(arms):
        vals = [d[K][k]["mean_abs_rel_err_pct"] for K in ("2", "3")]
        ax.bar(x + (j - 1) * 0.26, vals, 0.24, color=col, label=nm)
        for xi, v in zip(x + (j - 1) * 0.26, vals):
            ax.text(xi, v + 0.12, "%.1f" % v, ha="center", fontsize=10.5,
                    color=DARK, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(["2 candidates", "3 candidates"])
    frame(ax, "calibration roughly halves the plug-in error and overtakes "
              "Method 3 on average error",
          ylab="mean absolute error in FWER, per cent")
    ax.legend(frameon=False, fontsize=10.5)
    save(fig, "fig7_calibration.png")


if __name__ == "__main__":
    fig_loadings()
    fig_step()
    fig_scaling()
    fig_rules()
    fig_crossover()
    fig_decomp()
    fig_calibration()
    print("done")
