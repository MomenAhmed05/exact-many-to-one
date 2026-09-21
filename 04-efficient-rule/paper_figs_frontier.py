"""The seven new figures for the paper, drawn from the 04-efficient-rule and 03-certified-supremum
deposits.

THE FILE NAMES ARE NOT THE PAPER'S FIGURE NUMBERS.  These were written before
the document order was settled, and renaming them would break the deposits that
name them.  The mapping is fixed and `paper_verify_frontier.py` pins it by
hashing each embedded image against the file it should have come from:

    paper Figure 14 = fig20_bnbgap.png       paper Figure 18 = fig16_designmap.png
    paper Figure 15 = fig14_ledger.png       paper Figure 19 = fig18_serial.png
    paper Figure 16 = fig17_levels.png       paper Figure 20 = fig19_factor.png
    paper Figure 17 = fig15_frontier.png


Style, palette, frame() and save() follow 02-bounded-rule/paper_figs_w26.py, which follows
h12_paper_figs.py, so these sit beside Figures 1 to 13 without a visual seam.
Every series is read from a deposit; no figure carries a number that is not in
one.

  fig14_ledger      the four objects on one design set: what each costs in the
                    constant and in power, at equal variances
  fig15_frontier    the calibrated rule against the best rule of its class, so
                    the residual loss can be read as an information cost
  fig16_designmap   the planning surface over both groups' degrees of freedom,
                    with the symmetry that identifies which group binds
  fig17_levels      the internal level each certificate supports, by rule and
                    by candidate count
  fig18_serial      serial correlation, naive against whitened with a declared
                    coefficient
  fig19_factor      factor exposure by loading pattern, raw against residualised
  fig20_bnbgap      how tight the branch-and-bound certificate is, cell by cell

Outputs: 04-efficient-rule/figs/*.png
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEP = os.path.join(HERE, "deposits")
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)

GREY = "#9a9a9a"
RED = "#c24f34"
DARK = "#4d4d4d"
BLUE = "#4c6f9c"
GREEN = "#5c8a63"
AMBER = "#c89b3c"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.edgecolor": "#cccccc", "axes.labelcolor": DARK,
    "axes.labelsize": 12, "xtick.color": DARK, "ytick.color": DARK,
    "axes.grid": True, "grid.color": "#e8e8e8", "grid.linewidth": 0.9,
    "axes.axisbelow": True, "figure.dpi": 110,
})


def load(name):
    with open(os.path.join(DEP, name), encoding="utf-8") as fh:
        return json.load(fh)


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


# ------------------------------------------------------------- 14 the ledger
def fig14():
    d = load("wa_ledger.json")
    labs = ["oracle, knows the ratios", "plug-in at the full level",
            "region rule (Section 9)", "calibrated rule (Section 10)"]
    keys = ["oracle_conditional", "plugin_measured",
            "region_sup_conditional", "uncond_certified"]
    cols = [DARK, AMBER, RED, BLUE]
    des = [x for x in d["designs"] if x["design"] in ("df3", "df11", "df35")]
    x = np.arange(len(des))
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.9))

    hat = {1: "///"}          # the plug-in is not valid; say so on the bar
    handles = []
    for j, (k, lab, col) in enumerate(zip(keys, labs, cols)):
        r = [next(q for q in s["rows"] if abs(q["lam"] - 1.0) < 1e-9)[k]
             for s in des]
        h = axes[0].bar(x + (j - 1.5) * 0.2,
                        [e["ratio_to_oracle_mean"] for e in r], 0.19,
                        color=col, hatch=hat.get(j), edgecolor="white",
                        linewidth=0.6, label=lab)
        handles.append(h)
        axes[1].bar(x + (j - 1.5) * 0.2, [e["loss_pp_at_80"] for e in r],
                    0.19, color=col, hatch=hat.get(j), edgecolor="white",
                    linewidth=0.6)
    for ax, ttl, yl in ((axes[0], "what each rule asks of the data", None),
                        (axes[1], "what each rule gives up in power", None)):
        ax.set_xticks(x)
        ax.set_xticklabels(["3 df", "11 df", "35 df"])
    axes[0].axhline(1.0, color=DARK, lw=1.2)
    frame(axes[0], "threshold as a multiple of the oracle's",
          ylab="c / c oracle", xlab="degrees of freedom on the benchmark")
    frame(axes[1], "power given up at the effect where the oracle reaches 80%",
          ylab="percentage points lost",
          xlab="degrees of freedom on the benchmark")
    axes[1].text(0.02, 0.96, "the oracle is the zero line: it knows the ratios",
                 transform=axes[1].transAxes, fontsize=9.5, color=DARK,
                 va="top")
    axes[0].set_ylim(top=3.05)
    axes[0].text(0.02, 0.97, "hatched: exceeds the ceiling somewhere in the declared domain",
                 transform=axes[0].transAxes, fontsize=9.5, color=DARK,
                 va="top")
    fig.legend(handles=[h[0] for h in handles], labels=labs, frameon=False,
               fontsize=10, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.06))
    save(fig, "fig14_ledger.png")


# ---------------------------------------------------------- 15 the frontier
def fig15():
    d = load("wc_frontier.json")
    df = [int(x["design"][2:]) for x in d["designs"]]
    plug = [float(np.mean(x["plugin_loss_pp"])) for x in d["designs"]]
    opt = [float(np.mean(x["optimal_loss_pp"])) for x in d["designs"]]
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    ax.plot(df, plug, "o-", color=BLUE, lw=2.3, ms=8,
            label="the calibrated rule")
    ax.plot(df, opt, "s--", color=RED, lw=2.3, ms=8,
            label="the best rule of its class")
    # shade only where the optimum is actually below the calibrated rule; the
    # solver is local, so at eleven degrees of freedom it comes out above and
    # there is nothing there to recover
    gap = np.array(plug) - np.array(opt)
    ax.fill_between(df, opt, plug, where=gap > 0, color=RED, alpha=0.10,
                    interpolate=True)
    ax.set_xscale("log")
    ax.set_xticks(df)
    ax.set_xticklabels([str(v) for v in df])
    for xx, a, b in zip(df, plug, opt):
        if a - b > 0.3:
            ax.annotate("%.1f pp" % (a - b), xy=(xx, a), xytext=(xx, a + 1.6),
                        fontsize=9.5, color=DARK, ha="center")
    j = int(np.argmin(gap))
    ax.annotate("the search comes out %.2f pp above the rule here, which is\n"
                "its own resolution rather than anything left to gain"
                % abs(gap[j]),
                xy=(df[j], opt[j] + 0.3),
                xytext=(df[j] * 1.15, opt[j] + 9.0),
                fontsize=9.8, color=DARK, ha="center",
                arrowprops=dict(arrowstyle="->", color=DARK, lw=1.1))
    ax.set_ylim(0.0, max(plug) + 6.0)
    frame(ax, "the gap between the rule and the best rule of its class closes "
              "by eleven degrees of freedom",
          ylab="power lost against the oracle, percentage points",
          xlab="degrees of freedom on the benchmark")
    ax.legend(frameon=False, fontsize=10.5)
    save(fig, "fig15_frontier.png")


# --------------------------------------------------------- 16 the design map
def fig16():
    d = load("wd_designmap.json")
    nu0 = sorted({c["nu0"] for c in d["cells"]})
    nu1 = sorted({c["nu1"] for c in d["cells"]})
    M = np.full((len(nu1), len(nu0)), np.nan)
    for c in d["cells"]:
        M[nu1.index(c["nu1"]), nu0.index(c["nu0"])] = c["max_loss_pp"]
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    im = ax.imshow(M, origin="lower", cmap="YlOrRd", aspect="auto",
                   vmin=0, vmax=float(np.nanmax(M)))
    ax.set_xticks(range(len(nu0)))
    ax.set_xticklabels([str(v) for v in nu0])
    ax.set_yticks(range(len(nu1)))
    ax.set_yticklabels([str(v) for v in nu1])
    for i in range(len(nu1)):
        for j in range(len(nu0)):
            if np.isfinite(M[i, j]):
                ax.text(j, i, "%.0f" % M[i, j], ha="center", va="center",
                        fontsize=9.5,
                        color="white" if M[i, j] > 0.6 * np.nanmax(M)
                        else DARK)
    ax.grid(False)
    a = (nu0.index(3), nu1.index(35))
    b = (nu0.index(35), nu1.index(3))
    for (jx, iy) in (a, b):
        ax.add_patch(plt.Rectangle((jx - 0.5, iy - 0.5), 1, 1, fill=False,
                                   edgecolor=BLUE, lw=2.4))
    ax.annotate("", xy=(a[0] + 0.45, a[1]), xytext=(b[0] - 0.45, b[1]),
                arrowprops=dict(arrowstyle="<->", color=BLUE, lw=1.4,
                                connectionstyle="arc3,rad=-0.25"))
    ax.text(3.3, 2.62, "3 against 35 and 35 against 3 cost the same",
            fontsize=9.8, color=BLUE, ha="center")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("power lost, percentage points", color=DARK)
    frame(ax, "above seven degrees of freedom the loss follows the smaller "
              "group, not the benchmark",
          ylab="degrees of freedom on the candidate",
          xlab="degrees of freedom on the benchmark", grid_axis="none")
    save(fig, "fig16_designmap.png")


# ------------------------------------------------------ 17 certified levels
def fig17():
    ak = load("certify_anyk.json")
    jt = load("certify_joint.json")
    gen = load("wb_general.json")
    rc = load("rule_compare.json")
    names = ["K2_df35", "K3_df35", "K5_df35"]
    x = np.arange(len(names))

    def pick(rows, key, cond=lambda r: True):
        return [next((r[key] for r in rows
                      if r["design"] == nm and cond(r)), np.nan)
                for nm in names]

    meas = pick(gen["designs"], "alpha0")
    jo = pick(jt["designs"], "alpha0_certified_joint")
    pc = pick(ak["designs"], "alpha0", lambda r: not r["two_sided"])
    two = pick(ak["designs"], "alpha0", lambda r: r["two_sided"])

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0),
                             gridspec_kw={"width_ratios": [1.35, 1.0]})

    ax = axes[0]
    for j, (v, lab, col) in enumerate((
            (meas, "joint rule, measured", GREY),
            (jo, "joint rule, certified", BLUE),
            (pc, "per-candidate, certified", GREEN),
            (two, "per-candidate, certified, two-sided", AMBER))):
        ax.bar(x + (j - 1.5) * 0.2, v, 0.19, color=col, label=lab)
    ax.axhline(0.025, color=DARK, lw=1.3, ls="--")
    ax.text(-0.44, 0.0256, "ceiling 0.025", fontsize=9.5, color=DARK,
            ha="left")
    ax.text(x[-1] - 0.20, 0.0072, "the joint certificate" + chr(10) +
            "was not attempted at five", fontsize=8.8, color=DARK,
            ha="center", style="italic")
    ax.set_xticks(x)
    ax.set_xticklabels(["2 candidates", "3 candidates", "5 candidates"])
    ax.set_ylim(0, 0.0295)
    frame(ax, "internal level: comparable within a colour, not between them",
          ylab="internal level")

    # the quantity that does compare across rule types
    ax = axes[1]
    des = [d for d in rc["designs"] if d["design"] in ("K2_df35", "K3_df35")]
    xx = np.arange(len(des))
    jc = [d["points"][0]["c_joint_certified"] for d in des]
    jm = [d["points"][0]["c_joint_measured"] for d in des]
    pcc = [max(d["points"][0]["c_percand_certified"]) for d in des]
    for j, (v, lab, col) in enumerate(((jm, "joint rule, measured", GREY),
                                       (jc, "joint rule, certified", BLUE),
                                       (pcc, "per-candidate, certified",
                                        GREEN))):
        b = ax.bar(xx + (j - 1) * 0.24, v, 0.23, color=col, label=lab)
        for r, val in zip(b, v):
            ax.text(r.get_x() + r.get_width() / 2, val + 0.012,
                    "%.3f" % val, ha="center", fontsize=9, color=DARK)
    ax.set_xticks(xx)
    ax.set_xticklabels(["2 candidates", "3 candidates"])
    ax.set_ylim(2.0, max(jc) + 0.22)
    frame(ax, "threshold at equal ratios: lower is the better rule",
          ylab="critical constant")
    ax.annotate("the per-candidate certificate\nasks for less",
                xy=(xx[0] + 0.24, pcc[0] + 0.02),
                xytext=(xx[0] + 0.30, pcc[0] + 0.30),
                fontsize=9.5, color=GREEN, ha="center",
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.1))

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=9.8, ncol=4,
               loc="lower center", bbox_to_anchor=(0.5, -0.04))
    save(fig, "fig17_levels.png")


# ------------------------------------------------------- 18 serial correlation
def fig18():
    d = load("serial.json")
    phi = [r["phi"] for r in d["rows"]]
    naive = [r["naive"]["fwer"] for r in d["rows"]]
    wht = [r["whitened_correct"]["fwer"] for r in d["rows"]]
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    ax.plot(phi, naive, "o-", color=RED, lw=2.3, ms=8,
            label="the rule applied to the raw series")
    ax.plot(phi, wht, "s-", color=BLUE, lw=2.3, ms=8,
            label="the same rule after whitening at a declared coefficient")
    ax.axhline(0.025, color=DARK, lw=1.4, ls="--")
    ax.text(0.405, 0.0265, "ceiling 0.025", fontsize=9.5, color=DARK,
            ha="right")
    ax.annotate("%.1f times the ceiling by 0.1" % (naive[2] / 0.025),
                xy=(0.1, naive[2]), xytext=(0.118, naive[2] + 0.030),
                fontsize=10, color=DARK, ha="left",
                arrowprops=dict(arrowstyle="->", color=DARK, lw=1.1))
    ax.set_ylim(0.0, max(naive) * 1.12)
    frame(ax, "serial correlation stops being a boundary once the coefficient "
              "is declared",
          ylab="familywise error rate", xlab="lag-one correlation")
    ax.legend(frameon=False, fontsize=10.5, loc="upper left")
    save(fig, "fig18_serial.png")


# ------------------------------------------------------------- 19 the factor
def fig19():
    d = load("factor.json")
    pats, arms = [], ["raw returns", "declared loadings",
                      "estimated loadings, df corrected",
                      "estimated loadings, df NOT corrected"]
    for r in d["rows"]:
        if r["exposure"] not in pats:
            pats.append(r["exposure"])
    val = {(r["exposure"], r["arm"]): r["fwer"] for r in d["rows"]}
    x = np.arange(len(pats))
    cols = [RED, BLUE, GREEN, AMBER]
    short = ["raw returns", "declared loadings",
             "estimated, the degree of freedom returned",
             "estimated, not returned"]
    fig, ax = plt.subplots(figsize=(11.4, 5.2))
    for j, (a, col, lab) in enumerate(zip(arms, cols, short)):
        ax.bar(x + (j - 1.5) * 0.2, [val[(p, a)] for p in pats], 0.19,
               color=col, label=lab)
    ax.axhline(0.025, color=DARK, lw=1.4, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([p.replace(" (", "\n(") for p in pats], fontsize=9.5)
    ax.text(-0.44, 0.0262, "ceiling 0.025", fontsize=9.5, color=DARK,
            ha="left")
    frame(ax, "a factor moves the error rate in both directions, and declared "
              "loadings remove it", ylab="familywise error rate")
    ax.legend(frameon=False, fontsize=9.5, loc="upper left", ncol=2)
    save(fig, "fig19_factor.png")


# ------------------------------------------------- 20 the certificate's gap
def fig20():
    """How tight the branch-and-bound certificate is, cell by cell."""
    with open(os.path.join(HERE, "..", "03-certified-supremum", "deposits", "bnb.json"),
              encoding="utf-8") as fh:
        d = json.load(fh)
    cells = d["cells"]
    name = {"A": "family A", "B": "family B"}
    labs, gaps, cert = [], [], []
    for c in cells:
        k = c["cell"]
        if k.startswith(("A_", "B_")):
            lab = "%s, %s df" % (name[k[0]], k.split("df")[1])
        else:
            lab = {"K1_df3": "1 candidate, 3 df",
                   "K3_df6": "3 candidates, 6 df",
                   "K5_df11": "5 candidates, 11 df",
                   "K2_extreme": "extreme spread",
                   "K2_quiet": "quiet benchmark"}[k]
        labs.append(lab)
        gaps.append(c["gap"])
        cert.append(bool(c["certified_at_tol"]))

    order = np.argsort(gaps)
    labs = [labs[i] for i in order]
    gaps = [gaps[i] for i in order]
    cert = [cert[i] for i in order]
    y = np.arange(len(labs))

    fig, ax = plt.subplots(figsize=(10.6, 6.4))
    ax.barh(y, gaps, 0.62,
            color=[GREEN if c else AMBER for c in cert],
            edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=10)
    ax.set_xscale("log")
    ax.set_xlim(3e-7, 6e-2)
    ax.set_ylim(len(labs) + 0.3, -1.4)
    ax.axvline(d["tol"], color=DARK, lw=1.4, ls="--")
    ax.axvline(2.0 ** -20, color=GREY, lw=1.2, ls=":")
    ax.text(d["tol"] * 1.3, -1.0, "declared tolerance 1e-4", fontsize=9.5,
            color=DARK, ha="left", va="center")
    ax.text(2.0 ** -20 * 1.3, len(labs) - 0.1, "solver floor $2^{-20}$",
            fontsize=9.5, color=GREY, ha="left", va="center")
    for i, g in enumerate(gaps):
        ax.text(g * 1.35, i, "%.1e" % g, fontsize=9, color=DARK,
                va="center", ha="left")
    frame(ax, "the certified bound sits this far above a witness below it",
          xlab="gap in the critical constant", grid_axis="x")
    fig.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, color=GREEN,
                      label="certified at the declared tolerance"),
        plt.Rectangle((0, 0), 1, 1, color=AMBER,
                      label="split budget exhausted, gap reported")],
        frameon=False, fontsize=10, loc="lower center", ncol=2,
        bbox_to_anchor=(0.55, -0.045))
    save(fig, "fig20_bnbgap.png")


if __name__ == "__main__":
    for f in (fig14, fig15, fig16, fig17, fig18, fig19, fig20):
        f()
    print("done")
