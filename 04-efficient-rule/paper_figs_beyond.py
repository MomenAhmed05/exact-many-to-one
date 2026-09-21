"""The three new figures for the W1/W3 paper update.

Style, palette, frame() and save() are imported from `paper_figs_frontier`, so
these sit beside Figures 1 to 20 without a visual seam.  Every series is read
from a deposit; no figure carries a number that is not in one.

The file names are NOT the paper's figure numbers here either: inserting a
figure into Section 10.8 pushed the three that followed it down one, so the
floor figure is Figure 18 and the old 18, 19 and 20 became 19, 20 and 21.
`paper_verify_frontier.py` pins the whole mapping by content hash.

    paper Figure 18 = fig21_floor.png     the calibrated rule against the bound
                                          over every rule that sees only the
                                          estimated ratio
    paper Figure 22 = fig22_tiers.png     the three validity tiers for a serial
                                          nuisance, and what pooling recovers
    paper Figure 23 = fig23_sign.png      the aligned factor: which side of the
                                          ceiling the loading signs put you on

Outputs: 04-efficient-rule/figs/*.png
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from paper_figs_frontier import (AMBER, BLUE, DARK, GREEN, GREY,  # noqa: E402
                                 RED, frame, load, save)

HERE = os.path.dirname(os.path.abspath(__file__))

DF = [("df3", 3), ("df5", 5), ("df7", 7), ("df11", 11), ("df19", 19),
      ("df35", 35)]


# ------------------------------------------------------- 21 the floor
def fig21():
    """Left: three losses by benchmark length.  Right: the gap that remains."""
    fr = load("wc_frontier.json")
    rg = load("w1_region.json")
    fd = {d["design"]: d for d in fr["designs"]}
    rd = {d["design"]: d for d in rg["designs"]}

    x = [n for _k, n in DF]
    dep = [float(np.mean(fd[k]["plugin_loss_pp"])) for k, _n in DF]
    opt = [float(np.mean(fd[k]["optimal_loss_pp"])) for k, _n in DF]
    bound = [rd[k]["loss_pp_monotone"] for k, _n in DF]
    gap = [d - b for d, b in zip(dep, bound)]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.9))

    ax = axes[0]
    ax.fill_between(x, bound, dep, color=RED, alpha=0.13, zorder=1)
    ax.plot(x, dep, "o-", color=RED, lw=2.0, ms=6, label="the calibrated rule",
            zorder=3)
    ax.plot(x, opt, "s--", color=AMBER, lw=1.6, ms=5,
            label="Section 10.7's optimum", zorder=3)
    ax.plot(x, bound, "^-", color=BLUE, lw=2.0, ms=6,
            label="bound over every rule in the class", zorder=3)
    ax.set_xscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in x])
    ax.set_ylim(0, max(dep) * 1.08)
    frame(ax, "What the rule gives up against an oracle that knows the ratio",
          "power given up at 80% (points)", "benchmark degrees of freedom")
    ax.legend(frameon=False, fontsize=10, loc="upper right")

    ax = axes[1]
    cols = [RED if g > 0.5 else BLUE for g in gap]
    ax.bar(range(len(x)), gap, 0.55, color=cols, edgecolor="white",
           linewidth=0.8)
    for i, g in enumerate(gap):
        ax.annotate("%.2f" % g, (i, g), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=10, color=DARK)
    ax.set_xticks(range(len(x)))
    ax.set_xticklabels([str(n) for n in x])
    ax.axhline(0.5, color=GREY, lw=1.0, ls=":")
    frame(ax, "What is left on the table", "points above the bound",
          "benchmark degrees of freedom")
    save(fig, "fig21_floor.png")


# ------------------------------------------------------- 22 the tiers
def fig22():
    """Left: the tier ladder.  Right: pooling against the common mean."""
    ti = load("w3_tiers.json")
    rp = load("w3_serial_repair.json")
    abar = ti["abar"]

    labs = ["A, declared", "B, 36", "B, 120", "B, 360", "C, same period"]
    ser = [ti["serial_tierA_worst_fwer"],
           ti["serial_tierB_cal36_worst_fwer"],
           ti["serial_tierB_cal120_worst_fwer"],
           ti["serial_tierB_cal360_worst_fwer"],
           ti["serial_tierC_worst_fwer"]]
    fac = [ti["factor_tierA_worst_fwer"],
           ti["factor_tierB_cal36_worst_fwer"],
           ti["factor_tierB_cal120_worst_fwer"],
           ti["factor_tierB_cal360_worst_fwer"],
           ti["factor_tierC_worst_fwer"]]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.9))

    # plotted as a MULTIPLE of the ceiling: a log axis under bars that start
    # at the axis floor reads as though the floor were zero, which it is not
    ax = axes[0]
    i = np.arange(len(labs))
    rs = [v / abar for v in ser]
    rf = [v / abar for v in fac]
    ax.bar(i - 0.19, rs, 0.36, color=RED, edgecolor="white", linewidth=0.8,
           label="serial coefficient")
    ax.bar(i + 0.19, rf, 0.36, color=GREEN, edgecolor="white", linewidth=0.8,
           label="factor loadings")
    for k, (a_, b_) in enumerate(zip(rs, rf)):
        ax.annotate("%.1fx" % a_, (k - 0.19, a_), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=9, color=DARK)
        ax.annotate("%.1fx" % b_, (k + 0.19, b_), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=9, color=DARK)
    ax.axhline(1.0, color=DARK, lw=1.4, ls="--")
    ax.set_xticks(i)
    ax.set_xticklabels(labs, fontsize=10)
    ax.set_ylim(0, max(rs) * 1.16)
    frame(ax, "Worst realised error rate by validity tier",
          "multiple of the familywise ceiling", "tier")
    ax.legend(frameon=False, fontsize=10, loc="upper left")

    ax = axes[1]
    wc = max(rp["rows"], key=lambda r: r["fwer"])
    trip = [r for r in rp["rows"] if r["design"] == wc["design"]
            and r["phi"] == wc["phi"] and r["tier"] == wc["tier"]]
    mus = sorted(set(r["mu"] for r in trip))
    sep = [next(r["fwer"] for r in trip if r["mu"] == m and not r["pooled"])
           for m in mus]
    poo = [next(r["fwer"] for r in trip if r["mu"] == m and r["pooled"])
           for m in mus]
    ax.plot(mus, sep, "o-", color=RED, lw=2.0, ms=7,
            label="one coefficient per group")
    ax.plot(mus, poo, "s-", color=BLUE, lw=2.0, ms=7,
            label="one coefficient, pooled")
    ax.axhline(abar, color=DARK, lw=1.4, ls="--")
    ax.set_ylim(0, max(sep) * 1.12)
    frame(ax, "Why pooling matters, and how far it goes",
          "familywise error rate", "common mean of both groups")
    ax.legend(frameon=False, fontsize=10, loc="upper left")
    save(fig, "fig22_tiers.png")


# ------------------------------------------------------- 23 the loading sign
def fig23():
    """Left: realised rate by loading pair.  Right: the exact law's constant."""
    ce = load("w2_cert.json")
    fa = load("w2_factor.json")
    abar = ce["abar"]

    rows = ce["price_rows"]
    lams = sorted(set(r["lam"] for r in rows))
    pairs = []
    for r in rows:
        k = (r["b0"], r["b1"])
        if k not in pairs:
            pairs.append(k)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.9))

    ax = axes[0]
    w = 0.26
    for j, lam in enumerate(lams):
        y = [next(r["fwer"] for r in rows
                  if r["lam"] == lam and (r["b0"], r["b1"]) == p)
             for p in pairs]
        col = [RED if p[0] * p[1] < 0 else BLUE for p in pairs]
        ax.bar(np.arange(len(pairs)) + (j - 1) * w, y, w * 0.92,
               color=col, alpha=0.55 + 0.22 * j, edgecolor="white",
               linewidth=0.7,
               label="variance ratio %.2f" % lam)
    ax.axhline(abar, color=DARK, lw=1.4, ls="--")
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels(["%.1f, %.1f" % p for p in pairs], fontsize=10)
    frame(ax, "Blue: loadings share a sign.  Red: they do not.",
          "familywise error rate", "benchmark and candidate loading")
    ax.legend(frameon=False, fontsize=9, loc="upper left")

    ax = axes[1]
    labs = [t["design"].replace("_", " ") for t in fa["critical_table"]]
    vals = [t["c_alpha"] for t in fa["critical_table"]]
    cols = [GREY if "independent" in t["design"] else BLUE
            for t in fa["critical_table"]]
    ax.barh(range(len(labs)), vals, 0.55, color=cols, edgecolor="white",
            linewidth=0.8)
    for i, v in enumerate(vals):
        ax.annotate("%.3f" % v, (v, i), textcoords="offset points",
                    xytext=(5, -4), fontsize=10, color=DARK)
    ax.set_yticks(range(len(labs)))
    ax.set_yticklabels(labs, fontsize=10)
    ax.set_xlim(0, max(vals) * 1.18)
    frame(ax, "The exact constant the aligned law requires",
          None, "critical constant at the ceiling", grid_axis="x")
    save(fig, "fig23_sign.png")


def main():
    print("figures for the W1/W3 update")
    fig21()
    fig22()
    fig23()


if __name__ == "__main__":
    main()
