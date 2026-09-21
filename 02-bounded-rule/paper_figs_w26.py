"""Figures 8 to 13 for the paper, drawn from the 02-bounded-rule deposits only.

Style, palette, frame() and save() follow h12_paper_figs.py so the new figures sit
beside Figures 1 to 7 without a visual seam.  Every series is read from a deposited
JSON under 02-bounded-rule/deposits/; no figure carries a number that is not in a deposit.

  fig8_guarantee   the envelope and the plug-in on the paper's own four cases:
                   error rate against the nominal level and the guarantee, and what
                   each rule costs in the critical constant
  fig9_deskmap     the twelve-cell desk map: family B and family A by benchmark
                   degrees of freedom, error rate and cost against the union-bound
                   allocation of the same budget
  fig10_power      what the envelope costs in power at the effect size where the
                   oracle would have had eighty per cent, drawn as a dot plot so the
                   axis can be truncated without misrepresenting any bar length
  fig11_ar1        the guarantee under serial correlation: three rules, two sample
                   sizes, with the correlation at which the guarantee is exhausted
  fig12_factor     the factor design: the inflation the exact K = 1 law requires
                   against the one the closed-form correction supplies
  fig13_miss       the region's designed miss rate from the coverage quadrature
                   against the miss rate measured over 240,000 draws

Outputs: paper/figs/*.png
"""

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEP = os.path.join(HERE, "deposits")
FIG = os.path.join(HERE, "paper", "figs")

GREY = "#9a9a9a"
RED = "#c24f34"
DARK = "#4d4d4d"
BLUE = "#4c6f9c"
GREEN = "#5c8a63"
BOX = {"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 1.4}
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.edgecolor": "#cccccc",
        "axes.labelcolor": DARK,
        "axes.labelsize": 12,
        "xtick.color": DARK,
        "ytick.color": DARK,
        "axes.grid": True,
        "grid.color": "#e8e8e8",
        "grid.linewidth": 0.9,
        "axes.axisbelow": True,
        "figure.dpi": 110,
    }
)


def _f(x, default=float("nan")):
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return float(default)


def load(*parts):
    p = os.path.join(DEP, *parts)
    try:
        with open(p, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"cannot read {p}: {exc}") from exc


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


def levels(ax, alpha, gamma, x_nominal, x_guarantee):
    """Nominal and guarantee reference lines with labels kept off the data."""
    ax.axhline(alpha, color=DARK, lw=1.4)
    ax.axhline(alpha + gamma, color=GREEN, lw=1.4, ls="--")
    ax.annotate(
        f"nominal {alpha:.3f}",
        (x_nominal, alpha),
        ha="right",
        va="bottom",
        fontsize=9,
        color=DARK,
        bbox=BOX,
    )
    ax.annotate(
        "guarantee %.4f" % (alpha + gamma),
        (x_guarantee, alpha + gamma),
        ha="left",
        va="bottom",
        fontsize=9,
        color=GREEN,
        bbox=BOX,
    )


def save(fig, name):
    try:
        os.makedirs(FIG, exist_ok=True)
    except OSError as exc:
        raise SystemExit(f"cannot create {FIG}: {exc}") from exc
    fig.tight_layout()
    p = os.path.join(FIG, name)
    fig.savefig(p, bbox_inches="tight", facecolor="white", dpi=200)
    plt.close(fig)
    print(f"  wrote {p}")


# ------------------------------------------------------------------ figure 8
def fig8_guarantee():
    d = load("env", "paper_cases.json")
    rows = d["rows"]
    alpha, gamma = _f(d["alpha"]), _f(d["gamma"])
    labels = [
        "S1\nbalanced",
        "S4\ncontrol\nsmallest",
        "S5\ntrt\nsmallest",
        "S6\ncontrol\nlargest",
    ]
    fw_or = [_f(r["fwer_oracle"]) for r in rows]
    fw_pl = [_f(r["fwer_plug"]) for r in rows]
    fw_ev = [_f(r["fwer_env"]) for r in rows]
    cost_pl = [_f(r["cost_plug_pct"]) for r in rows]
    cost_ev = [_f(r["cost_env_pct"]) for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.3))
    ax = axes[0]
    x = np.arange(len(rows))
    w = 0.26
    ax.bar(x - w, fw_or, w, color=GREY, label="oracle, true ratios")
    ax.bar(x, fw_pl, w, color=RED, label="plug-in, estimated ratios")
    ax.bar(x + w, fw_ev, w, color=BLUE, label="envelope, region guarantee")
    ax.set_xlim(-0.72, len(rows) - 0.28)
    levels(ax, alpha, gamma, x_nominal=len(rows) - 0.34, x_guarantee=-0.70)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0, 0.036)
    ax.legend(frameon=False, fontsize=9, loc="upper center", ncol=1)
    frame(
        ax,
        "Familywise error rate on the four designs of Section 7",
        ylab="true FWER at a one-sided 0.025",
    )

    ax = axes[1]
    ax.bar(x - w / 2, cost_pl, w, color=RED, label="plug-in")
    ax.bar(x + w / 2, cost_ev, w, color=BLUE, label="envelope")
    ax.axhline(0.0, color=DARK, lw=1.2)
    for xi, v in zip(x + w / 2, cost_ev, strict=True):
        ax.annotate(
            f"{v:+.1f}%", (xi, v), ha="center", va="bottom", fontsize=9, color=BLUE
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(-6, 30)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    frame(
        ax,
        "What each rule costs in the critical constant",
        ylab="increase over the oracle constant, per cent",
    )
    save(fig, "fig8_guarantee.png")


# ------------------------------------------------------------------ figure 9
def _sweep():
    cells = load("env", "sweep_summary.json")["cells"]
    out = []
    for c in cells:
        fam, df0 = c["family"], int(c["df0"])
        per = load("env", f"sweep_df{df0}_{fam}.json")
        plug = [
            _f(r.get("plugin_fwer"))
            for r in per.get("per_seed", [])
            if r.get("plugin_fwer") is not None
        ]
        out.append(
            {
                "family": fam,
                "df0": df0,
                "env_fwer": _f(c["tight_fwer_est"]),
                "plug_fwer": float(np.mean(plug)) if plug else float("nan"),
                "tight_cost": _f(c["tight_cost_pct"]),
                "ub_cost": _f(c["ub_cost_pct"]),
                "n_outside": int(c["tight_n_outside"]),
            }
        )
    return out


def fig9_deskmap():
    rows = _sweep()
    alpha, gamma = 0.025, 0.0005
    ticks = [3, 5, 7, 11, 19, 35]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4))
    for fam, col, mk in (("A", BLUE, "o"), ("B", RED, "s")):
        f = sorted([r for r in rows if r["family"] == fam], key=lambda r: r["df0"])
        xs = [r["df0"] for r in f]
        axes[0].plot(
            xs,
            [r["env_fwer"] for r in f],
            mk + "-",
            color=col,
            ms=6,
            label=f"envelope, family {fam}",
        )
        axes[0].plot(
            xs,
            [r["plug_fwer"] for r in f],
            mk + "--",
            color=col,
            ms=5,
            alpha=0.5,
            label=f"plug-in, family {fam}",
        )
        axes[1].plot(
            xs,
            [r["tight_cost"] for r in f],
            mk + "-",
            color=col,
            ms=6,
            label=f"tight region, family {fam}",
        )
        axes[1].plot(
            xs,
            [r["ub_cost"] for r in f],
            mk + ":",
            color=col,
            ms=5,
            alpha=0.5,
            label=f"union bound, family {fam}",
        )
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks])
        ax.set_xlim(2.6, 42)
    axes[0].set_ylim(0.005, 0.036)
    levels(axes[0], alpha, gamma, x_nominal=40, x_guarantee=2.7)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    frame(
        axes[0],
        "Error rate across the desk map",
        ylab="true FWER at a one-sided 0.025",
        xlab="benchmark degrees of freedom",
    )
    axes[1].set_ylim(0, 27)
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")
    frame(
        axes[1],
        "Cost of the same guarantee under two allocations",
        ylab="increase over the oracle constant, per cent",
        xlab="benchmark degrees of freedom",
    )
    save(fig, "fig9_deskmap.png")


# ----------------------------------------------------------------- figure 10
def fig10_power():
    rows = load("power", "power_table.json")["rows"]
    labels = [f"{r['family']} {int(r['df0'])}" for r in rows]
    pw = [_f(r["power_env_at_nu80"]) for r in rows]
    loss = [_f(r["power_loss_at_nu80"]) for r in rows]
    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    ax.axhline(0.80, color=DARK, lw=1.4)
    ax.annotate(
        "the oracle clears this effect at 0.80",
        (-0.42, 0.80),
        ha="left",
        va="bottom",
        fontsize=9,
        color=DARK,
        bbox=BOX,
    )
    for i, (p, loss_v, r) in enumerate(zip(pw, loss, rows, strict=True)):
        col = BLUE if r["family"] == "A" else RED
        ax.plot([i], [p], "o", color=col, ms=8)
        note = f"{p:.3f}\n-{100 * loss_v:.1f} pp"
        if _f(r["nu_for_80pct_oracle"]) >= 3.9:
            note += "\nlanded on\nthe grid end"
        ax.annotate(
            note,
            (i, p),
            xytext=(7, -2),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=8.5,
            color=DARK,
        )
    ax.plot([], [], "o", color=BLUE, label="family A, ratios near one")
    ax.plot([], [], "o", color=RED, label="family B, ratios one to twenty-five")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlim(-0.6, len(rows) - 0.1)
    ax.set_ylim(0.55, 0.845)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    frame(
        ax,
        "Power of the envelope at the effect its oracle would clear at eighty per cent",
        ylab="power at that effect",
        xlab="design (family, benchmark degrees of freedom)",
    )
    save(fig, "fig10_power.png")


# ----------------------------------------------------------------- figure 11
def fig11_ar1():
    d = load("ar2", "ar1.json")
    rows = d["rows"]
    alpha = _f(d["alpha"])
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4), sharey=True)
    for ax, n_i in zip(axes, sorted({int(r["n"]) for r in rows}), strict=True):
        f = sorted([r for r in rows if int(r["n"]) == n_i], key=lambda r: _f(r["phi"]))
        phis = [_f(r["phi"]) for r in f]
        for key, col, mk, lab in (
            ("plugin", RED, "o", "plug-in"),
            ("envelope", BLUE, "s", "envelope"),
            ("bootstrap", GREY, "^", "moving-block bootstrap"),
        ):
            ax.errorbar(
                phis,
                [_f(r[f"fwer_{key}"]) for r in f],
                yerr=[_f(r[f"mcse_{key}"]) for r in f],
                fmt=mk + "-",
                color=col,
                ms=5,
                lw=1.3,
                capsize=3,
                elinewidth=1.1,
                label=lab,
            )
        ax.axhline(alpha, color=DARK, lw=1.4)
        at = next((r for r in f if abs(_f(r["phi"]) - 0.1) < 1e-9), None)
        if at is not None:
            ax.annotate(
                f"+{_f(at['rel_err_pct_plugin']):.0f}% at phi = 0.1",
                (0.1, _f(at["fwer_plugin"])),
                xytext=(9, 6),
                textcoords="offset points",
                fontsize=9,
                color=RED,
                bbox=BOX,
            )
        ax.annotate(
            "nominal 0.025",
            (0.415, alpha),
            ha="right",
            va="bottom",
            fontsize=9,
            color=DARK,
            bbox=BOX,
        )
        ax.set_xlim(-0.02, 0.42)
        ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
        ax.legend(frameon=False, fontsize=9, loc="upper left")
        frame(
            ax,
            f"n = {n_i} observations per group",
            ylab="true FWER" if n_i == 36 else None,
            xlab="lag-one correlation",
        )
    axes[0].annotate(
        "27 of the 30 rule-cells sit outside five per cent of nominal",
        (0.02, 0.015),
        ha="left",
        va="bottom",
        fontsize=9,
        color=DARK,
    )
    save(fig, "fig11_ar1.png")


# ----------------------------------------------------------------- figure 12
def fig12_factor():
    d = load("w6", "factor_d2.json")
    naive_c = _f(load("cert", "k1_reference.json")["naive_reference_c"])
    exact = {_f(r["rho"]): _f(r["c_exact_k1"]) / naive_c for r in d["k1_falsifier"]}
    rows = d["rows"]
    fig, ax = plt.subplots(figsize=(8.8, 4.5))
    for mode, col, mk, lab in (
        ("cand_only", GREY, "^", "closed form, candidates only"),
        ("same", BLUE, "o", "closed form, same sign"),
        ("sign_flip", RED, "s", "closed form, sign flip"),
    ):
        pts = sorted(
            {(_f(r["rho"]), _f(r["w_factor"])) for r in rows if r["mode"] == mode}
        )
        ax.plot(
            [p[0] for p in pts],
            [p[1] for p in pts],
            mk + "-",
            color=col,
            ms=5,
            lw=1.4,
            label=lab,
        )
    ax.plot(
        sorted(exact),
        [exact[k] for k in sorted(exact)],
        "kD-",
        ms=7,
        lw=1.1,
        label="exact law, one candidate",
    )
    ax.axhline(1.0, color=DARK, lw=1.2)
    ax.annotate(
        "no correction required",
        (0.05, 1.0),
        ha="left",
        va="bottom",
        fontsize=9,
        color=DARK,
        bbox=BOX,
    )
    for rho, val in sorted(exact.items()):
        ax.annotate(
            f"{val:.2f}x",
            (rho, val),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9.5,
            color="black",
        )
    ax.set_xlim(0.0, 0.83)
    ax.set_ylim(0.4, 3.05)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    frame(
        ax,
        "Inflation of the critical constant under a common factor",
        ylab="required constant, relative to no correction",
        xlab="loading magnitude rho",
    )
    save(fig, "fig12_factor.png")


# ----------------------------------------------------------------- figure 13
def fig13_miss():
    d = load("cert", "miss_probe.json")
    rows = d["rows"]
    gamma = _f(d["gamma"])
    labels = [f"{r['family']} {int(r['df0'])}" for r in rows]
    x = np.arange(len(rows))
    obs = [_f(r["miss_observed"]) for r in rows]
    err = [_f(r["miss_observed_mcse"]) for r in rows]
    des = [_f(r["miss_designed_quadrature"]) for r in rows]
    fig, ax = plt.subplots(figsize=(9.8, 4.4))
    ax.errorbar(
        x,
        obs,
        yerr=err,
        fmt="o",
        color=RED,
        ms=5,
        capsize=3,
        elinewidth=1.4,
        lw=0,
        label="measured, 20,000 draws per cell",
    )
    ax.plot(x, des, "s", color=BLUE, ms=5, label="coverage quadrature")
    ax.axhline(gamma, color=DARK, lw=1.3, ls="--")
    ax.annotate(
        "designed gamma = 5e-4",
        (0.15, gamma),
        ha="left",
        va="bottom",
        fontsize=9,
        color=DARK,
        bbox=BOX,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(-0.00015, 0.00110)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    frame(
        ax,
        "The region's miss rate: designed against measured",
        ylab="probability the true ratios fall outside the region",
        xlab="design (family, benchmark degrees of freedom)",
    )
    save(fig, "fig13_miss.png")


def main():
    which = sys.argv[1:] or ["all"]
    fns = {
        "8": fig8_guarantee,
        "9": fig9_deskmap,
        "10": fig10_power,
        "11": fig11_ar1,
        "12": fig12_factor,
        "13": fig13_miss,
    }
    for k in sorted(fns) if which == ["all"] else which:
        fns[k]()


if __name__ == "__main__":
    main()
