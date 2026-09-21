"""W5 analysis: read the cell deposits and locate the 2 per cent crossings.

Usage: python w5_analyse.py [rule]
Writes deposits/w5_summary_<rule>.json and prints a compact table.
"""

import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_util import dump_json  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "deposits", "w5")
NAN = math.nan


def load_cells(rule):
    cells = {}
    for p in sorted(glob.glob(os.path.join(OUTDIR, f"{rule}_*.json"))):
        try:
            with open(p, encoding="utf-8") as fh:
                d = json.load(fh)
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(f"bad deposit {p}") from exc
        key = (d["K"], d["shape"], d["vol"], d["alpha"])
        cells[key] = d
    return cells


def crossing(rows, alpha, band=2.0):
    """First df from which |rel err| stays inside the band."""
    dfs = [r["df"] for r in rows]
    err = [abs(r["rel_err_pct"]) for r in rows]
    keep = [df for df, e in zip(dfs, err, strict=False) if e <= band]
    return keep[0] if keep else None, dfs, err


def main(rule="plugin"):
    cells = load_cells(rule)
    out = {"rule": rule, "band_pct": 2.0, "cells": []}
    print(f"{rule}: {len(cells)} cells")
    print(
        f"{'K':>2} {'shape':<12} {'vol':<12} {'alpha':>6} "
        f"{'err@df3':>8} {'err@df5':>8} {'err@df7':>8} {'err@df35':>9} "
        f"{'first|e|<=2':>11} {'max power%':>10}"
    )
    for key in sorted(cells):
        d = cells[key]
        rows = d["rows"]
        first, dfs, err = crossing(rows, d["alpha"])
        pw = max(r["power_cost_pct"] for r in rows)
        e_by_df = {r["df"]: r["rel_err_pct"] for r in rows}
        out["cells"].append(
            {
                "K": d["K"],
                "shape": d["shape"],
                "vol": d["vol"],
                "alpha": d["alpha"],
                "lam": d["lam"],
                "dfs": dfs,
                "rel_err_pct": err,
                "seed_se_pct": [100 * r["seed_se"] / d["alpha"] for r in rows],
                "power_cost_pct": [r["power_cost_pct"] for r in rows],
                "first_df_inside_band": first,
                "worst_rel_err_pct": max(err),
                "worst_power_cost_pct": pw,
            }
        )
        print(
            f"{d['K']:>2} {d['shape']:<12} {d['vol']:<12} {d['alpha']:>6.3f} "
            f"{e_by_df.get(3, NAN):>8.2f} "
            f"{e_by_df.get(5, NAN):>8.2f} "
            f"{e_by_df.get(7, NAN):>8.2f} "
            f"{e_by_df.get(35, NAN):>9.2f} "
            f"{str(first):>11} {pw:>10.2f}"
        )
    p = os.path.join(HERE, "deposits", f"w5_summary_{rule}.json")
    dump_json(out, p)
    print("wrote", p)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "plugin")
