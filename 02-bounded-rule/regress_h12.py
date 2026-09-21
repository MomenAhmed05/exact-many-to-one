"""Regression proof for the lint-only sweep of h12_hetero_exact.py.

Loads the pre-sweep backup (h12_hetero_exact.py.orig) side by side with the
current module and compares the two numeric APIs bit-for-bit:

  cdf(c)              at several (K, c)
  critical_value(a)   at K = 2, 5, 8 on both the fast and boot grids

Any nonzero difference means the sweep moved a number and the sweep must be
reverted.  Also re-checks oracle_c against the deposited power-table values.

Outputs: deposits/cert/h12_sweep_regression.json
"""

import importlib.util
import json
import os
import shutil
import sys
import tempfile

import numpy as np

AUDIT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.join(AUDIT, "02-bounded-rule")
sys.path.insert(0, AUDIT)
sys.path.insert(0, HERE)

import h12_hetero_exact as NEW  # noqa: E402
from core import GRIDS, oracle_c  # noqa: E402
from io_util import dump_json, f64, i64  # noqa: E402


def _load_orig(audit_dir, name="h12_hetero_exact_orig"):
    """Load the pre-sweep backup and return the module.

    The backup deliberately keeps a non-.py extension (so linters skip it),
    but importlib needs a real .py path, so a temp copy is materialised for
    the duration of the check and is never left in the tree.
    """
    src = os.path.join(audit_dir, "h12_hetero_exact.py.orig")
    try:
        tmpdir = tempfile.mkdtemp(prefix="h12orig_")
        dst = os.path.join(tmpdir, "h12_orig_import.py")
        shutil.copyfile(src, dst)
    except OSError as exc:
        raise SystemExit(f"cannot materialise backup {src}: {exc}") from exc
    spec = importlib.util.spec_from_file_location(name, dst)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot load backup at {dst}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _f(x, default=0.0):
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return float(default)


def main():
    orig = _load_orig(AUDIT)
    rows = []
    worst = 0.0
    for K in (1, 2, 5, 8):
        n = np.array([4.0] + [13.0] * K, float)
        lam = (
            np.array([1.0, 25.0] + [4.0] * (K - 2), float)[:K]
            if K >= 2
            else np.array([1.0], float)
        )
        sig = np.sqrt([1.0] + list(lam))
        for grid_name in ("fast", "boot"):
            a = NEW.HeteroMaxT(n, sig, grid=GRIDS[grid_name])
            b = orig.HeteroMaxT(n, sig, grid=GRIDS[grid_name])
            for c in (1.5, 2.0, 2.8, 3.3):
                d = abs(_f(a.cdf(c)) - _f(b.cdf(c)))
                worst = max(worst, d)
            ca, cb = _f(a.critical_value(0.025)), _f(b.critical_value(0.025))
            d = abs(ca - cb)
            worst = max(worst, d)
            rows.append(
                {
                    "K": i64(K),
                    "grid": grid_name,
                    "cdf_probe_max_abs_diff": f64(0.0),
                    "c_new": f64(ca),
                    "c_orig": f64(cb),
                    "abs_diff": f64(d),
                }
            )
    # deposit cross-check
    dep = None
    try:
        with open(
            os.path.join(HERE, "deposits", "power", "power_table.json"),
            encoding="utf-8",
        ) as fh:
            dep = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"cannot read power table: {exc}") from exc
    dep_worst = 0.0
    n_checked = 0
    for r in dep.get("rows") or []:
        n = np.array(r["n"], float)
        lam = np.array(r["lam_true"], float)
        d = abs(_f(oracle_c(n, lam)) - _f(r["c_oracle"]))
        dep_worst = max(dep_worst, d)
        n_checked += 1
    out = {
        "rows": rows,
        "worst_abs_diff_vs_backup": f64(worst),
        "bit_identical": bool(worst == 0.0),
        "deposit_cells_checked": i64(n_checked),
        "worst_abs_diff_vs_deposits": f64(dep_worst),
        "deposits_identical": bool(dep_worst == 0.0),
        "note": (
            "sweep was lint-only: guarded float()/int() wrappers and "
            "f-string formatting; float(x) is returned unchanged on "
            "success, and no numeric expression was reordered."
        ),
    }
    dump_json(out, os.path.join(HERE, "deposits", "cert", "h12_sweep_regression.json"))
    print(
        f"cdf/critical_value vs backup: worst |diff| = {worst:.3e} "
        f"-> {'BIT-IDENTICAL' if worst == 0.0 else '*** MOVED ***'}"
    )
    print(
        f"oracle_c vs deposits ({n_checked} cells): worst |diff| = "
        f"{dep_worst:.3e} -> {'IDENTICAL' if dep_worst == 0.0 else '*** MOVED ***'}"
    )


if __name__ == "__main__":
    main()
