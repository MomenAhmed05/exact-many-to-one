"""W0: resolve the unresolved 3-df discrepancy and certify the fast solver.

The question: the paper's df-3 plug-in sweep column reads +13.8%, the independent
replication reads +8.9%/+9.8% on the two volatility orderings; which
ordering reproduces +13.8%?  The short answer this script tests: neither ordering
is at fault; the replication's vectorised per-draw solver ran only 8 bisection
steps, which leaves a ~1e-3 error in the critical value; and that error is worth
about 4 points of relative FWER at 3 df.

Outputs: 02-bounded-rule/deposits/w0_discrepancy.json
"""

import json
import os
import sys
import time
from math import sqrt

import numpy as np

sys.path.insert(
    0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)  # tools/audit
from io_util import dump_json, f64  # noqa: E402

sys.path.insert(0, os.environ.get("REPLICATION_DIR", "."))  # folder holding the replication modules

import exact_manytoone_replication as R  # noqa: E402
from h12_hetero_exact import HeteroMaxT  # noqa: E402

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "deposits", "w0_discrepancy.json"
)
GRID = {"u": (14, 14), "v": (14, 14), "s": (4, 14)}  # paper working grid
FAST = {"u": (8, 10), "v": (8, 10), "s": (3, 10)}  # paper fast grid
ALPHA = 0.025


def cond_noreject(cvec, s0, si, n, lam_true, npanels=40, order=10):
    """P(no rejection | drawn variances), TRUE lam; 1-D integral over U."""
    n0 = n[0]
    ni = np.asarray(n[1:], float)
    z, a = R.var_grid(npanels, order)
    d = len(cvec)
    prod = np.ones((d, len(z)))
    for i in range(len(ni)):
        b = cvec[:, None] * np.sqrt(si[:, i][:, None] / ni[i] + s0[:, None] / n0)
        arg = np.sqrt(ni[i] / lam_true[i]) * (z[None, :] / np.sqrt(n0) + b)
        prod = prod * R.ndtr(arg)
    return (prod * a[None, :]).sum(axis=1)


def run(
    design_label, n_l, sig_l, draws=400, seed=9001, solvers=("fast", "repl", "repl26")
):
    n = np.asarray(n_l, float)
    sig = np.asarray(sig_l, float)
    lam_true = list((sig / sig[0]) ** 2)[1:]
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    s0 = rng.chisquare(nu[0], draws) / nu[0]
    si = np.column_stack(
        [
            lam_true[i] * rng.chisquare(nu[i + 1], draws) / nu[i + 1]
            for i in range(len(lam_true))
        ]
    )
    lh = si / s0[:, None]

    c_or = HeteroMaxT(n, sig, grid=GRID).critical_value(ALPHA)
    oracle = cond_noreject(np.full(draws, c_or), s0, si, n, lam_true)

    out = {
        "design": design_label,
        "n": n_l,
        "sigma": sig_l,
        "lam_true": lam_true,
        "draws": draws,
        "seed": seed,
        "c_oracle_base_grid": c_or,
        "alpha": ALPHA,
        "solvers": {},
    }
    t0 = time.time()
    cs = {}
    if "fast" in solvers:
        cs["fast"] = np.array(
            [
                HeteroMaxT(n, np.sqrt([1.0] + list(lv)), grid=FAST).critical_value(
                    ALPHA
                )
                for lv in lh
            ]
        )
    if "base" in solvers:
        cs["base"] = np.array(
            [
                HeteroMaxT(n, np.sqrt([1.0] + list(lv)), grid=GRID).critical_value(
                    ALPHA
                )
                for lv in lh
            ]
        )
    if "repl" in solvers:
        cs["repl"] = R.exact_c_vec(n, lh, alpha=ALPHA, c0=c_or, half=0.45, iters=8)
    if "repl26" in solvers:
        cs["repl26"] = R.exact_c_vec(n, lh, alpha=ALPHA, c0=c_or, half=0.45, iters=26)
    for k, c in cs.items():
        nr = cond_noreject(c, s0, si, n, lam_true)
        d = oracle - nr
        out["solvers"][k] = {
            "mean_c": f64(c.mean()),
            "sd_c": f64(c.std(ddof=1)),
            "fwer": ALPHA + f64(d.mean()),
            "fwer_mcse": f64(d.std(ddof=1) / sqrt(draws)),
            "rel_err_pct": 100 * f64(d.mean()) / ALPHA,
        }
    # solver-vs-solver critical-value error on the identical draws
    if "base" in cs:
        for k in cs:
            out["solvers"][k]["max_abs_dc_vs_base"] = f64(
                np.abs(cs[k] - cs["base"]).max()
            )
            out["solvers"][k]["mean_abs_dc_vs_base"] = f64(
                np.abs(cs[k] - cs["base"]).mean()
            )
    if "repl26" in cs and "base" not in cs:
        for k in cs:
            out["solvers"][k]["max_abs_dc_vs_repl26"] = f64(
                np.abs(cs[k] - cs["repl26"]).max()
            )
    out["seconds"] = time.time() - t0
    return out


def grid_study(n_l=(4.0, 13.0, 13.0), lam=(1.0, 25.0), draws=12, seed=77):
    """Certify the vectorised solver's grid error against the base grid."""
    n = np.asarray(n_l, float)
    nu = n - 1.0
    rng = np.random.default_rng(seed)
    s0 = rng.chisquare(nu[0], draws) / nu[0]
    si = np.column_stack(
        [lam[i] * rng.chisquare(nu[i + 1], draws) / nu[i + 1] for i in range(len(lam))]
    )
    lh = si / s0[:, None]
    ref = np.array(
        [
            HeteroMaxT(n, np.sqrt([1.0] + list(lv)), grid=GRID).critical_value(ALPHA)
            for lv in lh
        ]
    )
    rows = []
    for g, order in (
        ((24, 6, 4), (5, 5, 8)),
        ((32, 10, 6), (6, 6, 10)),
        ((40, 12, 8), (6, 6, 10)),
        ((48, 16, 8), (8, 8, 12)),
    ):
        t0 = time.time()
        c = R.exact_c_vec(
            n,
            lh,
            alpha=ALPHA,
            c0=f64(ref.mean()),
            half=0.6,
            iters=26,
            g=g,
            orders=order,
        )
        rows.append(
            {
                "grid": list(g),
                "orders": list(order),
                "max_abs_dc": f64(np.abs(c - ref).max()),
                "mean_abs_dc": f64(np.abs(c - ref).mean()),
                "seconds": time.time() - t0,
            }
        )
    return rows


if __name__ == "__main__":
    res = {
        "alpha": ALPHA,
        "note": (
            "paper +13.8% is n=(4,13,13), sigma=(10,10,50), i.e. "
            "lam=(1,25); replication family B is the same design"
        ),
    }
    # the paper's df-3 design (family B)
    res["df3_familyB"] = run(
        "df3 family B n=(4,13,13) sd=(10,10,50)",
        [4, 13, 13],
        [10.0, 10.0, 50.0],
        draws=400,
        seed=9001,
        solvers=("fast", "base", "repl", "repl26"),
    )
    # the other ordering, for completeness
    res["df3_familyA"] = run(
        "df3 family A n=(4,13,13) sd=(50,10,10)",
        [4, 13, 13],
        [50.0, 10.0, 10.0],
        draws=400,
        seed=9002,
        solvers=("fast", "repl", "repl26"),
    )
    print(json.dumps(res, indent=1))
    dump_json(res, OUT)
    print("wrote", OUT)
