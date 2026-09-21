"""S3.2: selector v2 -- route a design to the plug-in rule or to the envelope.

The decision is per design, before any test is run: if the exact plug-in rule is
believed accurate to within BAND percent for this design's shape, use it (it is
the cheapest and most powerful option); otherwise route to the envelope rule
(which carries the alpha+gamma guarantee).

v1 used two features (candidate df, variance spread) and reached 0.752 held-out
accuracy on a proper random row split; that is documented as "weak but usable as
a default-with-fallback".  v2 adds the features the design actually exposes:

  K                 number of candidates
  bench_df          benchmark (control) degrees of freedom
  min_cand_df       smallest candidate df
  spread            max(lambda_hat)/min(lambda_hat)
  power_cost_pct    measured power cost of the envelope on that design
  alpha

A greedy depth<=3 classification tree is fitted on the TRAINING half of a fixed
random row split (seed 20260916, pre-registered) and evaluated on the held-out
half.  Bar: held-out accuracy >= 0.85, else document and stop.

Outputs: deposits/selector/selector2.json
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_util import dump_json, f64, i64  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "selector")
MAP = os.path.join(HERE, "deposits", "w5_summary_plugin.json")
SPLIT_SEED = 20260916
BAND = 2.0
BAR = 0.85
MAX_DEPTH = 3
MIN_LEAF = 20
FEATS = ("K", "bench_df", "min_cand_df", "spread", "power_cost_pct", "alpha")


try:
    NAN = float("nan")
except Exception:  # noqa: BLE001
    NAN = 0.0


def _f(x, default=None):
    """Guarded float conversion; a None default yields NaN."""
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return NAN if default is None else float(default)


def _i(x, default=0):
    try:
        return int(x)
    except Exception:  # noqa: BLE001
        return int(default)


def bench_cand_dfs(shape, m):
    """(benchmark df, smallest candidate df) implied by the shape and m = df."""
    if shape == "bench_small":
        return m, 13 - 1
    if shape == "short_cand":
        return 36 - 1, m
    return m, m


def load_rows(path=MAP):
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"cannot read map {path}: {exc}") from exc
    rows = []
    for cell in data["cells"]:
        try:
            lam = [_f(v) for v in cell["lam"]]
            spread = max(lam) / min(lam)
            shape = str(cell["shape"])
            K = _i(cell["K"])
            alpha = _f(cell["alpha"])
            pc = _f(cell.get("power_cost_pct"), 0.0)
            for df, err in zip(cell["dfs"], cell["rel_err_pct"], strict=True):
                b_df, c_df = bench_cand_dfs(shape, _i(df))
                rows.append(
                    {
                        "K": i64(K),
                        "shape": shape,
                        "alpha": f64(alpha),
                        "bench_df": i64(b_df),
                        "min_cand_df": i64(c_df),
                        "spread": f64(spread),
                        "power_cost_pct": f64(pc),
                        "err_pct": f64(_f(err)),
                        "inside": bool(abs(_f(err)) <= BAND),
                    }
                )
        except Exception as exc:  # noqa: BLE001
            raise SystemExit(f"malformed cell K={cell.get('K')}: {exc}") from exc
    return rows


def split(rows, seed=SPLIT_SEED):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(rows))
    half = len(rows) // 2
    return (
        [rows[i] for i in sorted(idx[:half])],
        [rows[i] for i in sorted(idx[half:])],
    )


def misclass(rows, keep_fn):
    return sum(1 for r in rows if keep_fn(r) != r["inside"])


def best_split(rows):
    """Greedy: best single (feature, threshold, direction) by misclassification."""
    best = None
    for feat in FEATS:
        vals = sorted({_f(r[feat]) for r in rows if np.isfinite(_f(r[feat]))})
        if len(vals) < 2:
            continue
        qs = np.quantile(vals, np.linspace(0.1, 0.9, 9))
        for thr in sorted({_f(v) for v in qs}):
            for le in (True, False):

                def keep_fn(r, _f_=feat, _t=thr, _le=le):
                    v = _f(r[_f_])
                    return bool(v <= _t) if _le else bool(v > _t)

                left = [r for r in rows if keep_fn(r)]
                right = [r for r in rows if not keep_fn(r)]
                if len(left) < MIN_LEAF or len(right) < MIN_LEAF:
                    continue
                err = misclass(rows, keep_fn)
                if best is None or err < best["err"]:
                    best = {
                        "feat": feat,
                        "thr": f64(thr),
                        "le": bool(le),
                        "err": i64(err),
                        "left": left,
                        "right": right,
                    }
    return best


def fit_tree(rows, depth=0):
    base = misclass(rows, lambda r: True)
    base0 = misclass(rows, lambda r: False)
    if depth >= MAX_DEPTH or min(base, base0) == 0:
        keep = base <= base0
        return {
            "leaf": True,
            "keep": bool(keep),
            "err": i64(min(base, base0)),
            "n": i64(len(rows)),
        }
    sp = best_split(rows)
    best_const = min(base, base0)
    if sp is None or sp["err"] >= best_const - 1:
        keep = base <= base0
        return {
            "leaf": True,
            "keep": bool(keep),
            "err": i64(best_const),
            "n": i64(len(rows)),
        }
    left = fit_tree(sp["left"], depth + 1)
    right = fit_tree(sp["right"], depth + 1)
    return {
        "leaf": False,
        "feat": sp["feat"],
        "thr": sp["thr"],
        "le": sp["le"],
        "left": left,
        "right": right,
        "err": i64(left["err"] + right["err"]),
        "n": i64(len(rows)),
    }


def predict(node, r):
    while not node["leaf"]:
        v = _f(r[node["feat"]])
        node = node["left"] if (v <= node["thr"]) == node["le"] else node["right"]
    return bool(node["keep"])


def evaluate(rows, tree):
    bad = [r for r in rows if predict(tree, r) != r["inside"]]
    kept = [r for r in rows if predict(tree, r)]
    errs = [abs(_f(r["err_pct"])) for r in kept]
    med = _f(np.median(errs)) if errs else NAN
    mx = _f(np.max(errs)) if errs else NAN
    return {
        "n": i64(len(rows)),
        "accuracy": f64(1.0 - len(bad) / max(1, len(rows))),
        "n_misclassified": i64(len(bad)),
        "frac_routed_to_envelope": f64(1.0 - len(kept) / max(1, len(rows))),
        "median_abs_err_kept_pct": f64(med),
        "max_abs_err_kept_pct": f64(mx),
        "missed": [
            {
                "K": r["K"],
                "shape": r["shape"],
                "df": r["min_cand_df"],
                "spread": r["spread"],
                "err_pct": r["err_pct"],
            }
            for r in bad[:10]
        ],
    }


def counters(tree):
    """Every threshold the tree actually uses (for the write-up)."""
    out = []
    stack = [tree]
    while stack:
        n = stack.pop()
        if n["leaf"]:
            continue
        out.append({"feat": n["feat"], "thr": n["thr"], "le": n["le"]})
        stack += [n["left"], n["right"]]
    return out


def main():
    rows = load_rows()
    train, test = split(rows)
    tree = fit_tree(train)
    tr = evaluate(train, tree)
    te = evaluate(test, tree)
    # v1 reference: candidate df + spread thresholds, fitted the same way
    best = None
    for t in (3, 5, 7, 9, 11, 13, 19, 35):
        for s in (1.5, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
            err = sum(
                1
                for r in train
                if ((r["min_cand_df"] >= t) and (r["spread"] <= s)) != r["inside"]
            )
            if best is None or err < best[0]:
                best = (err, t, s)
    _, t1, s1 = best
    v1_pred = [bool((r["min_cand_df"] >= t1) and (r["spread"] <= s1)) for r in test]
    v1_acc = 1.0 - sum(
        1 for p, r in zip(v1_pred, test, strict=True) if p != r["inside"]
    ) / len(test)
    out = {
        "band_pct": f64(BAND),
        "bar": f64(BAR),
        "split_seed": i64(SPLIT_SEED),
        "features": list(FEATS),
        "max_depth": i64(MAX_DEPTH),
        "n_rows": i64(len(rows)),
        "n_train": i64(len(train)),
        "n_test": i64(len(test)),
        "tree": tree,
        "thresholds_used": counters(tree),
        "train": tr,
        "held_out": te,
        "v1_reference": {
            "df_thr": i64(t1),
            "spread_thr": f64(s1),
            "held_out_accuracy": f64(v1_acc),
        },
        "bar_met": bool(te["accuracy"] >= BAR),
        "verdict": (
            "keep as a default-with-fallback"
            if te["accuracy"] >= BAR
            else "below bar: document and stop (fallback stays primary)"
        ),
    }
    dump_json(out, os.path.join(OUT, "selector2.json"))
    print(
        f"held-out accuracy {te['accuracy']:.3f} (v1 {v1_acc:.3f}) | bar "
        f"{BAR:.2f} -> {'MET' if out['bar_met'] else 'NOT MET'} | "
        f"routed {te['frac_routed_to_envelope'] * 100:.1f}% | median err kept "
        f"{te['median_abs_err_kept_pct']:.3f}%"
    )
    print("thresholds used:", counters(tree))


if __name__ == "__main__":
    main()
