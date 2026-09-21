"""Track C: desk selector -- which rule to use, decided from observable features.

Feature set (all observable at run time): K, benchmark df, the spread of the
estimated variance ratios, and alpha.  The label is whether the plug-in rule's
error is inside the +/-2% band, measured RELATIVE to each cell's own alpha (so
pooling alpha = 0.025 and 0.05 rows is legitimate; note the absolute tolerance is
twice as loose at alpha = 0.05, while the envelope's own bar is an ABSOLUTE
+0.0005 -- the two criteria are not the same).

Verdict history: the first version split the 84-cell map by a content hash of the
design key and reported held-out accuracy 0.646 as evidence that the two-feature
family "cannot predict".  An audit showed that split was GROUP-DISJOINT by
accident (every design group landed wholly on one side, so the halves were
different mixtures of designs, train prevalence 0.559 vs test 0.444), and that
another member of the SAME family (df >= 11) scores 0.803 on the same held-out
rows.  So 0.646 was a property of that split, not of the feature set.  This
version reports a proper random row-level split with a fixed seed, keeps the old
split's numbers for transparency, and reports both the train-fitted and the
test-fitted thresholds so the generalisation gap is visible.

Outputs: deposits/selector/rule.json and the one-page desk map
deposits/selector/desk_map.md.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_util import dump_json, f64, i64, load_json, write_text  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "deposits", "selector")
BAND = 2.0
SPLIT_SEED = 20260916
DF_GRID = (1, 3, 5, 7, 11, 19, 35)
SPREAD_GRID = (1.5, 4.0, 25.0, 100.0, 625.0)


def _median(xs, default=np.nan):
    try:
        return float(np.median(xs)) if xs else float(default)
    except Exception:  # noqa: BLE001
        return float(default)


def _max(xs, default=np.nan):
    try:
        return float(np.max(xs)) if xs else float(default)
    except Exception:  # noqa: BLE001
        return float(default)


def _key(cell):
    return f"{cell['K']}|{cell['shape']}|{cell['vol']}|{cell['alpha']}"


def _keep(df, spread, t, s):
    return (df >= t) and (spread <= s)


def rows_from_map(path):
    data = load_json(path)
    rows = []
    for cell in data["cells"]:
        try:
            lam = [float(v) for v in cell["lam"]]
            spread = max(lam) / min(lam)
            group = f"{cell['K']}|{cell['shape']}|{cell['vol']}"
            for df, err in zip(cell["dfs"], cell["rel_err_pct"], strict=False):
                rows.append(
                    {
                        "key": _key(cell),
                        "group": group,
                        "K": int(cell["K"]),
                        "shape": str(cell["shape"]),
                        "vol": str(cell["vol"]),
                        "alpha": float(cell["alpha"]),
                        "df": int(df),
                        "spread": float(spread),
                        "err_pct": float(err),
                        "inside": bool(abs(float(err)) <= BAND),
                    }
                )
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(f"malformed W5 cell {cell.get('K')}") from exc
    return rows


def split_rows(rows, seed=SPLIT_SEED):
    """Proper random row-level split with a fixed seed.

    Rows within a design group are dependent (same design, two alphas), so the
    accuracy carries a clustered standard error larger than binomial; the
    group-disjoint alternative is reported separately for transparency.
    """
    try:
        rng = np.random.default_rng(seed)
        idx = rng.permutation(len(rows))
        half = len(rows) // 2
        train = [rows[i] for i in sorted(idx[:half])]
        test = [rows[i] for i in sorted(idx[half:])]
        return train, test
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("row split failed") from exc


def split_group_disjoint(rows):
    """The earlier content-hash split (kept so the refuted number is traceable)."""
    tr = [r for r in rows if sum(ord(ch) for ch in r["key"]) % 2 == 0]
    te = [r for r in rows if sum(ord(ch) for ch in r["key"]) % 2 == 1]
    return tr, te


def fit(rows):
    """Pick (df_thr, spread_thr) minimising training misclassification."""
    best = None
    try:
        for t in DF_GRID:
            for s in SPREAD_GRID:
                err = sum(
                    int(_keep(r["df"], r["spread"], t, s) != r["inside"]) for r in rows
                )
                if best is None or err < best[0]:
                    best = (err, t, s)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("threshold fitting failed") from exc
    if best is None:
        raise RuntimeError("no threshold candidate evaluated")
    return best


def evaluate(rows, t, s):
    bad = [r for r in rows if _keep(r["df"], r["spread"], t, s) != r["inside"]]
    kept = [r for r in rows if _keep(r["df"], r["spread"], t, s)]
    routed = [r for r in rows if not _keep(r["df"], r["spread"], t, s)]
    errs = [abs(r["err_pct"]) for r in kept]
    try:
        acc = 1.0 - len(bad) / max(1, len(rows))
        frac_routed = len(routed) / max(1, len(rows))
    except Exception:  # noqa: BLE001
        acc = frac_routed = float("nan")
    return {
        "n": i64(len(rows)),
        "n_misclassified": i64(len(bad)),
        "accuracy": f64(acc),
        "n_routed_to_envelope": i64(len(routed)),
        "frac_routed": f64(frac_routed),
        "median_abs_err_plug_in_pct": f64(_median(errs)),
        "max_abs_err_plug_in_pct": f64(_max(errs)),
        "missed": [
            {
                "key": r["key"],
                "df": r["df"],
                "spread": f64(r["spread"]),
                "err_pct": f64(r["err_pct"]),
            }
            for r in bad[:12]
        ],
    }


def desk_families():
    """The two desk families from the envelope sweep (4 seeds each)."""
    rows = []
    d = os.path.join(HERE, "deposits", "env")
    for fam in ("A", "B"):
        for df in (3, 5, 7, 11, 19, 35):
            p = os.path.join(d, f"sweep_df{df}_{fam}.json")
            if not os.path.exists(p):
                continue
            cell = load_json(p)
            per = cell.get("per_seed", [])
            errs = []
            for r in per:
                try:
                    errs.append(abs(float(r["plugin_fwer"]) / 0.025 - 1.0) * 100.0)
                except Exception as exc:  # noqa: BLE001
                    print(f"  desk-family row unusable: {exc!r}")
            med = _median(errs)
            rows.append(
                {
                    "family": fam,
                    "df": i64(df),
                    "median_abs_err_pct": f64(med),
                    "n_seeds": i64(len(errs)),
                    "inside": bool(errs and med <= BAND),
                }
            )
    return rows


def main():
    rows = rows_from_map(os.path.join(HERE, "deposits", "w5_summary_plugin.json"))
    train, test = split_rows(rows)
    err_tr, t_tr, s_tr = fit(train)
    err_te, t_te, s_te = fit(test)
    tr_eval = evaluate(train, t_tr, s_tr)
    te_eval = evaluate(test, t_tr, s_tr)
    te_opt = evaluate(test, t_te, s_te)
    legacy_tr, legacy_te = split_group_disjoint(rows)
    l_err, l_t, l_s = fit(legacy_tr)
    legacy = {
        "thresholds": {"df": i64(l_t), "spread": f64(l_s)},
        "train": evaluate(legacy_tr, l_t, l_s),
        "test": evaluate(legacy_te, l_t, l_s),
    }
    fams = desk_families()
    out = {
        "rule_family": "plug-in iff bench_df >= t AND ratio_spread <= s",
        "band_pct": BAND,
        "split": {
            "kind": "random row-level, fixed seed",
            "seed": i64(SPLIT_SEED),
            "n_train": i64(len(train)),
            "n_test": i64(len(test)),
        },
        "train_fitted": {
            "thresholds": {"df": i64(t_tr), "spread": f64(s_tr)},
            "train_misclassified_at_fit": i64(err_tr),
            "train": tr_eval,
            "test": te_eval,
        },
        "test_fitted": {
            "thresholds": {"df": i64(t_te), "spread": f64(s_te)},
            "train_misclassified_at_fit": i64(err_te),
            "test": te_opt,
        },
        "legacy_group_disjoint_split": legacy,
        "desk_families": fams,
        "note": (
            "the previous report quoted the legacy split's 0.646 as a "
            "property of the feature family; the audit showed the split was "
            "group-disjoint and that another family member scores 0.803 on "
            "the same rows, so the honest statements are (a) the family "
            "contains a rule with test accuracy about 0.80, (b) the "
            "train-fitted thresholds generalise worse by the gap below, and "
            "(c) rows within a design group are dependent, so the clustered "
            "SE exceeds the binomial one."
        ),
        "desk_families_all_inside": bool(fams and all(r["inside"] for r in fams)),
    }
    dump_json(out, os.path.join(OUT, "rule.json"))
    lines = [
        "# Desk map: which rule to use",
        "",
        "Selector family: **plug-in iff benchmark df >= t and ratio spread <= s**, "
        "else the guaranteed envelope (K <= 3) or an unmeasured cell (K = 5/8).",
        "",
        f"- random row-level split (seed {SPLIT_SEED}), "
        f"{len(train)} train / {len(test)} test rows",
        f"- train-fitted thresholds (`df >= {t_tr}`, spread <= {s_tr:g}): "
        f"train accuracy **{tr_eval['accuracy']:.3f}**, "
        f"held-out accuracy **{te_eval['accuracy']:.3f}**",
        f"- test-fitted thresholds (`df >= {t_te}`, spread <= {s_te:g}): "
        f"held-out accuracy **{te_opt['accuracy']:.3f}**",
        f"- legacy group-disjoint split, same family: train "
        f"{legacy['train']['accuracy']:.3f} / test {legacy['test']['accuracy']:.3f} "
        f"(the number quoted in the first report; NOT a property of the family)",
        f"- desk families inside +/-2% under the plug-in: "
        f"**{out['desk_families_all_inside']}** (median over 4 seeds per cell)",
        "",
        "## Desk families",
        "",
        "| family | bench df | median |plug-in error| | seeds | inside +/-2% |",
        "|---|---|---|---|---|",
    ]
    for r in fams:
        lines.append(
            f"| {r['family']} | {r['df']} | {r['median_abs_err_pct']:.2f}% | "
            f"{r['n_seeds']} | {r['inside']} |"
        )
    write_text("\n".join(lines) + "\n", os.path.join(OUT, "desk_map.md"))
    print(f"split: random seed {SPLIT_SEED}, {len(train)}/{len(test)} rows")
    print(
        f"train-fitted df>={t_tr} spread<={s_tr:g}: train "
        f"{tr_eval['accuracy']:.3f} -> test {te_eval['accuracy']:.3f}"
    )
    print(f"test-fitted  df>={t_te} spread<={s_te:g}: test {te_opt['accuracy']:.3f}")
    print(
        f"legacy split: train {legacy['train']['accuracy']:.3f} / "
        f"test {legacy['test']['accuracy']:.3f}"
    )
    print(f"desk families all inside: {out['desk_families_all_inside']}")


if __name__ == "__main__":
    main()
