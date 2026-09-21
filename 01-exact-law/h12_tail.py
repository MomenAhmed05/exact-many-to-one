"""H12(b)+(c) — heavy-tail innovations and the deployment table.

(b) H12's exact estimated-scale max-t critical value was derived under Gaussian
    errors. Strategy panels are fat-tailed. This arm replaces the Gaussian errors
    with unit-variance t5 and measures the size of using the Gaussian-derived
    critical value.

    Design: K arms + one control, n obs each, equal variances; the per-arm
    statistic is the pooled-scale t
        t_j = (xbar_j - xbar_0) / (s_p * sqrt(2/n)),   s_p^2 = pooled within-arm
    (K+1)*(n-1) df. Under Gaussian errors this is exactly multivariate t; under
    t5 errors it is not. The 95% critical value of max|t| is measured by MC in
    both cases and the induced size of the Gaussian CV under t5 is reported.

(c) Deployment table: Gaussian-error critical values over K and n, alongside the
    known-scale (Dunnett-type) value, for the same design.

Guarded blocks re-raise with a stage label.
"""

import json
import math

import numpy as np

OUT = "tools/audit/series/h12_tail.json"
K_GRID = (3, 5, 10, 20)
N_GRID = (60, 120, 250, 500, 1000)
R_MC = 20000
SEEDS = (69001, 69002)


def max_t(rng, k_arms, n, family, r_draws):
    """Max |pooled-scale t| across K arms vs one control."""
    try:
        if family == "gauss":

            def draw(shape):
                return rng.standard_normal(shape)
        elif family == "t5":
            s = math.sqrt(5.0 / 3.0)

            def draw(shape):
                return rng.standard_t(5, size=shape) / s
        else:
            raise RuntimeError(f"unknown family {family}")
        out = np.empty(r_draws)
        for i in range(r_draws):
            x = draw((k_arms + 1, n))
            means = x.mean(axis=1)
            within = x - means[:, None]
            df = (k_arms + 1) * (n - 1)
            sp2 = float((within**2).sum() / df)
            t = (means[1:] - means[0]) / math.sqrt(sp2 * 2.0 / n)
            out[i] = np.abs(t).max()
        return out
    except RuntimeError:
        raise
    except Exception as exc:
        raise RuntimeError("max_t failed") from exc


def known_scale_cv(rng, k_arms, r_draws):
    """Max|Z| over K iid N(0,1) — the known-scale reference."""
    try:
        z = rng.standard_normal((r_draws, k_arms))
        return float(np.quantile(np.abs(z).max(axis=1), 0.95))
    except Exception as exc:
        raise RuntimeError("known_scale_cv failed") from exc


def run():
    res = {"tail_arm": {}, "table": {}, "falsified": []}
    # ---- (b) heavy tails ----
    for k_arms, n in ((5, 60), (5, 250), (10, 250)):
        try:
            rng = np.random.default_rng(SEEDS[0] + k_arms * 13 + n)
            g = max_t(rng, k_arms, n, "gauss", R_MC)
            rng2 = np.random.default_rng(SEEDS[1] + k_arms * 17 + n)
            t = max_t(rng2, k_arms, n, "t5", R_MC)
            cv_g = float(np.quantile(g, 0.95))
            cv_t = float(np.quantile(t, 0.95))
            size_g_cv_under_t = float((t > cv_g).mean())
            size_t_cv_under_t = float((t > cv_t).mean())
        except RuntimeError:
            raise
        except Exception as exc:
            raise RuntimeError(f"tail arm failed K={k_arms} n={n}") from exc
        key = f"K{k_arms}/n{n}"
        res["tail_arm"][key] = {
            "cv_gauss": round(cv_g, 4),
            "cv_t5": round(cv_t, 4),
            "size_of_gauss_cv_under_t5": round(size_g_cv_under_t, 4),
            "size_of_t5_cv_under_t5": round(size_t_cv_under_t, 4),
            "cv_gap": round(cv_t - cv_g, 4),
        }
        print(
            f"TAIL {key}: cv_gauss {cv_g:.4f} vs cv_t5 {cv_t:.4f} "
            f"(gap {cv_t - cv_g:+.4f}); using the gauss CV under t5 gives size "
            f"{size_g_cv_under_t:.4f}",
            flush=True,
        )
        if size_g_cv_under_t > 0.075:
            res["falsified"].append(
                f"{key}: gauss CV under t5 size {size_g_cv_under_t:.3f}"
            )
    # ---- (c) table ----
    for k_arms in K_GRID:
        row = {}
        for n in N_GRID:
            try:
                rng = np.random.default_rng(SEEDS[0] + k_arms * 31 + n)
                cv = float(np.quantile(max_t(rng, k_arms, n, "gauss", R_MC), 0.95))
            except RuntimeError:
                raise
            except Exception as exc:
                raise RuntimeError(f"table failed K={k_arms} n={n}") from exc
            row[f"n{n}"] = round(cv, 4)
        try:
            rngk = np.random.default_rng(SEEDS[1] + k_arms)
            row["known_scale"] = round(known_scale_cv(rngk, k_arms, R_MC), 4)
        except RuntimeError:
            raise
        res["table"][f"K{k_arms}"] = row
        print(f"TABLE K={k_arms}: {row}", flush=True)
    try:
        with open(OUT, "w") as f:
            json.dump(res, f, indent=1)
    except OSError as exc:
        raise RuntimeError(f"cannot write {OUT}") from exc
    print("falsified:", res["falsified"] or "none", flush=True)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    run()
