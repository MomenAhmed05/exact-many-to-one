"""High-precision independent check of the rate at a LOADING-DEPENDENT point.

It repeats an earlier 4e6-draw check under the two conditions that run did not meet:

  * evaluated at c = U (the earlier run used c0, which differs by ~6e-8);
  * with the standard error DECLARED IN ADVANCE (the earlier run computed it after the
    fact).

Why it exists. Every acceptance rule -- A-2..A-5, F-1, F-3, F-5 -- reuses the same rate
engine, so none of them can detect an error IN that engine. The independent evidence is a
direct Monte Carlo of the raw model's statistic T, which shares no code with the engine. The
earlier run left a blind spot of about 3e-4 in rate. This run is sized to close most of it.

PRE-DECLARATION, written before any draw is made and not revised afterwards
------------------------------------------------------------------------
  point            theta = (1/25, -3/2, 3/2)  (the argmax corner; b0,b1 both nonzero)
  critical value   c = U = 271740/100000
  estimator        binomial proportion of 1{T > c} under the raw aligned-factor model
  target SE        <= 1.5e-05  (the earlier run's SE was 7.8e-05)
  draws            2.5e8, in chunks of 1e6
  model            r_it = b_i*f_t + eps_it, i in {0,1}, t = 1..36,
                   f_t ~ N(0,1) iid, eps_0t ~ N(0,1) iid, eps_1t ~ N(0,lambda) iid,
                   all independent; the aligned (common-factor) null
  statistic        T = (ybar_1 - ybar_0) / sqrt((s1^2 + s0^2)/36), s^2 unbiased (ddof=1)
  decision rule    PASS iff |empirical - engine| <= 3 * declared SE target,
                   i.e. the two agree within 4.5e-5, and the engine value is <= alpha.
                   The declared target is the yardstick: it is not widened after seeing the
                   answer, and a FAIL here is a real finding, not a tolerance problem.
"""

from __future__ import annotations

import json
import pathlib
import time

import numpy as np

# ===================== PRE-DECLARED BEFORE ANY DRAW =========================
LAM, B0, B1, N_OBS = 1 / 25, -3 / 2, 3 / 2, 36
UC = 271740 / 100000
ALPHA = 1 / 40
CHUNK = 1_000_000
TOTAL_DRAWS = 250_000_000
SE_TARGET = 1.5e-5
DECISION_MULTIPLE = 3
SEED = 20260929
# ======================== the run ============================================


def main() -> int:
    rng = np.random.default_rng(SEED)
    hits = 0
    done = 0
    t0 = time.time()
    while done < TOTAL_DRAWS:
        k = min(CHUNK, TOTAL_DRAWS - done)
        f = rng.standard_normal((k, N_OBS))
        e0 = rng.standard_normal((k, N_OBS))
        e1 = np.sqrt(LAM) * rng.standard_normal((k, N_OBS))
        y0 = B0 * f + e0
        y1 = B1 * f + e1
        num = y1.mean(axis=1) - y0.mean(axis=1)
        s0 = y0.var(axis=1, ddof=1)
        s1 = y1.var(axis=1, ddof=1)
        stat = num / np.sqrt((s1 + s0) / N_OBS)
        hits += int((stat > UC).sum())
        done += k
        if (done // CHUNK) % 50 == 0:
            print(f"  {done:,} draws, {time.time() - t0:.0f}s", flush=True)
    p = hits / done
    se = float(np.sqrt(p * (1 - p) / done))
    # the engine's value at the same point, from the validated engine
    import importlib.util
    import sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    sp = importlib.util.spec_from_file_location("pe", "formal01_phi_engine.py")
    if sp is None or sp.loader is None:
        raise ImportError("engine")
    pe = importlib.util.module_from_spec(sp)
    sys.modules["pe"] = pe
    sp.loader.exec_module(pe)
    from fractions import Fraction as F
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import formal01_f8a_battery as bat
    R1, R2 = bat.r_enc((F(1, 25), F(-3, 2), F(3, 2)))
    eng_lo, eng_hi = pe.ball_bounds(pe.phi_r_beta(F(271740, 100000), R1.lo, R2.lo))
    eng = float(eng_lo)
    diff = p - eng
    tol = DECISION_MULTIPLE * SE_TARGET
    verdict = abs(diff) <= tol and eng <= ALPHA
    out = {
        "point": f"theta=({LAM},{B0},{B1})",
        "c": "Uc=271740/100000",
        "declared_se_target": SE_TARGET,
        "declared_draws": TOTAL_DRAWS,
        "draws": done,
        "seed": SEED,
        "empirical_P_T_gt_U": p,
        "achieved_se": se,
        "achieved_se_meets_target": se <= SE_TARGET,
        "engine_value": eng,
        "engine_ball": [str(eng_lo), str(eng_hi)],
        "difference": diff,
        "decision_tolerance": tol,
        "engine_below_alpha": eng <= ALPHA,
        "verdict": "PASS" if verdict else "FAIL",
        "seconds": round(time.time() - t0, 1),
    }
    (pathlib.Path(__file__).resolve().parent / "formal01_hp_result.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    lines = [f"CLAIM {k}: {k} = {v}" for k, v in out.items()]
    (pathlib.Path(__file__).resolve().parent
     / "deposits" / "formal01_highprec_claims.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
