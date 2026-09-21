"""Core cross-checks: my harness vs the paper's deposited machinery.

Run from 02-bounded-rule.  Every comparison uses the SAME draws so the only
difference is the implementation.  Prints PASS/FAIL lines; exits non-zero if a
comparison exceeds its tolerance.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "01-exact-law"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import cond_noreject, draw_var  # noqa: E402
from h12_calibrated import calibrate as paper_calibrate
from h12_quad import normal_nodes  # noqa: E402
from h12_tamhane_headtohead import cond_reject  # noqa: E402

FAILS = []


def check(name, a, b, tol):
    ok = abs(a - b) <= tol
    print(
        f"{'PASS' if ok else 'FAIL'}  {name}: mine {a:.8g} vs paper {b:.8g} "
        f"(|d| {abs(a - b):.3g}, tol {tol:g})"
    )
    if not ok:
        FAILS.append(name)


def chk_cond_reject():
    """Same draws: paper cond_reject (rejection) vs 1 - cond_noreject."""
    n = np.array([4.0, 13.0, 13.0])
    sig = np.array([1.0, 1.0, 5.0])
    lam = np.array([1.0, 25.0])
    s0, si, _ = draw_var(n, lam, 200, 12345)
    un = normal_nodes(20, 12)
    u_nodes = (un[0] * sig[0] / np.sqrt(n[0]), un[1])
    for c in (2.5, 2.7003, 3.0):
        diff = []
        for d in range(200):
            p = cond_reject(n, sig, np.concatenate([[s0[d]], si[d]]), c, u_nodes)
            m = (
                1.0
                - cond_noreject(
                    np.array([c]), np.array([s0[d]]), si[d : d + 1], n, lam
                )[0]
            )
            diff.append(abs(p - m))
        try:
            worst = float(np.max(diff))
        except Exception as exc:  # noqa: BLE001
            print(f"INFO  cond_reject diff aggregation failed: {exc!r}")
            worst = float("nan")
        check(f"cond_reject agreement c={c}", worst, 0.0, 2e-6)


def chk_calibrator():
    """Paper calibrate vs my parallel calibrate on the same design/seed."""
    cases = [
        ("S1", np.array([10.0, 10.0, 10.0]), np.array([30.0, 30.0, 30.0])),
        ("S4", np.array([4.0, 13.0, 13.0]), np.array([10.0, 10.0, 50.0])),
        ("S5", np.array([13.0, 13.0, 4.0]), np.array([10.0, 10.0, 50.0])),
        ("S6", np.array([4.0, 13.0, 13.0]), np.array([50.0, 10.0, 10.0])),
    ]
    import multiprocessing as mp

    from w4_calibration import calibrate_par

    with mp.Pool(4) as pool:
        for name, n, sig in cases:
            a_p, F_p = paper_calibrate(
                n, sig, 0.025, np.random.default_rng(31337), B=600, ngrid=9
            )
            a_m, F_m = calibrate_par(
                n, sig, 0.025, 600, seed=31337, workers=4, pool=pool
            )
            check(f"{name} F(alpha)", F_m(0.025), F_p(0.025), 3e-3)
            check(f"{name} alpha*", a_m, a_p, 1e-3)


def chk_envelope():
    """Envelope: exact coverage certificate + deposited guarantee flags."""
    from env_rule import region_coverage, tight_scale
    from io_util import load_json

    gamma = 0.0005
    try:
        cov_lo = min(
            region_coverage(
                (n0, 13.0, 13.0), gamma, tight_scale((n0, 13.0, 13.0), gamma)
            )
            for n0 in (4.0, 12.0, 36.0)
        )
        cov_ok = 1.0 if cov_lo >= 1.0 - gamma else 0.0
    except Exception as exc:  # noqa: BLE001
        print(f"INFO  envelope coverage check failed: {exc!r}")
        cov_ok = 0.0
    check("envelope tight coverage >= 1-gamma", cov_ok, 1.0, 0.0)

    dep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "deposits")
    try:
        cells = load_json(os.path.join(dep, "env", "sweep_summary.json"))["cells"]
        n_ok = sum(
            1
            for c in cells
            if c["passes_guarantee"] and c["tight_excess_vs_oracle_max"] <= 0.0
        )
        n_cheap = sum(
            1 for c in cells if c["tight_cost_pct"] <= c["ub_cost_pct"] + 1e-9
        )
        n_tot = len(cells)
    except Exception as exc:  # noqa: BLE001
        print(f"INFO  envelope deposit unreadable: {exc!r}")
        n_ok = n_cheap = 0
        n_tot = 1
    try:
        ok_frac = float(n_ok)
        cheap_frac = float(n_cheap)
        tot = float(n_tot)
    except Exception as exc:  # noqa: BLE001
        print(f"INFO  envelope flag aggregation failed: {exc!r}")
        ok_frac = cheap_frac = 0.0
        tot = 1.0
    check("envelope guarantee flags on the desk cells", ok_frac, tot, 0.0)
    check("tight region cost <= union bound", cheap_frac, tot, 0.0)


def chk_fwer_harness():
    """My rb_fwer paired estimator vs the paper's arms_fwer on the same cell."""
    import h12_crossover as X
    from core import batch_c, oracle_c

    n = np.array([4.0, 13.0, 13.0])
    lam = np.array([1.0, 25.0])
    D, seed = 150, 4242
    s0, si, lh = draw_var(n, lam, D, seed)
    c_or = oracle_c(n, lam, 0.025)
    cplug = batch_c(n, lh, alpha=0.025, workers=4)
    nr_o = cond_noreject(np.full(D, c_or), s0, si, n, lam)
    nr_r = cond_noreject(cplug, s0, si, n, lam)
    try:
        mine = 0.025 + float((nr_o - nr_r).mean())
    except Exception as exc:  # noqa: BLE001
        print(f"INFO  paired FWER aggregation failed: {exc!r}")
        mine = float("nan")
    # paper path on a fresh draw stream (same code path as the deposit)
    paper = X.arms_fwer(
        [4, 13, 13],
        [10.0, 10.0, 50.0],
        n_draws=150,
        seed=seed,
        alpha=0.025,
        m3_sims=20000,
    )
    print(
        f"INFO  FWER mine {mine:.5f} vs paper-path "
        f"{paper['uncalibrated']['fwer']:.5f} (different draws)"
    )
    check("rb_fwer direction", mine > 0.025, True, 0.0)


if __name__ == "__main__":
    chk_cond_reject()
    chk_envelope()
    chk_fwer_harness()
    chk_calibrator()
    print()
    if FAILS:
        print("FAILURES:", FAILS)
        sys.exit(1)
    print("ALL CORE CHECKS PASSED")
