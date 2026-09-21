"""P3 step 17: an INDEPENDENT high-precision cross-check of certificate points (plan section 3.5).

This is the third implementation class of the plan: it shares no code with the enclosure or with
`core.cond_noreject`.  For each sample point it evaluates

    F(c; lambda) = int prod_i Phi(s_i (z/sqrt(n0) + c d_i)) phi(z) dz

at high precision with `mpmath`, splitting the range at every KINK `z*_i = -c d_i sqrt(n0)` (where
the integrand is only piecewise smooth) and using a dense Gauss-Legendre rule on each smooth piece.
The resulting value is then checked for CONTAINMENT inside the Arb enclosure from `enclose_gl`.

It does not certify anything (the plan is explicit: the stochastic/numeric falsifier "does not
certify a result, but it catches shared conceptual bugs"); it is a falsifier.  A value outside the
enclosure is a bug and is reported as one.

Run:  python enclose_iv.py         (~1-2 min)
Out:  deposits/enclose_iv.json     (tier: region_sup_certified for the enclosure it checks)
"""

from __future__ import annotations

import importlib
import os
import time

import enclose_gl
import numpy as np
import paths
import tiers
from paths import NAN, _f, _i

OUT = os.path.join(paths.ensure_dirs(), "enclose_iv.json")
CORE = importlib.import_module("core")
ENV = importlib.import_module("env_rule")
draw_var = CORE.draw_var

POINTS = (
    ("desk_B_df35", [36.0, 13.0, 13.0], [1.0, 25.0], 2.6),
    ("desk_A_df3", [4.0, 13.0, 13.0], [0.04, 0.04], 2.4),
    ("sharp_A_df3_lo", [4.0, 13.0, 13.0], [0.04, 0.04], 5.0),
    ("K3_mixed", [8.0, 8.0, 10.0, 12.0], [0.25, 1.0, 4.0], 3.0),
    ("K5_balanced", [13.0] * 6, [1.0, 0.5, 0.25, 2.0, 4.0], 3.0),
)
DPS = (30, 40, 50)


def mpmath_F(c, d_row, n0, s_pt, dps=40, nodes=600, z_max=12.0):
    """High-precision F by kink-split adaptive quadrature (independent of the programme).

    `mp.quad` is used on each piece between consecutive kinks, so every piece is smooth.  (The first
    version called `mp.legendre(n)`, which is the Legendre POLYNOMIAL, not quadrature nodes - hence
    the NaN values.)
    """
    from mpmath import mp

    mp.dps = _i(dps, 40)
    c = mp.mpf(_f(c))
    sq = mp.sqrt(mp.mpf(_f(n0, 1.0)))
    s_mp = [mp.mpf(_f(v)) for v in np.asarray(s_pt, float)]
    d_mp = [mp.mpf(_f(v)) for v in np.asarray(d_row, float).reshape(-1)]
    kinks = sorted(
        _f(-_f(c) * _f(d_row[i]) * _f(np.sqrt(_f(n0, 1.0)))) for i in range(len(s_pt))
    )
    edges = (
        [-_f(z_max, 12.0)]
        + [k for k in kinks if -_f(z_max) < k < _f(z_max)]
        + [_f(z_max, 12.0)]
    )

    def integrand(z):
        acc = mp.mpf(1)
        for i in range(len(s_mp)):
            arg = s_mp[i] * (z / sq + c * d_mp[i])
            acc *= (1 + mp.erf(arg / mp.sqrt(2))) / 2
        return acc * mp.e ** (-(z * z) / 2) / mp.sqrt(2 * mp.pi)

    total = mp.mpf(0)
    for a, b in zip(edges[:-1], edges[1:], strict=False):
        total += mp.quad(integrand, [mp.mpf(a), mp.mpf(b)])
    return _f(total)


def point(name, n, lam, c, seed=20260918):
    n = np.asarray(n, float)
    lam = np.asarray(lam, float)
    s0, si, _lh = draw_var(n, lam, 1, seed)
    d_row = np.sqrt(si[0] / n[1:] + s0[0] / n[0])
    s_pt = np.sqrt(n[1:] / lam)
    lo, hi = enclose_gl.bounds_gl_arb(c, d_row, n[0], s_pt, s_pt)
    refs = []
    for dps in DPS:
        try:
            refs.append({"dps": _i(dps), "F": mpmath_F(c, d_row, n[0], s_pt, dps=dps)})
        except Exception as exc:  # noqa: BLE001
            refs.append({"dps": _i(dps), "F": NAN, "error": repr(exc)})
    vals = [_f(r["F"]) for r in refs if _f(r["F"]) == _f(r["F"])]
    conv = abs(vals[-1] - vals[0]) if len(vals) >= 2 else NAN
    inside = bool(vals and lo <= vals[-1] <= hi)
    return {
        "point": name,
        "n": list(n),
        "lam": list(lam),
        "c": _f(c),
        "K": _i(lam.size),
        "enclosure": [lo, hi],
        "width": _f(hi) - _f(lo),
        "mpmath_refs": refs,
        "convergence_abs": _f(conv),
        "contains": inside,
        "margin_lo": _f(vals[-1]) - _f(lo) if vals else NAN,
        "margin_hi": _f(hi) - _f(vals[-1]) if vals else NAN,
    }


def main():
    t0 = time.time()
    rows = [point(*spec) for spec in POINTS]
    ok = all(r["contains"] for r in rows)
    out = {
        "n_points": _i(len(rows)),
        "dps_set": [_i(d) for d in DPS],
        "all_contained": bool(ok),
        "max_width": _f(max(_f(r["width"]) for r in rows)),
        "max_convergence": _f(max(_f(r["convergence_abs"]) for r in rows)),
        "note": (
            "independent mpmath evaluation (kink-split dense Gauss-Legendre at 30/40/50 digits) "
            "checked for containment inside the Arb enclosure.  This is a FALSIFIER, not a "
            "certificate: it catches shared conceptual bugs and would fail loudly on one."
        ),
        "seconds": _f(time.time() - t0),
        "rows": rows,
    }
    for r in rows:
        print(
            f"{r['point']:<16} enclosure [{_f(r['enclosure'][0]):.12f}, "
            f"{_f(r['enclosure'][1]):.12f}] width {_f(r['width']):.2e} | mpmath "
            f"{_f(r['mpmath_refs'][-1]['F']):.12f} | contains {r['contains']} "
            f"| margins {_f(r['margin_lo']):+.2e}/{_f(r['margin_hi']):+.2e}"
        )
    print(f"\nall contained: {ok} | max convergence {out['max_convergence']:.2e}")
    tiers.write_deposit(
        out,
        OUT,
        "region_sup_certified",
        notes="independent mpmath cross-check of certificate points (falsifier)",
    )
    return out


if __name__ == "__main__":
    import sys

    try:
        res = main()
        sys.exit(0 if res.get("all_contained") else 1)
    except Exception as exc:  # noqa: BLE001
        print("IV CHECK ERROR:", repr(exc))
        sys.exit(2)
