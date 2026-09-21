"""Corrected implicit-function derivatives of the conditional critical value.

For the conditional law `F_cond(c; lambda; draw) = E_U[ prod_i Phi(s_i (U + c d_i)) ]` with
`s_i = sqrt(n_i/lambda_i)` and `U ~ N(0, 1/n_0)`, the implicit-function theorem gives

    d c* / d log lambda_i = - (dF/d log lambda_i) / (dF/dc)

and differentiating under the integral sign (writing `a_i = s_i (U + c d_i)`,
`P = prod_j Phi(a_j)`, `ratio_i = phi(a_i)/Phi(a_i)`):

    dF/dc        = E_U[ P * sum_i ratio_i * s_i d_i ]
    dF/dloglam_i = E_U[ P * ratio_i * (U + c d_i) * (-s_i/2) ]

**Why this module exists.**  The 02-bounded-rule implementation (`cert_sup.dcdloglam`) evaluated the
integrand on the 201-node grid `cert_sup._z_grid(m)` but paired it with the 200 weights
`cert_sup._weights(z)`, so its `np.dot` always raised `ValueError: shapes (201,) and (200,)
not aligned`; a bare `except` then returned `[]`.  Every row of
`02-bounded-rule/deposits/cert/rigorous.json` therefore carries an EMPTY `dcdloglam_at_oracle` column
(0 of 12 cells non-empty), and the Stage 1 report's derivative diagnostic was silently
empty.  No published number is affected (no 02-bounded-rule report quotes a derivative value), but the
diagnostic existed only in name.  This is recorded as verifier-known defect #11.

This implementation uses the probability-scale composite Gauss-Legendre nodes of
`02-bounded-rule/h12_quad.normal_nodes`, the same quadrature the exact-law rejection probabilities use,
and returns plain lists.  `deriv_diag.py` checks it against central finite differences of
the conditional critical value solved by bisection on the independent solver.
"""

from __future__ import annotations

import numpy as np
import paths  # noqa: F401  (wires 02-bounded-rule/ onto sys.path)


def normal_nodes(npanels=40, order=10):
    """Probability-scale composite Gauss-Legendre nodes for `int f dPhi`."""
    from h12_quad import normal_nodes as _nn

    return _nn(npanels, order)


def dcdloglam(n, d_row, lam, c, npanels=40, order=10):
    """[d c*/d log lambda_i] at a draw; empty list only on a genuine failure."""
    try:
        z, w = normal_nodes(npanels, order)
        n = np.asarray(n, float)
        lam = np.asarray(lam, float)
        d_row = np.asarray(d_row, float).reshape(-1)
        u = z / np.sqrt(n[0])
        P = np.ones(z.size)
        terms_dc = np.zeros(z.size)
        terms_dlam = np.zeros((lam.size, z.size))
        for i in range(lam.size):
            s = np.sqrt(n[i + 1] / lam[i])
            a = s * (u + c * d_row[i])
            phi = np.exp(-0.5 * a * a) / np.sqrt(2.0 * np.pi)
            ph = _ndtr(a)
            ratio = np.where(ph > 0.0, phi / np.maximum(ph, 1e-300), 0.0)
            P = P * ph
            terms_dc += ratio * s * d_row[i]
            terms_dlam[i] = ratio * (u + c * d_row[i]) * (-s / 2.0)
        dFdc = float(np.dot(P * terms_dc, w))
        if not np.isfinite(dFdc) or abs(dFdc) < 1e-300:
            return []
        return [float(-np.dot(P * terms_dlam[i], w) / dFdc) for i in range(lam.size)]
    except Exception:  # noqa: BLE001
        return []


def _ndtr(a):
    from scipy.special import ndtr

    return ndtr(a)


def w26_version_defect(n=None, d_row=None, lam=None, c=2.6):
    """Evidence for defect #11: the 02-bounded-rule call returns [] where this one does not."""
    import importlib

    out: dict[str, object] = {
        "w26_empty": None,
        "fixed_nonempty": None,
        "w26_error_text": None,
    }
    try:
        cert = importlib.import_module("cert_sup")
        n = np.array([13.0, 13.0, 13.0]) if n is None else np.asarray(n, float)
        lam = np.array([1.0, 25.0]) if lam is None else np.asarray(lam, float)
        d = np.array([0.32, 1.69]) if d_row is None else np.asarray(d_row, float)
        old = cert.dcdloglam(n, d, lam, c)
        out["w26_empty"] = len(old) == 0
        fixed = dcdloglam(n, d, lam, c)
        out["fixed_nonempty"] = len(fixed) > 0
        out["w26_error_text"] = (
            "np.dot shapes (201,) vs (200,) -> ValueError caught by a "
            "bare except; see the module docstring"
        )
        out["fixed_values"] = [round(v, 6) for v in fixed]
    except Exception as exc:  # noqa: BLE001
        out["w26_error_text"] = repr(exc)
    return out


if __name__ == "__main__":
    import json

    print(json.dumps(w26_version_defect(), indent=2))
