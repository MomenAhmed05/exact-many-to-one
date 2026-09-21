"""Shared paths, constants and guarded converters for the 03-certified-supremum programme.

`02-bounded-rule/` is read-only input; this module wires it onto `sys.path` once and pins the
programme-wide constants declared before any result was produced:

    ABAR    = 0.025     hard one-sided FWER ceiling (desk default)
    GAMMA   = 0.0005    ratio-region miss budget
    ALPHA0  = ABAR - GAMMA = 0.0245   internal level used by the certified tier
    TOL     = 1e-4      predeclared absolute target on c_U - c_L
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
W26 = os.path.join(os.path.dirname(HERE), "02-bounded-rule")
EXACT_LAW = os.path.join(os.path.dirname(HERE), "01-exact-law")
DEPOSITS = os.path.join(HERE, "deposits")
VERIFY = os.path.join(HERE, "verify")
RELEASE = os.path.join(HERE, "release")

for _p in (W26, HERE, EXACT_LAW):
    if _p not in sys.path:
        sys.path.insert(0, _p)

ABAR = 0.025
GAMMA = 0.0005
ALPHA0 = ABAR - GAMMA
TOL = 1e-4


CERT_DP_TARGET = 1e-6  # predeclared Delta p target at certificate points
CERT_MIN_PER_FAMILY = 60.0  # predeclared runtime budget, seconds per comparison family


def _f(x, default=None):
    try:
        return float(x)
    except Exception:  # noqa: BLE001
        return float("nan") if default is None else default


def _i(x, default=0):
    try:
        return int(x)
    except Exception:  # noqa: BLE001
        return default


def _env_int(name, default):
    try:
        return int(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return int(default)


def _env_float(name, default):
    try:
        return float(os.environ.get(name) or default)
    except Exception:  # noqa: BLE001
        return float(default)


ERRLOG: list[str] = []


def ensure_dirs(*sub):
    """Create deposit subdirectories on demand and return the path.

    A failure is recorded in `ERRLOG` rather than raised (and rather than silently
    swallowed): the caller can keep working and the run reports the problem.
    """
    p = os.path.join(DEPOSITS, *sub) if sub else DEPOSITS
    try:
        os.makedirs(p, exist_ok=True)
    except OSError as exc:
        ERRLOG.append(f"{p}: {exc!r}")
    return p


NAN = _f("nan")

# Reference-solver resolution for the falsifiers.  Finding #13: the 02-bounded-rule solver
# core.cond_noreject with its DEFAULT 40 panels is inaccurate by ~2e-4 at sharp
# configurations (large s = sqrt(n/lambda)), where the integrand is a step of width
# 1/s in z; scipy.quad and npanels=4000 agree with each other and disagree with the
# default.  The falsifiers therefore use a refined reference.
REF_PANELS = _env_int("WA2H_REF_PANELS", 8000)
