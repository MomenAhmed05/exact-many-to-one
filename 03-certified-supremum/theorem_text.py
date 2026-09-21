"""P4 step 22: the required theorem statement and the model exclusions.

Plan section 4.4 requires every theorem-grade report to print the statement below and to list the
model exclusions explicitly.  This module is the single source for that text, so a report cannot
paraphrase it by accident.

    python theorem_text.py            # print the statement and the exclusions
"""
from __future__ import annotations

import paths
from paths import _f

STATEMENT = (
    "Under independent Gaussian samples across groups, with the stated sample sizes and the "
    "specified separate-variance Welch statistics, this decision uses internal level "
    "{a0:g} = {a:g} - {g:g} and a ratio-region miss budget {g:g}. The certified region coverage is "
    "at least 1 - {g:g}, and the reported threshold is a certified upper bound on the oracle "
    "critical value throughout that region. Therefore the familywise error rate is at most "
    "{a0:g} + {g:g} = {a:g}."
)

EXCLUSIONS = (
    "serial dependence",
    "cross-group dependence beyond the model",
    "factor exposure",
    "non-Gaussianity",
    "rolling re-estimation",
    "research-process multiplicity",
)

NUMBERS = {
    "abar": paths.ABAR,
    "gamma": paths.GAMMA,
    "alpha0": paths.ALPHA0,
    "tol": paths.TOL,
    "dp_target": paths.CERT_DP_TARGET,
}


def statement(abar=None, gamma=None):
    """The theorem statement with the levels substituted."""
    a = paths.ABAR if abar is None else _f(abar, paths.ABAR)
    g = paths.GAMMA if gamma is None else _f(gamma, paths.GAMMA)
    return STATEMENT.format(a=a, g=g, a0=a - g)


def exclusions():
    return list(EXCLUSIONS)


def full_text(abar=None, gamma=None):
    """Statement + exclusions, as printed by a theorem-grade report."""
    a = paths.ABAR if abar is None else _f(abar, paths.ABAR)
    g = paths.GAMMA if gamma is None else _f(gamma, paths.GAMMA)
    lines = ["THEOREM-GRADE STATEMENT", "", statement(a, g), "", "MODEL EXCLUSIONS", ""]
    lines += [f"  - {x}" for x in EXCLUSIONS]
    lines += ["", f"levels: abar = {a:g}, gamma = {g:g}, alpha_0 = {a - g:g}; "
                  f"tol = {paths.TOL:g}, Delta p target = {paths.CERT_DP_TARGET:g}"]
    return "\n".join(lines)


if __name__ == "__main__":
    print(full_text())
