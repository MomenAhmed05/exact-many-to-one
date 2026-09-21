"""Check HUNT-STEIN.md against its deposits, and check the argument's shape.

Two kinds of check.

  numeric     every quantity the document quotes is formatted from a deposit
              and looked for in the text, the same discipline as the rest of
              the folder;
  structural  the claims the argument is NOT allowed to make.  The trap this
              document exists to avoid is quoting a weighted-average-power
              number under a minimax theorem, so the checks below fail if the
              two criteria are ever conflated, if condition (C4) stops being
              flagged as taken from the literature, or if the falsifier's
              sanity arm stops showing that the mesh can express a
              non-invariant rule.

    python verify_huntstein.py
"""

from __future__ import annotations

import io
import os
import re
import sys

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.join(HERE, "HUNT-STEIN.md")

ok = True
npass = nfail = 0
FAILED = []


def chk(label, cond, detail=""):
    global ok, npass, nfail
    cond = bool(cond)
    ok &= cond
    npass += cond
    nfail += not cond
    if not cond:
        FAILED.append(label)
    print(("  OK   " if cond else "  FAIL ") + label
          + ("   " + detail if detail else ""))


def main():
    if not os.path.exists(DOC):
        print("HUNT-STEIN.md not found")
        return 1
    TXT = re.sub(r"\s+", " ", io.open(DOC, encoding="utf-8").read())

    def has(label, value, fmt="%s"):
        v = fmt % value
        chk("%s -> %s" % (label, v), v in TXT)

    print("LEMMA 1, AT EVERY K")
    mk = F.read_deposit("w1_maximal_k.json")
    has("invariance drift", mk["worst_invariance_drift"], "%.2e")
    has("recovery error", mk["worst_recovery_error"], "%.2e")
    has("pairs", "{:,}".format(mk["pairs_per_design"]), "%s independently")
    has("max K", mk["max_K"], "K = 1, 2, 3, %d")
    chk("the dimension counts agree at every K", mk["all_counts_agree"])
    chk("the count is stated", "2K + 2` coordinates" in TXT
        or "2K + 2 coordinates" in TXT)
    chk("the element is fitted on the benchmark alone",
        "fitted on the benchmark pair alone" in TXT)

    print("\nTHE CRITERION")
    chk("the weighted-average criterion is excluded by name",
        "does not cover that criterion" in TXT)
    chk("the reason is given, not asserted",
        "orbit mean can sit below the value at the representative" in TXT)
    chk("the loss is shown invariant, which is the load-bearing step",
        "the loss is invariant" in TXT)
    chk("minimax regret is named as the criterion that is covered",
        "worst-case regret against the known-ratio" in TXT)

    print("\nLEMMA 2 AND THE CONDITIONS")
    chk("solvable is shown, not asserted",
        "abelian" in TXT and "solvable" in TXT)
    chk("the counterexample is named",
        "GL(p)" in TXT and "not amenable" in TXT)
    for c in ("(C1)", "(C2)", "(C3)", "(C4)"):
        chk("condition %s present" % c, c in TXT)
    chk("C4 is flagged as taken from the literature",
        "taken from the literature rather than verified from" in TXT)
    chk("the references are named",
        "Bondar and Milnes (1981)" in TXT and "Lehmann and Romano" in TXT)

    print("\nWHAT IS NOT CLAIMED")
    for s in ("Not weighted average power",
              "Not a claim about randomised versus deterministic",
              "Not a statement outside the declared domain",
              "Not an existence claim for a minimax test"):
        chk("disclaimer: %s" % s[:42], s in TXT)

    print("\nTHE MINIMAX FLOOR")
    try:
        mm = F.read_deposit("w1_minimax.json")
    except Exception:
        chk("deposit present: w1_minimax.json", False)
        mm = None
    if mm:
        has("gap low", mm["gap_above_floor_low_pp"], "%.2f")
        has("gap high", mm["gap_above_floor_high_pp"], "%.2f")
        chk("the two criteria are reported as different numbers",
            "it is a different number" in TXT)

    print("\nTHE FALSIFIER")
    try:
        nv = F.read_deposit("w1_noninvariant.json")
    except Exception:
        chk("deposit present: w1_noninvariant.json", False)
        nv = None
    if nv:
        # this deposit records a SUPERSEDED design, so the check is that it
        # says so, not that its numbers came out a particular way
        chk("the superseded falsifier is marked invalid",
            not nv.get("VALID", True))
        chk("it says why", bool(nv.get("invalid_because")))
        chk("its verdict cannot be read as evidence",
            "claim_survives" not in nv)
        has("superseded gain", nv.get("superseded_first_version_gain_pp", 0.0),
            "%.1f percentage points")
        chk("the first version's error is recorded, not hidden",
            "weighted average power at ONE value" in TXT
            or "weighted average power at one value" in TXT)
        chk("the objective is a supremum over the orbit",
            "objective is a supremum over the orbit" in TXT)
        chk("the mesh can express a non-invariant rule",
            nv["sanity_shows_the_mesh_can_cheat"])
        chk("the dependence of Xbar_0 on T is handled",
            "is **not** independent of `T`" in TXT
            or "is not independent of" in TXT)
        chk("the sanity arm is described",
            "sanity arm" in TXT)
        chk("both generators are probed",
            "`s_0^2` is the scale" in TXT and "`Xbar_0` the location" in TXT)

    print(chr(10) + "THE FALSIFIER IS UNBUILT, AND MUST SAY SO")
    try:
        rf = F.read_deposit("w1_noninv_refine.json")
    except Exception:
        rf = None
    if rf is not None:
        # write_deposit stores booleans as 0 and 1, so compare truthiness
        chk("the refinement run is marked invalid",
            not rf.get("VALID", True))
        chk("it says why", bool(rf.get("invalid_because")))
        chk("its shrinking flag cannot be read as evidence",
            "monotonically_shrinking" not in rf)
    print(chr(10) + "THE REPLACEMENT INSTRUMENT")
    try:
        hy = F.read_deposit("w1_hypotheses.json")
    except Exception:
        hy = None
    if hy is None:
        chk("the hypothesis check has a deposit", False)
    else:
        # the three hypotheses that are specific to this problem
        chk("Lemma 1 closes in exact algebra, not at a tolerance",
            bool(hy.get("lemma1_closed_symbolically")))
        chk("the algebra was closed at more than one K",
            len(hy.get("lemma1_symbolic", [])) >= 2)
        chk("the null and the alternative are invariant",
            bool(hy.get("H2_null_and_alternative_invariant")))
        chk("the loss is invariant", bool(hy.get("H3_loss_invariant")))
        # H3 is vacuous if beta* sits at an endpoint, so the draws have to
        # be interior for the check to have been capable of failing
        chk("the loss check was made where it could have failed",
            int(hy.get("hypotheses", {}).get("n_beta_star_interior", 0)) >= 40)
        # a rule that is conservative everywhere loses by construction and
        # probes nothing, so the family has to contain liberal shapes
        chk("the probe family can be liberal somewhere",
            int(hy.get("n_rules_liberal_somewhere", 0)) >= 12)
        chk("the probes were not optimised over a mesh",
            "no mesh sits between the two"
            in hy.get("note", "").lower())
        # the outcome is recorded either way, and the document has to agree
        # with it rather than with what was hoped for
        survives = bool(hy.get("conclusion_survives"))
        chk("the outcome is recorded", "conclusion_survives" in hy)
        chk("the document agrees with the measured outcome",
            ("Zero of 36 beat the floor" in TXT) == survives,
            "deposit says survives=%s" % survives)

    chk("the write-up says the instrument is built",
        "The instrument is built and the conclusion survives it" in TXT)
    chk("it still records that the linear program was abandoned",
        "ABANDONED after three" in TXT)

    print("\nTHE SWEEP, AND THE FIXED BASELINE")
    try:
        sr = F.read_deposit("w1_sigma_range.json")
    except Exception:
        sr = None
    if sr is None:
        chk("the sweep convergence run has a deposit", False)
    else:
        cv = sr.get("convergence", {})
        # the gap between two guesses MUST shrink; it is zero in exact
        # arithmetic, so a gap that does not shrink means widening is not the
        # fix and the instrument would need rethinking
        chk("the guess gap shrinks at every gamma",
            bool(cv) and all(c.get("shrank") for c in cv.values()))
        chk("and does so monotonically",
            bool(cv) and all(c.get("monotone") for c in cv.values()))
        # convergence is NOT claimed, because one arm did not reach it
        chk("non-convergence is recorded rather than smoothed",
            not sr.get("sup_has_converged", True))
        chk("the document says which arm is unquotable",
            "has not converged, so its regrets are directional" in TXT)
        chk("the document states which way the bias runs",
            "FLATTERS" in TXT)

    try:
        bf = F.read_deposit("w1_baseline_fixed.json")
    except Exception:
        bf = None
    if bf is None:
        chk("the fixed-baseline check has a deposit", False)
    else:
        # the rising losses are only a comparison if the baseline holds still,
        # and it should hold EXACTLY, not approximately
        chk("the baseline is exactly fixed across the sweeps",
            bool(bf.get("exactly_fixed")))
        chk("measured at more than two sweep widths",
            len(bf.get("rows", [])) >= 4)
        has("the fixed baseline", bf.get("worst_regret_pp", 0.0), "%.6f")

    print("\nREGISTER")
    chk("no em dash", "—" not in TXT)
    chk("no spelled-out per cent", "per cent" not in TXT.lower())
    for w in ("caveat", "the twist", "uncomfortable truth"):
        chk("no filler: %s" % w, w not in TXT.lower())

    print("\n%d checks passed, %d failed" % (npass, nfail))
    if FAILED:
        print("failures:")
        for f in FAILED[:25]:
            print("   " + f)
    print("HUNT-STEIN VERIFIED" if ok else "HUNT-STEIN HAS FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
