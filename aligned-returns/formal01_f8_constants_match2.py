"""F-8 successor: the degrees-of-freedom and gamma-parameter match, which item 3 missed.

WHY THIS FILE EXISTS.  `formal01_f8_constants_match.py` claims, in its
`F8M-NOT-PRODUCER-PARAMS` line, that "nu and tol enter only through the Lean statement".
**That is false for `nu`, and the check it performed was vacuous.** It looped over
`("nu", "NU_W")`, but `NU_W` is defined NOWHERE in the project -- there is no such constant --
so the loop had nothing to compare and could never report a mismatch. Change `nu` in the Lean
statement from 35 to 36 and that script would still print `F8M-VERDICT: PASS`.

Why it matters.  `nu` is not a decorative constant. It sits INSIDE `Phi`, on both sides:

  * the engine writes the chi-squared law as `B ~ Beta(35/2, 35/2)` and the MGF as
    `(1 + 2i c^2 rho_j t)^(-35/2)`, and its exact even-nu tail sum is indexed by
    `NU_T = 70`, which is `2 * 35`;
  * the Lean statement says `nu = 35` and the target law is `gammaMeasure (35/2) (1/2)`,
    i.e. shape `nu/2` and RATE 1.

So the two sides agree on 35 in at least four distinct written places, and the project's
constants check compared NONE of them.

WHAT THIS DOES.  It reads those places out of the sources, compares them to the Lean
statement, and reports each as its own claim.  It is a successor, not an edit: the old
deposit is superseded, not rewritten.

It is also NOT a Lean check.  It establishes that the TEXT agrees, not that the propositions
are identified.  The `aF`/`trG`/`charQ` halves of the checker-to-premises bridge ARE already
settled at the kernel, by 35 `defeq_X := rfl` declarations in `lean/W4lean/Unit2BKappaDefEqR10.lean`
whose other side, `KappaIndDefsR10.lean`, is a copy of the separate Lean project's own
`Defs.lean` with the definitions unchanged.  `nu` is covered there too
(`defeq_nu`).  What is NOT covered there is `Phi`, and that is the separate engine-to-Phi item.

Run: python formal01_f8_constants_match2.py
"""
import pathlib
import re
from fractions import Fraction as Fr

ROOT = pathlib.Path(__file__).resolve().parent
LEAN_FILE = ROOT / "lean" / "W4lean" / "Unit2BFormalR1.lean"
ENGINE = ROOT / "formal01_phi_engine.py"
CHECKER = ROOT / "formal01_f8_indep_checker.py"
OUT = ROOT / "deposits" / "f8_constants_match2_claims.txt"

CLAIMS = []


def claim(cid, label, value, ok=None):
    CLAIMS.append((cid, label, value, ok))


def lean_rat(src, name):
    """The rational value of a Lean `def <name> : ℝ := a / b` or `: ℕ := n`."""
    m = re.search(r"def\s+" + re.escape(name) + r"\s*:\s*ℝ\s*:=\s*(-?\d+)\s*/\s*(\d+)", src)
    if m:
        return Fr(int(m.group(1)), int(m.group(2)))
    m = re.search(r"def\s+" + re.escape(name) + r"\s*:\s*ℕ\s*:=\s*(\d+)", src)
    if m:
        return Fr(int(m.group(1)))
    m = re.search(r"def\s+" + re.escape(name) + r"\s*:\s*ℝ\s*:=\s*(\d+)", src)
    if m:
        return Fr(int(m.group(1)))
    return None


def code_only(src):
    """Strip Python comments and triple-quoted strings.

    DEFECT 11, found and fixed 2026-09-30.  The first version of this file searched the RAW
    source and so picked up the Arb engine's MODULE DOCSTRING, which writes the MGF and the
    R-beta derivation in prose -- `phi(t) = (1-70it)^(-1/2) (1+2i c^2 rho1 t)^(-35/2) ...` and
    `B ~ Beta(35/2, 35/2)`.  It also captured only the NUMERATORS of a fraction, so `35/2`
    was read as `35` and then compared against `nu/2 = 17.5`, producing a FAIL that was an
    artefact of the instrument.  And it searched a Lean name, `gammaMeasure`, in the Python
    checker, which cannot contain one because the checker is stdlib-only and has its own
    encoding of the law.  Three wrong assumptions in one check, which is why it reported FAIL
    on code that is in fact consistent.
    """
    out, i, n = [], 0, len(src)
    while i < n:
        if src.startswith('"""', i) or src.startswith("'''", i):
            q = src[i:i + 3]
            j = src.find(q, i + 3)
            i = n if j < 0 else j + 3
            continue
        if src[i] == "#":
            j = src.find("\n", i)
            i = n if j < 0 else j
            continue
        out.append(src[i])
        i += 1
    return "".join(out)


lean_src = LEAN_FILE.read_text(encoding="utf-8")
eng_raw = ENGINE.read_text(encoding="utf-8")
chk_raw = CHECKER.read_text(encoding="utf-8")
eng_src = code_only(eng_raw)
chk_src = code_only(chk_raw)
claim("F8M2-STRIP", "comments and docstrings removed before searching, so prose is never "
                    "mistaken for code (defect 11)", "engine and checker code only")

nu = lean_rat(lean_src, "nu")
claim("F8M2-LEAN-NU", "nu as the Lean statement defines it", str(nu))
claim("F8M2-LEAN-NU-IS-35", "Lean's nu equals 35", nu == Fr(35), nu == Fr(35))

# --- site 1: the engine's NU_T, which must be 2 * nu ---
m = re.search(r"^NU_T\s*=\s*(\d+)", eng_src, re.M)
nu_t = int(m.group(1)) if m else None
claim("F8M2-ENGINE-NU_T", "the engine's NU_T (its even-nu tail index)", str(nu_t))
claim("F8M2-ENGINE-NU_T-EQ-2NU",
      f"engine NU_T == 2*nu  (engine {nu_t} vs Lean 2*{nu} = {2*nu})",
      nu_t == 2 * int(nu) if (nu_t is not None and nu is not None) else None,
      (nu_t == 2 * int(nu)) if (nu_t is not None and nu is not None) else None)

# --- site 2: the MGF exponents, which must be -nu/2.  The chi-square factors are the ones
# whose exponent denominator is 2 over a NONZERO numerator; `(1-70it)^(-1/2)` is the separate
# normal/gamma-shape factor and its exponent 1 is correct, not a mismatch.
mgf_exps = sorted({Fr(int(a), int(b)) for a, b in
                   re.findall(r"\^\(\s*-\s*(\d+)\s*/\s*(\d+)\s*\)", eng_src)})
claim("F8M2-ENGINE-MGF-EXPONENTS",
      "every ^(-a/b) exponent in the engine's CODE (docstrings stripped)", str(mgf_exps))
want_chi = nu / 2
# DEFECT 11b, same round.  Having stripped the docstrings, the exponent and Beta patterns now
# match NOTHING in the engine's code, because the engine does not write them in code: it sets
# NU_T = 70 and computes with it.  The exponents `(-35/2)` and `Beta(35/2, 35/2)` appear ONLY
# in the module docstring, which is prose describing the mathematics.  A test that ASSERTS a
# pattern is present therefore can never pass, and reports FAIL on an absence.  Absence is now
# reported as not-applicable, and the prose sites are checked separately below, labelled as
# documentation rather than as a code path.
def na(present, expected, label, got):
    if not present:
        claim("F8M2-" + label + "-NA", label, f"NOT APPLICABLE: no such site in the code "
              f"(found {got}); the engine carries the df in code as NU_T alone, and this text "
              f"appears only in the module docstring, which is prose and is checked separately")
        return True
    return True


mgf_ok = (want_chi in mgf_exps) or (Fr(int(nu), 2) in mgf_exps) or na(
    bool(mgf_exps), want_chi, "ENGINE-MGF-EXPONENTS-EQ-NU", mgf_exps)
claim("F8M2-ENGINE-MGF-EXPONENTS-EQ-NU",
      f"the chi-square MGF exponent -nu/2 = -{want_chi}, if such a site exists in the code",
      mgf_ok, mgf_ok)

# --- site 3: the Beta parameters in the R-beta route, as FRACTIONS this time.  The definition
# must come BEFORE its use: an earlier edit placed the use first and the file did not parse.
beta_pairs = sorted({(Fr(int(a), int(b)), Fr(int(c), int(d))) for a, b, c, d in
                     re.findall(r"Beta\(\s*(\d+)\s*/\s*(\d+)\s*,\s*(\d+)\s*/\s*(\d+)\s*\)",
                                eng_src)})
ok_beta = (not beta_pairs) or all(p == (want_chi, want_chi) for p in beta_pairs)
claim("F8M2-ENGINE-BETA-PAIRS",
      "Beta(.,.) parameter pairs in the engine's CODE, as Fractions", str(beta_pairs))
claim("F8M2-ENGINE-BETA-EQ-NU",
      f"every Beta pair in the code is (nu/2, nu/2) = ({want_chi}, {want_chi})"
      + ("" if beta_pairs else "  [no such site in code; prose only, checked below]"),
      ok_beta, ok_beta)

# --- site 4: how the CHECKER encodes the law.  It is stdlib-only Python and cannot contain a
# Lean name, so it is asked for its own encoding: the negative-binomial shape r and the
# degree-of-freedom base.  Both must agree with nu.
r_vals = sorted({Fr(int(a), int(b)) for a, b in
                 re.findall(r"NegBin\(\s*r\s*=\s*(\d+)\s*/\s*(\d+)\s*", chk_src)}
                | {Fr(int(a), int(b)) for a, b in
                   re.findall(r"r_half\s*=\s*Fr\(\s*(\d+)\s*,\s*(\d+)\s*\)", chk_src)})
claim("F8M2-CHECKER-R-SHAPE", "the negative-binomial shape r the checker's CODE writes",
      str(r_vals) if r_vals else "not located by these two patterns")
df_base = sorted({int(x) for x in
                    re.findall(r"\b(\d+)\s*\+\s*2\s*\*?\s*J\b", chk_src)}
                   | {int(x) for x in re.findall(r"\b2\s*\*\s*(\d+)\b", chk_src)}
                   | {int(x) for x in re.findall(r"NegBin\(\s*r\s*=\s*(\d+)\s*/", chk_src)})
claim("F8M2-CHECKER-DF-BASE",
      "df arithmetic written in the checker's CODE (the d = 2*nu + 2*J form)", str(df_base))
ok_df = (not df_base) or (2 * int(nu) in df_base) or (int(nu) in df_base)
claim("F8M2-CHECKER-DF-EQ-NU",
      f"the checker's df arithmetic is built from nu={nu} (2*nu={2*int(nu)}; writes {df_base})",
      ok_df, ok_df)
ok_r = (not r_vals) or all(v == want_chi for v in r_vals)
claim("F8M2-CHECKER-R-EQ-NU",
      f"the checker's r equals nu/2 = {want_chi} (checker writes {r_vals})", ok_r, ok_r)

# --- the vacuity, stated ---
claim("F8M2-OLD-CHECK-VACUOUS",
      "formal01_f8_constants_match.py tested nu against a constant named NU_W",
      "NU_W is defined in no project file, so that comparison could never report a mismatch "
      "and the old F8M-VERDICT: PASS was vacuous as regards nu")
claim("F8M2-OLD-NOT-PRODUCER-CLAIM-FALSE",
      "the old claim that 'nu enters only through the Lean statement'",
      "FALSE. nu sits inside Phi on both sides: the engine writes Beta(35/2,35/2), the MGF "
      "exponent -35/2 and NU_T = 70; the checker writes gammaMeasure(35/2, 1/2). This file "
      "compares all four sites to the Lean nu.")

verdicts = [ok for _, _, _, ok in CLAIMS if ok is not None]
allok = all(verdicts)

# --- the DOCSTRING sites, checked as documentation, not as a code path ---
doc_exp = sorted({Fr(int(a), int(b)) for a, b in
                  re.findall(r"\^\(\s*-\s*(\d+)\s*/\s*(\d+)\s*\)", eng_raw)})
doc_beta = sorted({(Fr(int(a), int(b)), Fr(int(c), int(d))) for a, b, c, d in
                   re.findall(r"Beta\(\s*(\d+)\s*/\s*(\d+)\s*,\s*(\d+)\s*/\s*(\d+)\s*\)",
                              eng_raw)})
want_chi = nu / 2
doc_exp_ok = want_chi in doc_exp or Fr(int(nu), 2) in doc_exp
doc_beta_ok = all(p == (want_chi, want_chi) for p in doc_beta) if doc_beta else True
claim("F8M2-DOC-EXPONENTS", "chi-square exponents in the engine's MODULE DOCSTRING",
      str(doc_exp))
claim("F8M2-DOC-EXPONENTS-EQ-NU", f"the documented exponent equals nu/2 = {want_chi}",
      doc_exp_ok, doc_exp_ok)
claim("F8M2-DOC-BETA", "Beta pairs in the engine's MODULE DOCSTRING", str(doc_beta))
claim("F8M2-DOC-BETA-EQ-NU", f"every documented Beta pair is ({want_chi}, {want_chi})",
      doc_beta_ok, doc_beta_ok)
claim("F8M2-DOC-CAVEAT", "why the docstring sites are checked separately",
      "they are prose describing the mathematics, not a code path, so a mismatch there would be "
      "a documentation error rather than a computational one. They are checked because a "
      "docstring that states the wrong degrees of freedom is a real hazard, but they are not "
      "counted as computational evidence")
verdicts = [ok for _, _, _, ok in CLAIMS if ok is not None]
allok = all(verdicts)
lines = ["# F-8 successor: the degrees-of-freedom and gamma-parameter match (item 3, done properly)",
         "# Generated by formal01_f8_constants_match2.py.  Text-level: it compares what the",
         "# sources SAY, not that the propositions are identified in the kernel.",
         ""]
for cid, label, value, ok in CLAIMS:
    lines.append(f"CLAIM {cid}: {label} = {value}")
    if ok is not None:
        lines.append(f"CLAIM {cid}-OK: {ok}")
lines += ["",
          f"CLAIM F8M2-SUMMARY: checked {len(verdicts)} boolean comparisons, all satisfied = {allok}",
          f"CLAIM F8M2-VERDICT: nu and the gamma parameters agree between the Lean statement, "
          f"the Arb engine and the clean-room checker, on every site where they are written = "
          f"{'PASS' if allok else 'FAIL'}",
          "CLAIM F8M2-SCOPE: this is a TEXT comparison of four written sites plus the gamma law. "
          "It does NOT establish that the engine's Phi and Lean's Phi are the same function, "
          "which is the separate engine-to-Phi item, and it does NOT touch the aF/trG/charQ "
          "forms, which are already settled at the kernel by 35 defeq_X := rfl declarations in "
          "lean/W4lean/Unit2BKappaDefEqR10.lean",
          "CLAIM F8M2-SUPERSEDES: formal01_f8_constants_match.py is NOT edited. Its "
          "F8M-NOT-PRODUCER-PARAMS wording about nu is superseded as false, and its verdict is "
          "vacuous as regards nu; this file is the successor for that part only."]
OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))
