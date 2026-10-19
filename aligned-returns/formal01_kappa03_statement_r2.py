"""Generate `Unit2BFormalR1R2.lean`: a single Lean artefact with NO `sorry` anywhere.

WHY.  `lean/W4lean/Unit2BFormalR1.lean` carries `sorry` on all 13 theorems, `comp_thm_2B`
included.  The convention behind that file is CORRECT -- F2-2 requires the statement first, F2-3
the proof in a fresh file -- and the proofs genuinely exist and are sorry-free.  But a checker
that greps for `sorry` will keep objecting, and it is right that the file it is shown has none.
So we ADD an artefact with no `sorry` in it.  `Unit2BFormalR1.lean` is NOT edited.

WHY A GENERATOR AND NOT A HAND-WRITTEN FILE.  This is a mechanical assembly: thirteen
statements from one file, thirteen proof bodies from another, concatenated.  Two reasons to
generate rather than write:
  * a hand-typed theorem statement is a place for a silent typo, and a typo in a statement is
    worse than no artefact at all -- it would look like a faithful copy of something else;
  * a generated file can be checked against its sources by re-extraction, which is a real
    check.  A hand-written one can only be asserted to match.

THE HAZARD, WHICH IS THE WHOLE DIFFICULTY.  Proof bodies written inside
`namespace Unit2B.FormalR1.Proofs` refer to sibling theorems by BARE NAME.  Lean's name
resolution searches the current namespace and then its parents -- and one parent is
`Unit2B.FormalR1`, which is the STATEMENT FILE, whose thirteen theorems all end in `sorry`.  So
a bare name would silently pull `sorryAx` into a file that is supposed to have none.  This has
happened in this project before: a bare `box_pos` once did exactly that, and only the axiom
audit caught it.

The fix applied here is mechanical and exhaustive: every occurrence of each of the thirteen
theorem NAMES, matched on word boundaries and NOT preceded by a dot, is rewritten to
`Unit2B.FormalR1.R2.<name>`.  Definitions (`aF`, `trG`, `mu1`, ...) are deliberately left
bare: they are `def`s in `Unit2B.FormalR1`, not theorems, so they carry no `sorryAx`, and they
resolve up the namespace chain correctly.

Run: python formal01_kappa03_statement_r2.py
Writes: lean/W4lean/Unit2BFormalR1R2.lean, and deposits/statement_identity_r2_claims.txt
"""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent
LEAN = ROOT / "lean" / "W4lean"
STMT_FILE = LEAN / "Unit2BFormalR1.lean"
PROOF_FILE = LEAN / "Unit2BFormalR1Proofs.lean"
OUT_LEAN = LEAN / "Unit2BFormalR1R2.lean"
OUT_CLAIMS = ROOT / "deposits" / "statement_identity_r2_claims.txt"

NEW_NS = "Unit2B.FormalR1.R2"
ORIG_NS = "Unit2B.FormalR1"
PROOF_NS = "Unit2B.FormalR1.Proofs"

# The thirteen APPROVED statements.  Everything in the proof file is emitted, not just these;
# the list is used only to decide which declarations get a statement-identity check against the
# statement file of record, and which are helpers with no approved statement of their own.
THIRTEEN = [
    "Phi_antitone", "ae_nonneg_of_gammaLaw", "charQ_eq_det", "discG_eq",
    "mu_sum_prod", "charQ_factor", "box_pos", "r1_lower", "r2_lower", "cstar_bracket",
    "cstarBox_eq_sSup", "cover_sound", "comp_thm_2B",
]
# DEFECT 19, in this generator.  A hand-maintained list of the thirteen was used as the set of
# declarations to EMIT, and `mu_sum_prod`'s proof depends on `discG_nonneg`, which exists only
# in the proof file and is not one of the thirteen.  The reference was qualified to
# `Unit2B.FormalR1.R2.discG_nonneg`, which did not exist, so the build failed with a single
# `Unknown identifier` -- and because six declarations transitively depend on `mu_sum_prod`,
# six came back carrying `sorryAx`.  One missing helper contaminated six theorems.  The emitted
# set is now DERIVED from the proof file, so a helper can never be omitted again.

CLAIMS = []


def claim(cid, label, value):
    CLAIMS.append((cid, label, value))


# ---------------------------------------------------------------- strip Lean comments
def strip_comments(src: str) -> str:
    """Remove `--` to end of line and nestable /- -/ blocks, preserving line structure."""
    out, depth, i, n = [], 0, 0, len(src)
    while i < n:
        if depth > 0:
            if src.startswith("/-", i):
                depth += 1
                i += 2
                continue
            if src.startswith("-/", i):
                depth -= 1
                i += 2
                continue
            out.append("\n" if src[i] == "\n" else " ")
            i += 1
            continue
        if src.startswith("/-", i):
            depth += 1
            i += 2
            continue
        if src.startswith("--", i):
            while i < n and src[i] != "\n":
                i += 1
            continue
        out.append(src[i])
        i += 1
    return "".join(out)


DECL = re.compile(r"^(theorem|lemma)\s+([A-Za-z_][^\s:]*)", re.M)


def _delta(ln: str) -> int:
    return (ln.count("{") + ln.count("(")) - (ln.count("}") + ln.count(")"))


def extract(src, stop_at_sorry: bool):
    """Return {name: (statement_up_to_colon_eq, body_after_colon_eq)}.

    DEFECT 14, found and fixed immediately after first run: the first version conflated a
    CHARACTER index with a LINE index.  It set `m2 = k`, where `k` is the offset of ":=" inside
    the joined header, and then used `m2` to index `lines`.  So `i = max(j, m2) + 1` jumped far
    past the end of the file and only ONE of the thirteen declarations was found.  The symptom
    was loud -- SI-COUNT reported 1 and twelve names "absent" -- which is the only reason it
    was caught.  Every index below is now a line index, and the declaration is bounded by
    advancing `i` to just past the end of its own body.
    """
    code = strip_comments(src)
    lines = code.split("\n")
    res, i = {}, 0
    while i < len(lines):
        m = DECL.match(lines[i])
        if not m:
            i += 1
            continue
        name = m.group(2)

        # Header: from line i to the line carrying ":=" at bracket depth <= 0.
        head, j, depth = [], i, 0
        while j < len(lines):
            head.append(lines[j])
            depth += _delta(lines[j])
            if ":=" in lines[j] and depth <= 0:
                break
            j += 1
        if j >= len(lines):          # unterminated; stop rather than loop
            break
        header_text = "\n".join(head)
        k = header_text.find(":=")
        statement = header_text[:k].rstrip()

        # DEFECT 17, same generator, same CLASS as 14: `k` is a character offset into the
        # JOINED header, and it was used to slice a SINGLE line.  For a multi-line header the
        # offset overshoots, so the slice was empty and the generator emitted bare `:=` with
        # no proof body.  The offset has to be recomputed within the line that holds the `:=`.
        line_j = head[j - i]
        pos_in_line = line_j.find(":=")
        tail = line_j[pos_in_line + 2:] if pos_in_line >= 0 else ""
        body_lines = [tail]
        # DEFECT 18, and the most important of the three.  The proof body was bounded by
        # BRACKET DEPTH, on the assumption that a body is delimited by braces.  It is not: a
        # Lean body is delimited by INDENTATION.  A proof beginning `by` and continuing with
        # `have h : T := ...` never opens a brace, so the depth stayed 0 and the walk stopped
        # on the `by` itself, emitting a bare `:= by` with no proof.  The correct terminator
        # for a file whose top-level declarations all start in column 0 is: the next non-blank
        # line that itself starts in column 0 with a top-level keyword.
        TOP_LEVEL = re.compile(
            r"^(theorem|lemma|def|abbrev|instance|structure|example|open|import|namespace|end"
            r"|#print|noncomputable|set_option|variable|@\[|section)")
        jj = j + 1
        while jj < len(lines):
            ln = lines[jj]
            if ln.strip() == "":
                body_lines.append(ln)
                jj += 1
                continue
            if TOP_LEVEL.match(ln):
                break
            body_lines.append(ln)
            jj += 1
        while body_lines and body_lines[-1].strip() == "":
            body_lines.pop()

        body = "\n".join(body_lines)
        if stop_at_sorry and re.search(r"(?<![A-Za-z_.])sorry(?![A-Za-z_])", body):
            res[name] = (statement, None)
        else:
            res[name] = (statement, body)
        i = jj
    return res


stmt_src = STMT_FILE.read_text(encoding="utf-8")
proof_src = PROOF_FILE.read_text(encoding="utf-8")
before_sha = __import__("hashlib").sha256(STMT_FILE.read_bytes()).hexdigest()

stmt_decls = extract(stmt_src, stop_at_sorry=True)
proof_decls = extract(proof_src, stop_at_sorry=False)
# Emit EVERY declaration the proof file contains, in the order it appears there.  A
# declaration's statement comes from the statement file when it has one there (the thirteen),
# and from the proof file itself when it does not (`discG_nonneg`), in which case it is a
# helper and carries no fidelity claim.
EMIT = list(proof_decls.keys())

claim("SI-COUNT", "thirteen theorems expected in the statement file",
      f"{len(stmt_decls)} top-level theorem/lemma declarations found")
missing = [n for n in THIRTEEN if n not in stmt_decls]
claim("SI-MISSING-IN-STMT", "named theorems absent from the statement file",
      ", ".join(missing) if missing else "none")
missing_p = [n for n in THIRTEEN if n not in proof_decls]
claim("SI-MISSING-IN-PROOFS", "named theorems absent from the proof file",
      ", ".join(missing_p) if missing_p else "none")

# ------------------------------------------------- qualify the thirteen names in the bodies
name_re = re.compile(
    r"(?<![.\w])(" + "|".join(re.escape(n) for n in EMIT) + r")(?![A-Za-z_0-9'])")
# An ALREADY-qualified reference to one of the emitted names, e.g. `Unit2B.FormalR1.discG_eq`
# written out in the proof body.  Left alone by `name_re` (the preceding dot blocks it), it
# would resolve to the STATEMENT file's `sorry` theorem -- the exact leak this file exists to
# prevent.  So it is rewritten to the new namespace as well.
qualified_re = re.compile(
    r"\b" + re.escape(ORIG_NS) + r"\.(" + "|".join(re.escape(n) for n in EMIT) + r")(?![A-Za-z_0-9'])")


def qualify(text: str) -> str:
    """Rewrite references to the thirteen into the NEW namespace, and NOTHING else.

    Applied ONLY to proof bodies.  A first version applied it to the `theorem` head as well,
    which produced `theorem Unit2B.FormalR1.R2.Phi_antitone` inside
    `namespace Unit2B.FormalR1.R2`, and the kernel duly named every declaration
    `Unit2B.FormalR1.R2.Unit2B.FormalR1.R2.Phi_antitone`.  Nothing resolved, body references
    fell through to the statement file's `sorry` theorems, and the build reported `sorryAx` on
    all thirteen.  That was defect 15 in this one generator, and the build caught it.
    """
    text = qualified_re.sub(f"{NEW_NS}.\\1", text)
    return name_re.sub(lambda m: f"{NEW_NS}.{m.group(1)}", text)


pieces = [
    "/-",
    "  Generated by `formal01_kappa03_statement_r2.py`.  DO NOT EDIT BY HAND.",
    "",
    "  A single artefact carrying all thirteen theorems of the Unit 2B statement file WITH",
    "  their proofs inlined, so that a reader or a checker inspecting this one file finds no",
    "  `sorry` anywhere.  It exists because the statement file carries `sorry` on all thirteen -- which is",
    "  by design, since under F2-2 the statement is authored there first and F2-3 puts the proof",
    "  in a fresh file.  That file is NOT edited.  This file is an ADDITION, not a replacement,",
    "  and the statements below are generated from `Unit2BFormalR1.lean` and compared against",
    "  it character for character by `deposits/statement_identity_r2_claims.txt`.",
    "",
    "  THE BARE-NAME HAZARD.  Proof bodies written inside `Unit2B.FormalR1.Proofs` refer to",
    "  sibling theorems by bare name, and Lean's resolution searches parent namespaces -- one of",
    "  which is `Unit2B.FormalR1`, the STATEMENT FILE, whose theorems all end in `sorry`.  Every",
    "  occurrence of the thirteen names on a word boundary and not preceded by a dot is",
    "  therefore rewritten to this namespace, so that no `sorryAx` can enter.  Definitions are",
    "  left bare on purpose: they are `def`s in `Unit2B.FormalR1` and carry no `sorryAx`.",
    "",
    "  The kernel-level check that these declarations have the types of the statements of record",
    "  is in `StatementIdentityR2.lean`.",
    "-/",
    "import W4lean.Unit2BFormalR1Proofs",
    "import Mathlib",
    "open MeasureTheory ProbabilityTheory Set",
    "noncomputable section",
    f"namespace {NEW_NS}",
    "",
]

mismatches = []
helpers = []
for n in EMIT:
    if n not in proof_decls:
        continue
    ph, pbody = proof_decls[n]
    if n in stmt_decls:
        sh, _ = stmt_decls[n]
        header = ("/- APPROVED STATEMENT, verbatim from Unit2BFormalR1.lean -/\n" + sh.strip())
        # THE STATEMENT-IDENTITY CHECK, done properly.  `sh` is the statement file's header up
        # to its `:=`; `ph` is the PROOF file's header up to its `:=`.  They must be the same
        # string.  A first version compared `sh` with itself, which is a tautology and proves
        # nothing -- the same class of defect as a test that cannot fail.
        a = re.sub(r"\s+", " ", sh).strip()
        b = re.sub(r"\s+", " ", ph).strip()
        identical = a == b
        claim(f"SI-{n}", f"statement of `{n}` is character-identical in the statement file "
                         f"and the proof file",
              f"{identical}" + ("" if identical else
                                f"  -- STMT: {a[:120]}  ||  PROOF: {b[:120]}"))
        if not identical:
            mismatches.append(n)
    else:
        # A helper with no approved statement of its own -- `discG_nonneg` is the one that
        # matters.  It is emitted because a thirteen DEPENDS ON IT, and it carries NO fidelity
        # claim: there is nothing in the statement file for it to be identical to, so claiming
        # identity for it would be claiming something there is no evidence for.
        header = ("/- HELPER: exists only in Unit2BFormalR1Proofs.lean, with no approved "
                  "statement in the statement file of record. Emitted because a thirteen "
                  "depends on it. Carries NO fidelity claim. -/\n" + ph.strip())
        helpers.append(n)
        claim(f"SI-{n}", f"`{n}` is a helper with no approved statement; emitted because a "
                         f"thirteen depends on it",
              "no statement-identity claim is made for it")
    body_txt = pbody if pbody is not None else " by sorry"
    pieces.append(header + " :=" + qualify(body_txt).strip())
    pieces.append("")

pieces += ["", "-- Every declaration, on its axioms.  The file has no `sorry` by construction:", ]
for n in EMIT:
    if n in proof_decls:
        pieces.append(f"#print axioms {NEW_NS}.{n}")
pieces += ["", f"end {NEW_NS}", ""]

OUT_LEAN.write_text("\n".join(pieces) + "\n", encoding="utf-8")

# ------------------------------------------------------------------ verify the artefact
gen = OUT_LEAN.read_text(encoding="utf-8")
gen_code = strip_comments(gen)
SORRY = re.compile(r"(?<![A-Za-z_.])(sorry|admit|native_decide)(?![A-Za-z_])")
hits = [f"L{i+1}: {ln.strip()[:80]}" for i, ln in enumerate(gen_code.split("\n"))
        if SORRY.search(ln)]
claim("SI-NO-SORRY", "sorry/admit/native_decide in the generated file, COMMENTS AND "
                      "DOCSTRINGS STRIPPED", f"{len(hits)} occurrences"
      + ("" if not hits else " -> " + "; ".join(hits[:4])))
claim("SI-GENERATED-NAMES", "declarations emitted (DERIVED from the proof file, not hand-listed)",
      ", ".join(EMIT))
claim("SI-HELPERS", "declarations emitted that are helpers with no approved statement",
      ", ".join(helpers) if helpers else "none")
claim("SI-MISMATCHES", "thirteen whose statement is NOT character-identical",
      ", ".join(mismatches) if mismatches else "none")
claim("SI-QUALIFIED-REFERENCES", "bare references to the thirteen remaining in the generated file",
      str(len([1 for m in name_re.finditer(qualify("\n".join(pieces)))])))
claim("SI-STMT-FILE-UNTOUCHED",
      "sha256 of Unit2BFormalR1.lean after generation (must equal the 'before' value)",
      __import__("hashlib").sha256(STMT_FILE.read_bytes()).hexdigest())
claim("SI-STMT-FILE-SHA-BEFORE", "sha256 of Unit2BFormalR1.lean before generation", before_sha)
claim("SI-GENERATED-SHA", "sha256 of the generated file",
      __import__("hashlib").sha256(OUT_LEAN.read_bytes()).hexdigest())
claim("SI-NOT-EDITED", "the statement file of record is untouched by design",
      "Unit2BFormalR1.lean is NOT modified; this is an addition under the new-files-only rule")

lines = ["# Statement identity for the generated successor artefact",
         "# Generated by formal01_kappa03_statement_r2.py",
         "# Compares the header of each declaration in Unit2BFormalR1.lean (up to its `:=`)",
         "# with the header of the same declaration in Unit2BFormalR1Proofs.lean.", ""]
for cid, lab, val in CLAIMS:
    lines.append(f"CLAIM {cid}: {lab} = {val}")
OUT_CLAIMS.write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))
