"""Statement identity and forbidden-token scan for the FORMAL-01 successor Unit2BFormalR1R2.lean.

Part A.  For each of the thirteen theorems, extract the statement from the statement file of
record (`Unit2BFormalR1.lean`) and from the successor (`Unit2BFormalR1R2.lean`): the text from
`theorem <name>` up to (not including) the first `:=` that follows it.  Normalise whitespace
(every run of whitespace -> one space, ends stripped) and compare character for character.
Also reported: whether the RAW span is byte-identical before normalisation.

Part B.  Scan both new Lean files for `sorry`, `admit`, `native_decide`:
  (i)  on the RAW text, comments and docstrings included, case-insensitively, as substrings
       -- the strongest form: nothing anywhere in the file;
  (ii) on CODE ONLY, after a real Lean comment stripper (nested `/- -/` blocks, which covers
       `/--` and `/-!` docstrings, and `--` line comments, with string literals respected so a
       `--` inside a string does not start a comment).
The stripper is self-tested on a synthetic sample first, so a stripper that silently removes
everything (which would make (ii) vacuously zero) is caught.

Writes deposits/statement_identity_scan_claims.txt.  It does not invoke Lean.
"""
import hashlib
import pathlib
import re
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = pathlib.Path(__file__).resolve().parent
W = ROOT / "lean" / "W4lean"
STMT = W / "Unit2BFormalR1.lean"
SUCC = W / "Unit2BFormalR1R2.lean"
IDF = W / "StatementIdentityR2.lean"

NAMES = ["Phi_antitone", "ae_nonneg_of_gammaLaw", "charQ_eq_det", "discG_eq", "charQ_factor",
         "mu_sum_prod", "box_pos", "r1_lower", "r2_lower", "cstar_bracket", "cstarBox_eq_sSup",
         "cover_sound", "comp_thm_2B"]
TOKENS = ["sorry", "admit", "native_decide"]


def strip_comments(src):
    """Remove Lean comments.  Nested block comments, line comments, strings respected."""
    out, i, n, depth = [], 0, len(src), 0
    in_str = False
    while i < n:
        if depth > 0:
            if src.startswith("/-", i):
                depth += 1
                i += 2
            elif src.startswith("-/", i):
                depth -= 1
                i += 2
            else:
                if src[i] == "\n":
                    out.append("\n")      # keep line structure
                i += 1
            continue
        if in_str:
            out.append(src[i])
            if src[i] == "\\" and i + 1 < n:
                out.append(src[i + 1])
                i += 2
                continue
            if src[i] == '"':
                in_str = False
            i += 1
            continue
        if src.startswith("/-", i):
            depth = 1
            i += 2
            continue
        if src.startswith("--", i):
            j = src.find("\n", i)
            i = n if j < 0 else j
            continue
        if src[i] == '"':
            in_str = True
        out.append(src[i])
        i += 1
    if depth != 0 or in_str:
        raise ValueError(f"unterminated comment/string (depth={depth}, in_str={in_str})")
    return "".join(out)


def code_hits(code, tok):
    return len(re.findall(rf"(?<![\w.]){tok}(?![\w])", code))


def extract(text, name):
    ms = list(re.finditer(rf"^theorem\s+{re.escape(name)}(?=[\s(\[{{:])", text, re.M))
    if len(ms) != 1:
        return None, len(ms)
    start = ms[0].start()
    end = text.find(":=", start)
    return text[start:end], 1


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


out = []

# ---- stripper self-test -------------------------------------------------------------------
sample = ('/-- doc sorry -/\n/-! mod admit -/\n/- outer /- nested sorry -/ still -/\n'
          'theorem t : True := by sorry -- trailing sorry\n'
          'def s := "-- not a comment" ++ "x"\n-- line native_decide\nexample := admit\n')
sc = strip_comments(sample)
st_ok = (code_hits(sc, "sorry") == 1 and code_hits(sc, "admit") == 1
         and code_hits(sc, "native_decide") == 0 and '"-- not a comment"' in sc
         and "theorem t" in sc and "nested" not in sc and "doc" not in sc)
out.append(f"CLAIM SI-STRIPPER-SELFTEST: synthetic sample with 5 commented and 2 code "
           f"occurrences yields code hits sorry=1 admit=1 native_decide=0, docstrings and nested "
           f"blocks removed, string literal preserved = {st_ok}")

# ---- Part A: statement identity ------------------------------------------------------------
ts, tr = STMT.read_text(encoding="utf-8"), SUCC.read_text(encoding="utf-8")
all_ok = True
for nm in NAMES:
    a, ca = extract(ts, nm)
    b, cb = extract(tr, nm)
    if a is None or b is None:
        all_ok = False
        out.append(f"CLAIM SI-{nm}: identical = False  (EXTRACTION FAILED: occurrences "
                   f"statement={ca}, successor={cb})")
        continue
    same = norm(a) == norm(b)
    raw = a == b
    all_ok &= same
    out.append(f"CLAIM SI-{nm}: identical = {same}  (raw span byte-identical = {raw}; "
               f"normalised length {len(norm(a))}/{len(norm(b))})")
    if not same:
        na, nb = norm(a), norm(b)
        k = next((i for i in range(min(len(na), len(nb))) if na[i] != nb[i]), min(len(na), len(nb)))
        out.append(f"  !!! STATEMENT DIFFERS at char {k}: statement=...{na[max(0,k-40):k+40]!r}"
                   f" successor=...{nb[max(0,k-40):k+40]!r}")
out.append(f"CLAIM SI-ALL: all thirteen statements identical after whitespace normalisation = "
           f"{all_ok}")
extra = sorted(set(re.findall(r"^theorem\s+([A-Za-z_][\w']*)", tr, re.M)) - set(NAMES))
out.append(f"CLAIM SI-EXTRA-DECLS: theorems in the successor that are not among the thirteen = "
           f"{', '.join(extra) or 'none'}")

# ---- Part B: forbidden-token scan ----------------------------------------------------------
for p in (SUCC, IDF):
    raw = p.read_text(encoding="utf-8")
    code = strip_comments(raw)
    rawc = {t: raw.lower().count(t) for t in TOKENS}
    codec = {t: code_hits(code, t) for t in TOKENS}
    out.append(f"CLAIM SCAN-RAW-{p.name}: case-insensitive substring counts over the FULL text, "
               f"comments and docstrings included = {rawc}")
    out.append(f"CLAIM SCAN-CODE-{p.name}: word counts over code with comments and docstrings "
               f"stripped = {codec} (code chars kept {len(code)} of {len(raw)})")
    out.append(f"CLAIM SCAN-ZERO-{p.name}: = {sum(rawc.values()) == 0 and sum(codec.values()) == 0}")

out.append(f"CLAIM SI-STMT-FILE-SHA256: {STMT.name} = {sha(STMT)}")
out.append(f"CLAIM SI-SUCC-SHA256: {SUCC.name} = {sha(SUCC)}")
out.append(f"CLAIM SI-IDF-SHA256: {IDF.name} = {sha(IDF)}")

(ROOT / "deposits" / "statement_identity_scan_claims.txt").write_text("\n".join(out) + "\n",
                                                                     encoding="utf-8")
print("\n".join(out))
