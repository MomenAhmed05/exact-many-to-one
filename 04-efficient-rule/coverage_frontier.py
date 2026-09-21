"""How much of the new material is actually covered by a check.

`paper_verify_frontier.py` reports how many checks pass.  That is not the same
question as how many of the paper's numbers a check ever looks at: a claim no
check formats passes unconditionally, and counting passes cannot see it.

This reads every numeric token out of the sections the update added, and asks
of each whether it appears inside a string some check looked for.  What it
prints is the uncovered list, which is the only part worth reading.

    DUMP_EXPECTED=_expected.txt python paper_verify_frontier.py
    python coverage_frontier.py
"""

from __future__ import annotations

import io
import os
import re
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DOCX = os.path.join(HERE, "paper",
                    "Exact-Many-to-One-Under-Estimated-Volatilities.docx")
EXPECTED = os.path.join(HERE, "_expected.txt")

# The spans the update owns.  Everything outside them belongs to earlier
# programmes and is checked by their own verifiers.
SPANS = [("9.11 A certified supremum", "12. Discussion"),
         ("That bound is then made both proven and efficient", "1. Introduction"),
         ("(9) A certified supremum over the continuum", "2. Background"),
         ("Five of the seventeen branch-and-bound cells", "14. Conclusion"),
         ("That second rule is not the one to run", "References")]

# A number as a reader meets it: 0.025, 72.1 pp, 3.17e-03, 44.9%, 2,400, 1e-4.
NUM = re.compile(r"\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?%?")


def tokens(text):
    """Numeric tokens, with list punctuation stripped off the ends."""
    return [t.rstrip(",.") for t in NUM.findall(text) if t.rstrip(",.")]

# Quoted numbers no deposit in this programme can carry, with the reason.
# They are still claims; they are claims sourced elsewhere, and the text says so.
KNOWN = {
    "18.5": "Section 9.7's own figure, from the preceding programme",
    "33%": "Table 12's own figure, attributed in the text",
    "23%": "Section 9's own figure for the region rule's price in the constant",
    "324%": "Table 12's own figure, attributed in the text",
    "1e-6": "a rounded characterisation of the gap, not a deposited value",
    "80%": "the definition of the effect size power is measured at",
}

# Tokens that are structural rather than quantitative claims.
IGNORE = {
    # section, table, figure and equation numbers, and the years in citations
    "9.11", "10.1", "10.2", "10.3", "10.4", "10.5", "10.6", "10.7", "10.8",
    "10.9", "10.10", "11.1", "11.2", "11.3", "12.1", "12.2", "12.3",
    "9.3", "9.6", "9.7", "9.8", "9.10", "3", "4", "9", "10", "11", "12",
    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
    "25", "26", "27", "37", "38", "39", "40", "41", "42", "43", "44", "45",
    "46", "1", "2", "5", "6", "7", "8",
}


def doc_paragraphs(path):
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    out = []
    for para in re.findall(r"<w:p[ >].*?</w:p>", xml, re.S):
        txt = "".join(re.findall(r"<w:[tT][^>]*>([^<]*)</w:t>", para))
        out.append(re.sub(r"\s+", " ", txt).strip())
    return [p for p in out if p]


def main():
    paras = doc_paragraphs(DOCX)
    if not os.path.exists(EXPECTED):
        sys.exit("run: DUMP_EXPECTED=_expected.txt python "
                 "paper_verify_frontier.py")
    expected = io.open(EXPECTED, encoding="utf-8").read()

    body = []
    for start, stop in SPANS:
        i = next((k for k, p in enumerate(paras) if p.startswith(start)), None)
        if i is None:
            print("  span not found: %r" % start)
            continue
        j = next((k for k in range(i + 1, len(paras))
                  if paras[k].startswith(stop)), len(paras))
        body.extend(paras[i:j])

    seen, covered, missing, explained = set(), 0, [], 0
    for para in body:
        for tok in tokens(para):
            if tok in IGNORE or tok in seen:
                continue
            seen.add(tok)
            if tok in expected:
                covered += 1
            elif tok in KNOWN:
                explained += 1
            else:
                where = para[:70]
                missing.append((tok, where))

    total = covered + len(missing)
    print("numeric tokens in the new material : %d" % total)
    print("covered by at least one check      : %d  (%.1f%%)"
          % (covered, 100.0 * covered / total if total else 0.0))
    print("uncovered, sourced elsewhere       : %d" % explained)
    print("uncovered and unexplained          : %d" % len(missing))
    for tok, why in sorted(KNOWN.items()):
        print("  sourced elsewhere: %-6s %s" % (tok, why))
    if missing:
        print()
    for tok, where in sorted(missing):
        print("  %-12s  %s" % (tok, where))


if __name__ == "__main__":
    main()
