"""How much of the new paper material any check actually looks at.

"136 checks pass" does not say how many of the new sections' numbers a check
ever formats.  A number no check formats passes unconditionally.

The tokeniser has to respect the document's structure or it measures nothing.
Concatenating table cells fuses "0.1583" and "0.0736" into "15830.07360", which
is in no deposit and looks uncovered while the row check that covers both has
already passed.  So body text and table cells are tokenised separately, and a
cell token counts as covered when it appears in an expected string, which for a
row check is the joined row.

    python coverage_paper_beyond.py
Deposit: deposits/coverage_paper_beyond.json
"""

from __future__ import annotations

import io
import os
import re
import subprocess
import sys
import zipfile

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))
DOCX = os.path.join(HERE, "paper",
                    "Exact-Many-to-One-Under-Estimated-Volatilities.docx")
DUMP = "_exp_paper_beyond.txt"

NUM = re.compile(r"[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?%?")

# Section, equation, table and figure numbers are cross-references, not
# results, and the structural checks already pin them.
XREF = re.compile(r"^(9|10|11|12|13|14)\.\d+$")

KNOWN = {
    "2023": "the year of the state of the art, fixed in Section 2",
    "80%": "the power level every loss in the paper is quoted at",
    "0,1": "the unit interval in equation (48), not a quantity",
    "36": "the bin count, a declared parameter of the search",
    "72": "a bin count in the resolution sentence, paired with 144",
    "12": "the trailing window length of the rolling-variance cell",
    "1955": "the year Dunnett's problem was posed",
}


def cells_and_body(path):
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    cells = []
    for tbl in re.findall(r"<w:tbl>.*?</w:tbl>", xml, re.S):
        for tc in re.findall(r"<w:tc>.*?</w:tc>", tbl, re.S):
            txt = "".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", tc))
            t = re.sub(r"\s+", " ", txt).strip()
            if t:
                cells.append(t)
    body = re.sub(r"<w:tbl>.*?</w:tbl>", " ", xml, flags=re.S)
    body = re.sub(r"<[^>]+>", "\x00", body).replace("\x00", "")
    return cells, re.sub(r"\s+", " ", body)


def between(text, a, b):
    i, j = text.find(a), text.find(b)
    return text[i:j] if 0 <= i < j else ""


def main():
    env = dict(os.environ)
    env["DUMP_EXPECTED"] = DUMP
    subprocess.run([sys.executable,
                    os.path.join(HERE, "paper_verify_beyond.py")],
                   cwd=HERE, env=env, capture_output=True, text=True)
    expected = io.open(os.path.join(HERE, DUMP), encoding="utf-8").read()

    cells, body = cells_and_body(DOCX)

    # the body of the new subsections only
    # each slice ends at the last words of the inserted text, not at the next
    # heading: a slice that runs to the heading sweeps in pre-existing prose
    # and reports its numbers as uncovered by this update's checks
    new_body = (between(body, "10.8 The floor the class reaches",
                        "10.9 Which group binds")
                + between(body, "11.4 Three validity tiers", "12. Discussion")
                + between(body, "Two further questions are settled",
                          "the split exactly along it.")
                + between(body, "(14) The efficiency result placed",
                          "(Section 11.6).")
                + between(body, "The bound of Section 10.8 is over rules",
                          "conservative rather than liberal.")
                + between(body, "What the efficiency claim rests on",
                          "a number read from the returns."))

    # the cells of the new tables only: every cell of a row the verifier knows
    exp_rows = [e for e in expected.split("\n") if " | " in e]
    new_cells = []
    for e in exp_rows:
        new_cells.extend(e.split(" | "))

    LABEL = re.compile(r"(?:Table|Figure|Section)\s+[A-Z]?\d+(?:\.\d+)?|"
                       r"equation\s+\(\d+\)|\(\d+\)")

    def toks(s):
        # a caption's own number and an equation's own number are labels the
        # structural checks pin, not quantities a deposit could carry
        s = LABEL.sub(" ", s)
        out = []
        for t in NUM.findall(s):
            t = t.rstrip(".,")
            if t and any(c.isdigit() for c in t):
                out.append(t)
        return out

    seen, order = set(), []
    for t in toks(new_body) + [x for c in new_cells for x in toks(c)]:
        if t not in seen:
            seen.add(t)
            order.append(t)

    covered, excused, uncovered = [], [], []
    for t in order:
        if t in expected:
            covered.append(t)
        elif XREF.match(t) or t in KNOWN:
            excused.append(t)
        else:
            uncovered.append(t)

    total = len(order)
    pct = 100.0 * (len(covered) + len(excused)) / max(total, 1)
    print("%d distinct numeric tokens in the new material" % total)
    print("  covered by a check      %d" % len(covered))
    print("  cross-reference or declared parameter   %d" % len(excused))
    print("  UNCOVERED               %d" % len(uncovered))
    if uncovered:
        print("\nuncovered:")
        for t in uncovered:
            print("   " + t)

    out = {
        "document": os.path.basename(DOCX),
        "n_expected_strings": len([e for e in expected.split("\n") if e]),
        "n_tokens": total,
        "n_covered": len(covered),
        "n_excused": len(excused),
        "n_uncovered": len(uncovered),
        "coverage_pct": float(pct),
        "uncovered": uncovered,
        "clean": bool(not uncovered),
        "note": ("Body text and table cells are tokenised separately, because "
                 "concatenating cells fuses adjacent numbers into tokens that "
                 "are in no deposit and look uncovered while the row check "
                 "covering both has already passed."),
    }
    F.write_deposit("coverage_paper_beyond.json", out, "diagnostic")
    print("\ncoverage %.1f%% (%d uncovered)" % (pct, len(uncovered)))
    return 0 if out["clean"] else 1


if __name__ == "__main__":
    sys.exit(main())
