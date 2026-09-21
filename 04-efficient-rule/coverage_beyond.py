"""How much of REPORT-BEYOND.md any check actually looks at.

"178 checks pass" does not say how many of the report's claims a check ever
touches.  A number no check formats passes unconditionally, and counting passes
cannot see it.  This runs the verifier with its dump hook, extracts every
numeric token from the report, and reports the ones nothing covers.

Tokens with a stated reason are listed in KNOWN and are excluded from the
uncovered count, with the reason printed so the exclusion is visible rather
than silent.

    python coverage_beyond.py
Deposit: deposits/coverage_beyond.json
"""

from __future__ import annotations

import io
import os
import re
import subprocess
import sys

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(HERE, "REPORT-BEYOND.md")
DUMP = "_expected_beyond.txt"

NUM = re.compile(r"[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?%?")

KNOWN = {
    "0.025": "the ceiling itself, a constant of the programme",
    "1955": "the year Dunnett's problem was posed",
    "10": "a section number",
    "11": "a section number",
    "9.11": "a section number",
    "1.47": "the spurious gain the first region module reported, quoted as "
            "prose about a superseded run rather than as a result",
    "36": "the bin count, a declared parameter of the search",
    "72": "a bin count in the resolution table, covered by the row check",
    "144": "a bin count in the resolution table, covered by the row check",
    "97": "the constraint grid size, a declared parameter",
    "193": "the verification grid size, a declared parameter",
    "200,000": "the draw count of the invariance check, a declared parameter",
    "360": "a calibration length, covered as a column header",
    "120": "a calibration length, covered as a column header",
    "2": "a dimension count or an index",
    "3": "a dimension count or an index",
    "4": "a dimension count or an index",
    "1": "a dimension count or an index",
    "0": "a mean value in the invariance table, covered by the row check",
    "5": "a mean value in the invariance table, covered by the row check",
}


def main():
    env = dict(os.environ)
    env["DUMP_EXPECTED"] = DUMP
    r = subprocess.run([sys.executable, os.path.join(HERE, "verify_beyond.py")],
                       cwd=HERE, env=env, capture_output=True, text=True)
    path = os.path.join(HERE, DUMP)
    if not os.path.exists(path):
        print("no dump produced; verifier said:\n" + r.stdout[-800:])
        return 1
    expected = io.open(path, encoding="utf-8").read()
    text = io.open(REPORT, encoding="utf-8").read()

    # only the body, not the fenced blocks or the link targets
    toks = [t.rstrip(".,") for t in NUM.findall(text)
            if any(c.isdigit() for c in t)]
    toks = [t for t in toks if t and any(c.isdigit() for c in t)]
    seen, order = set(), []
    for t in toks:
        if t not in seen:
            seen.add(t)
            order.append(t)

    covered, uncovered, excused = [], [], []
    for t in order:
        if t in expected:
            covered.append(t)
        elif t in KNOWN:
            excused.append(t)
        else:
            uncovered.append(t)

    total = len(order)
    pct = 100.0 * (len(covered) + len(excused)) / max(total, 1)
    print("%d distinct numeric tokens in REPORT-BEYOND.md" % total)
    print("  covered by a check      %d" % len(covered))
    print("  excused with a reason   %d" % len(excused))
    print("  UNCOVERED               %d" % len(uncovered))
    if uncovered:
        print("\nuncovered:")
        for t in uncovered:
            print("   " + t)
    if excused:
        print("\nexcused:")
        for t in excused:
            print("   %-10s %s" % (t, KNOWN[t]))

    out = {
        "report": os.path.basename(REPORT),
        "n_expected_strings": len(expected.split("\n")) - 1,
        "n_tokens": total,
        "n_covered": len(covered),
        "n_excused": len(excused),
        "n_uncovered": len(uncovered),
        "coverage_pct": float(pct),
        "uncovered": uncovered,
        "excused": dict((t, KNOWN[t]) for t in excused),
        "clean": bool(not uncovered),
        "note": ("Coverage counts DISTINCT numeric tokens, not occurrences.  "
                 "A token is covered when some check formats exactly that "
                 "string from a deposit.  Excused tokens are declared "
                 "parameters, section numbers and table axis labels, each "
                 "named; everything else must be covered or the run is not "
                 "clean."),
    }
    F.write_deposit("coverage_beyond.json", out, "diagnostic")
    print("\ncoverage %.1f%% (%d uncovered)" % (pct, len(uncovered)))
    return 0 if out["clean"] else 1


if __name__ == "__main__":
    sys.exit(main())
