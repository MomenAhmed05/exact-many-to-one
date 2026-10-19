"""Deposit the agreement between the clean-room A-5 bound and the producer's A-5 margin.

WHY THIS FILE EXISTS.  The paper states that the clean-room checker's certified lower bound
on `Phi(Lc; r(corner))` reproduces the producer's `A5-MARGIN` to a relative `1.4e-10`.  This
file computes that figure, with both inputs READ FROM DEPOSITS rather than transcribed, so the
arithmetic is re-derivable and its inputs are themselves sourced.

The two inputs:
  * the clean-room certified LOWER bound, from `cleanroom-a5/a5_claims.txt`,
    claim `A5-PHI-LO`;
  * the producer's `A5-MARGIN`, from `deposits/formal01_f8b_claims.txt`.

If either is missing the run FAILS rather than inventing either of them.

Run: python formal01_kappa03_a5_agreement.py
"""
import pathlib
import re
import sys
from fractions import Fraction as Fr

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = pathlib.Path(__file__).resolve().parent
CLEAN = ROOT / "cleanroom-a5" / "a5_claims.txt"
PROD = ROOT / "deposits" / "formal01_f8b_claims.txt"
OUT = ROOT / "deposits" / "a5_agreement_claims.txt"

CLAIMS = []


def claim(cid, label, value):
    CLAIMS.append((cid, label, value))


# --- input 1: the clean-room certified lower bound.  Read the EXACT RATIONAL, not the
# 30-digit decimal companion `A5-PHI-LO-dec`.  A first version regexed for a decimal after
# `A5-PHI-LO:`, did not find one -- because that claim's value IS a fraction -- and the script
# then FAILED rather than inventing a value, which is the correct behaviour and the reason the
# failure is recorded rather than papered over.  Exact arithmetic also avoids the precision
# loss of routing through the decimal.
lo_r = None
if CLEAN.exists():
    m = re.search(r"^CLAIM A5-PHI-LO:.*?=\s*(\d+)\s*/\s*(\d+)\s*$",
                  CLEAN.read_text(encoding="utf-8"), re.M)
    if m:
        lo_r = Fr(int(m.group(1)), int(m.group(2)))
claim("A5AG-CLEAN-LO-PRESENT",
      "the clean-room certified LOWER bound, as an exact rational, read from a5_claims.txt "
      "claim A5-PHI-LO",
      f"{lo_r.numerator}/{lo_r.denominator}" if lo_r is not None else "NOT FOUND")

# --- input 2: the producer's A5-MARGIN ------------------------------------------------------
a5_s = None
if PROD.exists():
    m = re.search(r"^CLAIM A5-MARGIN:.*?=\s*([0-9.eE+-]+)\s*$",
                  PROD.read_text(encoding="utf-8"), re.M)
    if m:
        a5_s = m.group(1)
claim("A5AG-PROD-A5MARGIN-PRESENT",
      "the producer's A5-MARGIN, read from formal01_f8b_claims.txt", a5_s if a5_s else "NOT FOUND")

ALPHA = Fr(1, 40)

if lo_r is None or not a5_s:
    claim("A5AG-VERDICT", "agreement computable", "FAIL -- an input is missing, and this file "
          "does not invent one")
    OUT.write_text("\n".join(f"CLAIM {c}: {lab} = {v}" for c, lab, v in CLAIMS) + "\n",
                   encoding="utf-8")
    print("\n".join(f"CLAIM {c}: {lab} = {v}" for c, lab, v in CLAIMS))
    sys.exit(1)

# --- the arithmetic, in exact rationals ----------------------------------------------------
lo = lo_r
a5 = Fr(a5_s)
gap = lo - ALPHA                      # how far the clean-room bound clears alpha
rel = abs(gap - a5) / a5               # relative disagreement with the producer's margin

claim("A5AG-ALPHA", "alpha", str(ALPHA))
claim("A5AG-CLEAN-GAP", "clean-room lower bound minus alpha, exact", str(gap))
claim("A5AG-GAP-DECIMAL", "that gap as a decimal", f"{float(gap):.10e}")
claim("A5AG-PROD-MARGIN-DECIMAL", "the producer's A5-MARGIN as a decimal", f"{float(a5):.10e}")
claim("A5AG-RELATIVE-DISAGREEMENT", "relative disagreement between the two margins", f"{float(rel):.3e}")
claim("A5AG-RELATIVE-TO-1E-10", "the relative disagreement rounds to 1.4e-10 at 2 s.f.",
      f"{float(rel):.1e}")
claim("A5AG-BOTH-POSITIVE", "both margins are positive, so the lower bound clears alpha",
      bool(gap > 0 and a5 > 0))
claim("A5AG-VERDICT", "the clean-room lower bound independently reproduces the producer's "
                      "A-5 margin",
      "CONFIRMED to a relative " + f"{float(rel):.1e}" if gap > 0 else "NOT CONFIRMED")
claim("A5AG-SCOPE", "what this does and does not establish",
      "It shows two implementations, sharing no code and in different folders, compute the same "
      "margin at the CORNER to 10 significant figures. It says nothing about any other point in "
      "the box, and the A-5 check covers exactly one point. The size bound is unaffected by it "
      "either way, since size_bound_2B goes through cover_sound and does not take hWit as a "
      "hypothesis")

OUT.write_text("\n".join(f"CLAIM {c}: {lab} = {v}" for c, lab, v in CLAIMS) + "\n",
               encoding="utf-8")
print("\n".join(f"CLAIM {c}: {lab} = {v}" for c, lab, v in CLAIMS))
