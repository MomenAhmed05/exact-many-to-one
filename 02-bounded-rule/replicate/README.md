# Independent replication kit

**What this is.** Everything needed to reproduce the headline numbers of the desk-guarantee work
from scratch, by a party that has not seen the original codebase. `METHOD.md` is the complete method
specification; `inputs/designs.json` fixes the designs, seeds and protocol constants.

**What is withheld.** The expected numbers. Do not look for them in this kit. Compute yours, report
them, and they will be compared against the sealed originals. The comparison is the point: if this
kit leaks the answers, the replication proves nothing.

## How to run it

1. Read `METHOD.md` end to end. Implement §1 (the exact conditional law) and validate it against the
   F1 identity check before building anything on top of it — if F1 does not reproduce the Student t
   quantile to ~1e-6, the rest is meaningless.
2. Implement §2 (the region and its exact coverage) and confirm that your `cov(s)` is monotone
   decreasing in `s`, with the reported scale giving `cov >= 1 - gamma + 1e-6`.
3. Implement §3 (the two rules) and run §4 over the 12 cells in `inputs/designs.json`.
4. Run the four robustness checks in §5.

Use any language and any library. Nothing about the result depends on the implementation choices,
except that §1's quadrature must converge under refinement (F3).

## What to report

One row per desk cell, with the fields listed in `inputs/designs.json -> report_template.per_cell`,
plus the four falsifier numbers. Then answer three questions in your own words:

1. **Does the envelope rule dominate the oracle per draw?** (paired dominance `<= 0` up to Monte
   Carlo error). If it does not, say on which cells and whether the true ratios lay inside the
   admitted region on those draws.
2. **What did the guarantee cost in power?** Report `cost_pct` per cell and say whether the cost is
   concentrated in one family.
3. **Did any falsifier fail?** In particular F3: if the quadrature does not converge as panels are
   added, the exactness claim is void and that is a reportable failure, not a nuisance.

## Ground rules

- Do not read `../*.py`, `../deposits/**`, `../PRIOR-ART.md`, `../STAGE*.md` or
  `../REPORT-PROVEN-BOUNDARY.md` before reporting your numbers — several carry the expected answers.
- If a step in `METHOD.md` is ambiguous, record the ambiguity and your resolution rather than
  guessing silently.
- Report negative results as prominently as positive ones. The value of this kit is falsification,
  not confirmation.
