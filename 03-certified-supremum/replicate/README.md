# Blind replication kit — the certified Welch–Dunnett rule

**What this is.** A self-contained guide for an independent party to reproduce the certified
threshold, the branch-and-bound certificate, and the falsifiers, from the specification alone.

**Order of operations (do not open the targets first).**

1. Read `../SPEC.md` in full. It is the method: the reduced form (section 3), the region and its
   coverage identity (section 4), the enclosure and the certificate quadrature (section 5), the
   certified critical values and their directions (section 6), the branch-and-bound and its
   `c_U = max over ALL leaf boxes` rule (section 7), the tiers (section 8), the error budget
   (section 9), the theorem text and exclusions (sections 10-11), the deposit conventions
   (section 12) and the known failures (section 13).
2. Implement it. `inputs/designs.json` has the designs, the draws' seeds, the levels, and the
   certificate points. No target value appears in it.
3. Run the four falsifiers in `FALSIFIERS.md` and record your outcomes.
4. **Only then** open `targets/expected.json`, which contains the SHA-256 of each deposit plus the
   canonical fields. Compare field by field and report every difference, including any you cannot
   explain.

**What a correct replica should see** (stated qualitatively here so the guide is usable before the
targets are opened): the region's exact quadrature coverage at or above `1 - gamma`; a certified
bracket whose width is far below `1e-4`; the branch-and-bound gap at `1e-6`-scale on narrow regions
and up to a few `1e-3` on the widest ones; and every falsifier passing.
