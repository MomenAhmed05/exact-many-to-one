# Blind replication kit — the certified adaptive Welch rule (K = 1)

Read `METHOD.md` and `inputs/designs.json`.  Do not read the programme source.  Do not open
`targets/expected.json` until your numbers exist.

    replicate/
      METHOD.md              the specification: model, law, rule, constraint, certificate
      FALSIFIERS.md          the four checks you must run before comparing
      inputs/designs.json    designs, domain, bin rule, quadrature and enclosure settings
      targets/expected.json  SEALED: canonical values plus a SHA-256 per deposit

## Why this kit is shaped this way

The preceding programme's blind replication reproduced its draw-independent half exactly and failed
its draw-dependent half, because the specification never pinned down the draw recipe, the definition
of `d`, or the semantics of a "certificate point".  Eight ambiguities were logged.

This construction has **no draws in the decision path at all**.  The rule is a deterministic function
of the sample variances, the constraint is a deterministic integral, and the certificate is a
deterministic enclosure.  There is therefore no draw recipe for a replica to guess, which removes the
entire class of defect that broke the last attempt.  Monte Carlo appears only in falsifier R3, where
it is a check on the bound and its own seed does not need to match.

## Status

**NOT EXECUTED.**  This kit has been written and sealed but no independent party has run it.  Until
one has, the correct description of this work is "cross-checked by an independent numerical route in
the same session" (see `deposits/crosscheck.json`), which is weaker, and the report says so.
