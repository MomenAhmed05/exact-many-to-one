# Exact Many-to-One Comparison Under Estimated Volatilities

Code, deposits and paper for an exact null distribution of Dunnett's
many-to-one comparison when the group variances are unequal **and estimated**.

Dunnett solved the equal-variance case in 1955. The unequal-variance case has
been treated as intractable since, and the current state of the art
approximates the joint law of the denominators in four places and then
simulates it. This work computes that law exactly, by conditioning on both of
the benchmark's sufficient statistics rather than on its mean alone, which
makes the comparisons exactly independent and turns a `(2K+2)`-dimensional
integral into a two-dimensional outer rule over `K` one-dimensional inner
rules. The outer dimension does not grow with the number of candidates.

The papers are in [`paper/`](paper/): the working paper, and a short companion
with the statements, proofs and evidence alone. The mathematics is checked in
Lean 4 in [`Lean/`](Lean/).

The second paper, *Many-to-One Strategy Comparison on Aligned Returns*, treats
candidates and a benchmark observed on the same dates. Its code, results and
Lean sources are in [`aligned-returns/`](aligned-returns/), which has its own
[README](aligned-returns/README.md) and uses this paper's `04-efficient-rule/` code.

## What is here

Every number printed in the paper is computed by a script in this repository
and written to a JSON **deposit** beside it. Nothing is typed into the document
by hand: the build reads the deposits, and a verifier re-reads the document and
checks each figure against the deposit it came from. That is why the layout
below is organised by which part of the paper a directory produces.

| directory | paper sections | appendix table | what it produces |
| --- | --- | --- | --- |
| [`01-exact-law/`](01-exact-law/) | 3 to 8 | A1 | the exact law, its quadrature, the six falsifiers, the true error rates of the rules in use, the head-to-head against the 2023 method, and the calibrated rule |
| [`02-bounded-rule/`](02-bounded-rule/) | 9 | A2 | the admission region for the variance ratios, its exact coverage, the envelope rule and its guarantee, the power it gives up, and the two boundaries where the guarantee stops |
| [`03-certified-supremum/`](03-certified-supremum/) | 9.11 | A3 | the branch-and-bound certificate over the continuum of ratio vectors, with interval arithmetic and a Bernstein-ellipse quadrature remainder |
| [`04-efficient-rule/`](04-efficient-rule/) | 10 and 11 | A3 | the calibrated rule that meets the same ceiling far more cheaply, its optimality floor over the invariant class, and serial correlation and factor exposure handled as declared nuisances |
| [`Lean/`](Lean/) | 3, 9, 10, 11 | — | the Lean 4 formalisation: every theorem the paper proves, checked against Mathlib, and every deposited certificate's arithmetic checked exactly |
| [`paper/`](paper/) | all | — | the working paper; the companion `Exact-Many-to-One-Proofs.docx`, nine pages of statements, proofs and evidence; and the gate that checks the companion |
| [`aligned-returns/`](aligned-returns/) | paper 2 | — | the second paper: the certified constant for the Welch rule on aligned returns, the exact paired test, and their Lean development |

Each directory keeps its own `deposits/` next to its scripts, because the
scripts resolve deposits relative to their own location. Moving a script away
from its deposits will break it.

The machinery that **builds** the document from the deposits is not here. It is
authoring tooling rather than a result, and the paper it produced is in
`paper/` already. What is kept is the machinery that **checks** the document,
because that is what makes the deposits worth anything.

## Reproducing

The scripts are plain Python and use `numpy`, `scipy` and `python-docx`. They
write deposits; they need no database and no network. Run them from inside
their own directory, which is how they resolve both their deposits and their
siblings.

```bash
cd 04-efficient-rule
python paper_verify_frontier.py ../paper/Exact-Many-to-One-Under-Estimated-Volatilities.docx
python paper_verify_beyond.py   ../paper/Exact-Many-to-One-Under-Estimated-Volatilities.docx
python verify_huntstein.py
python tests_frontier.py
```

### The gates

Every gate checks the document as it stands. Run each from its own directory.

| gate | covers | result |
| --- | --- | --- |
| `01-exact-law/h12_paper_verify.py` | Sections 3 to 8 | 229 checks |
| `02-bounded-rule/paper_verify_w26.py` | Section 9 | 295 checks |
| `03-certified-supremum/tests_gates.py` | the certificate | 9 of 9 gates |
| `04-efficient-rule/paper_verify_frontier.py` | Sections 9.11, 10, 11 | 558 checks |
| `04-efficient-rule/paper_verify_beyond.py` | the invariance and dependence material | 185 checks |
| `04-efficient-rule/verify_huntstein.py` | the Section 10.8 argument | 63 checks |
| `04-efficient-rule/tests_frontier.py` | the deposited gates | 11 of 11 |
| `Lean/scripts/audit.py` | the Lean development | 4 of 4: clean build, no `sorry` or axioms, certificates match the deposits, 41 headline theorems on the standard axioms only |
| `paper/verify_proofs_paper.py` | the companion paper | 8 checks: its numbers appear in the working paper, its Lean names resolve, its numbering is consistent |

All pass, none failing.

A numeric check formats its expected value from a deposit and looks for that
string in the document, so editing the paper cannot make a check pass and
editing a deposit cannot make one pass either. Structural checks run beside
them: the figure and table counts, the heading numbers, and the
cross-references that have to move when a section is inserted.

## The Lean development

[`Lean/`](Lean/) proves in Lean 4 what the paper proves:
the reduction behind the exact law, the envelope guarantee and the coverage
integral, the corner enclosure against the Gaussian law, the branch-and-bound
bracket, the union step, the two-sided identity, the joint projection, the
maximal invariant and the whitening and factor results. The interval-arithmetic
enclosures stay in Python. Lean checks exactly that each deposited certificate
supports the claim the paper makes from it. The Hunt–Stein theorem is quoted,
as in the paper, and enters one theorem as a named hypothesis. Its
[README](Lean/README.md) maps each result to its Lean theorem and
lists every assumption.

```bash
cd Lean
lake exe cache get
python scripts/audit.py
```

## Reading the record

The working notes are kept rather than tidied away, because several results in
the paper exist only because a check failed first. Each module has its own:

- `02-bounded-rule/STAGE1-CERTIFICATE.md` to `STAGE5-PUBLICATION.md`, and
  `PRIOR-ART.md`.
- `03-certified-supremum/PHASE1-ENCLOSURE.md` to `PHASE6-EXTENSIONS.md`, with
  `SPEC.md` and `GOVERNANCE.md`.
- `04-efficient-rule/REPORT-FRONTIER.md`, `REPORT-BEYOND.md`,
  `DESK-PROCEDURE.md`, and `HUNT-STEIN.md`.

Two are worth naming. `DESK-PROCEDURE.md` is the practitioner's version: what
to fix before seeing the data, how to read the threshold off, what is
guaranteed and where the guarantee stops. `HUNT-STEIN.md` carries the
invariance argument of Section 10.8 in full, including the condition it takes
from the literature rather than proving, and three attempts at a falsifier that
failed before one worked.

## Known rough edges

- The `.log` files written by long runs are not tracked.
- A few scripts read deposits that a later script supersedes. The appendix
  tables name the source for each table and figure.

## Licence

MIT. See [LICENSE](LICENSE).
