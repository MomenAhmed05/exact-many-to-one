# PHASE 4 — the shipped theorem-grade procedure (plan M5: sections 1.3, 4, 4.4, 4.2)

tier: region_sup_certified

**Status: complete; P4 verifier round fired after the blind replication returns.** Every number is read
from a strict-JSON deposit under `deposits/`.

---

## 1. What P4 delivered

| file | role |
|---|---|
| `rule_certified.py` | the eight-step algorithm at `alpha0 = abar - gamma`, with a per-application decision certificate |
| `region_shapes_cert.py` | shapes/allocations compared at **equal certified coverage** using the certified supremum |
| `theorem_text.py` | the required theorem statement and the model exclusions, in one place |
| `desk_map_alpha0.py` | the desk map re-run at `alpha0` with the tier boundary and the price of the proof |
| `tests_gates.py` | the E1-E8 verification gates, including the tier/reserved-word gate |
| deposits | `rule_certified.json`, `region_shapes_cert.json`, `desk_map_alpha0.json`, `gates.json` |

## 2. The shipped rule — `deposits/rule_certified.json`

`alpha0 = abar - gamma = 0.0245` (the hard-ceiling parameterisation: the guarantee reads
`FWER <= alpha0 + gamma = abar`, never "`alpha + gamma` described as a hard `abar`").  For four
comparison families the certificate records the inputs, the ratio estimates, the region with its exact
coverage, the certified threshold with its bracket and gap, the statistics, the decisions, the
branch-and-bound certificate identity, the error budget, the software version and the tier:

| family | `c_hat` | gap | certified at tol | promoted | seconds (target 60) |
|---|---|---|---|---|---|
| desk_B_df3 | 31.512092 | 9.54e-07 | yes | 0 | 16 |
| desk_A_df3 | 12.500415 | 3.81e-05 | yes | 0 | 26 |
| desk_B_df35 | 7.480059 | 9.54e-07 | yes | 0 | 10 |
| K5_df11 | 12.874008 | 9.54e-07 | yes | 0 | 36 |

**4/4 certified at `tol`, 4/4 regions coverage-certified, all inside the runtime target** (max 36 s).
The decisions on synthetic null data promote nothing, which is the expected behaviour of a valid rule
at these designs.

## 3. Region shapes and allocations at equal certified coverage — `deposits/region_shapes_cert.json`

Every member is evaluated at its **own** tight scale, so all members attain the same certified coverage
(`all_coverage_equal = True`).  The declared family: the shipped equal-tail split (`w_frac = 0.5`),
a candidate-heavy split (`0.25`) and a benchmark-heavy split (`0.75`).

| cell | equal-tail (`0.5`) | candidate-heavy (`0.25`) | benchmark-heavy (`0.75`) | best in family | shipped-vs-best |
|---|---|---|---|---|---|
| A_df3 | 12.500415 | **11.749724** | 13.359213 | candidate-heavy | **+6.39%** |
| B_df3 | 31.512092 | **29.082952** | 34.243819 | candidate-heavy | **+8.35%** |
| B_df11 | 10.400811 | **10.072907** | 10.834241 | candidate-heavy | +3.26% |
| K3_df6 | 20.375705 | **19.046060** | 22.198781 | candidate-heavy | +6.98% |

**The candidate-heavy allocation is best in the declared family in every cell tested, and the shipped
equal-tail split is 3.3-8.4% worse than it.**  That is a substantive finding, and it is reported as
"best in the declared family" — never as globally optimal.  A wider search over `w_frac` (and over
non-box shapes) is the natural next step and is recorded as an open item rather than claimed.

**Reconciliation with the 02-bounded-rule shape floor.** The 02-bounded-rule result that "every shape changes `c` by at
most 0.39%" concerned the *population / vertex-envelope* object under its declared constructions.
The comparison above concerns the **certified continuous-supremum threshold**, which includes the
geometry-induced enclosure slack, and that slack penalises different region geometries differently.
The two results are therefore not in conflict, and together they are a design finding: *the
allocation that minimises a population or vertex-envelope proxy need not minimise the certified
threshold.*  In the tested three-allocation family the candidate-heavy split reduces the certified
threshold by 3.3-8.4% relative to the equal-tail split, plausibly because it reduces the
enclosure slack as well as changing the region itself.  The claim is limited to that family.

## 4. The required theorem statement and the exclusions — `theorem_text.py`

    Under independent Gaussian samples across groups, with the stated sample sizes and the specified
    separate-variance Welch statistics, this decision uses internal level 0.0245 = 0.025 - 0.0005 and a
    ratio-region miss budget 0.0005.  The certified region coverage is at least 1 - 0.0005, and the
    reported threshold is a certified upper bound on the oracle critical value throughout that region.
    Therefore the familywise error rate is at most 0.0245 + 0.0005 = 0.025.

Model exclusions, always printed with it: serial dependence; cross-group dependence beyond the model;
factor exposure; non-Gaussianity; rolling re-estimation; research-process multiplicity.

## 5. The tier gate and the other gates — `deposits/gates.json`

**8/8 gates pass**, including **E4** (the tier/reserved-word gate: every deposit carries a valid tier,
no measured-tier text overclaims, no report violations) and **E2** (certificate-tree integrity: 68 leaf
proofs re-checked).  The gate set is the enforcement arm of the claim discipline: a measured-tier
deposit cannot be described as guaranteed/certified/controlled, and the certified tier is the only one
that carries the `alpha0 + gamma` theorem.

## 6. The desk map at `alpha0` and the boundary between the tiers — `deposits/desk_map_alpha0.json`

Per cell the deposit records four objects, each labelled with its tier: the population oracle
(`oracle_exact_law`, a reference), the 02-bounded-rule vertex envelope (`envelope_measured`, the previously
shipped rule), the certified conditional supremum (`region_sup_certified`, the shipped rule) and its
gap.

**The like-for-like price of the proof is the branch-and-bound gap** (`c_U - c_L`, where `c_L` is the
certified witness lower bound): `9.5e-07` on the narrow cells and at most `3.2e-03` on the widest
family-A cells.  Comparing the conditional `c_U` against the *population* envelope or the *population*
vertex maximum is an object mismatch — three separate instances of that mistake were found and removed
during verification — and the deposit keeps such comparisons only as clearly-labelled informational
columns.

The boundary this establishes:

| tier | what it is | price |
|---|---|---|
| `region_sup_certified` | the shipped threshold: a certified upper bound on the continuous supremum of the conditional critical value over the coverage-certified region, at `alpha0` | the gap: `1e-06` … `3.2e-03` in `c` |
| `envelope_measured` | the 02-bounded-rule vertex envelope (population critical values at the vertices) | measured-only; its sup step is empirically exact (16/16 cells: the grid and multi-start local maxima find nothing above the vertices) |
| `oracle_exact_law` | the population oracle at the true ratios | a reference, not a rule |

## 7. Verifier round P4

Fired after the blind replication returns (one verifier at a time): an independent audit of
`rule_certified.py`'s eight steps against the plan's §4.3 text, of the equal-coverage claim in
`region_shapes_cert.py` (are the members really at the same certified coverage?), and of the tier
labelling in `desk_map_alpha0.py`.

## 8. Open items carried forward

1. The allocation search covers three `w_frac` values; the candidate-heavy optimum is a *family*
   result, and a finer search (plus non-box shapes) is the next step.
2. The desk map is complete: **7 of the 12 desk cells** certify at `tol = 1e-4` (all six family-B
   cells and A_df3), the other five report gaps of `7.1e-04 ... 3.17e-03`, and the like-for-like
   price against the conditional vertex maximum is **median `0.0003%`, max `0.0691%`**
   (`deposits/desk_map_alpha0.json`).
3. The runtime target is met with margin on these families; the widest family-A cells cost ~2 min per
   family, which is above the 1-minute certification target and is recorded as a cost, not hidden.

## 9. Reproduction

    python rule_certified.py        # the shipped rule + decision certificates -> rule_certified.json
    python region_shapes_cert.py    # shapes at equal certified coverage     -> region_shapes_cert.json
    python theorem_text.py          # the statement and the exclusions
    python desk_map_alpha0.py       # the desk map at alpha0                 -> desk_map_alpha0.json
    python tests_gates.py           # the E1-E8 gates                       -> gates.json
