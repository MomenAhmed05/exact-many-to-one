# The promotion protocol

`DESK-PROCEDURE.md` says how to compute the threshold.  This says what has to be
true about the *process* for that threshold to mean what it claims.  The two are
not separable: a familywise error rate is a statement about a family, and a
family that was chosen after looking at the outcome is not the family the
mathematics is about.

Nothing here can be checked by a statistic.  All of it can be checked by a file.
That file is the manifest, its shape is `manifest.schema.json`, and
`w4_manifest.py` refuses runs that do not satisfy it.

---

## The three roles data can play

Every observation belongs to exactly one, declared in advance.

| Role | What happens to it | What it costs you |
|---|---|---|
| **Research** | Generate hypotheses, tune, sweep, inspect, discard | Nothing is guaranteed about anything done here, and nothing needs to be |
| **Calibration** | Estimate declared dependence nuisances: an AR(1) coefficient, factor loadings | The transform becomes random but stays independent of the evaluation data |
| **Evaluation** | Apply one locked family and one locked threshold table, once | This is the only window the guarantee describes |

Calibration may be empty.  Research may be empty.  Evaluation may not, and it
may not overlap either of the others.

**The verifier refuses** when the research window ends after the evaluation
window opens, and when a calibration window overlaps the evaluation window.

---

## The lock

Before any evaluation-period performance is inspected, write the manifest and
record its hash.  It carries:

* every candidate's identifier and code hash, and the benchmark's;
* the return construction: frequency, commission, slippage, financing,
  rebalance timing, and what happens to a missing period;
* the three windows;
* the hypothesis direction and the familywise target;
* the dependence model and its validity tier, for serial and cross-sectional
  dependence separately;
* the variance-ratio domain **and where it came from**;
* the threshold table's version and hash, and the certificate deposit it
  came from;
* the primary specification and any sensitivity specifications.

Amendments are allowed and are recorded with a timestamp and a reason.  An
amendment dated inside the evaluation window does not void the run, but it turns
a `GUARANTEED` verdict into a `CONDITIONAL` one and the condition is printed.

---

## The three verdicts

`w4_manifest.py` returns one of three labels, and the report must carry it.

**GUARANTEED.**  Every condition the certificate needs is recorded and
internally consistent.  The familywise error rate is at most the declared target
over the declared nuisance domain.

**CONDITIONAL.**  The procedure is valid given something the manifest declares
but cannot evidence.  The verdict names it.  The common cases:

* a dependence nuisance estimated on a calibration period rather than declared,
  so the transform's own error is outside the certificate;
* a nuisance estimated on the evaluation period, which makes the result a
  plug-in and not a certificate;
* `independent_groups` asserted for a panel observed on common dates, which is
  usually false and is the subject of `w2_factor.py`;
* the family locked inside the evaluation window, so only the part after the
  lock is covered.

**SCREEN.**  At least one condition is violated.  The test may still be worth
running and its number may still be informative, but it is a screen and the
report says so.  The refusals, each on its own ground:

| Refusal | Why it is fatal |
|---|---|
| family locked after the evaluation window closed | the candidate set could have been chosen on the outcome |
| research window ends after the evaluation window opens | hypotheses were generated on evaluation data |
| calibration window overlaps the evaluation window | the transform is a function of the data under test |
| variance-ratio domain inferred from the returns under test | the domain is then a function of the data, and a uniform bound over a data-dependent set is not a bound |
| Prais-Winsten whitening declared | it rescales the first observation, so the whitened mean depends on the group's length and the rule drifts off its ceiling at unequal lengths |
| internal level above the familywise target | the arithmetic is backwards |
| the benchmark also listed as a candidate | the comparison is with itself |
| a required field missing | the condition cannot be checked, which is the same as it not holding |
| threshold table hash mismatch | the table applied is not the table certified |

`same_returns_being_tested` is a permitted *value* for the ratio domain's source
precisely so that the verifier can refuse it.  Removing it from the schema would
only push the author towards writing something else.

---

## What this protocol does not do

It controls multiplicity **within a declared family**.  It does not repair an
adaptive search, and no procedure does: White's Reality Check, Hansen's SPA and
Romano-Wolf all require the family to be declared before the data speak.  Once
two hundred signals have been inspected and five retained, the object the
guarantee describes no longer exists.  The answer to that is the three-role
split above, not a threshold.

It is therefore a **last-stage promotion gate**, complementary to and not a
replacement for:

* White's Reality Check and Hansen's SPA, which control a different target under
  bootstrap resampling assumptions;
* Romano-Wolf stepdown, which gains power over a single-step rule on a declared
  family at the cost of a resampling approximation;
* false-discovery-rate procedures, which control a different quantity entirely
  and are the right choice when a fixed small number of false promotions is
  acceptable;
* deflated Sharpe ratio and probability-of-backtest-overfitting methods, which
  price search directly and are the right tool for the research window;
* walk-forward and embargoed cross-validation, which are design devices rather
  than tests.

---

## The claim taxonomy

Every number in a report carries exactly one label, and the deposit it came from
carries the same one.

| Label | Meaning |
|---|---|
| **Exact in model** | deterministic evaluation of the correct finite-sample quantity, conditional on stated nuisance values |
| **Certified bound** | a proved upper or lower bound over a stated nuisance domain |
| **Measured** | a simulation estimate, reported with its Monte Carlo standard error |
| **Asymptotic** | valid only in a stated limiting regime |
| **Diagnostic** | a stress-test result, not a guarantee |

A grid maximum is **measured**, not **certified**.  The distinction is the one
the branch-and-bound work of Section 9.11 exists to close, and calling a grid
maximum a certificate is the specific failure this programme is set up to
prevent.

---

## Report template

```
PROMOTION REPORT
  manifest            <path>  sha256 <hash>
  verdict             GUARANTEED | CONDITIONAL | SCREEN
  conditions          <one line each, verbatim from the verifier>

  family              K = <n> candidates against <benchmark>
  window              <evaluation_start> to <evaluation_end>, <n> observations
  target              FWER <= <abar>, <direction>

  dependence          serial <tier> <transform> phi <value|bound>
                      cross-sectional <tier> <model>
  ratio domain        [<lo>, <hi>]  source <source>

  threshold table     <version>  alpha0 <value>   [certified bound]
  certificate         <deposit>  sup <value> over <domain>

  RESULTS
    candidate    T       lambda_hat   threshold   decision
    ...

  LABELS
    threshold           certified bound
    realised FWER       measured (mcse <value>)
    power at <delta>    exact in model
    stress results      diagnostic
```

A report whose verdict is `SCREEN` prints the refusals immediately under the
verdict and does not print the word *guarantee* anywhere.
