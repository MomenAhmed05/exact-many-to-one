# Aligned returns: code, results and Lean sources for paper 2

This folder holds everything behind the second paper, *Many-to-One Strategy Comparison on Aligned Returns: a certified constant for the Welch rule, and the exact paired test that removes the need for it*. The paper is [`../paper/Aligned-Returns-Paper-2.docx`](../paper/Aligned-Returns-Paper-2.docx); its Appendix A names the script and deposit behind each section. It sits in the same repository as paper 1 because several scripts here use paper 1's code in [`../04-efficient-rule/`](../04-efficient-rule/).

Each script writes its results to a text or JSON file under [`deposits/`](deposits/). The paper calls such a file a deposit. Most results files are lists of lines of the form `CLAIM <id>: <description> = <value>`, and every number in the paper is one of those values, a Lean declaration, or a cited figure.

## Requirements

- Python 3.12 with `numpy`, `scipy`, `sympy`, `mpmath`, `python-flint` (ball arithmetic) and `matplotlib` (figures only).
- For the Lean part: [elan](https://github.com/leanprover/elan). The toolchain (Lean 4 v4.34.0) and the Mathlib version are pinned in [`lean/`](lean/).

Run every script from this folder, for example `python formal01_f8c_certified.py`. Scripts find their inputs relative to their own location.

## Layout

| path | contents |
| --- | --- |
| `*.py` | the scripts, one per check or measurement |
| `deposits/` | their results files |
| `cleanroom-a5/` | the independent A-5 checker of Section 5.5, which shares no code with the rest |
| `declarations/` | the declaration of the Section 4.5 kill test, fixed before its grid was run |
| `lean/` | the Lean 4 project |
| `w2a_kill_*_rows.jsonl`, `w2a_kill_anchors_*.json` | the grid rows and anchors of the Section 4.5 kill test |
| `make_figs.py` | draws Figures 1 to 3 from the deposits into `figs/` |

## Which script produces what

| section | script | results file in `deposits/` |
| --- | --- | --- |
| 3.3, Table 2 | `cstar_spread_probe.py` | `cstar_spread_claims.txt` |
| 4 | the Lean files in `lean/`; `hkappa_check.py` (`t_helmert_check.py` prints the same Helmert check) | `formal01_hkappa_claims.txt` |
| 4.5 | `w2a_kill.py`; `formal01_t1_schur.py` | `w2a_kill.json`; `formal01_t1_claims.txt` |
| 5.1 to 5.4 | `formal01_f8c_certified.py`; `formal01_phi_engine.py`; `formal01_f8c_witness_dump.py`; `formal01_f8a_battery.py` | `formal01_f8c_claims-R3.txt`; `formal01_phi_claims.txt`; `f8c_cover_witness-R1.txt`; `formal01_p3a_claims.txt` |
| 5.5 | `formal01_f8b_ratecheck.py`; `cleanroom-a5/a5_run.py`; `formal01_kappa03_a5_agreement.py` | `formal01_f8b_claims.txt`; `cleanroom-a5/a5_claims.txt`; `a5_agreement_claims.txt` |
| 6.1 | `formal01_f8_indep_checker.py` | `f8_indep_claims-R1.txt` |
| 6.2 | `formal01_highprec.py`; `corner_check_20m.py` | `formal01_highprec_claims.txt`; `corner_check_20m_claims.txt` |
| 6.3 | `formal01_f8a_battery.py`; `formal01_f6_gridcheck.py` | `formal01_p3a_claims.txt`; `formal01_f6_claims.txt` |
| 6.4 | `l2_challenge_v5.py` | `l2_challenge_v5_claims.txt` |
| 6.5 | `formal01_ax_all.py`; `formal01_kappa03_statement_r2.py`; `statement_identity_r2.py`; `formal01_f8_constants_match2.py` | `ax_all_claims.txt`; `statement_identity_r2_claims.txt`; `f8_constants_match2_claims.txt` |
| 7.1, 7.2 | `paired_vs_welch_probe.py` | `paired_vs_welch_probe_claims.txt` |
| 7.3 | `u4_crossterm_probe.py` | `u4_crossterm_probe_claims.txt` |
| 7.4 | `l3_killtest.py`, then `l3_claims.py` | `l3_killtest.json`; `l3_killtest_claims.txt`; `l3_boundary_claims.txt` |
| 7.5 | `paired_vs_welch_probe.py` | `paired_vs_welch_probe_claims.txt` |
| 8.2 | `lean/W4lean/UnitPReduction.lean` | the Lean file |
| 8.3 | `p_law.py gates`; `p_verify.py C`; `p_verify_triage_f2.py`; `p_killround.py` | `p_law_gates_claims.txt`; `p_verify_C_claims.txt`; `p_verify_triage_F2_claims.txt`; `p_killround_claims.txt` |
| 8.4 to 8.6 | `p_killround.py`; `p_killround_k5.py`; `p_verify.py A`; `p_verify_t4k5.py` | `p_killround_claims.txt`; `p_killround_k5_claims.txt`; `p_verify_A_claims.txt`; `p_verify_A2_*_claims.txt` |
| 8.7 | `p_killround.py`; `p_verify.py B` | `p_killround_claims.txt`; `p_verify_B_claims.txt` |
| 8.8 | `p_uniform_bar.py`; `lean/W4lean/UnitPUniformBar.lean` | `p_uniform_bar_claims.txt`; the Lean file |
| 8.9 | `p_calib_k2.py gates`, `p_calib_k2.py sweep 12`, `p_calib_k2.py sweep 36`; `pb_certify.py` (uses `pb_cert.py`, `pb_joint.py`); `pb_checks.py`; `pb_facts.py`; `lean/W4lean/UnitPBStructure.lean` | `p_calib_k2_*_claims.txt`; `pb_certify_n*_claims.txt` and `pb_certify_n*_intervals.json`; `pb_checks_n*_claims.txt`; `pb_facts_claims.txt` |
| 9, Appendix D | `l1_selfnorm.py uq`, `large`, `anchors`, `table` | `l1_selfnorm_*_claims.txt` |
| Figures 1 to 3 | `make_figs.py` | reads the deposits above |

`w2a_kill.py` has three commands: `smoke` writes the anchors, `run --level primary` and `run --level h12` compute the grid, and `verdict` turns the two grids into `w2a_kill.json`.

## Run times, and what was not re-run

Every script in the first two rows below was run again for this repository, from this folder, after its comments were edited for publication and its file paths made relative. Each results file it wrote was then compared with the one from the original run, and they agree in every value except for the two cases noted after the table. Times are from that re-run on a 12-core desktop.

| time | scripts |
| --- | --- |
| under 1 minute | `cstar_spread_probe.py`, `formal01_t1_schur.py`, `formal01_phi_engine.py`, `cleanroom-a5/a5_run.py`, `formal01_kappa03_a5_agreement.py`, `formal01_f8_constants_match2.py`, `formal01_kappa03_statement_r2.py`, `statement_identity_r2.py`, `hkappa_check.py`, `corner_check_20m.py`, `paired_vs_welch_probe.py`, `u4_crossterm_probe.py`, `l3_killtest.py`, `l3_claims.py`, `p_verify.py C`, `p_verify_triage_f2.py`, `w2a_kill.py verdict` |
| 1 to 25 minutes | `formal01_f8b_ratecheck.py` (2), `formal01_f8c_certified.py` (7), `formal01_f8c_witness_dump.py` (2), `formal01_f8_indep_checker.py` (2), `formal01_f6_gridcheck.py` (10), `formal01_highprec.py` (6), `l2_challenge_v5.py` (3), `p_law.py gates` (2), `w2a_kill.py smoke` (4), `l1_selfnorm.py uq` (18) |
| Lean | `lake build` about 10 minutes once the Mathlib cache is in place; `formal01_ax_all.py` 7 minutes |

The two exceptions:

- `l3_killtest.py` seeds each design from Python's built-in `hash()` of a string, which Python randomises in every new process, so a re-run draws different samples. The results files from the original run, which the paper quotes, are kept. The re-run is in `l3_killtest_rerun_claims.txt` and `l3_boundary_rerun_claims.txt`: it finds the same eight breaches and the same boundary (only identical windows survive), with z-values that differ by Monte Carlo noise.
- `formal01_f8b_ratecheck.py` writes the A-5 enclosures of r1 and r2 with different exact endpoints from one run to the next. The intervals are about 1e-62 wide and overlap, and every value derived from them is identical.

The following were not run again, because each takes from twenty minutes to over an hour: `formal01_f8a_battery.py` (its F-1 sweep alone took 63 minutes), `p_killround.py` (20 minutes), `p_killround_k5.py` (32), `p_verify.py A` and `B` (B took 64), `p_verify_t4k5.py` (about 31 per correlation structure) and the grid of `w2a_kill.py run`. The `large`, `anchors` and `table` stages of `l1_selfnorm.py`, about four minutes together, were also not run again. Their results files, and the grid rows, are those of the original runs. Their comments were edited for publication afterwards; apart from file paths, their code was not. `formal01_p3a_claims.txt` records the sha256 of four scripts as they were when the battery ran, so those four hashes do not match the files here. The battery also cannot be run from this repository as published, because it reads inputs from an earlier stage of the work that are not included.

`p_uniform_bar.py`, added for Section 8.8, was written for this folder and run once from it, in about 4 minutes; its results file is that run.

The scripts of Section 8.9 were also written for this folder and run once from it, and their results files are those runs. `pb_certify.py N LEVEL [TAG ROUNDS CELLS DEPTH]` certifies one table: the 12-date rule (`pb_certify.py 12 0.0247`) took 22 minutes, and the 36-date rule of record (`pb_certify.py 36 0.02468 _r2a02468 20 12000000 4`) 56 minutes, on a 12-core desktop. The 36-date runs that did not certify are kept beside it: `_a02458` (certified, with the default budget), `_a02468` and the untagged run (not certified), and the runs at 0.02478 and 0.02488. `p_calib_k2.py sweep` takes about 40 minutes per sample size.

`statement_identity_r2.py` writes `statement_identity_scan_claims.txt`, a second, independent check of the same statement identity that `formal01_kappa03_statement_r2.py` records in `statement_identity_r2_claims.txt`.

## Lean

The Lean project is in [`lean/`](lean/). It is built against Mathlib, and Table B1 of the paper lists the declarations the paper relies on.

```bash
cd lean
lake exe cache get
lake build
cd ..
python formal01_ax_all.py
```

`formal01_ax_all.py` generates `lean/W4lean/AxAll.lean`, which prints the axioms of every theorem and lemma in the proof files, builds it, and checks that each depends only on `propext`, `Classical.choice` and `Quot.sound`. `Unit2BFormalR1R2.lean` is generated by `formal01_kappa03_statement_r2.py` from the statement file and the proof files; it carries all thirteen statements of the certificate layer with their proofs, and contains no `sorry`. `UnitPUniformBar.lean` and `UnitPBStructure.lean`, the Lean files of Sections 8.8 and 8.9, print the axioms of their results when they are built, and each depends only on those three.

`ChiSquaredSpikeR1.lean` is a modified copy of a file from the [StatLean](https://github.com/StatLean/Stat-Lean) project (Apache-2.0). Its header lists the changes, and the unmodified source files and their licence are in `deposits/statlean_ref_R1/`.
