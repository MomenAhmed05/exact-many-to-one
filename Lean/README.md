# Lean: the paper's results in Lean 4

A machine-checked version of the mathematics in
*Exact Many-to-One Comparison Under Estimated Volatilities*. Every theorem here is
checked by Lean's kernel against Mathlib. None uses `sorry`, and none depends on any
axiom beyond Lean's standard three (`propext`, `Classical.choice`, `Quot.sound`).

## Building and auditing

With [elan](https://github.com/leanprover/elan) installed:

```
lake exe cache get      # Mathlib's prebuilt files, once
lake build
python scripts/audit.py
```

The audit checks four things and fails if any of them does not hold:

1. `Certificates.lean` regenerates from the JSON deposits byte for byte.
2. `lake build` finishes with no errors and no warnings.
3. No source file contains `sorry`, `admit`, `axiom` or `native_decide`.
4. `#print axioms` shows only the standard three for every headline theorem.

Toolchain: Lean `v4.35.0-rc3`, Mathlib at the matching tag, both pinned in
`lean-toolchain` and `lake-manifest.json`.

## What each paper result rests on

| Paper | Claim | Lean | File |
| --- | --- | --- | --- |
| (16) | The exact law: the joint probability is an average over the benchmark pair of a product of `K` inner functions | `exact_law` | Reduction |
| §3.3–3.4 | The outer rule is over `ℝ × ℝ` for every `K`, and it factors into a normal against a scaled chi-square | `exact_law_iterated` | Reduction |
| (15) | The inner rule is a one-dimensional integral of the normal distribution function | `inner_gaussian` | Reduction |
| — | The general factorisation behind (16), (29) and (42), for any laws | `prob_forall_eq_lintegral_prod` | Reduction |
| (32) | Envelope rule: `FWER ≤ α + P(escape) ≤ α + γ` | `envelope_guarantee` | Envelope |
| (27)–(29) | Coverage of the box is an exact one-dimensional integral | `mem_box_iff`, `coverage_eq` | Envelope |
| (37) | Corner substitution and the product enclosure | `corner_le`, `prod_enclosure` | Enclosure |
| (35), (37) | The certified constant holds at every variance vector in the box, against the Gaussian law itself | `certified_constant_valid` | Gaussian |
| §9.11 | The direction: a lower bound on the distribution function bounds the constant from above | `reject_le_of_accept_ge`, `cpop_le_crig` | Enclosure |
| (38) | Branch and bound: the leaf maximum dominates the region; a certified gap pins the supremum | `tiling_bound`, `certificate_brackets` | Enclosure |
| (40) | The calibrated rule is valid once its level is feasible; maximality is not needed | `calibrated_valid` | Calibrated |
| (41) | Partition bound with a tail term, and the tail is at most `4 p_min` | `partition_bound`, `tail_mass_le` | Enclosure |
| §10.4 | The three substitutions each move the rejection probability up | `three_substitutions` | Enclosure |
| (42) | Union step: `FWER ≤ K · sup FWER⁽¹⁾` | `union_step`, `one_sub_prod_le_sum` | Calibrated |
| (43) | Two-sided at one candidate is exactly twice one-sided | `two_sided_eq_two_mul`, `tstat_symm`, `gaussian_pair_symm` | Calibrated |
| (43) | Two-sided certificate: `FWER₂ₛ ≤ 2K · b` | `two_sided_union` | Calibrated |
| (44) | The joint rule is certified by projecting its table | `joint_projection_bound` | Calibrated |
| (47) | The affine group is solvable | `Aff.instIsSolvable` | Invariance |
| §10.8 | `(t, λ̂)` is a maximal invariant | `Suff.maximal_invariant_one` | Invariance |
| §10.8 | The invariant class is exactly the functions of `(t, λ̂)` | `Suff.invariant_factors` | Invariance |
| §10.8 | The floor extends from invariant tests to all tests, given Hunt–Stein | `Suff.floor_all_tests` | Invariance |
| (11) | Numerator correlation is the product of the loadings | `corr_eq_loading_mul` | Dependence |
| (45), §11.1 | Whitening with the first observation dropped is exact | `whiten_eq`, `whiten_iIndep` | Dependence |
| §11.1 | Rescaling the first observation biases unequal lengths | `pw_mean_diff_ne_zero` | Dependence |
| §11.3 | Group-specific declared coefficients need each group divided by `1 − ϕ`; the dropped-first transform alone needs a shared coefficient | `whiten_rescaled` | Dependence |
| (49) | Per-group coefficients leave `μ(ϕ̂₀ − ϕ̂ᵢ)`; pooling removes it | `estimated_mean_term`, `pooled_mean_term` | Dependence |
| (50) | The factor sign decides liberal or conservative | `den_sub_num`, `opposite_sign_liberal`, `same_sign_conservative` | Dependence |
| Tables 15, 17, 19, 20 | Every deposited certificate supports the claim made from it | `table15_*`, `table17_*`, `table19_*`, `table20_*` | Certificates |

`Headlines.lean` lists the same results as `H1` to `H10` in the paper's order.

## One result goes further than the paper

The paper proves the maximal invariant at one candidate and records the extension
above one candidate as unwritten. `Suff.maximal_invariant` proves it at any number
of candidates `K`: two positive sufficient statistics that share all `K` statistics
and all `K` estimated ratios are related by the group element `b = √(q₀'/q₀)`,
`a = x₀' − b x₀`. The one-candidate statement is the corollary
`Suff.maximal_invariant_one`. This settles the invariance half of the extension. The
Hunt–Stein step above one candidate also needs the loss to be invariant there, and
that is not formalised.

## What is assumed, and where it enters

Nothing is assumed silently. Each assumption is a hypothesis of the theorem that
needs it.

- **The Gaussian model.** The means are normal, and within each sample the mean is
  independent of the variance estimate. These are hypotheses (`hlaw`, `hind`,
  `hBY`, `hYY`), not theorems. The reduction itself (`exact_law`) needs only the
  independence, not normality.
- **The numerical enclosures.** The interval-arithmetic certificates are computed in
  Python and enter as the hypothesis `hcert`. Lean checks exactly that each deposited
  number supports the claim the paper makes from it, and that each pass or fail label
  is correct. It does not recompute the enclosures.
- **The marginal identity in the union step.** Each candidate's marginal rejection
  probability is the one-candidate error rate at that candidate's ratio (`hmarg`).
- **Hunt–Stein.** Quoted from the literature as the hypothesis `HuntSteinConclusion`,
  in `floor_all_tests` only, as the paper does. Amenability of solvable groups is also
  quoted. Mathlib has neither result.

`calibrated_valid` restates its hypothesis on purpose. The paper's point in (40) is
that validity needs a feasible level and nothing more, and the theorem records
exactly that.

## Layout

```
ExactManyToOne/      the library, imported as ExactManyToOne
  Reduction.lean     Section 3
  Envelope.lean      Sections 9.2-9.3
  Enclosure.lean     Sections 9.8, 9.11, 10.4
  Gaussian.lean      the enclosure against the normal law
  Calibrated.lean    Sections 10.2, 10.5, 10.6
  Invariance.lean    Section 10.8
  Dependence.lean    Sections 3.2, 9.10, 11
  Certificates.lean  generated from the deposits
  Headlines.lean     the index
Audit.lean           #print axioms for every headline
scripts/
  gen_certificates.py
  audit.py
```
