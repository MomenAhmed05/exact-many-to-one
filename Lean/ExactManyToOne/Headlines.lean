/-
Copyright (c) 2026 Momen Ahmed. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Momen Ahmed
-/
import ExactManyToOne.Reduction
import ExactManyToOne.Envelope
import ExactManyToOne.Enclosure
import ExactManyToOne.Gaussian
import ExactManyToOne.Calibrated
import ExactManyToOne.Invariance
import ExactManyToOne.Dependence
import ExactManyToOne.Certificates

/-!
# The headline results, in the order the paper states them

Each entry names a claim of the paper and the theorem that proves it. The proofs
live in the section files; this file is the index, and `Audit.lean` checks that
every entry depends on no axiom beyond Lean's three standard ones.

What is proved and what is assumed:

* **Proved** — the reduction and its dimension (Section 3), the envelope bound and
  the coverage integral (9.2–9.3), the corner enclosure against the Gaussian law and
  the certified constant's validity (9.8, 9.11), the branch-and-bound bracket (38),
  the partition bound and its tail term (41), the three substitutions (10.4), the
  union step (42), the two-sided identity (43), the joint projection (44), the
  maximal invariant at any number of candidates and the group's solvability (10.8),
  the whitening identities and the factor sign (11).
* **Checked exactly from the deposits** — every certificate's arithmetic and label
  (Tables 15, 17, 19 and 20), in `Certificates.lean`.
* **Assumed, and named where it enters** — the Gaussian facts the model supplies
  (normal means, and each mean independent of its variance estimate), which are
  hypotheses of the theorems that use them; the numerical enclosures computed in
  Python, which enter as `hcert`; and the Hunt–Stein theorem, which enters as
  `HuntSteinConclusion` in one theorem only.
-/

namespace ExactManyToOne

/-! ## Section 3 -/

/-- **H1. The exact law (equation (16)).** Conditioning on both of the benchmark's
sufficient statistics makes the `K` comparisons exactly independent, so the joint
law is an average over the benchmark pair of a product of `K` inner functions. -/
alias H1_exact_law := exact_law

/-- **H1a. The outer rule is two-dimensional and a product (Section 3.3–3.4).** The
outer integral is over `ℝ × ℝ` for every `K`, and it factors into a normal against
a scaled chi-square. -/
alias H1a_outer_product := exact_law_iterated

/-- **H1b. The inner rule is one-dimensional (equation (15)).** -/
alias H1b_inner_one_dimensional := inner_gaussian

/-- **H1c. The general factorisation behind all of it.** -/
alias H1c_factorisation := prob_forall_eq_lintegral_prod

/-! ## Section 9 -/

/-- **H2. The envelope guarantee (equation (32)).** `FWER ≤ α + P(escape) ≤ α + γ`. -/
alias H2_envelope := envelope_guarantee

/-- **H2a. The escape probability is an exact one-dimensional integral (equation (29)).** -/
alias H2a_coverage_integral := coverage_eq

/-- **H2b. Box membership in the chi-square factors (equation (27) to (29)).** -/
alias H2b_box_membership := mem_box_iff

/-- **H3. The certified constant is valid at every ratio vector in the box (equations
(35) and (37)).** Proved against the Gaussian law itself. -/
alias H3_certified_constant := certified_constant_valid

/-- **H3a. The corner enclosure (equation (37)).** -/
alias H3a_corner_enclosure := prod_enclosure

/-- **H3b. Branch and bound brackets the supremum over the continuum (equation (38)).** -/
alias H3b_tiling := tiling_bound

/-- **H3c. A certified gap pins the supremum to the tolerance (Table 15).** -/
alias H3c_bracket := certificate_brackets

/-! ## Section 10 -/

/-- **H4. The calibrated rule is valid as soon as its level is feasible (equation (40)).** -/
alias H4_calibrated := calibrated_valid

/-- **H4a. The internal-level certificate: partition bound with tail (equation (41)).** -/
alias H4a_partition := partition_bound

/-- **H4b. The tail term is `4 p_min`.** -/
alias H4b_tail := tail_mass_le

/-- **H4c. The three substitutions of Section 10.4 move the rejection probability up.** -/
alias H4c_substitutions := three_substitutions

/-- **H5. Any number of candidates by the union step (equation (42)).** -/
alias H5_union_step := union_step

/-- **H5a. The Weierstrass inequality behind (42).** -/
alias H5a_weierstrass := one_sub_prod_le_sum

/-- **H6. Two-sided at one candidate is exactly twice one-sided (equation (43)).** -/
alias H6_two_sided_identity := two_sided_eq_two_mul

/-- **H6a. The null law of one comparison is symmetric.** -/
alias H6a_symmetry := tstat_symm

/-- **H6b. Two-sided certificate at any `K`, `FWER₂ₛ ≤ 2K · b`.** -/
alias H6b_two_sided_union := two_sided_union

/-- **H7. The joint rule is certified by projecting the deployed table (equation (44)).** -/
alias H7_joint_projection := joint_projection_bound

/-- **H8. The maximal invariant (Section 10.8), at any number of candidates.** -/
alias H8_maximal_invariant := Suff.maximal_invariant

/-- **H8a. The paper's one-candidate statement.** -/
alias H8a_maximal_invariant_one := Suff.maximal_invariant_one

/-- **H8b. The invariant class is the class of functions of `(t, λ̂)`.** -/
alias H8b_invariant_factors := Suff.invariant_factors

/-- **H8c. The floor extends to all tests, given Hunt–Stein (Section 10.8).** -/
alias H8c_floor_all_tests := Suff.floor_all_tests

/-! ## Section 11 -/

/-- **H9. Whitening a declared AR(1) with the first observation dropped is exact.** -/
alias H9_whitening_exact := whiten_iIndep

/-- **H9a. Rescaling the first observation biases unequal lengths (the falsifier of 11.1).** -/
alias H9a_rescaling_fails := pw_mean_diff_ne_zero

/-- **H9c. Group-specific declared coefficients need the `1/(1 − ϕ)` rescaling.** -/
alias H9c_whiten_rescaled := whiten_rescaled

/-- **H9b. Estimating the coefficient per group leaves `μ(ϕ̂₀ − ϕ̂ᵢ)` (equation (49)).** -/
alias H9b_estimated_term := estimated_mean_term

/-- **H10. Opposite-sign loadings make the independent-groups rule liberal (equation (50)).** -/
alias H10_factor_liberal := opposite_sign_liberal

/-- **H10a. Same-sign loadings make it conservative.** -/
alias H10a_factor_conservative := same_sign_conservative

end ExactManyToOne
