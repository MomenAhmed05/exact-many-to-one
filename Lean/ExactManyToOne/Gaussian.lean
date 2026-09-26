/-
Copyright (c) 2026 Momen Ahmed. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Momen Ahmed
-/
import ExactManyToOne.Enclosure

/-!
# Sections 9.8 and 9.11 against the Gaussian law itself

`Enclosure.lean` proves the corner substitution for any monotone `Φ`. This file
evaluates the conditional law of Section 9.11 with the actual normal distribution
function and closes the loop: conditional on the variance estimates, the
probability that no candidate is rejected is the integral of equation (37), and
the lower corner envelope bounds it from below at **every** variance vector in the
box. A constant at which that envelope reaches `1 - α` therefore has conditional
rejection probability at most `α` everywhere in the box, which is the statement the
proven tier (35) and the continuum certificate (38) supply.

Scale: the benchmark variance is set to one, which loses nothing because the
statistic is scale invariant (`Suff.t_smul`). Candidate `i`'s mean then has variance
`wᵢ = λᵢ/nᵢ` and the benchmark's `w₀ = 1/n₀`, so `1/√wᵢ = √(nᵢ/λᵢ)` is the paper's
`sᵢ` and `√w₀ z + eᵢ` is its `xᵢ(z)` with `eᵢ = c dᵢ`.
-/

open MeasureTheory ProbabilityTheory Set
open scoped ENNReal NNReal

namespace ExactManyToOne

/-- The standard normal distribution function. -/
noncomputable def stdNormalCdf : ℝ → ℝ := cdf (gaussianReal 0 1)

lemma stdNormalCdf_mono : Monotone stdNormalCdf := monotone_cdf _

lemma stdNormalCdf_nonneg (x : ℝ) : 0 ≤ stdNormalCdf x := cdf_nonneg _ x

/-- A centred normal with variance `v` is the image of the standard one under
`x ↦ √v x`. -/
lemma gaussianReal_eq_map (v : ℝ≥0) :
    gaussianReal 0 v = (gaussianReal 0 1).map (fun x => Real.sqrt v * x) := by
  rw [gaussianReal_map_const_mul, mul_zero]
  congr 1
  ext
  simp

/-- **The normal distribution function in standard form**: `P(N(0, v) ≤ x) = Φ(x/√v)`. -/
lemma gaussian_measure_Iic (v : ℝ≥0) (hv : 0 < v) (x : ℝ) :
    gaussianReal 0 v (Iic x) = ENNReal.ofReal (stdNormalCdf (x / Real.sqrt v)) := by
  have hs : 0 < Real.sqrt v := Real.sqrt_pos.2 (by exact_mod_cast hv)
  rw [gaussianReal_eq_map v, Measure.map_apply (by fun_prop) measurableSet_Iic, stdNormalCdf,
    ofReal_cdf]
  congr 1
  ext y
  simp only [mem_preimage, mem_Iic]
  rw [le_div_iff₀ hs, mul_comm]

/-- **The conditional law of Section 9.11.** Conditional on the variance estimates,
the only randomness left is in the centred means, which are independent normals.
The probability that every candidate is accepted is then the integral, against the
standard normal, of the product of `K` normal distribution functions. -/
theorem conditional_law {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsProbabilityMeasure P] {K : ℕ}
    (Z : Fin K → Ω → ℝ) (U : Ω → ℝ) (hZ : ∀ i, Measurable (Z i)) (hU : Measurable U)
    (hUZ : IndepFun U (fun ω i => Z i ω) P) (hZZ : iIndepFun Z P)
    (w : Fin K → ℝ≥0) (hw : ∀ i, 0 < w i) (w0 : ℝ≥0)
    (hZlaw : ∀ i, P.map (Z i) = gaussianReal 0 (w i)) (hUlaw : P.map U = gaussianReal 0 w0)
    (e : Fin K → ℝ) :
    P {ω | ∀ i, Z i ω ≤ U ω + e i}
      = ∫⁻ z, ENNReal.ofReal
          (∏ i, stdNormalCdf ((Real.sqrt w0 * z + e i) / Real.sqrt (w i))) ∂(gaussianReal 0 1) := by
  have hSm : ∀ i, MeasurableSet {p : ℝ × ℝ | p.1 ≤ p.2 + e i} :=
    fun i => measurableSet_le (by fun_prop) (by fun_prop)
  have h := prob_forall_eq_lintegral_prod U Z hU hZ hUZ hZZ
    (fun i => {p : ℝ × ℝ | p.1 ≤ p.2 + e i}) hSm
  have hf : Measurable fun b : ℝ =>
      ∏ i, (P.map (Z i)) {y | (y, b) ∈ {p : ℝ × ℝ | p.1 ≤ p.2 + e i}} :=
    Finset.measurable_prod _ fun i _ => measurable_measure_prodMk_right (hSm i)
  rw [show {ω | ∀ i, Z i ω ≤ U ω + e i}
      = {ω | ∀ i, (Z i ω, U ω) ∈ {p : ℝ × ℝ | p.1 ≤ p.2 + e i}} from rfl, h,
    hUlaw, gaussianReal_eq_map w0, lintegral_map hf (by fun_prop)]
  refine lintegral_congr fun z => ?_
  rw [ENNReal.ofReal_prod_of_nonneg (fun i _ => stdNormalCdf_nonneg _)]
  refine Finset.prod_congr rfl fun i _ => ?_
  rw [hZlaw i, ← gaussian_measure_Iic (w i) (hw i)]
  rfl

/-- **The corner envelope bounds the conditional law at every variance vector in the
box.** `wm ≤ w ≤ wp` componentwise with `wm > 0`. The scales `1/√wᵢ` then lie in
`[1/√wpᵢ, 1/√wmᵢ]`, and the lower corner of (37) evaluated there is a lower bound on
the acceptance probability that does not depend on `w`. -/
theorem conditional_lower_envelope {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsProbabilityMeasure P] {K : ℕ}
    (Z : Fin K → Ω → ℝ) (U : Ω → ℝ) (hZ : ∀ i, Measurable (Z i)) (hU : Measurable U)
    (hUZ : IndepFun U (fun ω i => Z i ω) P) (hZZ : iIndepFun Z P)
    (w wm wp : Fin K → ℝ≥0) (hwm : ∀ i, 0 < wm i) (hbox : ∀ i, w i ∈ Icc (wm i) (wp i))
    (w0 : ℝ≥0) (hZlaw : ∀ i, P.map (Z i) = gaussianReal 0 (w i))
    (hUlaw : P.map U = gaussianReal 0 w0) (e : Fin K → ℝ) :
    ∫⁻ z, ENNReal.ofReal (∏ i, stdNormalCdf
        (lowCorner (1 / Real.sqrt (wp i)) (1 / Real.sqrt (wm i)) (Real.sqrt w0 * z + e i)
          * (Real.sqrt w0 * z + e i))) ∂(gaussianReal 0 1)
      ≤ P {ω | ∀ i, Z i ω ≤ U ω + e i} := by
  have hw : ∀ i, 0 < w i := fun i => lt_of_lt_of_le (hwm i) (hbox i).1
  rw [conditional_law Z U hZ hU hUZ hZZ w hw w0 hZlaw hUlaw e]
  have hs : ∀ i, 1 / Real.sqrt (w i) ∈ Icc (1 / Real.sqrt (wp i)) (1 / Real.sqrt (wm i)) := by
    intro i
    have h1 : 0 < Real.sqrt (wm i) := Real.sqrt_pos.2 (by exact_mod_cast hwm i)
    have h2 : 0 < Real.sqrt (w i) := Real.sqrt_pos.2 (by exact_mod_cast hw i)
    constructor
    · exact one_div_le_one_div_of_le h2 (Real.sqrt_le_sqrt (by exact_mod_cast (hbox i).2))
    · exact one_div_le_one_div_of_le h1 (Real.sqrt_le_sqrt (by exact_mod_cast (hbox i).1))
  have := lintegral_lower_envelope stdNormalCdf_mono stdNormalCdf_nonneg (gaussianReal 0 1)
    (fun i => 1 / Real.sqrt (w i)) (fun i => 1 / Real.sqrt (wp i)) (fun i => 1 / Real.sqrt (wm i))
    hs (fun i z => Real.sqrt w0 * z + e i)
  refine this.trans (le_of_eq (lintegral_congr fun z => ?_))
  congr 1
  refine Finset.prod_congr rfl fun i _ => ?_
  rw [one_div_mul_eq_div]

/-- **Theorem (the proven tier and the continuum certificate, (35) and (38)).** If the
corner envelope reaches `1 - α` at the constant, then at every variance vector in the
box the conditional probability that some candidate is rejected is at most `α`. -/
theorem certified_constant_valid {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsProbabilityMeasure P] {K : ℕ}
    (Z : Fin K → Ω → ℝ) (U : Ω → ℝ) (hZ : ∀ i, Measurable (Z i)) (hU : Measurable U)
    (hUZ : IndepFun U (fun ω i => Z i ω) P) (hZZ : iIndepFun Z P)
    (w wm wp : Fin K → ℝ≥0) (hwm : ∀ i, 0 < wm i) (hbox : ∀ i, w i ∈ Icc (wm i) (wp i))
    (w0 : ℝ≥0) (hZlaw : ∀ i, P.map (Z i) = gaussianReal 0 (w i))
    (hUlaw : P.map U = gaussianReal 0 w0) (e : Fin K → ℝ) {α : ℝ≥0∞} (hα : α ≤ 1)
    (hcert : 1 - α ≤ ∫⁻ z, ENNReal.ofReal (∏ i, stdNormalCdf
        (lowCorner (1 / Real.sqrt (wp i)) (1 / Real.sqrt (wm i)) (Real.sqrt w0 * z + e i)
          * (Real.sqrt w0 * z + e i))) ∂(gaussianReal 0 1)) :
    P {ω | ∃ i, U ω + e i < Z i ω} ≤ α := by
  have hmeas : MeasurableSet {ω | ∀ i, Z i ω ≤ U ω + e i} := by
    simp only [ofPred_forall]
    exact MeasurableSet.iInter fun i => measurableSet_le (hZ i) (hU.add measurable_const)
  have hcompl : {ω | ∃ i, U ω + e i < Z i ω} = {ω | ∀ i, Z i ω ≤ U ω + e i}ᶜ := by
    ext ω; simp [not_le]
  rw [hcompl]
  exact reject_le_of_accept_ge P hmeas hα
    (hcert.trans (conditional_lower_envelope Z U hZ hU hUZ hZZ w wm wp hwm hbox w0 hZlaw hUlaw e))

end ExactManyToOne
