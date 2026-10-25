import Mathlib

/-!
# Unit P: a fixed critical value valid at every correlation is Bonferroni's (two candidates)

Two paired statistics `t₁, t₂` share a family of laws indexed by the correlation `ρ` of their
difference series. Under the null each statistic has the same marginal law at every `ρ` (Student's
t on `n − 1` degrees of freedom; this enters as the hypotheses `hmarg₁`, `hmarg₂`, it is not proved
here). Two facts then settle which fixed critical value is valid at every correlation.

* Boole's inequality: at every `ρ`, `P(max(t₁, t₂) > c) ≤ P(t₁ > c) + P(t₂ > c)`.
* At `ρ = −1` the second difference series is a negative multiple of the first, so `t₂ = −t₁` on
  every sample (`pairedT_neg_smul`). For `c ≥ 0` the two exceedance events are then disjoint and
  Boole's bound holds with equality.

Hence a fixed `c ≥ 0` keeps the family-wise error rate at most `α` for every correlation exactly
when twice the marginal tail is at most `α` (`uniform_bar_iff`), which is Bonferroni's critical
value.
-/

open MeasureTheory Set
open scoped ENNReal

namespace UnitPUniform

/-! ## Algebra: a negatively scaled series has the negated paired statistic -/

/-- Sample mean of `n` observations. -/
noncomputable def mean {n : ℕ} (d : Fin n → ℝ) : ℝ := (∑ t, d t) / n

/-- Sample standard deviation, divisor `n − 1`. -/
noncomputable def sd {n : ℕ} (d : Fin n → ℝ) : ℝ :=
  Real.sqrt ((∑ t, (d t - mean d) ^ 2) / ((n : ℝ) - 1))

/-- The paired (one-sample) t statistic of a difference series. -/
noncomputable def pairedT {n : ℕ} (d : Fin n → ℝ) : ℝ :=
  Real.sqrt n * mean d / sd d

lemma mean_smul {n : ℕ} (k : ℝ) (d : Fin n → ℝ) :
    mean (fun t => k * d t) = k * mean d := by
  unfold mean
  rw [← Finset.mul_sum]
  ring

lemma sd_smul {n : ℕ} (k : ℝ) (d : Fin n → ℝ) :
    sd (fun t => k * d t) = |k| * sd d := by
  unfold sd
  rw [mean_smul]
  have h : ∀ t, (k * d t - k * mean d) ^ 2 = k ^ 2 * (d t - mean d) ^ 2 := fun t => by ring
  simp_rw [h, ← Finset.mul_sum, mul_div_assoc]
  rw [Real.sqrt_mul (sq_nonneg k), Real.sqrt_sq_eq_abs]

/-- A difference series that is a negative multiple of another has the negated paired statistic.
This is the pathwise content of correlation `−1`. -/
theorem pairedT_neg_smul {n : ℕ} {k : ℝ} (hk : 0 < k) (d : Fin n → ℝ) :
    pairedT (fun t => -k * d t) = -pairedT d := by
  unfold pairedT
  rw [mean_smul, sd_smul, abs_neg, abs_of_pos hk]
  by_cases hs : sd d = 0
  · simp [hs]
  · field_simp

/-! ## Probability: Boole's bound, and when it is attained -/

variable {Ω : Type*} [MeasurableSpace Ω]

omit [MeasurableSpace Ω] in
lemma max_gt_eq_union (T₁ T₂ : Ω → ℝ) (c : ℝ) :
    {ω | c < max (T₁ ω) (T₂ ω)} = {ω | c < T₁ ω} ∪ {ω | c < T₂ ω} := by
  ext ω
  simp [lt_max_iff]

/-- Boole's inequality for the larger of two statistics. -/
theorem fwer_le_sum (μ : Measure Ω) (T₁ T₂ : Ω → ℝ) (c : ℝ) :
    μ {ω | c < max (T₁ ω) (T₂ ω)} ≤ μ {ω | c < T₁ ω} + μ {ω | c < T₂ ω} := by
  rw [max_gt_eq_union]
  exact measure_union_le _ _

/-- If `T₂ = −T₁` almost surely and `c ≥ 0`, the two exceedance events cannot both occur, so
Boole's bound holds with equality. -/
theorem fwer_eq_sum_of_antithetic (μ : Measure Ω) {T₁ T₂ : Ω → ℝ} (hT₂ : Measurable T₂)
    {c : ℝ} (hc : 0 ≤ c) (hanti : ∀ᵐ ω ∂μ, T₂ ω = -T₁ ω) :
    μ {ω | c < max (T₁ ω) (T₂ ω)} = μ {ω | c < T₁ ω} + μ {ω | c < T₂ ω} := by
  rw [max_gt_eq_union, ← measure_union_add_inter _ (measurableSet_lt measurable_const hT₂)]
  have hnull : μ ({ω | c < T₁ ω} ∩ {ω | c < T₂ ω}) = 0 := by
    refine measure_mono_null ?_ (ae_iff.mp hanti)
    intro ω hω heq
    simp only [mem_inter_iff, mem_ofPred_eq] at hω
    rw [heq] at hω
    linarith [hω.1, hω.2]
  rw [hnull, add_zero]

/-- **The uniform critical value.** Let a family of laws `μ ρ` share the marginal upper tail `G`
of each statistic at `c`, and let one member `ρ₀` make the statistics antithetic. Then a fixed
`c ≥ 0` keeps the family-wise error rate at most `α` for every member exactly when `2·G ≤ α`. -/
theorem uniform_bar_iff {ι : Type*} (μ : ι → Measure Ω) {T₁ T₂ : Ω → ℝ} (hT₂ : Measurable T₂)
    {c : ℝ} (hc : 0 ≤ c) (G α : ℝ≥0∞)
    (hmarg₁ : ∀ ρ, μ ρ {ω | c < T₁ ω} = G) (hmarg₂ : ∀ ρ, μ ρ {ω | c < T₂ ω} = G)
    (ρ₀ : ι) (hanti : ∀ᵐ ω ∂(μ ρ₀), T₂ ω = -T₁ ω) :
    (∀ ρ, μ ρ {ω | c < max (T₁ ω) (T₂ ω)} ≤ α) ↔ 2 * G ≤ α := by
  constructor
  · intro h
    have h0 := h ρ₀
    rwa [fwer_eq_sum_of_antithetic (μ ρ₀) hT₂ hc hanti, hmarg₁, hmarg₂, ← two_mul] at h0
  · intro h ρ
    calc μ ρ {ω | c < max (T₁ ω) (T₂ ω)}
        ≤ μ ρ {ω | c < T₁ ω} + μ ρ {ω | c < T₂ ω} := fwer_le_sum (μ ρ) T₁ T₂ c
      _ = 2 * G := by rw [hmarg₁, hmarg₂, two_mul]
      _ ≤ α := h

/-- At correlation `−1` the second difference series is `−k` times the first on every sample, so
the two paired statistics are antithetic: the hypothesis `hanti` of `uniform_bar_iff` for the
paired statistics themselves. -/
theorem antithetic_of_scaled {n : ℕ} (μ : Measure Ω) (D₁ D₂ : Ω → (Fin n → ℝ)) {k : ℝ}
    (hk : 0 < k) (hD : ∀ᵐ ω ∂μ, D₂ ω = fun t => -k * D₁ ω t) :
    ∀ᵐ ω ∂μ, pairedT (D₂ ω) = -pairedT (D₁ ω) := by
  filter_upwards [hD] with ω hω
  rw [hω, pairedT_neg_smul hk]

#print axioms UnitPUniform.pairedT_neg_smul
#print axioms UnitPUniform.fwer_le_sum
#print axioms UnitPUniform.fwer_eq_sum_of_antithetic
#print axioms UnitPUniform.uniform_bar_iff
#print axioms UnitPUniform.antithetic_of_scaled

end UnitPUniform
