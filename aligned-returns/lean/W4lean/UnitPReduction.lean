import Mathlib

/-!
# Unit P, T1: the reduction for the paired many-to-one max-t

Component T1 of the paired verification (Section 8). On aligned dates the K difference series give
`Z = √n · D̄` and the diagonal `a` of `A = Σ (Dₜ − D̄)(Dₜ − D̄)ᵀ`; each paired statistic is
`tᵢ = Zᵢ / sᵢ` with its OWN `sᵢ = √(aᵢ / ν)`. If `Z` is independent of `a` (ASM-2: the sample mean
and the sample covariance of a normal sample are independent), then the probability that every
statistic is at most `c` is an average over the law of `a` alone of the probability that `Z`
lies below a box whose corner depends on `a`. No property of the laws is used beyond the
independence.
-/

open MeasureTheory ProbabilityTheory Set
open scoped ENNReal

namespace UnitP

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P] {K : ℕ}

/-- The event of a statistic at most `c`, rewritten in the centred mean (equation (5) of the
paired setting): for a positive standard deviation, `z / s ≤ c ↔ z ≤ c · s`. -/
lemma div_le_iff_mul {z s c : ℝ} (hs : 0 < s) : z / s ≤ c ↔ z ≤ c * s := by
  rw [div_le_iff₀ hs]

/-- The box whose corner the variance estimates set. -/
def box (ν c : ℝ) (a : Fin K → ℝ) : Set (Fin K → ℝ) :=
  {z | ∀ i, z i ≤ c * Real.sqrt (a i / ν)}

lemma measurableSet_pairs (ν c : ℝ) :
    MeasurableSet {p : (Fin K → ℝ) × (Fin K → ℝ) | ∀ i, p.1 i ≤ c * Real.sqrt (p.2 i / ν)} := by
  have : {p : (Fin K → ℝ) × (Fin K → ℝ) | ∀ i, p.1 i ≤ c * Real.sqrt (p.2 i / ν)}
      = ⋂ i, {p | p.1 i ≤ c * Real.sqrt (p.2 i / ν)} := by
    ext p; simp
  rw [this]
  refine MeasurableSet.iInter fun i => measurableSet_le ?_ ?_
  · exact (measurable_pi_apply i).comp measurable_fst
  · exact measurable_const.mul (Real.continuous_sqrt.measurable.comp (by fun_prop))

/-- **T1 (the reduction).** With `Z` independent of the variance diagonal `a`, and every
standard deviation positive, the joint probability that all K paired statistics are at most
`c` equals the `a`-average of the probability that `Z` lies in the box set by `a`. -/
theorem paired_reduction (Z a : Ω → (Fin K → ℝ)) (hZ : Measurable Z) (ha : Measurable a)
    (hind : IndepFun Z a P) (ν c : ℝ) (hpos : ∀ ω i, 0 < Real.sqrt (a ω i / ν)) :
    P {ω | ∀ i, Z ω i / Real.sqrt (a ω i / ν) ≤ c}
      = ∫⁻ x, (P.map Z) (box ν c x) ∂(P.map a) := by
  have hev : {ω | ∀ i, Z ω i / Real.sqrt (a ω i / ν) ≤ c}
      = (fun ω => (Z ω, a ω)) ⁻¹'
          {p : (Fin K → ℝ) × (Fin K → ℝ) | ∀ i, p.1 i ≤ c * Real.sqrt (p.2 i / ν)} := by
    ext ω
    simp only [mem_ofPred_eq, mem_preimage]
    exact forall_congr' fun i => div_le_iff_mul (hpos ω i)
  have hjoint : P.map (fun ω => (Z ω, a ω)) = (P.map Z).prod (P.map a) :=
    (indepFun_iff_map_prod_eq_prod_map_map hZ.aemeasurable ha.aemeasurable).1 hind
  rw [hev, ← Measure.map_apply (hZ.prodMk ha) (measurableSet_pairs ν c), hjoint,
    Measure.prod_apply_symm (measurableSet_pairs ν c)]
  rfl

#print axioms UnitP.paired_reduction

end UnitP
