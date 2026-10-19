import W4lean.Unit2BSizeConsR11
open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.SizeConsAsmR11

/-!
# 2B size bound, part 2b: the Phi step and `size_bound_conservative`

A size bound for the raw rate with NO numerical premise, NO engine and NO bridge.  It rests
only on the proved `kappa_reduction`, `Phi_antitone`, `r1 ≥ 1/2`, `r2 > 0`, and the existence
and attainment of `cCons`.  `cCons` is a STRICTLY WEAKER critical value than the certified
`Uc`; this file does not replace `size_bound_2B` and does not touch it.

Every project LEMMA is used fully qualified.
-/

open Unit2B.FormalR1.Nonvac

/-- The Phi step: on the Box, the rate at the true `(r1 θ, r2 θ)` is at most the rate at the
conservative corner `(1/2, 0)`. -/
theorem Phi_le_cons {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsFiniteMeasure P]
    {N X1 X2 : Ω → ℝ} (hX1 : ∀ᵐ ω ∂P, 0 ≤ X1 ω) (hX2 : ∀ᵐ ω ∂P, 0 ≤ X2 ω)
    {c : ℝ} (hc : 0 ≤ c) {θ : Param} (hθ : θ ∈ Box) :
    Phi P N X1 X2 c (r1 θ) (r2 θ) ≤ Phi P N X1 X2 c (1 / 2) 0 :=
  Unit2B.FormalR1.Proofs.Phi_antitone hX1 hX2 hc le_rfl
    (Unit2B.FormalR1.SizeConsR11.r1_ge_half hθ) (Unit2B.FormalR1.SizeConsR11.r2_pos hθ).le

/-- Part 2: the conservative size bound.  For every mean `m` and every `θ ∈ Box`, the raw
rejection rate at `cCons = sInf {c | 0 < c ∧ Phi P0 N0 X10 X20 c (1/2) 0 ≤ alpha}` is at
most `alpha`. -/
theorem size_bound_conservative (m : ℝ) :
    ∀ θ ∈ Box, Unit2B.FormalR1.Kappa.RawRate m Unit2B.FormalR1.SizeConsR11.cCons θ ≤ alpha := by
  obtain ⟨hNm, hX1m, hX2m, hN, hX1, hX2, hInd⟩ := Unit2B.FormalR1.SizeR11.canonical_model
  have hX1ae : ∀ᵐ ω ∂P0, 0 ≤ X10 ω :=
    Unit2B.FormalR1.Proofs.ae_nonneg_of_gammaLaw hX1m (by norm_num) (by norm_num) hX1
  have hX2ae : ∀ᵐ ω ∂P0, 0 ≤ X20 ω :=
    Unit2B.FormalR1.Proofs.ae_nonneg_of_gammaLaw hX2m (by norm_num) (by norm_num) hX2
  obtain ⟨hpos, hle⟩ := Unit2B.FormalR1.SizeConsR11.cCons_spec
  intro θ hθ
  rw [Unit2B.FormalR1.Kappa.ProofsR10.kappa_reduction hNm hX1m hX2m hN hX1 hX2 hInd m
    Unit2B.FormalR1.SizeConsR11.cCons hpos θ hθ]
  exact (Unit2B.FormalR1.SizeConsAsmR11.Phi_le_cons hX1ae hX2ae hpos.le hθ).trans hle

#print axioms Unit2B.FormalR1.SizeConsAsmR11.Phi_le_cons
#print axioms Unit2B.FormalR1.SizeConsAsmR11.size_bound_conservative

end Unit2B.FormalR1.SizeConsAsmR11
