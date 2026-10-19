import W4lean.Unit2BFormalR1Proofs
import W4lean.Unit2BFormalNonvacR1
import W4lean.Unit2BKappaProofsR10
open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.SizeR11

/-!
# 2B size bound, part 1: `kappa_reduction` composed with the cover

Every project LEMMA below is used fully qualified: `Unit2B.FormalR1.Proofs.*` for the FORMAL-01
proofs and `Unit2B.FormalR1.Kappa.ProofsR10.*` for the kappa port.  A bare lemma name inside
`Unit2B.FormalR1.*` resolves to the statement file, whose bodies are not proofs.

The composition goes through `cover_sound`, NOT `comp_thm_2B`, so the lower witness `hWit`
(A-5) is not a hypothesis of anything in this file.
-/

open Unit2B.FormalR1.Nonvac

/-- Non-vacuity (i), PROMOTED to a named lemma from the first `example` of
`Unit2BFormalNonvacR1.lean` (its lines 37-50); the proof text is that example's, verbatim.
The canonical product model `P0 = Measure.pi ![N(0,1), Γ(35/2, 1/2), Γ(35/2, 1/2)]` with the
coordinate maps satisfies every hypothesis of `kappa_reduction`; `iIndepFun_pi` supplies the
three-way independence. -/
theorem canonical_model : Measurable N0 ∧ Measurable X10 ∧ Measurable X20 ∧
    P0.map N0 = gaussianReal 0 1 ∧ P0.map X10 = gammaMeasure (35 / 2) (1 / 2) ∧
    P0.map X20 = gammaMeasure (35 / 2) (1 / 2) ∧ iIndepFun ![N0, X10, X20] P0 := by
  have hmap : ∀ i : Fin 3, P0.map (Function.eval i) = laws i := fun i => by
    show (Measure.pi laws).map (Function.eval i) = laws i
    rw [Measure.pi_map_eval]
    simp only [measure_univ, Finset.prod_const_one, one_smul]
  have hind := iIndepFun_pi (Ω := fun _ => ℝ) (μ := laws) (𝓧 := fun _ => ℝ)
    (X := fun _ => id) (fun _ => aemeasurable_id)
  have hfun : (![N0, X10, X20] : Fin 3 → (Fin 3 → ℝ) → ℝ) = fun i ω => id (ω i) := by
    funext i; fin_cases i <;> rfl
  refine ⟨measurable_pi_apply 0, measurable_pi_apply 1, measurable_pi_apply 2,
    hmap 0, hmap 1, hmap 2, ?_⟩
  rw [hfun]; exact hind

/-- Non-vacuity (ii), RE-TARGETED at `kappa_reduction`: on the canonical model
the raw rejection rate IS `Phi` at `(r1 θ, r2 θ)`.  This is the non-reflexive instance of
`hKappa` with `R := RawRate m`, discharged by the proved reduction. -/
theorem kappa_reduction_canonical (m : ℝ) :
    ∀ c, 0 < c → ∀ θ ∈ Box,
      Unit2B.FormalR1.Kappa.RawRate m c θ = Phi P0 N0 X10 X20 c (r1 θ) (r2 θ) := by
  obtain ⟨hNm, hX1m, hX2m, hN, hX1, hX2, hInd⟩ := Unit2B.FormalR1.SizeR11.canonical_model
  exact Unit2B.FormalR1.Kappa.ProofsR10.kappa_reduction hNm hX1m hX2m hN hX1 hX2 hInd m

/-- Non-vacuity (ii), existential form: the hypotheses of `kappa_reduction` are jointly
satisfiable, and on a model satisfying them the reduction's conclusion holds for every `m`. -/
theorem kappa_reduction_nonvacuous :
    ∃ (Ω : Type) (_ : MeasurableSpace Ω) (P : Measure Ω) (_ : IsProbabilityMeasure P)
      (N X1 X2 : Ω → ℝ), Measurable N ∧ Measurable X1 ∧ Measurable X2 ∧
      P.map N = gaussianReal 0 1 ∧ P.map X1 = gammaMeasure (35 / 2) (1 / 2) ∧
      P.map X2 = gammaMeasure (35 / 2) (1 / 2) ∧ iIndepFun ![N, X1, X2] P ∧
      ∀ m c, 0 < c → ∀ θ ∈ Box,
        Unit2B.FormalR1.Kappa.RawRate m c θ = Phi P N X1 X2 c (r1 θ) (r2 θ) := by
  obtain ⟨hNm, hX1m, hX2m, hN, hX1, hX2, hInd⟩ := Unit2B.FormalR1.SizeR11.canonical_model
  exact ⟨Fin 3 → ℝ, inferInstance, P0, inferInstance, N0, X10, X20,
    hNm, hX1m, hX2m, hN, hX1, hX2, hInd,
    fun m => Unit2B.FormalR1.SizeR11.kappa_reduction_canonical m⟩

/-- Part 1, model-generic form: the size bound for `RawRate` at `c = Uc`, from the cover
certificate (`hCover`, `hLow1`, `hLow2`, `hNum`) and the proved `kappa_reduction`.
No lower witness is assumed. -/
theorem size_bound_2B_of {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    [IsProbabilityMeasure P]
    {N X1 X2 : Ω → ℝ} (hNm : Measurable N) (hX1m : Measurable X1) (hX2m : Measurable X2)
    (hN : P.map N = gaussianReal 0 1)
    (hX1 : P.map X1 = gammaMeasure (35 / 2) (1 / 2))
    (hX2 : P.map X2 = gammaMeasure (35 / 2) (1 / 2))
    (hInd : iIndepFun ![N, X1, X2] P) (m : ℝ)
    (cells : List Cell)
    (hCover : ∀ θ ∈ Box, ∃ k ∈ cells, k.mem θ)
    (hLow1 : ∀ k ∈ cells, ∀ θ ∈ Box, k.mem θ →
      (k.ρ1 * aF θ ≤ trG θ / 2 ∨ charQ (k.ρ1 * aF θ) θ ≤ 0))
    (hLow2 : ∀ k ∈ cells, ∀ θ ∈ Box, k.mem θ →
      (k.ρ2 * aF θ ≤ trG θ / 2 ∧ 0 ≤ charQ (k.ρ2 * aF θ) θ))
    (hNum : ∀ k ∈ cells, Phi P N X1 X2 Uc k.ρ1 k.ρ2 ≤ alpha) :
    ∀ θ ∈ Box, Unit2B.FormalR1.Kappa.RawRate m Uc θ ≤ alpha := by
  have hX1ae : ∀ᵐ ω ∂P, 0 ≤ X1 ω :=
    Unit2B.FormalR1.Proofs.ae_nonneg_of_gammaLaw hX1m (by norm_num) (by norm_num) hX1
  have hX2ae : ∀ᵐ ω ∂P, 0 ≤ X2 ω :=
    Unit2B.FormalR1.Proofs.ae_nonneg_of_gammaLaw hX2m (by norm_num) (by norm_num) hX2
  have hUc : (0 : ℝ) < Uc := by norm_num [Uc]
  intro θ hθ
  rw [Unit2B.FormalR1.Kappa.ProofsR10.kappa_reduction hNm hX1m hX2m hN hX1 hX2 hInd m
    Uc hUc θ hθ]
  exact Unit2B.FormalR1.Proofs.cover_sound hX1ae hX2ae cells hUc.le hCover hLow1 hLow2
    hNum θ hθ

/-- Part 1, `size_bound_2B`, on the canonical product model: the size bound for `RawRate` at
`c = Uc`, for every `θ ∈ Box` and every `m`, given the cover certificate. -/
theorem size_bound_2B (m : ℝ) (cells : List Cell)
    (hCover : ∀ θ ∈ Box, ∃ k ∈ cells, k.mem θ)
    (hLow1 : ∀ k ∈ cells, ∀ θ ∈ Box, k.mem θ →
      (k.ρ1 * aF θ ≤ trG θ / 2 ∨ charQ (k.ρ1 * aF θ) θ ≤ 0))
    (hLow2 : ∀ k ∈ cells, ∀ θ ∈ Box, k.mem θ →
      (k.ρ2 * aF θ ≤ trG θ / 2 ∧ 0 ≤ charQ (k.ρ2 * aF θ) θ))
    (hNum : ∀ k ∈ cells, Phi P0 N0 X10 X20 Uc k.ρ1 k.ρ2 ≤ alpha) :
    ∀ θ ∈ Box, Unit2B.FormalR1.Kappa.RawRate m Uc θ ≤ alpha := by
  obtain ⟨hNm, hX1m, hX2m, hN, hX1, hX2, hInd⟩ := Unit2B.FormalR1.SizeR11.canonical_model
  exact Unit2B.FormalR1.SizeR11.size_bound_2B_of hNm hX1m hX2m hN hX1 hX2 hInd m cells
    hCover hLow1 hLow2 hNum

#print axioms Unit2B.FormalR1.SizeR11.canonical_model
#print axioms Unit2B.FormalR1.SizeR11.kappa_reduction_canonical
#print axioms Unit2B.FormalR1.SizeR11.kappa_reduction_nonvacuous
#print axioms Unit2B.FormalR1.SizeR11.size_bound_2B_of
#print axioms Unit2B.FormalR1.SizeR11.size_bound_2B

end Unit2B.FormalR1.SizeR11
