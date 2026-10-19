import W4lean.Unit2BFormalR1Proofs
open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.Nonvac

open Unit2B.FormalR1

/-! ## (i) H-law witness: Ω = Fin 3 → ℝ, P0 = Measure.pi of the three marginals, coordinates. -/

def laws : Fin 3 → Measure ℝ :=
  ![gaussianReal 0 1, gammaMeasure (35 / 2) (1 / 2), gammaMeasure (35 / 2) (1 / 2)]

instance laws_isProb (i : Fin 3) : IsProbabilityMeasure (laws i) := by
  have hg : IsProbabilityMeasure (gammaMeasure (35 / 2 : ℝ) (1 / 2)) :=
    isProbabilityMeasure_gammaMeasure (by norm_num) (by norm_num)
  fin_cases i
  · show IsProbabilityMeasure (gaussianReal 0 1); infer_instance
  · exact hg
  · exact hg

def P0 : Measure (Fin 3 → ℝ) := Measure.pi laws
instance : IsProbabilityMeasure P0 := by unfold P0; infer_instance
def N0 : (Fin 3 → ℝ) → ℝ := fun ω => ω 0
def X10 : (Fin 3 → ℝ) → ℝ := fun ω => ω 1
def X20 : (Fin 3 → ℝ) → ℝ := fun ω => ω 2

lemma map_eval (i : Fin 3) : P0.map (Function.eval i) = laws i := by
  show (Measure.pi laws).map (Function.eval i) = laws i
  rw [Measure.pi_map_eval]
  simp only [measure_univ, Finset.prod_const_one, one_smul]

lemma hN0 : P0.map N0 = gaussianReal 0 1 := by
  have : N0 = Function.eval 0 := rfl
  rw [this, map_eval]
  rfl

example : Measurable N0 ∧ Measurable X10 ∧ Measurable X20 ∧
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

example : ∃ (Ω : Type) (_ : MeasurableSpace Ω) (P : Measure Ω) (_ : IsProbabilityMeasure P)
    (N X1 X2 : Ω → ℝ), Measurable N ∧ Measurable X1 ∧ Measurable X2 ∧
    P.map N = gaussianReal 0 1 ∧ P.map X1 = gammaMeasure (35 / 2) (1 / 2) ∧
    P.map X2 = gammaMeasure (35 / 2) (1 / 2) ∧ iIndepFun ![N, X1, X2] P := by
  have hmap : ∀ i : Fin 3, P0.map (Function.eval i) = laws i := fun i => by
    show (Measure.pi laws).map (Function.eval i) = laws i
    rw [Measure.pi_map_eval]
    simp only [measure_univ, Finset.prod_const_one, one_smul]
  have hind := iIndepFun_pi (Ω := fun _ => ℝ) (μ := laws) (𝓧 := fun _ => ℝ)
    (X := fun _ => id) (fun _ => aemeasurable_id)
  have hfun : (![N0, X10, X20] : Fin 3 → (Fin 3 → ℝ) → ℝ) = fun i ω => id (ω i) := by
    funext i; fin_cases i <;> rfl
  refine ⟨Fin 3 → ℝ, inferInstance, P0, inferInstance, N0, X10, X20,
    measurable_pi_apply 0, measurable_pi_apply 1, measurable_pi_apply 2,
    hmap 0, hmap 1, hmap 2, ?_⟩
  rw [hfun]; exact hind

/-! ## (ii) H-κ.  WARNING — READ BEFORE CITING.
`hKappa` with `R c θ := Phi P N X1 X2 c (r1 θ) (r2 θ)` holds by `rfl`. This is the REFLEXIVE
instance: it witnesses only that `hKappa` is consistent, NOT that any model rate equals `Phi`.
What is NOT witnessed here: instantiating `R` with the RAW model's rejection rate and proving
`hKappa` for it is exactly the S-κ reduction, which is STATEMENT-ONLY (`sorry` bodies) in
`lean/W4lean/Unit2BKappaStatementR1.lean`. Nothing in this file discharges S-κ. -/

example : ∃ R : ℝ → Param → ℝ, ∀ c, 0 < c → ∀ θ ∈ Box,
    R c θ = Phi P0 N0 X10 X20 c (r1 θ) (r2 θ) :=
  ⟨fun c θ => Phi P0 N0 X10 X20 c (r1 θ) (r2 θ), fun _ _ _ _ => rfl⟩

-- hKappa pins R on (0,∞) × Box: it is a binding, not a constraint.
example {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω) (N X1 X2 : Ω → ℝ)
    (R R' : ℝ → Param → ℝ)
    (hR : ∀ c, 0 < c → ∀ θ ∈ Box, R c θ = Phi P N X1 X2 c (r1 θ) (r2 θ))
    (hR' : ∀ c, 0 < c → ∀ θ ∈ Box, R' c θ = Phi P N X1 X2 c (r1 θ) (r2 θ)) :
    ∀ c, 0 < c → ∀ θ ∈ Box, R c θ = R' c θ :=
  fun c hc θ hθ => (hR c hc θ hθ).trans (hR' c hc θ hθ).symm

/-! ## (iii) H-cert with one cell = whole box, `ρ₁ = ρ₂ = 0`. -/

def cell0 : Cell := ⟨(1 / 25, -(3 / 2), -(3 / 2)), (25, 3 / 2, 3 / 2), 0, 0⟩

example : ∀ θ ∈ Box, ∃ k ∈ [cell0], k.mem θ :=
  fun θ hθ => ⟨cell0, List.mem_singleton_self _, hθ⟩

example : ∀ k ∈ [cell0], ∀ θ ∈ Box, k.mem θ →
    (k.ρ1 * aF θ ≤ trG θ / 2 ∨ charQ (k.ρ1 * aF θ) θ ≤ 0) := by
  intro k hk θ hθ _
  rw [List.mem_singleton] at hk
  subst hk
  simp only [Box, Set.mem_setOf_eq] at hθ
  obtain ⟨hl, -⟩ := hθ
  left
  show (0 : ℝ) * aF θ ≤ trG θ / 2
  unfold trG g00 g11
  nlinarith [sq_nonneg θ.2.1, sq_nonneg θ.2.2]

example : ∀ k ∈ [cell0], ∀ θ ∈ Box, k.mem θ →
    (k.ρ2 * aF θ ≤ trG θ / 2 ∧ 0 ≤ charQ (k.ρ2 * aF θ) θ) := by
  intro k hk θ hθ _
  rw [List.mem_singleton] at hk
  subst hk
  have hdet := (Unit2B.FormalR1.Proofs.box_pos hθ).2.1
  refine ⟨?_, ?_⟩
  · show (0 : ℝ) * aF θ ≤ trG θ / 2
    simp only [Box, Set.mem_setOf_eq] at hθ
    obtain ⟨hl, -⟩ := hθ
    unfold trG g00 g11
    nlinarith [sq_nonneg θ.2.1, sq_nonneg θ.2.2]
  · show 0 ≤ charQ ((0 : ℝ) * aF θ) θ
    have hq : charQ ((0 : ℝ) * aF θ) θ = detG θ := by unfold charQ; ring
    rw [hq]; exact hdet.le

lemma gauss_half {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {N X1 X2 : Ω → ℝ} (hNm : Measurable N) (hN : P.map N = gaussianReal 0 1) (c : ℝ) :
    Phi P N X1 X2 c 0 0 = 1 / 2 := by
  have hset : {ω | c * Real.sqrt ((0 * X1 ω + 0 * X2 ω) / nu) < N ω} = N ⁻¹' Ioi 0 := by
    ext ω; simp
  haveI : NullSingletonClass (gaussianReal 0 1) := nullSingletonClass_gaussianReal one_ne_zero
  have hP : P (N ⁻¹' Ioi 0) = gaussianReal 0 1 (Ioi 0) := by
    rw [← Measure.map_apply hNm measurableSet_Ioi, hN]
  have hsym : gaussianReal 0 1 (Iio 0) = gaussianReal 0 1 (Ioi 0) := by
    have h2 : (fun x : ℝ => -x) ⁻¹' Ioi 0 = Iio 0 := by ext x; simp
    rw [← h2, ← Measure.map_apply (f := fun x : ℝ => -x) measurable_neg measurableSet_Ioi,
      gaussianReal_map_neg, neg_zero]
  have hIic : gaussianReal 0 1 (Iic 0) = gaussianReal 0 1 (Iio 0) :=
    (measure_congr Iio_ae_eq_Iic).symm
  have hfin : gaussianReal 0 1 (Ioi 0) ≠ ⊤ := measure_ne_top _ _
  have hunion : gaussianReal 0 1 (Iic 0) + gaussianReal 0 1 (Ioi 0) = 1 := by
    rw [← measure_union (Iic_disjoint_Ioi le_rfl) measurableSet_Ioi, Iic_union_Ioi,
      measure_univ]
  rw [hIic, hsym] at hunion
  unfold Phi
  rw [hset, hP]
  have h := congrArg ENNReal.toReal hunion
  rw [ENNReal.toReal_add hfin hfin, ENNReal.toReal_one] at h
  linarith

-- hNum FAILS for cell0 under ANY P with N ~ N(0,1): Phi Uc 0 0 = 1/2 > 1/40.
lemma hhalf : Phi P0 N0 X10 X20 Uc 0 0 = 1 / 2 :=
  gauss_half (Ω := Fin 3 → ℝ) (P := P0) (N := N0) (X1 := X10) (X2 := X20)
    (measurable_pi_apply 0) hN0 Uc

example {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {N X1 X2 : Ω → ℝ}
    (hNm : Measurable N) (hN : P.map N = gaussianReal 0 1) :
    ¬ (∀ k ∈ [cell0], Phi P N X1 X2 Uc k.ρ1 k.ρ2 ≤ alpha) := by
  intro h
  have h1 := h cell0 (List.mem_singleton_self _)
  have h2 : Phi P N X1 X2 Uc cell0.ρ1 cell0.ρ2 = 1 / 2 := gauss_half hNm hN Uc
  rw [h2] at h1
  norm_num [alpha] at h1

example : Phi P0 N0 X10 X20 Uc 0 0 = 1 / 2 := hhalf

example : alpha < Phi P0 N0 X10 X20 Uc 0 0 := by
  rw [hhalf]
  norm_num [alpha]

-- Joint: on (P0, N0, X10, X20, reflexive R, [cell0]) hKappa AND hCover AND hLow1 AND hLow2
-- hold, while hNum is FALSE. So H-cert's geometric half is satisfiable and its numeric
-- half is genuinely load-bearing.
example :
    let R : ℝ → Param → ℝ := fun c θ => Phi P0 N0 X10 X20 c (r1 θ) (r2 θ)
    (∀ c, 0 < c → ∀ θ ∈ Box, R c θ = Phi P0 N0 X10 X20 c (r1 θ) (r2 θ)) ∧
    (∀ θ ∈ Box, ∃ k ∈ [cell0], k.mem θ) ∧
    (∀ k ∈ [cell0], ∀ θ ∈ Box, k.mem θ →
      (k.ρ1 * aF θ ≤ trG θ / 2 ∨ charQ (k.ρ1 * aF θ) θ ≤ 0)) ∧
    (∀ k ∈ [cell0], ∀ θ ∈ Box, k.mem θ →
      (k.ρ2 * aF θ ≤ trG θ / 2 ∧ 0 ≤ charQ (k.ρ2 * aF θ) θ)) ∧
    ¬ (∀ k ∈ [cell0], Phi P0 N0 X10 X20 Uc k.ρ1 k.ρ2 ≤ alpha) := by
  intro R
  have hhalf' : Phi P0 N0 X10 X20 Uc 0 0 = 1 / 2 := hhalf
  have htr : ∀ θ ∈ Box, (0 : ℝ) * aF θ ≤ trG θ / 2 := by
    intro θ hθ
    simp only [Box, Set.mem_setOf_eq] at hθ
    obtain ⟨hl, -⟩ := hθ
    unfold trG g00 g11
    nlinarith [sq_nonneg θ.2.1, sq_nonneg θ.2.2]
  refine ⟨fun _ _ _ _ => rfl, fun θ hθ => ⟨cell0, List.mem_singleton_self _, hθ⟩, ?_, ?_, ?_⟩
  · intro k hk θ hθ _
    rw [List.mem_singleton] at hk
    subst hk
    exact Or.inl (htr θ hθ)
  · intro k hk θ hθ _
    rw [List.mem_singleton] at hk
    subst hk
    refine ⟨htr θ hθ, ?_⟩
    show 0 ≤ charQ ((0 : ℝ) * aF θ) θ
    have hq : charQ ((0 : ℝ) * aF θ) θ = detG θ := by unfold charQ; ring
    rw [hq]; exact (Unit2B.FormalR1.Proofs.box_pos hθ).2.1.le
  · intro h
    have h1 := h cell0 (List.mem_singleton_self _)
    have h2 : Phi P0 N0 X10 X20 Uc cell0.ρ1 cell0.ρ2 = 1 / 2 := hhalf'
    rw [h2] at h1
    norm_num [alpha] at h1

/-! ## (iv) -/
example : corner ∈ Box := by norm_num [Box, corner]

end Unit2B.FormalR1.Nonvac
