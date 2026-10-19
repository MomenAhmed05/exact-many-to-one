import Mathlib
-- Stage P′, sub-lemma P-4 `hInd_split` + `pair_transfer`.
-- Abstract measure theory: arbitrary measures and measurable functions.  No dependence on
-- D-6, on P-1's block sums, or on any panel geometry.
-- Signatures confirmed by source inspection at Mathlib 5ed2965256430c3649e86755f9576b54eca72435:
--   `iIndepFun` (Probability/Independence/Basic.lean:136, `[m : ∀ x, MeasurableSpace (β x)]`),
--   `iIndepFun.indepFun` (:449), `IndepFun.symm` (:710), `IndepFun.comp` (:756),
--   `iIndepFun.indepFun_prodMk` (:809; kernel form Kernel/IndepFun.lean:439, no instance
--   hypotheses), `Matrix.cons_val_zero` (Data/Fin/VecNotation.lean:128, `rfl`),
--   `Matrix.cons_val_one` (:271, `rfl`), `Matrix.cons_val_two` (:274, `rfl`),
--   `indepFun_iff_map_prod_eq_prod_map_map` (Independence/Basic.lean:703, `[IsFiniteMeasure μ]`),
--   `Measure.map_apply` (MeasureTheory/Measure/Map.lean:170),
--   `Measurable.prodMk` (MeasureTheory/MeasurableSpace/Constructions.lean:409).
-- H-12.  The split is stated once for an explicit triple `F : Fin 3 → Ω → ℝ` and bridged to
-- `![N, X1, X2]` once, in `hInd_split`, by defeq only: the three `cons_val_*` above are `rfl`,
-- so no `rw` or `simp` with them is used anywhere.

open MeasureTheory ProbabilityTheory
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR5

/-- The two pairwise consequences of three-way independence, for an explicit triple `F`:
`F 1 ⟂ F 2`, and the pair `(F 1, F 2)` against `F 0`. -/
theorem hInd_split_fun {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {F : Fin 3 → Ω → ℝ} (hInd : iIndepFun F P) (hFm : ∀ i, Measurable (F i)) :
    IndepFun (F 1) (F 2) P ∧ IndepFun (fun ω => (F 1 ω, F 2 ω)) (F 0) P :=
  ⟨hInd.indepFun (i := 1) (j := 2) (by decide),
    hInd.indepFun_prodMk hFm 1 2 0 (by decide) (by decide)⟩

/-- **P-4 `hInd_split`.**  From `hInd : iIndepFun ![N, X1, X2] P`, the pairwise facts that
`kappa_reduction` and K-7 consume: `X1 ⟂ X2`, and `N` against the pair `(X1, X2)` in both
orders. -/
theorem hInd_split {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {N X1 X2 : Ω → ℝ}
    (hNm : Measurable N) (hX1m : Measurable X1) (hX2m : Measurable X2)
    (hInd : iIndepFun ![N, X1, X2] P) :
    IndepFun X1 X2 P ∧ IndepFun (fun ω => (X1 ω, X2 ω)) N P ∧
      IndepFun N (fun ω => (X1 ω, X2 ω)) P := by
  have hFm : ∀ i, Measurable (![N, X1, X2] i) := by
    intro i
    fin_cases i
    · exact hNm
    · exact hX1m
    · exact hX2m
  have h := hInd_split_fun hInd hFm
  exact ⟨h.1, h.2, h.2.symm⟩

/-- **P-4 `hInd_split`, composed.**  `N` is independent of `μ₁ X1 + μ₂ X2`, in both orders. -/
theorem hInd_split_lin {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} {N X1 X2 : Ω → ℝ}
    (hNm : Measurable N) (hX1m : Measurable X1) (hX2m : Measurable X2)
    (hInd : iIndepFun ![N, X1, X2] P) (μ₁ μ₂ : ℝ) :
    IndepFun (fun ω => μ₁ * X1 ω + μ₂ * X2 ω) N P ∧
      IndepFun N (fun ω => μ₁ * X1 ω + μ₂ * X2 ω) P := by
  have hpair := (hInd_split hNm hX1m hX2m hInd).2.1
  -- typed first, so that `IndepFun.comp`'s implicit `_mγ`, `_mγ'` are fixed (as in P-1)
  have hφ : Measurable (fun p : ℝ × ℝ => μ₁ * p.1 + μ₂ * p.2) := by fun_prop
  have hψ : Measurable (id : ℝ → ℝ) := measurable_id
  have h := hpair.comp hφ hψ
  exact ⟨h, h.symm⟩

/-- **P-4 `pair_transfer`.**  Two independent pairs, on possibly different spaces, whose
marginals agree have the same joint law: both joint laws are the product of the marginals. -/
theorem pair_transfer {Ω Ω' β γ : Type*} [MeasurableSpace Ω] [MeasurableSpace Ω']
    [MeasurableSpace β] [MeasurableSpace γ] {μ : Measure Ω} {ν : Measure Ω'}
    [IsFiniteMeasure μ] [IsFiniteMeasure ν] {f : Ω → β} {g : Ω → γ} {f' : Ω' → β}
    {g' : Ω' → γ} (hf : AEMeasurable f μ) (hg : AEMeasurable g μ) (hf' : AEMeasurable f' ν)
    (hg' : AEMeasurable g' ν) (hfg : IndepFun f g μ) (hfg' : IndepFun f' g' ν)
    (h1 : μ.map f = ν.map f') (h2 : μ.map g = ν.map g') :
    μ.map (fun ω => (f ω, g ω)) = ν.map (fun ω => (f' ω, g' ω)) := by
  have e := (indepFun_iff_map_prod_eq_prod_map_map hf hg).1 hfg
  have e' := (indepFun_iff_map_prod_eq_prod_map_map hf' hg').1 hfg'
  rw [e, e', h1, h2]

/-- **P-4 `pair_transfer`, on events.**  The form `kappa_reduction` consumes through
`map_apply`: the two pairs put the same mass on every measurable set. -/
theorem pair_transfer_set {Ω Ω' β γ : Type*} [MeasurableSpace Ω] [MeasurableSpace Ω']
    [MeasurableSpace β] [MeasurableSpace γ] {μ : Measure Ω} {ν : Measure Ω'}
    [IsFiniteMeasure μ] [IsFiniteMeasure ν] {f : Ω → β} {g : Ω → γ} {f' : Ω' → β}
    {g' : Ω' → γ} (hf : Measurable f) (hg : Measurable g) (hf' : Measurable f')
    (hg' : Measurable g') (hfg : IndepFun f g μ) (hfg' : IndepFun f' g' ν)
    (h1 : μ.map f = ν.map f') (h2 : μ.map g = ν.map g') {A : Set (β × γ)}
    (hA : MeasurableSet A) :
    μ ((fun ω => (f ω, g ω)) ⁻¹' A) = ν ((fun ω => (f' ω, g' ω)) ⁻¹' A) := by
  have e := pair_transfer hf.aemeasurable hg.aemeasurable hf'.aemeasurable hg'.aemeasurable
    hfg hfg' h1 h2
  calc μ ((fun ω => (f ω, g ω)) ⁻¹' A)
      = μ.map (fun ω => (f ω, g ω)) A := Eq.symm (Measure.map_apply (hf.prodMk hg) hA)
    _ = ν.map (fun ω => (f' ω, g' ω)) A := by rw [e]
    _ = ν ((fun ω => (f' ω, g' ω)) ⁻¹' A) := Measure.map_apply (hf'.prodMk hg') hA

/-- **P-4, the P side.**  From `hInd`, the joint law of `N` and `μ₁ X1 + μ₂ X2` is
`N(0,1) ⊗ law(μ₁ X1 + μ₂ X2)`. -/
theorem pair_law_lin {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsFiniteMeasure P]
    {N X1 X2 : Ω → ℝ} (hNm : Measurable N) (hX1m : Measurable X1) (hX2m : Measurable X2)
    (hN : P.map N = gaussianReal 0 1) (hInd : iIndepFun ![N, X1, X2] P) (μ₁ μ₂ : ℝ) :
    P.map (fun ω => (N ω, μ₁ * X1 ω + μ₂ * X2 ω))
      = (gaussianReal 0 1).prod (P.map (fun ω => μ₁ * X1 ω + μ₂ * X2 ω)) := by
  have hL : Measurable (fun ω => μ₁ * X1 ω + μ₂ * X2 ω) := by fun_prop
  have hNL := (hInd_split_lin hNm hX1m hX2m hInd μ₁ μ₂).2
  have e := (indepFun_iff_map_prod_eq_prod_map_map hNm.aemeasurable hL.aemeasurable).1 hNL
  rw [hN] at e
  exact e

end Unit2B.FormalR1.Kappa.ProofsR5

#print axioms Unit2B.FormalR1.Kappa.ProofsR5.hInd_split_fun
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.hInd_split
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.hInd_split_lin
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.pair_transfer
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.pair_transfer_set
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.pair_law_lin
