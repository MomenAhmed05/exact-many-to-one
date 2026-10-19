import Mathlib
-- Stage P′, sub-lemma P-1 `pi_coord_indep`.
-- Stated for an arbitrary `Fintype ι`, so it does not depend on D-6.
-- Signatures confirmed by source inspection at Mathlib 5ed2965256430c3649e86755f9576b54eca72435:
--   `iIndepFun_pi` (Probability/Independence/Basic.lean:888),
--   `iIndepFun.indepFun_finset` (:796), `IndepFun.comp` (:756),
--   `Finset.sum_coe_sort` (to_additive of Algebra/BigOperators/Group/Finset/Basic.lean:461),
--   `Finset.measurable_sum` (to_additive of MeasureTheory/Group/Arithmetic.lean:831).

open MeasureTheory ProbabilityTheory
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR5

/-- The coordinates of a product of standard Gaussians are mutually independent. -/
theorem pi_coord_iIndepFun (ι : Type*) [Fintype ι] :
    iIndepFun (fun i (x : ι → ℝ) => x i) (Measure.pi fun _ : ι => gaussianReal 0 1) :=
  iIndepFun_pi (μ := fun _ : ι => gaussianReal 0 1) (X := fun _ => (id : ℝ → ℝ))
    (fun _ => measurable_id.aemeasurable)

/-- **P-1 (a).**  The coordinate at `a` is independent of any measurable function of the
coordinates in a set `T` with `a ∉ T`. -/
theorem pi_coord_indep_single {ι : Type*} [Fintype ι] {γ : Type*} [MeasurableSpace γ]
    (a : ι) (T : Finset ι) (haT : a ∉ T) {F : (T → ℝ) → γ} (hF : Measurable F) :
    IndepFun (fun x : ι → ℝ => x a) (fun x : ι → ℝ => F (fun i : T => x i))
      (Measure.pi fun _ : ι => gaussianReal 0 1) := by
  have hsep : Disjoint ({a} : Finset ι) T := Finset.disjoint_singleton_left.mpr haT
  have ha : a ∈ ({a} : Finset ι) := Finset.mem_singleton_self a
  have h := (pi_coord_iIndepFun ι).indepFun_finset {a} T hsep (fun i => measurable_pi_apply i)
  -- typed first, so that `IndepFun.comp`'s implicit `_mγ` is fixed to the instance on `ℝ`;
  -- inline, it was left as an mvar that `?X ?i` could not be unified against (an earlier build)
  have hφ : Measurable (fun v : ({a} : Finset ι) → ℝ => v ⟨a, ha⟩) := measurable_pi_apply _
  have h2 := h.comp hφ hF
  exact h2

/-- **P-1 (b).**  Sums of squares over disjoint blocks of coordinates are independent. -/
theorem pi_coord_indep_blocks {ι : Type*} [Fintype ι] (B₁ B₂ : Finset ι)
    (hB : Disjoint B₁ B₂) :
    IndepFun (fun x : ι → ℝ => ∑ i ∈ B₁, x i ^ 2) (fun x : ι → ℝ => ∑ i ∈ B₂, x i ^ 2)
      (Measure.pi fun _ : ι => gaussianReal 0 1) := by
  have h := (pi_coord_iIndepFun ι).indepFun_finset B₁ B₂ hB (fun i => measurable_pi_apply i)
  have h2 := h.comp (φ := fun v : B₁ → ℝ => ∑ i, v i ^ 2) (ψ := fun v : B₂ → ℝ => ∑ i, v i ^ 2)
    (Finset.measurable_sum _ (fun i _ => (measurable_pi_apply i).pow_const 2))
    (Finset.measurable_sum _ (fun i _ => (measurable_pi_apply i).pow_const 2))
  have e1 : (fun v : B₁ → ℝ => ∑ i, v i ^ 2) ∘ (fun (x : ι → ℝ) (i : B₁) => x i)
      = (fun x : ι → ℝ => ∑ i ∈ B₁, x i ^ 2) := by
    funext x
    exact Finset.sum_coe_sort (s := B₁) (f := fun i => x i ^ 2)
  have e2 : (fun v : B₂ → ℝ => ∑ i, v i ^ 2) ∘ (fun (x : ι → ℝ) (i : B₂) => x i)
      = (fun x : ι → ℝ => ∑ i ∈ B₂, x i ^ 2) := by
    funext x
    exact Finset.sum_coe_sort (s := B₂) (f := fun i => x i ^ 2)
  -- transport along the two function equalities without `rw`, by `subst` on local variables
  have key : ∀ {f f' g g' : (ι → ℝ) → ℝ}, f = f' → g = g' →
      IndepFun f g (Measure.pi fun _ : ι => gaussianReal 0 1) →
      IndepFun f' g' (Measure.pi fun _ : ι => gaussianReal 0 1) := by
    intro f f' g g' hf hg hfg
    subst hf
    subst hg
    exact hfg
  exact key e1 e2 h2

/-- **P-1 `pi_coord_indep`.**  The two pairwise facts, (a) and (b), for a product of i.i.d.
standard normal coordinates over an arbitrary finite index type `ι`. -/
theorem pi_coord_indep {ι : Type*} [Fintype ι] {γ : Type*} [MeasurableSpace γ] :
    (∀ (a : ι) (T : Finset ι), a ∉ T → ∀ F : (T → ℝ) → γ, Measurable F →
      IndepFun (fun x : ι → ℝ => x a) (fun x : ι → ℝ => F (fun i : T => x i))
        (Measure.pi fun _ : ι => gaussianReal 0 1)) ∧
    (∀ B₁ B₂ : Finset ι, Disjoint B₁ B₂ →
      IndepFun (fun x : ι → ℝ => ∑ i ∈ B₁, x i ^ 2) (fun x : ι → ℝ => ∑ i ∈ B₂, x i ^ 2)
        (Measure.pi fun _ : ι => gaussianReal 0 1)) :=
  ⟨fun a T haT _ hF => pi_coord_indep_single a T haT hF,
    fun B₁ B₂ hB => pi_coord_indep_blocks B₁ B₂ hB⟩

end Unit2B.FormalR1.Kappa.ProofsR5

#print axioms Unit2B.FormalR1.Kappa.ProofsR5.pi_coord_indep
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.pi_coord_indep_single
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.pi_coord_indep_blocks
