import Mathlib
import W4lean.ChiSquaredSpikeR1
import W4lean.ProductMarginalR1
import W4lean.Unit2BKappaProofsR5
-- Stage P′, sub-lemma P-2 `block_chiSq`.
-- Stated for an arbitrary `Fintype ι`, so it does not depend on D-6.
-- Primary route of R2 §2: `iIndepFun.precomp`, then `iIndepFun_iff_map_fun_eq_pi_map` on
-- `Fin 35`, then the banked marginal, then the banked `chiSq_eq_gamma 35`, then the H-13 cast.
-- Signatures confirmed by source inspection at Mathlib 5ed2965256430c3649e86755f9576b54eca72435:
--   `iIndepFun.precomp` (Probability/Independence/Basic.lean:324),
--   `iIndepFun_iff_map_fun_eq_pi_map` (:859), `Measure.map_map` (MeasureTheory/Measure/Map.lean:208),
--   `Measurable.of_eval` (MeasureTheory/MeasurableSpace/Constructions.lean:596),
--   `measurable_pi_apply` (:584), `Finset.measurable_sum` (to_additive of
--   MeasureTheory/Group/Arithmetic.lean:831), `Measurable.pow_const` (:192).
-- Banked in this project: `pi_coord_iIndepFun` (Unit2BKappaProofsR5.lean:15),
--   `gaussianReal_pi_marginal'` (ProductMarginalR1.lean:31, namespace `ProofsR4`),
--   `chiSq_eq_gamma` (ChiSquaredSpikeR1.lean:464).

open MeasureTheory ProbabilityTheory
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR5

/-- **H-13.**  The one cast between the `ℕ`-indexed degrees of freedom of `chiSq_eq_gamma`
and the real literal of the statement.  Closed once, here, and reused. -/
theorem chiSq35_cast : ((35 : ℕ) : ℝ) / 2 = (35 : ℝ) / 2 := by norm_num

/-- The coordinates of a product of standard Gaussians along an injective
`g : Fin 35 → ι` have law the product of 35 standard Gaussians. -/
theorem block_coords_law {ι : Type*} [Fintype ι] {g : Fin 35 → ι}
    (hg : Function.Injective g) :
    (Measure.pi fun _ : ι => gaussianReal 0 1).map (fun (x : ι → ℝ) (i : Fin 35) => x (g i))
      = Measure.pi fun _ : Fin 35 => gaussianReal 0 1 := by
  classical
  have hind : iIndepFun (fun (i : Fin 35) (x : ι → ℝ) => x (g i))
      (Measure.pi fun _ : ι => gaussianReal 0 1) :=
    (pi_coord_iIndepFun ι).precomp hg
  have hf : ∀ i : Fin 35, AEMeasurable (fun x : ι → ℝ => x (g i))
      (Measure.pi fun _ : ι => gaussianReal 0 1) :=
    fun i => (measurable_pi_apply (g i)).aemeasurable
  have h : (Measure.pi fun _ : ι => gaussianReal 0 1).map
        (fun (x : ι → ℝ) (i : Fin 35) => x (g i))
      = Measure.pi (fun i : Fin 35 =>
          (Measure.pi fun _ : ι => gaussianReal 0 1).map (fun x : ι → ℝ => x (g i))) :=
    (iIndepFun_iff_map_fun_eq_pi_map (f := fun (i : Fin 35) (x : ι → ℝ) => x (g i)) hf).1 hind
  have hm : Measure.pi (fun i : Fin 35 =>
        (Measure.pi fun _ : ι => gaussianReal 0 1).map (fun x : ι → ℝ => x (g i)))
      = Measure.pi fun _ : Fin 35 => gaussianReal 0 1 := by
    congr 1
    funext i
    exact Unit2B.FormalR1.Kappa.ProofsR4.gaussianReal_pi_marginal' (g i)
  exact h.trans hm

/-- **P-2 `block_chiSq`.**  For an injective `g : Fin 35 → ι`, the block sum of squares of a
product of i.i.d. standard normal coordinates has law `χ²₃₅`, written in Mathlib's
shape-and-rate parameterisation as `gammaMeasure (35 / 2) (1 / 2)`. -/
theorem block_chiSq {ι : Type*} [Fintype ι] {g : Fin 35 → ι} (hg : Function.Injective g) :
    (Measure.pi fun _ : ι => gaussianReal 0 1).map (fun x : ι → ℝ => ∑ i, x (g i) ^ 2)
      = gammaMeasure (35 / 2) (1 / 2) := by
  have hR : Measurable (fun (x : ι → ℝ) (i : Fin 35) => x (g i)) :=
    Measurable.of_eval (fun i => measurable_pi_apply (g i))
  have hS : Measurable (fun y : Fin 35 → ℝ => ∑ i, y i ^ 2) :=
    Finset.measurable_sum _ (fun i _ => (measurable_pi_apply i).pow_const 2)
  show (Measure.pi fun _ : ι => gaussianReal 0 1).map
      ((fun y : Fin 35 → ℝ => ∑ i, y i ^ 2) ∘ (fun (x : ι → ℝ) (i : Fin 35) => x (g i))) = _
  rw [← Measure.map_map hS hR, block_coords_law hg, chiSq_eq_gamma 35 (by norm_num),
    chiSq35_cast]

end Unit2B.FormalR1.Kappa.ProofsR5

#print axioms Unit2B.FormalR1.Kappa.ProofsR5.chiSq35_cast
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.block_coords_law
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.block_chiSq
