import W4lean.Unit2BKappaProofsR3
-- Stage D needs the multivariate Gaussian API, which is stated on `EuclideanSpace`.
-- Located by source inspection at the pinned revision: `map_pi_eq_stdGaussian`
-- (Multivariate.lean:135), `stdGaussian_map` (:126), `stdGaussian_eq_map_pi_orthonormalBasis`
-- (:144), and `iIndepFun_iff_map_fun_eq_pi_map` (Independence/Basic.lean:862).
import Mathlib.Probability.Distributions.Gaussian.Multivariate

open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR4

open Unit2B.FormalR1.Kappa
open Unit2B.FormalR1.Kappa.ProofsR3

/-!
Stage D, under this method:

  * do NOT build 71 explicit vectors and do NOT flatten to R^108;
  * prove ONE general tensor-product inner-product lemma, indexed by `Fin 3 x Fin n`;
  * use ONE carrier, `EuclideanSpace R PanelIdx`, with a single measurable equivalence
    to `Panel`, because the `stdGaussian` API lives there.

Sub-lemmas were declared BEFORE the first build.

**D-5a (the orthonormal extension of the mean direction) is NOT COMPLETE here and is not
claimed in this file.**  The Mathlib `Orthonormal` API needed to discharge it was not
located precisely enough within the attempts allowed, so the file contains only the
sub-lemmas that closed cleanly: D-0a (the carrier) and D-7a (the tensor inner-product
identity, plus its orthogonal corollary).
-/

/-- **D-0a.**  The panel law, pushed along the Euclidean identification, *is* the
standard Gaussian on `EuclideanSpace R PanelIdx`.  This is the carrier step the whole
stage rests on, and it is exactly `map_pi_eq_stdGaussian` with `PanelIdx` as the index. -/
theorem panelLaw_map_toLp :
    panelLaw.map (WithLp.toLp 2) = stdGaussian (EuclideanSpace ℝ PanelIdx) := by
  exact map_pi_eq_stdGaussian

/-- The tensor used throughout Stage D.  `tens a x` is the product vector
`a t * x i` at the index `(i, t)`, i.e. the elementary tensor of a vector in the time
space with a vector in the three-dimensional loading space. -/
def tens {ι κ : Type*} (a : ι → ℝ) (x : κ → ℝ) : κ × ι → ℝ := fun p => x p.1 * a p.2

/-- **D-7a.**  The single lemma that replaces building 71 explicit vectors: the inner
product of two elementary tensors is the product of the inner products.

Stated as a finite-sum identity rather than via the `⟪,⟫` notation, so that it is a
direct computation depending on no inner-product API.  The inner-product form
`<a (x) x, b (x) y> = <a,b> * <x,y>` is this equation read through
`EuclideanSpace.inner_eq_sum_mul`, which is what the rest of Stage D will use.

`Finset.mul_sum` puts the scalar on the LEFT, so every factoring step below is written
`sum_congr` (to commute) followed by `Finset.mul_sum ... |>.symm`.  Getting that
direction wrong is the single most common error in this file. -/
theorem tens_inner {ι κ : Type*} [Fintype ι] [Fintype κ] (a b : ι → ℝ) (x y : κ → ℝ) :
    (∑ p : κ × ι, (tens a x) p * (tens b y) p)
      = (∑ t : ι, a t * b t) * ∑ i : κ, x i * y i := by
  simp only [tens, Fintype.sum_prod_type]
  calc ∑ i : κ, ∑ t : ι, (x i * a t) * (y i * b t)
      = ∑ i : κ, ∑ t : ι, (a t * b t) * (x i * y i) := by
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro t _
          ring
    _ = ∑ i : κ, (x i * y i) * ∑ t : ι, (a t * b t) := by
          apply Finset.sum_congr rfl
          intro i _
          calc ∑ t : ι, (a t * b t) * (x i * y i)
              = ∑ t : ι, (x i * y i) * (a t * b t) := by
                  apply Finset.sum_congr rfl
                  intro t _
                  rw [mul_comm]
            _ = (x i * y i) * ∑ t : ι, (a t * b t) := (Finset.mul_sum _ _ _).symm
    _ = (∑ t : ι, (a t * b t)) * ∑ i : κ, (x i * y i) := by
          calc ∑ i : κ, (x i * y i) * ∑ t : ι, (a t * b t)
              = ∑ i : κ, (∑ t : ι, (a t * b t)) * (x i * y i) := by
                  apply Finset.sum_congr rfl
                  intro i _
                  rw [mul_comm]
            _ = (∑ t : ι, (a t * b t)) * ∑ i : κ, (x i * y i) := (Finset.mul_sum _ _ _).symm

/-- **D-7a, immediate consequence.**  The two elementary tensors of two orthogonal pairs
are orthogonal.  This is the whole reason 71 explicit vectors are unnecessary. -/
theorem tens_inner_eq_zero {ι κ : Type*} [Fintype ι] [Fintype κ] {a b : ι → ℝ} {x y : κ → ℝ}
    (ha : ∑ t, a t * b t = 0) (hx : ∑ i, x i * y i = 0) :
    (∑ p : κ × ι, (tens a x) p * (tens b y) p) = 0 := by
  rw [tens_inner, ha, hx]
  ring

end Unit2B.FormalR1.Kappa.ProofsR4
