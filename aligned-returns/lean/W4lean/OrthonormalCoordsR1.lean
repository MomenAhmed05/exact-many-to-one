import W4lean.ProductMarginalR1
import Mathlib.Probability.Distributions.Gaussian.Multivariate
import Mathlib.Tactic.FunProp

/-!
Stage D′, sub-lemma **D-7c**: the i.i.d. conclusion.

The marginal it needs is already in `ProductMarginalR1.lean`, so this file holds only
the two statements that compose it with the standard Gaussian on the carrier.  Both are for an
ARBITRARY orthonormal basis, so neither waits on D-6 or D-7b.

There is **no `Orthonormal -> iIndepFun` theorem in Mathlib** — established by search and
recorded.  So the i.i.d. conclusion has to come from the joint law instead, and the content of
D-7c is that the joint law is available: an orthonormal change of coordinates carries the
standard Gaussian to the product of standard Gaussians.
-/

open MeasureTheory ProbabilityTheory Real Set
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR4

variable {ι : Type*} [Fintype ι] [DecidableEq ι] [Zero ι] [AddCommMonoid ι]

/-- The marginal, on the carrier `b.repr` lands in. -/
theorem product_marginal (i : ι) :
    ((Measure.pi (fun _ : ι => gaussianReal 0 1)).map (WithLp.toLp 2)).map
      (fun p : EuclideanSpace ℝ ι => p i) = gaussianReal 0 1 := by
  have h := gaussianReal_pi_marginal' i
  have e : (fun p : EuclideanSpace ℝ ι => p i) ∘ (WithLp.toLp 2)
      = (fun x : ι → ℝ => x i) := by
    funext x
    rfl
  have hm := Measure.map_map (μ := (Measure.pi (fun _ : ι => gaussianReal 0 1) : Measure (ι → ℝ)))
    (f := WithLp.toLp 2) (g := fun p : EuclideanSpace ℝ ι => p i)
    (by fun_prop) (by fun_prop)
  rw [hm, e]
  exact h

/-- **D-7c, the heart.**  Pushing the standard Gaussian on the carrier through the
representation map of an ORTHONORMAL basis returns the product of standard Gaussians. -/
theorem stdGaussian_repr (b : OrthonormalBasis ι ℝ (EuclideanSpace ℝ ι)) :
    Measure.map b.repr (stdGaussian (EuclideanSpace ℝ ι))
      = (Measure.pi (fun _ : ι => gaussianReal 0 1)).map (WithLp.toLp 2) := by
  exact stdGaussian_map b.repr |>.trans map_pi_eq_stdGaussian.symm

/-- **D-7c.**  Each basis coordinate of a standard Gaussian is standard normal. -/
theorem basis_coords_law (b : OrthonormalBasis ι ℝ (EuclideanSpace ℝ ι)) (i : ι) :
    Measure.map (fun p : EuclideanSpace ℝ ι => (b.repr p) i)
      (stdGaussian (EuclideanSpace ℝ ι)) = gaussianReal 0 1 := by
  have h1 := stdGaussian_repr b
  have e : (fun p : EuclideanSpace ℝ ι => p i) ∘ b.repr
      = (fun p : EuclideanSpace ℝ ι => (b.repr p) i) := by
    funext p
    rfl
  have hm := Measure.map_map (μ := (stdGaussian (EuclideanSpace ℝ ι) : Measure (EuclideanSpace ℝ ι)))
    (f := b.repr) (g := fun p : EuclideanSpace ℝ ι => p i) (by fun_prop) (by fun_prop)
  calc Measure.map (fun p : EuclideanSpace ℝ ι => (b.repr p) i) (stdGaussian (EuclideanSpace ℝ ι))
      = Measure.map ((fun p : EuclideanSpace ℝ ι => p i) ∘ b.repr) (stdGaussian (EuclideanSpace ℝ ι)) := by
        rw [e]
      _ = Measure.map (fun p : EuclideanSpace ℝ ι => p i)
          (Measure.map b.repr (stdGaussian (EuclideanSpace ℝ ι))) := hm.symm
      _ = Measure.map (fun p : EuclideanSpace ℝ ι => p i)
          ((Measure.pi (fun _ : ι => gaussianReal 0 1)).map (WithLp.toLp 2)) := congrArg _ h1
      _ = gaussianReal 0 1 := product_marginal i

end Unit2B.FormalR1.Kappa.ProofsR4
