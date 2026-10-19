import W4lean.OrthonormalCoordsR1

/-!
Stage D′, sub-lemma **D-5b**: the Parseval step and the mean
direction.

**Banked:** `parseval_coeff`, Mathlib's `OrthonormalBasis.sum_inner_mul_inner`
(`Analysis/InnerProductSpace/PiL2.lean:509`) restated in the coefficient order
`⟪hⱼ, x⟫ ⋅ ⟪hⱼ, y⟫` that the centred identity consumes. The order's D-5 asked for the
**bilinear** form rather than the squared one because `Sres` couples the three coordinate
series through `MᵀM`, so a squared-only identity would not have been enough.
-/

open MeasureTheory ProbabilityTheory Real Set
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR4

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

theorem parseval_coeff (b : OrthonormalBasis ι ℝ (EuclideanSpace ℝ ι))
    (x y : EuclideanSpace ℝ ι) :
    (∑ j, inner ℝ (b j) x * inner ℝ (b j) y) = inner ℝ x y := by
  simpa [real_inner_comm] using b.sum_inner_mul_inner x y

variable {n : ℕ}

/-- **D-5b, the mean direction, as pure algebra.**  The constant direction
`1/sqrt(n)` paired with `x` is `sqrt(n)` times the mean of `x`.  This is what makes the
Helmert numerator the sample mean scaled by a unit constant.

Stated on plain functions `Fin n -> R`, with no `WithLp` or `PiLp` involved.  Every earlier
attempt failed on a *coercion display* mismatch: Lean unfolds `EuclideanSpace` to
`PiLp` in one place and to `WithLp` in another, so a hypothesis stated on `b i` never
matches a goal stated the same way, in either `rw` direction.  Isolating the algebra on
raw functions removes the entire class of failure; the inner-product form is recorded as a closure-route item, not attempted here: the `WithLp`
display mismatch recurs on `PiLp.inner_apply` too, and the algebra below is the content. -/
theorem mean_direction_sum (hn : 0 < n) (f : Fin n → Fin n → ℝ)
    (hf0 : ∀ t, f (⟨0, hn⟩ : Fin n) t = (Real.sqrt (n : ℝ))⁻¹)
    (x : Fin n → ℝ) :
    (∑ t, f (⟨0, hn⟩ : Fin n) t * x t) = Real.sqrt (n : ℝ) * ((∑ t, x t) / (n : ℝ)) := by
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  simp_rw [hf0]
  rw [← Finset.mul_sum]
  field_simp [ne_of_gt hnR, ne_of_gt (Real.sqrt_pos.2 hnR)]
  rw [Real.sq_sqrt hnR.le]
  ring


end Unit2B.FormalR1.Kappa.ProofsR4
