import W4lean.ChiSquaredSpikeR1

/-!
Kept from the attempt at **D-7c** (`OrthonormalCoordsR1.lean`), which succeeded.

`gaussianReal_pi_marginal'` compiled in every build in which it was reached, and is retained
here so the work is not discarded with the sub-lemma that did not close.  It is stated for an
arbitrary finite index type, which is the generality D-7c needs; the `Fin k` instance used by
the D-3 statement is the special case.

**What did not close, recorded precisely.** In `OrthonormalCoordsR1.lean` the two statements
that compose this marginal with the standard Gaussian on the carrier — `stdGaussian_repr`
and `basis_coords_law` — did not close in the attempts allowed for them.  The
obstruction is a `WithLp` coercion mismatch, not a mathematical one: `Measure.map_map`
infers its second function as `fun p ↦ p.ofLp i`, while the statement is written `p i`, and
the two do not unify syntactically at the `Measure.map` level.  Both are definitionally
equal, and Mathlib's own `map_pi_eq_stdGaussian` is stated over the same `toLp` layer, so
this is API friction of the same class met in A′.  It is recorded as a
closure route item, not as a mathematical obstruction.
-/

open MeasureTheory ProbabilityTheory Real Set
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR4

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The marginal of a product of i.i.d. standard Gaussians is again a standard Gaussian.
Proved from `Measure.pi_pi` on a cylinder set, so no characteristic functions are involved —
they were needed for the `Fin k` case and the same route generalises. -/
theorem gaussianReal_pi_marginal' (i : ι) :
    (Measure.pi fun _ : ι => gaussianReal 0 1).map (fun x => x i) = gaussianReal 0 1 := by
  have hmeas : Measurable (fun x : ι → ℝ => x i) := measurable_pi_apply i
  refine Measure.ext fun A hA => ?_
  have hpull : (fun x : ι → ℝ => x i) ⁻¹' A
      = Set.pi Set.univ (fun j => if j = i then A else Set.univ) := by
    ext x
    simp only [Set.mem_preimage, Set.mem_pi]
    constructor
    · intro hx j hj
      by_cases hji : j = i
      · simpa [hji] using hx
      · simp [hji]
    · intro hx
      have := hx i (Set.mem_univ i)
      simpa using this
  rw [Measure.map_apply hmeas hA, hpull, Measure.pi_pi]
  have hval : (∏ j, (gaussianReal 0 1) (if j = i then A else Set.univ)) = (gaussianReal 0 1) A := by
    have h1 : (gaussianReal 0 1) (Set.univ : Set ℝ) = 1 :=
      (inferInstance : IsProbabilityMeasure (gaussianReal 0 1)).measure_univ
    calc (∏ j, (gaussianReal 0 1) (if j = i then A else Set.univ))
        = (∏ j, if j = i then (gaussianReal 0 1) A else 1) := by
            apply Finset.prod_congr rfl
            intro j _
            by_cases hji : j = i
            · simp [hji, h1]
            · simp [hji, h1]
      _ = (gaussianReal 0 1) A := by
          rw [Finset.prod_ite_eq']
          simp
  rw [hval]

end Unit2B.FormalR1.Kappa.ProofsR4
