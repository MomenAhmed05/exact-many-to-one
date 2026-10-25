import Mathlib

/-!
# Option B: the structure the calibrated-rule certificate rests on (two candidates)

The certificate bounds the family-wise error rate of a table rule: reject candidate `i` when its
paired statistic exceeds the constant `c(r̂)` of the bin that the estimated correlation `r̂`
falls in.
Three facts turn that into one-dimensional bounds; this file proves them.

* `fwer_exchangeable`: if the two rejection events have the same probability (the two candidates are
  exchangeable), the family-wise error rate is twice one candidate's rate minus the joint rate.
* `abel_bins`: one candidate's rejection rate, written as a sum over bins of the bin's rejection
  rate times the bin's probability, equals the tail form `g 0 + ∑ (g (j+1) - g j) * U (j+1)`, where
  `U j` is the probability that `r̂` lies in bin `j` or above.
* `tail_form_mono`: if the bin rates `g j` are nondecreasing in `j` (the constants are
  nonincreasing) and each tail `U (j+1)` is nondecreasing in a parameter, the tail form is
  nondecreasing in that parameter. With the parameter the correlation, this is why the certificate
  may evaluate one candidate's rate at the upper end of each correlation interval.

That the tails are nondecreasing in the correlation is a property of the noncentral t law of `r̂`
given the first variance, established in the paper; it enters `tail_form_mono` as the
hypothesis `hU`.
-/

open MeasureTheory Finset

namespace UnitPB

/-- Two exchangeable rejection events: `P(A ∪ B) = 2 P(A) - P(A ∩ B)`. -/
theorem fwer_exchangeable {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsFiniteMeasure μ]
    {A B : Set Ω} (hB : MeasurableSet B) (hAB : μ A = μ B) :
    μ.real (A ∪ B) = 2 * μ.real A - μ.real (A ∩ B) := by
  have h := congrArg ENNReal.toReal (measure_union_add_inter A hB (μ := μ))
  rw [ENNReal.toReal_add (measure_ne_top _ _) (measure_ne_top _ _),
    ENNReal.toReal_add (measure_ne_top _ _) (measure_ne_top _ _)] at h
  have hAB' : μ.real A = μ.real B := by simp [Measure.real, hAB]
  simp only [Measure.real] at hAB' ⊢
  linarith

/-- Summation by parts over the bins. `U j` is the probability of bin `j` or above (`U 0 = 1` in
use), so bin `j < N` has probability `U j - U (j+1)` and the last bin `N` has probability `U N`. -/
theorem abel_bins (g U : ℕ → ℝ) (N : ℕ) :
    (∑ j ∈ range N, g j * (U j - U (j + 1))) + g N * U N
      = g 0 * U 0 + ∑ j ∈ range N, (g (j + 1) - g j) * U (j + 1) := by
  induction N with
  | zero => simp
  | succ N ih =>
    rw [sum_range_succ, sum_range_succ]
    linarith [ih]

/-- The tail form is nondecreasing in the parameter when the bin rates are nondecreasing in the bin
index and every tail above the first bin is nondecreasing in the parameter. -/
theorem tail_form_mono (g : ℕ → ℝ) (hg : Monotone g) (U : ℕ → ℝ → ℝ)
    (hU0 : ∀ ρ, U 0 ρ = 1) (hU : ∀ j, Monotone (U (j + 1))) (N : ℕ) :
    Monotone (fun ρ => g 0 * U 0 ρ + ∑ j ∈ range N, (g (j + 1) - g j) * U (j + 1) ρ) := by
  intro a b hab
  simp only [hU0]
  gcongr with j hj
  · exact sub_nonneg.mpr (hg (Nat.le_succ j))
  · exact hU j hab

#print axioms UnitPB.fwer_exchangeable
#print axioms UnitPB.abel_bins
#print axioms UnitPB.tail_form_mono

end UnitPB
