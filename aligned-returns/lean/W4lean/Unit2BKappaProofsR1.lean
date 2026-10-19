import W4lean.Unit2BKappaStatementR1
open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.Kappa.Proofs

/-- K-0: bridge constants. `nu = 35` by definition and `nObs = 36`, so this is arithmetic. -/
theorem nu_eq : nu = (nObs : ℝ) - 1 := by
  norm_num [nu, nObs]

/-- K-1: the model's Gram matrix is §6's `Gmat`.

  `Mmat = !![b₀, 1, 0; b₁, 0, √λ]`, so `Mmat * Mmatᵀ` is
  `![b₀²+1, b₀b₁; b₀b₁, b₁²+λ]`, which is `Gmat` entry by entry. No hypothesis on θ is
  needed, so `hθ` is unused. -/
theorem gram_eq {θ : Param} (hθ : θ ∈ Box) : Mmat θ * (Mmat θ).transpose = Gmat θ := by
  have h0 : (0 : ℝ) ≤ θ.1 := by
    have hx := hθ
    simp only [Box, Set.mem_ofPred_eq] at hx
    have h1 : (1 / 25 : ℝ) ≤ θ.1 := hx.1
    linarith
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Mmat, Gmat, g00, g01, g11, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_succ, Real.mul_self_sqrt h0]
  all_goals ring

end Unit2B.FormalR1.Kappa.Proofs
