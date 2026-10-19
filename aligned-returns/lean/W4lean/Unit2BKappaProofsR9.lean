import W4lean.Unit2BFormalR1Proofs
-- D″, item 6, sub-lemma D-6a `eig2`.
-- A unit eigenvector `(c, s)` of `Gmat θ` for the eigenvalue `mu1 θ`.  One split, on
-- `g01 θ = 0 ∧ g11 θ ≤ g00 θ`, exactly where `v = (g01 θ, mu1 θ - g00 θ)` vanishes.
-- Only the `mu1` eigen-equation is proved here; the `mu2` one is D-6b's.
-- Project lemmas are used FULLY QUALIFIED as `Unit2B.FormalR1.Proofs.*`: the unqualified
-- names resolve, from this namespace, to the `sorry` statements in `Unit2BFormalR1.lean`.
-- Signatures confirmed by source inspection:
--   project: `Unit2B.FormalR1.Proofs.discG_eq` (Unit2BFormalR1Proofs.lean:64),
--   `Unit2B.FormalR1.Proofs.mu_sum_prod` (:72), `Unit2B.FormalR1.Proofs.charQ_factor` (:82);
--   Mathlib 5ed2965256430c3649e86755f9576b54eca72435: `Real.sqrt_sq` (Analysis/Real/Sqrt.lean:181,
--   `(h : 0 ≤ x) : √(x ^ 2) = x`), `Real.sq_sqrt` (:178, `(h : 0 ≤ x) : √x ^ 2 = x`),
--   `Real.sqrt_pos` (:286, `0 < √x ↔ 0 < x`), `Real.sqrt_nonneg` (:144),
--   `Matrix.mulVec` (Data/Matrix/Mul.lean:698), `dotProduct` (:72),
--   `Fin.sum_univ_two` (to_additive of Algebra/BigOperators/Fin.lean:111),
--   `div_pow` (Algebra/Group/Basic.lean:598), `add_div` (Algebra/Field/Basic.lean:36),
--   `div_self` (Algebra/GroupWithZero/Units/Basic.lean:294), `pow_ne_zero`
--   (Algebra/GroupWithZero/Basic.lean:259), `not_le` (Order/Defs/LinearOrder.lean:119),
--   `lt_of_le_of_ne` (Order/Defs/PartialOrder.lean:200).
-- `*ᵥ` is scoped (Mul.lean:702), so `Matrix.mulVec` is written out.

open MeasureTheory ProbabilityTheory
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR5

/-- The two scalar eigen-equations for `mu1 θ` give the vector equation. -/
theorem gmat_mulVec_of {θ : Param} {c s : ℝ}
    (h1 : g00 θ * c + g01 θ * s = mu1 θ * c) (h2 : g01 θ * c + g11 θ * s = mu1 θ * s) :
    Matrix.mulVec (Gmat θ) ![c, s] = mu1 θ • ![c, s] := by
  ext i
  fin_cases i <;> simp [Gmat, Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;>
    first | linear_combination h1 | linear_combination h2

/-- INSIDE the split, `mu1 θ = g00 θ`, via `discG_eq` and `Real.sqrt_sq`. -/
theorem mu1_eq_g00_of_inside {θ : Param} (h0 : g01 θ = 0) (hle : g11 θ ≤ g00 θ) :
    mu1 θ = g00 θ := by
  have hD : discG θ = (g00 θ - g11 θ) ^ 2 := by
    rw [Unit2B.FormalR1.Proofs.discG_eq, h0]
    ring
  have hr : Real.sqrt (discG θ) = g00 θ - g11 θ := by
    rw [hD]
    exact Real.sqrt_sq (by linarith)
  unfold mu1 trG
  rw [hr]
  ring

/-- INSIDE the split, consequently `mu2 θ = g11 θ`. -/
theorem mu2_eq_g11_of_inside {θ : Param} (h0 : g01 θ = 0) (hle : g11 θ ≤ g00 θ) :
    mu2 θ = g11 θ := by
  have hs := (Unit2B.FormalR1.Proofs.mu_sum_prod θ).1
  have hm := mu1_eq_g00_of_inside h0 hle
  unfold trG at hs
  linarith

/-- **D-6a, INSIDE branch.**  On `g01 θ = 0 ∧ g11 θ ≤ g00 θ`, `e₀ = (1, 0)` is a unit
eigenvector of `Gmat θ` for `mu1 θ`. -/
theorem eig2_inside {θ : Param} (h0 : g01 θ = 0) (hle : g11 θ ≤ g00 θ) :
    ∃ c s : ℝ, c ^ 2 + s ^ 2 = 1 ∧ Matrix.mulVec (Gmat θ) ![c, s] = mu1 θ • ![c, s] := by
  have hm := mu1_eq_g00_of_inside h0 hle
  refine ⟨1, 0, by norm_num, gmat_mulVec_of ?_ ?_⟩
  · linear_combination (-1 : ℝ) * hm
  · linear_combination h0

end Unit2B.FormalR1.Kappa.ProofsR5

#print axioms Unit2B.FormalR1.Kappa.ProofsR5.gmat_mulVec_of
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.mu1_eq_g00_of_inside
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.mu2_eq_g11_of_inside
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.eig2_inside

-- OUTSIDE region, added below the banked INSIDE region (lines 1-69, sha256 9cff031b...).

namespace Unit2B.FormalR1.Kappa.ProofsR5

/-- `mu1 θ` is a root of the characteristic polynomial, from `charQ_factor`. -/
theorem mu1_charQ (θ : Param) : (mu1 θ - g00 θ) * (mu1 θ - g11 θ) = g01 θ ^ 2 := by
  have h := Unit2B.FormalR1.Proofs.charQ_factor (mu1 θ) θ
  unfold charQ trG detG at h
  linear_combination h

/-- OUTSIDE the split, `v = (g01 θ, mu1 θ - g00 θ) ≠ 0`, stated as `0 < ‖v‖²`.  Either
`g01 θ ≠ 0`, or `g01 θ = 0` and `g00 θ < g11 θ`, whence `mu1 θ > g00 θ`. -/
theorem v_sq_pos_of_outside {θ : Param} (hout : ¬ (g01 θ = 0 ∧ g11 θ ≤ g00 θ)) :
    0 < g01 θ ^ 2 + (mu1 θ - g00 θ) ^ 2 := by
  by_cases h0 : g01 θ = 0
  · have hlt : g00 θ < g11 θ := not_le.mp (fun hle => hout ⟨h0, hle⟩)
    have hs := Real.sqrt_nonneg (discG θ)
    have hv : 0 < mu1 θ - g00 θ := by
      unfold mu1 trG
      linarith
    have hv2 := pow_pos hv 2
    have hb := sq_nonneg (g01 θ)
    linarith
  · have hsq : 0 < g01 θ ^ 2 := lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 h0))
    have hb := sq_nonneg (mu1 θ - g00 θ)
    linarith

/-- **D-6a, OUTSIDE branch.**  Off `g01 θ = 0 ∧ g11 θ ≤ g00 θ`, `e₀ = v / ‖v‖` with
`v = (g01 θ, mu1 θ - g00 θ)` is a unit eigenvector of `Gmat θ` for `mu1 θ`. -/
theorem eig2_outside {θ : Param} (hout : ¬ (g01 θ = 0 ∧ g11 θ ≤ g00 θ)) :
    ∃ c s : ℝ, c ^ 2 + s ^ 2 = 1 ∧ Matrix.mulVec (Gmat θ) ![c, s] = mu1 θ • ![c, s] := by
  have hpos := v_sq_pos_of_outside hout
  have hq := mu1_charQ θ
  obtain ⟨n, hn, hn2⟩ : ∃ n : ℝ, 0 < n ∧ n ^ 2 = g01 θ ^ 2 + (mu1 θ - g00 θ) ^ 2 :=
    ⟨Real.sqrt (g01 θ ^ 2 + (mu1 θ - g00 θ) ^ 2), Real.sqrt_pos.mpr hpos, Real.sq_sqrt hpos.le⟩
  refine ⟨g01 θ / n, (mu1 θ - g00 θ) / n, ?_, gmat_mulVec_of ?_ ?_⟩
  · rw [div_pow, div_pow, ← add_div, ← hn2]
    exact div_self (pow_ne_zero 2 hn.ne')
  · ring
  · linear_combination (-1 / n) * hq

/-- **D-6a `eig2`.**  A unit eigenvector `(c, s)` of `Gmat θ` for the eigenvalue `mu1 θ`.
One split, on `g01 θ = 0 ∧ g11 θ ≤ g00 θ`.  No `Box` hypothesis is needed. -/
theorem eig2 (θ : Param) :
    ∃ c s : ℝ, c ^ 2 + s ^ 2 = 1 ∧ Matrix.mulVec (Gmat θ) ![c, s] = mu1 θ • ![c, s] := by
  by_cases hin : g01 θ = 0 ∧ g11 θ ≤ g00 θ
  · exact eig2_inside hin.1 hin.2
  · exact eig2_outside hin

end Unit2B.FormalR1.Kappa.ProofsR5

#print axioms Unit2B.FormalR1.Kappa.ProofsR5.mu1_charQ
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.v_sq_pos_of_outside
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.eig2_outside
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.eig2
