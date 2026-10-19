import W4lean.Unit2BKappaProofsR2
-- Fully qualified for the same reason as in R2: `Unit2BFormalR1.lean` declares
-- sorry-backed twins of ten lemmas in namespace `Unit2B.FormalR1`, and the bare names
-- resolve to those.  See deposits/kappa01_sorry_leak.txt.
import W4lean.Unit2BFormalR1Proofs

open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR3

open Unit2B.FormalR1.Kappa
open Unit2B.FormalR1.Kappa.ProofsR2

/-!
Stage C: the pathwise identity at general `n`, then the original K-2.

**Stage C in its own file.**  It is placed in a new file rather than in
`Unit2BKappaProofsR2.lean`, so that a failure here cannot regress the K-4 proof already
in R2, which is exactly the failure mode that forced earlier K-2 and K-4 attempts to be
reverted.

D-1 (general n) is followed: nothing below mentions `nObs` until the bridge in C-4.
-/

section StageC

variable {n : ℕ}

def sMeanN (n : ℕ) (x : Fin n → ℝ) : ℝ := (∑ t, x t) / (n : ℝ)

def sVarN (n : ℕ) (x : Fin n → ℝ) : ℝ :=
  (∑ t, (x t - sMeanN n x) ^ 2) / ((n : ℝ) - 1)

def ret0N (m : ℝ) (θ : Param) (z : (Fin 3 × Fin n) → ℝ) (t : Fin n) : ℝ :=
  m + θ.2.1 * z (0, t) + z (1, t)

def ret1N (m : ℝ) (θ : Param) (z : (Fin 3 × Fin n) → ℝ) (t : Fin n) : ℝ :=
  m + θ.2.2 * z (0, t) + Real.sqrt θ.1 * z (2, t)

def zcN (z : (Fin 3 × Fin n) → ℝ) (t : Fin n) : Fin 3 → ℝ :=
  fun k => z (k, t) - sMeanN n (fun s => z (k, s))

theorem nR_ne_zero (hn : 0 < n) : (n : ℝ) ≠ 0 := by
  have h1 : (0 : ℝ) < (n : ℝ) := by exact_mod_cast hn
  exact ne_of_gt h1

theorem sum_const_n (m : ℝ) : (∑ _s : Fin n, m) = ((n : ℝ)) * m := by
  rw [Finset.sum_const, Finset.card_fin n]
  simp

theorem sum_mul_congr (a : ℝ) (v : Fin n → ℝ) : (∑ s, a * v s) = a * (∑ s, v s) := by
  calc (∑ s, a * v s) = ∑ s, v s * a :=
        Finset.sum_congr rfl fun s _ => mul_comm _ _
    _ = (∑ s, v s) * a := (Finset.sum_mul _ _ _).symm
    _ = a * (∑ s, v s) := mul_comm _ _

/-- **C-1.**  A centred vector sums to zero.  This is the lemma that kills the cross
terms in C-2; it is the only reason the pathwise identity is easier than it looks. -/
theorem sum_centered_zero (hn : 0 < n) (x : Fin n → ℝ) :
    ∑ t, (x t - sMeanN n x) = 0 := by
  rw [Finset.sum_sub_distrib, sum_const_n (sMeanN n x)]
  simp only [sMeanN]
  field_simp [nR_ne_zero hn] <;> ring

/-- **C-2, first half.**  Centring removes the common mean and the linear part, leaving
exactly the two rows of `Mmat θ` acting on the centred coordinates. -/
theorem ret0N_centered (m : ℝ) (θ : Param) (z : (Fin 3 × Fin n) → ℝ) (hn : 0 < n) (t : Fin n) :
    ret0N m θ z t - sMeanN n (ret0N m θ z) = θ.2.1 * zcN z t 0 + zcN z t 1 := by
  have hs : sMeanN n (ret0N m θ z)
      = m + θ.2.1 * sMeanN n (fun s => z (0, s)) + sMeanN n (fun s => z (1, s)) := by
    simp only [ret0N, sMeanN, Finset.sum_add_distrib]
    rw [sum_const_n, sum_mul_congr θ.2.1 (fun s => z (0, s))]
    field_simp [nR_ne_zero hn]
  rw [hs]
  simp only [ret0N, zcN, sMeanN]
  ring

/-- **C-2, second half.** -/
theorem ret1N_centered (m : ℝ) (θ : Param) (z : (Fin 3 × Fin n) → ℝ) (hn : 0 < n) (t : Fin n) :
    ret1N m θ z t - sMeanN n (ret1N m θ z) = θ.2.2 * zcN z t 0 + Real.sqrt θ.1 * zcN z t 2 := by
  have hs : sMeanN n (ret1N m θ z)
      = m + θ.2.2 * sMeanN n (fun s => z (0, s))
        + Real.sqrt θ.1 * sMeanN n (fun s => z (2, s)) := by
    simp only [ret1N, sMeanN, Finset.sum_add_distrib]
    rw [sum_const_n, sum_mul_congr θ.2.2 (fun s => z (0, s)),
        sum_mul_congr (Real.sqrt θ.1) (fun s => z (2, s))]
    field_simp [nR_ne_zero hn]
  rw [hs]
  simp only [ret1N, zcN, sMeanN]
  ring

/-- The two rows of `Mmat θ` acting on a vector, read off the matrix literal. -/
theorem mulVec_Mmat_row0 (θ : Param) (v : Fin 3 → ℝ) :
    Matrix.mulVec (Mmat θ) v 0 = θ.2.1 * v 0 + v 1 := by
  simp [Mmat, Matrix.mulVec, Matrix.vecHead, Matrix.vecTail, Fin.isValue]

theorem mulVec_Mmat_row1 (θ : Param) (v : Fin 3 → ℝ) :
    Matrix.mulVec (Mmat θ) v 1 = θ.2.2 * v 0 + Real.sqrt θ.1 * v 2 := by
  simp [Mmat, Matrix.mulVec, Matrix.vecHead, Matrix.vecTail, Fin.isValue]

/-- `Mmat θ : Matrix (Fin 2) (Fin 3) ℝ` has two rows, so the sum here is over
`Fin 2`, and `Fin.sum_univ_two` collapses it to the two row readings. -/
theorem sum_sq_mulVec (θ : Param) (v : Fin 3 → ℝ) :
    (∑ i, (Matrix.mulVec (Mmat θ) v i) ^ 2)
      = (θ.2.1 * v 0 + v 1) ^ 2 + (θ.2.2 * v 0 + Real.sqrt θ.1 * v 2) ^ 2 := by
  have e0 := mulVec_Mmat_row0 θ v
  have e1 := mulVec_Mmat_row1 θ v
  rw [Fin.sum_univ_two, e0, e1]

/-- **C-3.**  The pathwise identity at general `n`. -/
theorem SresN_eq_gram_form (m : ℝ) (θ : Param) (z : (Fin 3 × Fin n) → ℝ) (hn2 : 2 ≤ n) :
    ((n : ℝ) - 1) * (sVarN n (ret0N m θ z) + sVarN n (ret1N m θ z))
      = ∑ t, ∑ i, (Matrix.mulVec (Mmat θ) (zcN z t) i) ^ 2 := by
  have hn1 : n ≠ 1 := by omega
  have hn : 0 < n := by omega
  have hden : ((n : ℝ) - 1) ≠ 0 := by
    intro h
    have h2 : (n : ℝ) = 1 := by linarith
    have h3 : n = 1 := by exact_mod_cast h2
    exact hn1 h3
  -- The left side, with the (n-1) factors cancelled.
  have hlhs : ((n : ℝ) - 1) * (sVarN n (ret0N m θ z) + sVarN n (ret1N m θ z))
      = (∑ t, (ret0N m θ z t - sMeanN n (ret0N m θ z)) ^ 2)
        + ∑ t, (ret1N m θ z t - sMeanN n (ret1N m θ z)) ^ 2 := by
    simp only [sVarN]
    field_simp
  rw [hlhs]
  -- Replace each summand by its row reading, then split the sum.
  calc (∑ t, (ret0N m θ z t - sMeanN n (ret0N m θ z)) ^ 2)
        + ∑ t, (ret1N m θ z t - sMeanN n (ret1N m θ z)) ^ 2
      = (∑ t, (θ.2.1 * zcN z t 0 + zcN z t 1) ^ 2)
        + ∑ t, (θ.2.2 * zcN z t 0 + Real.sqrt θ.1 * zcN z t 2) ^ 2 := by
        congr 1
        · exact Finset.sum_congr rfl fun t _ => by rw [ret0N_centered m θ z hn t]
        · exact Finset.sum_congr rfl fun t _ => by rw [ret1N_centered m θ z hn t]
    _ = ∑ t, ((θ.2.1 * zcN z t 0 + zcN z t 1) ^ 2
        + (θ.2.2 * zcN z t 0 + Real.sqrt θ.1 * zcN z t 2) ^ 2) := by
        rw [Finset.sum_add_distrib]
    _ = ∑ t, ∑ i, (Matrix.mulVec (Mmat θ) (zcN z t) i) ^ 2 := by
        exact (Finset.sum_congr rfl fun t _ => sum_sq_mulVec θ (zcN z t)).symm

/-- **C-4.**  The original K-2 statement, byte-identical to `Unit2BKappaStatementR1.lean`
lines 38-39, obtained from C-3 by the `n := nObs` bridge.  Every bridge step is `rfl`
except the `nu = nObs - 1` step, which is K-0. -/
theorem Sres_eq_gram_form (m : ℝ) (θ : Param) (z : Panel) :
    Sres m θ z = ∑ t, ∑ i, (Matrix.mulVec (Mmat θ) (zc z t) i) ^ 2 := by
  have hnu : nu = ((nObs : ℝ) - 1) :=
    Unit2B.FormalR1.Kappa.Proofs.nu_eq
  have hmain := SresN_eq_gram_form (n := nObs) m θ z (by norm_num [nObs])
  simp only [Sres, hnu, sVar, sVarN, sMean, sMeanN, ret0, ret0N, ret1, ret1N, zc, zcN]
  exact hmain

end StageC

end Unit2B.FormalR1.Kappa.ProofsR3
