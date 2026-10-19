import W4lean.Unit2BKappaProofsR1
-- REQUIRED, and not merely convenient.  `Unit2BFormalR1.lean:65` declares a
-- sorry-backed `box_pos` in namespace `Unit2B.FormalR1`, while
-- `Unit2BFormalR1Proofs.lean:88` declares the proved one in
-- `Unit2B.FormalR1.Proofs`.  The kappa chain does not otherwise import the proofs
-- file, so a bare `box_pos` silently resolves to the sorry version and poisons the
-- axiom closure with `sorryAx`.  Importing the proofs file makes the namespace
-- hierarchy prefer the proved declaration.  See deposits kappa01_sorry_leak.txt.
import W4lean.Unit2BFormalR1Proofs

open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR2

open Unit2B.FormalR1.Kappa

/-- **Stage B, sub-lemma B-1.**  The K-4 event identity abstracted to reals.

Both sides reduce to the same inequality.  Writing `V = s0 + s1`, the left side is
`c < d / √(V / n)`, i.e. after `Real.sqrt_div` and clearing the positive `√V`,
`c * √V < d * √n`.  The right side is `c * √(ν V / (ν a))`, i.e. after `field_simp`,
`c * √(V / a)`, and after `Real.sqrt_div` and clearing the positive `√a` it is again
`c * √V < √n * d`.  No hypothesis is used beyond the stated nonnegativities. -/
theorem event_iff_abstract {d s0 s1 a c n nu : ℝ}
    (hn : 0 < n) (ha : 0 < a) (hc : 0 < c) (hnu : 0 < nu)
    (_hs0 : 0 ≤ s0) (_hs1 : 0 ≤ s1) (hV : 0 < s0 + s1) :
    c < d / Real.sqrt (s1 / n + s0 / n) ↔
      c * Real.sqrt (nu * (s0 + s1) / (nu * a)) < Real.sqrt n * d / Real.sqrt a := by
  have hsV : 0 ≤ s0 + s1 := le_of_lt hV
  have hsqrtV : 0 < Real.sqrt (s0 + s1) := Real.sqrt_pos.2 hV
  have hsqrtN : 0 < Real.sqrt n := Real.sqrt_pos.2 hn
  have hsqrtA : 0 < Real.sqrt a := Real.sqrt_pos.2 ha
  have hX : 0 < Real.sqrt (s0 + s1) / Real.sqrt n := div_pos hsqrtV hsqrtN
  rw [show s1 / n + s0 / n = (s0 + s1) / n from by ring]
  rw [show nu * (s0 + s1) / (nu * a) = (s0 + s1) / a from by field_simp]
  rw [Real.sqrt_div hsV n, Real.sqrt_div hsV a]
  constructor
  · intro h
    have h1 : c * (Real.sqrt (s0 + s1) / Real.sqrt n) < d := (lt_div_iff₀ hX).mp h
    have h2 : c * Real.sqrt (s0 + s1) < d * Real.sqrt n := by
      calc c * Real.sqrt (s0 + s1)
          = c * (Real.sqrt (s0 + s1) / Real.sqrt n) * Real.sqrt n := by
            field_simp
        _ < d * Real.sqrt n := mul_lt_mul_of_pos_right h1 hsqrtN
    rw [lt_div_iff₀ hsqrtA]
    calc c * (Real.sqrt (s0 + s1) / Real.sqrt a) * Real.sqrt a
        = c * Real.sqrt (s0 + s1) := by
          field_simp
      _ < d * Real.sqrt n := h2
      _ = Real.sqrt n * d := by ring
  · intro h
    have h := (lt_div_iff₀ hsqrtA).mp h
    have h2 : c * Real.sqrt (s0 + s1) < d * Real.sqrt n := by
      calc c * Real.sqrt (s0 + s1)
          = c * (Real.sqrt (s0 + s1) / Real.sqrt a) * Real.sqrt a := by
            field_simp
        _ < Real.sqrt n * d := h
        _ = d * Real.sqrt n := by ring
    have h1 : c * (Real.sqrt (s0 + s1) / Real.sqrt n) < d := by
      calc c * (Real.sqrt (s0 + s1) / Real.sqrt n)
          = c * Real.sqrt (s0 + s1) / Real.sqrt n := by ring
        _ < d * Real.sqrt n / Real.sqrt n := div_lt_div_of_pos_right h2 hsqrtN
        _ = d := by field_simp
    exact (lt_div_iff₀ hX).mpr h1

/-- **Stage B, sub-lemma B-2.**  An unbiased sample variance is nonnegative. -/
theorem sVar_nonneg (x : Fin nObs → ℝ) : 0 ≤ sVar x := by
  have hd : 0 < ((nObs : ℝ) - 1) := by norm_num [nObs]
  have hs : 0 ≤ (∑ t, (x t - sMean x) ^ 2) :=
    Finset.sum_nonneg fun t _ => sq_nonneg _
  exact div_nonneg hs (le_of_lt hd)

/-- **Stage B, sub-lemma B-3.**  The original K-4 statement, obtained by instantiating
B-1.  This file's version carries a proof; the statement-file version at line 44 of
`Unit2BKappaStatementR1.lean` keeps its `sorry` and is never edited. -/
theorem kappa_event {m c : ℝ} {θ : Param} (hθ : θ ∈ Box) (hc : 0 < c) {z : Panel}
    (hS : 0 < Sres m θ z) :
    c < Tstat m θ z ↔ c * Real.sqrt (Sres m θ z / (nu * aF θ)) < Nstd m θ z := by
  have ha : 0 < aF θ := (Unit2B.FormalR1.Proofs.box_pos hθ).1
  have hnu : 0 < nu := by
    rw [Unit2B.FormalR1.Kappa.Proofs.nu_eq]; norm_num [nObs]
  have hs0 : 0 ≤ sVar (ret0 m θ z) := sVar_nonneg _
  have hs1 : 0 ≤ sVar (ret1 m θ z) := sVar_nonneg _
  have hV : 0 < sVar (ret0 m θ z) + sVar (ret1 m θ z) := by
    have h : 0 < nu * (sVar (ret0 m θ z) + sVar (ret1 m θ z)) := hS
    exact (mul_pos_iff_of_pos_left hnu).mp h
  show c < (sMean (ret1 m θ z) - sMean (ret0 m θ z)) /
      Real.sqrt (sVar (ret1 m θ z) / (nObs : ℝ) + sVar (ret0 m θ z) / (nObs : ℝ)) ↔
      c * Real.sqrt (nu * (sVar (ret0 m θ z) + sVar (ret1 m θ z)) / (nu * aF θ)) <
        Real.sqrt (nObs : ℝ) * (sMean (ret1 m θ z) - sMean (ret0 m θ z)) / Real.sqrt (aF θ)
  exact event_iff_abstract (by norm_num [nObs]) ha hc hnu hs0 hs1 hV

end Unit2B.FormalR1.Kappa.ProofsR2
