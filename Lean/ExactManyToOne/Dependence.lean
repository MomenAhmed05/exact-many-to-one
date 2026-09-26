/-
Copyright (c) 2026 Momen Ahmed. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Momen Ahmed
-/
import Mathlib

/-!
# Sections 3.2, 9.10 and 11: loadings, whitening and the factor sign

* The numerators' correlation is the product of the loadings (equation (11)).
* Whitening an AR(1) series by `yₜ = xₜ − ϕ xₜ₋₁` and **dropping** the first
  observation gives a constant mean `μ(1 − ϕ)` and exactly the innovations as noise,
  so the whitened coordinates are independent with a common law whenever the
  innovations are (equation (45) and the repair of Section 11.1).
* Rescaling the first observation instead gives a mean vector that is not constant,
  and two groups of different lengths then have different whitened means, so the
  numerator is biased by a term proportional to the common mean. This is the failure
  the falsifier of Section 11.1 found.
* With the coefficient estimated per group, the numerator acquires `μ(ϕ̂₀ − ϕ̂ᵢ)`
  (equation (49)), and pooling one estimate across both groups removes it.
* Under a common factor, the numerator's variance and the denominator's expectation
  differ by `2βᵢβ₀σ_f²/n` (equation (50)), so same-sign loadings make the
  independent-groups rule conservative and opposite-sign loadings make it liberal.
-/

open MeasureTheory ProbabilityTheory

namespace ExactManyToOne

/-! ## Equation (11): the one-factor structure of the numerators -/

/-- The loading `hᵢ = √(w₀/(wᵢ + w₀))`. -/
noncomputable def loading (w0 w : ℝ) : ℝ := Real.sqrt (w0 / (w + w0))

/-- **The correlation of two numerators is the product of their loadings.** The
covariance of `Zᵢ − U` and `Zⱼ − U` is `Var U = w₀` and their variances are
`wᵢ + w₀` and `wⱼ + w₀`. -/
theorem corr_eq_loading_mul {w0 wi wj : ℝ} (h0 : 0 < w0) (hi : 0 ≤ wi) (hj : 0 ≤ wj) :
    w0 / Real.sqrt ((wi + w0) * (wj + w0)) = loading w0 wi * loading w0 wj := by
  have hi' : 0 < wi + w0 := by linarith
  have hj' : 0 < wj + w0 := by linarith
  unfold loading
  rw [← Real.sqrt_mul (div_nonneg h0.le hi'.le), div_mul_div_comm,
    Real.sqrt_div' _ (mul_pos hi' hj').le, Real.sqrt_mul_self h0.le]

/-- The loading lies in `(0, 1]` and equals one exactly for a candidate with no noise
of its own. -/
theorem loading_mem {w0 w : ℝ} (h0 : 0 < w0) (hw : 0 ≤ w) :
    0 < loading w0 w ∧ loading w0 w ≤ 1 := by
  unfold loading
  refine ⟨Real.sqrt_pos.2 (div_pos h0 (by linarith)), ?_⟩
  rw [Real.sqrt_le_one]
  exact (div_le_one (by linarith)).2 (by linarith)

/-! ## Section 11.1: whitening a declared AR(1) -/

/-- **Whitening identity.** If `xₜ = μ + eₜ` with `eₜ₊₁ = ϕ eₜ + ηₜ₊₁`, the whitened value
`xₜ₊₁ − ϕ xₜ` is the constant `μ(1 − ϕ)` plus the innovation, with no remainder. -/
theorem whiten_eq (μ ϕ : ℝ) (x e η : ℕ → ℝ) (hx : ∀ t, x t = μ + e t)
    (hAR : ∀ t, e (t + 1) = ϕ * e t + η (t + 1)) (t : ℕ) :
    x (t + 1) - ϕ * x t = μ * (1 - ϕ) + η (t + 1) := by
  rw [hx, hx, hAR]; ring

/-- **The whitened sample is exactly independent with a common law.** With the first
observation dropped, the whitened coordinates are `μ(1 − ϕ) + ηₜ`, a fixed measurable
function of the innovations alone, so they inherit the innovations' mutual
independence exactly, not asymptotically. -/
theorem whiten_iIndep {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} (μ ϕ : ℝ)
    (x e η : ℕ → Ω → ℝ) (hx : ∀ t ω, x t ω = μ + e t ω)
    (hAR : ∀ t ω, e (t + 1) ω = ϕ * e t ω + η (t + 1) ω)
    (hη : iIndepFun (fun t : ℕ => η (t + 1)) P) :
    iIndepFun (fun t : ℕ => fun ω => x (t + 1) ω - ϕ * x t ω) P := by
  have : (fun t : ℕ => fun ω => x (t + 1) ω - ϕ * x t ω)
      = fun t : ℕ => (fun y => μ * (1 - ϕ) + y) ∘ η (t + 1) := by
    funext t ω
    exact whiten_eq μ ϕ (fun s => x s ω) (fun s => e s ω) (fun s => η s ω)
      (fun s => hx s ω) (fun s => hAR s ω) t
  rw [this]
  exact hη.comp _ fun _ => by fun_prop

/-- **Group-specific coefficients: rescale by `1 − ϕ`.** With different declared
coefficients the whitened means `μ(1 − ϕ₀)` and `μ(1 − ϕᵢ)` differ under the null by
`μ(ϕ₀ − ϕᵢ)` (the same term as `estimated_mean_term`), so the dropped-first transform
alone needs one coefficient shared by every group. Dividing each group's whitened values
by its own `1 − ϕ` restores the common mean `μ` for any coefficients, leaving
independent normals with group-specific variances: the unequal-variance location model
the rest of the development handles. -/
theorem whiten_rescaled (μ ϕ : ℝ) (hϕ : ϕ ≠ 1) (x e η : ℕ → ℝ) (hx : ∀ t, x t = μ + e t)
    (hAR : ∀ t, e (t + 1) = ϕ * e t + η (t + 1)) (t : ℕ) :
    (x (t + 1) - ϕ * x t) / (1 - ϕ) = μ + η (t + 1) / (1 - ϕ) := by
  have h : (1 - ϕ) ≠ 0 := sub_ne_zero.2 (Ne.symm hϕ)
  rw [whiten_eq μ ϕ x e η hx hAR t, add_div, mul_div_assoc, div_self h, mul_one]

/-- **The dropped-first transform has a constant mean.** Every whitened coordinate has
mean `μ(1 − ϕ)` whatever the group's length, so two groups share it and the
difference of whitened means is zero under the null. -/
theorem dropfirst_mean_diff (μ ϕ : ℝ) (n0 n1 : ℕ) (hn0 : 0 < n0) (hn1 : 0 < n1) :
    (∑ _t ∈ Finset.range n1, μ * (1 - ϕ)) / n1
      - (∑ _t ∈ Finset.range n0, μ * (1 - ϕ)) / n0 = 0 := by
  simp only [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  field_simp
  ring

/-- Mean of a Prais–Winsten whitened sample of length `n`: the rescaled first
coordinate has mean `√(1 − ϕ²) μ` and the other `n − 1` have mean `(1 − ϕ) μ`. -/
noncomputable def pwMean (μ ϕ : ℝ) (n : ℕ) : ℝ :=
  μ * (Real.sqrt (1 - ϕ ^ 2) + (n - 1) * (1 - ϕ)) / n

/-- **The rescaled first coordinate breaks the constant mean.** For `|ϕ| < 1`, the
first coordinate's multiplier `√(1 − ϕ²)` equals the others' `1 − ϕ` only at `ϕ = 0`. -/
theorem pw_first_eq_iff {ϕ : ℝ} (h1 : -1 < ϕ) (h2 : ϕ < 1) :
    Real.sqrt (1 - ϕ ^ 2) = 1 - ϕ ↔ ϕ = 0 := by
  constructor
  · intro h
    have hsq : 1 - ϕ ^ 2 = (1 - ϕ) ^ 2 := by
      rw [← h, Real.sq_sqrt]; nlinarith
    have : ϕ * (ϕ - 1) = 0 := by nlinarith
    rcases mul_eq_zero.1 this with h0 | h0
    · exact h0
    · exact absurd (by linarith : ϕ = 1) (ne_of_lt h2)
  · rintro rfl; simp

/-- **Unequal lengths bias the Prais–Winsten numerator.** The difference of the two
groups' whitened means under the null is `μ (√(1 − ϕ²) − (1 − ϕ)) (1/n₁ − 1/n₀)`,
which is proportional to the common mean. -/
theorem pw_mean_diff (μ ϕ : ℝ) (n0 n1 : ℕ) (hn0 : 0 < n0) (hn1 : 0 < n1) :
    pwMean μ ϕ n1 - pwMean μ ϕ n0
      = μ * (Real.sqrt (1 - ϕ ^ 2) - (1 - ϕ)) * (1 / n1 - 1 / n0) := by
  unfold pwMean
  have h0 : (n0 : ℝ) ≠ 0 := by exact_mod_cast hn0.ne'
  have h1 : (n1 : ℝ) ≠ 0 := by exact_mod_cast hn1.ne'
  field_simp
  ring

/-- ... and that bias is nonzero exactly when the common mean is nonzero, the
coefficient is nonzero, and the lengths differ. -/
theorem pw_mean_diff_ne_zero {μ ϕ : ℝ} {n0 n1 : ℕ} (hn0 : 0 < n0) (hn1 : 0 < n1)
    (hμ : μ ≠ 0) (h1 : -1 < ϕ) (h2 : ϕ < 1) (hϕ : ϕ ≠ 0) (hn : n0 ≠ n1) :
    pwMean μ ϕ n1 - pwMean μ ϕ n0 ≠ 0 := by
  rw [pw_mean_diff μ ϕ n0 n1 hn0 hn1]
  refine mul_ne_zero (mul_ne_zero hμ ?_) ?_
  · exact sub_ne_zero.2 fun h => hϕ ((pw_first_eq_iff h1 h2).1 h)
  · have h0 : (n0 : ℝ) ≠ 0 := by exact_mod_cast hn0.ne'
    have h1' : (n1 : ℝ) ≠ 0 := by exact_mod_cast hn1.ne'
    rw [sub_ne_zero]
    intro h
    apply hn
    have : (n1 : ℝ) = n0 := by
      field_simp at h; linarith
    exact_mod_cast this.symm

/-- **Equation (49).** Whitening each group at its own estimated coefficient splits
the transformed value into `μ(1 − ϕ̂)` and a part free of `μ`, so the difference of
the two groups' transformed means carries `μ(ϕ̂₀ − ϕ̂ᵢ)`. -/
theorem whiten_estimated (μ ϕh : ℝ) (x e : ℕ → ℝ) (hx : ∀ t, x t = μ + e t) (t : ℕ) :
    x (t + 1) - ϕh * x t = μ * (1 - ϕh) + (e (t + 1) - ϕh * e t) := by
  rw [hx, hx]; ring

theorem estimated_mean_term (μ ϕ0 ϕi : ℝ) : μ * (1 - ϕi) - μ * (1 - ϕ0) = μ * (ϕ0 - ϕi) := by
  ring

/-- **Pooling removes the common mean.** One coefficient applied to both groups gives
both the same multiplier, so the mean term cancels whatever the estimate is. -/
theorem pooled_mean_term (μ ϕ : ℝ) : μ * (1 - ϕ) - μ * (1 - ϕ) = 0 := sub_self _

/-! ## Section 11.2 and 11.6: a common factor -/

/-- Residualising on a declared loading returns the factor alpha plus the idiosyncratic
noise: the transformed series is independent of the factor by construction. -/
theorem residualise (α β f ε : ℝ) : (α + β * f + ε) - β * f = α + ε := by ring

/-- Variance of the numerator under one factor, `((βᵢ − β₀)²σ_f² + σᵢ² + σ₀²)/n`. -/
noncomputable def numVar (bi b0 sf2 si2 s02 n : ℝ) : ℝ := ((bi - b0) ^ 2 * sf2 + si2 + s02) / n

/-- Expectation of the denominator's square, `((βᵢ² + β₀²)σ_f² + σᵢ² + σ₀²)/n`. -/
noncomputable def denExp (bi b0 sf2 si2 s02 n : ℝ) : ℝ := ((bi ^ 2 + b0 ^ 2) * sf2 + si2 + s02) / n

/-- **Equation (50).** The two differ by exactly `2βᵢβ₀σ_f²/n`. -/
theorem den_sub_num (bi b0 sf2 si2 s02 n : ℝ) :
    denExp bi b0 sf2 si2 s02 n - numVar bi b0 sf2 si2 s02 n = 2 * bi * b0 * sf2 / n := by
  unfold denExp numVar; ring

/-- **Same-sign loadings: the rule is conservative.** The denominator overstates the
numerator's variance. -/
theorem same_sign_conservative {bi b0 sf2 n : ℝ} (si2 s02 : ℝ) (hb : 0 < bi * b0)
    (hf : 0 < sf2) (hn : 0 < n) : numVar bi b0 sf2 si2 s02 n < denExp bi b0 sf2 si2 s02 n := by
  have := den_sub_num bi b0 sf2 si2 s02 n
  have hpos : 0 < 2 * bi * b0 * sf2 / n := by
    apply div_pos _ hn; nlinarith
  linarith

/-- **Opposite-sign loadings: the rule is liberal.** The denominator understates it. -/
theorem opposite_sign_liberal {bi b0 sf2 n : ℝ} (si2 s02 : ℝ) (hb : bi * b0 < 0)
    (hf : 0 < sf2) (hn : 0 < n) : denExp bi b0 sf2 si2 s02 n < numVar bi b0 sf2 si2 s02 n := by
  have := den_sub_num bi b0 sf2 si2 s02 n
  have hneg : 2 * bi * b0 * sf2 / n < 0 := by
    apply div_neg_of_neg_of_pos _ hn; nlinarith
  linarith

/-- **An unloaded benchmark: the two agree.** This is why the closed-form scale factor
of (36) is exact when the benchmark carries no loading. -/
theorem unloaded_benchmark (bi sf2 si2 s02 n : ℝ) :
    numVar bi 0 sf2 si2 s02 n = denExp bi 0 sf2 si2 s02 n := by
  unfold numVar denExp; ring

end ExactManyToOne
