/-
SPDX-License-Identifier: Apache-2.0

# ChiSquaredSpikeR1 — the chi-square law, K-7a

**PROVENANCE AND ATTRIBUTION (Apache-2.0 §4(b), §4(c)).** This file is a MODIFIED copy of
work from the StatLean project:

  repository : https://github.com/StatLean/Stat-Lean
  commit     : 86b3823f77fc2b0f5a9923f85b588a4f8a65f94c   (2026-09-18)
  source file: StatLean/MultipleTesting/ForMathlib/GaussianSquareMGF.lean
  licence    : Apache License 2.0 (https://www.apache.org/licenses/LICENSE-2.0)

**This file was modified from the source**, in the following respects, and those changes are
material:

  * the `StatLean.MultipleTesting` namespace is replaced by
    `Unit2B.FormalR1.Kappa.ProofsR5`;
  * Mathlib's `ProbabilityTheory` namespace is `open`ed rather than fully qualified;
  * the call site of `gaussianPDF_lt_top` is adapted to this Mathlib revision's argument
    convention — the lemma is fully implicit here, and Mathlib's own usage at
    `Probability/Distributions/Gaussian/Real.lean:280` is `ae_of_all _ fun _ ↦ gaussianPDF_lt_top`,
    whereas the source passes `0 1` explicitly.

No other change has been made to this file. Changes belonging to later sub-lemmas — routing `isProbabilityMeasure_map` through `inferInstance` in A-3,
and re-checking `integrableExpSet_eq_of_mgf` in A-2 — are **not** present here and are not
claimed.

The mathematical content is unchanged. The reference bytes are vendored unmodified at
`deposits/statlean_ref_R1/` with recorded sha256 values, so the copy can be diffed against its
source.

Author of the original work: Tingzhou Wei, under Apache-2.0.
-/
import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.Analysis.SpecialFunctions.Gaussian.GaussianIntegral
import Mathlib.Probability.Distributions.Gamma
import Mathlib.Probability.Moments.Basic
import Mathlib.Probability.Moments.ComplexMGF
import Mathlib.Probability.Moments.IntegrableExpMul
import Mathlib.MeasureTheory.Measure.CharacteristicFunction.Basic
import Mathlib.Probability.Independence.Basic
import Mathlib.Tactic.FunProp

/-!
Sub-lemma **A-1**: the MGF of the square of a standard
Gaussian, together with its two integrability halves. Ported from StatLean's
`GaussianSquareMGF.lean`; see the file header for provenance and licence.

This is the chi-squared-1 ingredient. A′ is gated: D′ does not start unless A′ closes clean on
`[propext, Classical.choice, Quot.sound]`.
-/

open MeasureTheory ProbabilityTheory Real
open scoped NNReal

namespace Unit2B.FormalR1.Kappa.ProofsR5

/-- Pointwise identity collapsing the standard-Gaussian density times `exp(l x²)` into a
single Gaussian kernel `(√(2π))⁻¹ · exp(−(1/2 − l) x²)`. -/
private lemma gaussianPDF_mul_exp_sq (l x : ℝ) :
    gaussianPDFReal 0 1 x * Real.exp (l * x ^ 2)
      = (Real.sqrt (2 * π))⁻¹ * Real.exp (-(1 / 2 - l) * x ^ 2) := by
  have hpdf : gaussianPDFReal 0 1 x = (Real.sqrt (2 * π))⁻¹ * Real.exp (-x ^ 2 / 2) := by
    rw [gaussianPDFReal]
    simp only [NNReal.coe_one, mul_one, sub_zero]
  rw [hpdf, mul_assoc, ← Real.exp_add,
    show -x ^ 2 / 2 + l * x ^ 2 = -(1 / 2 - l) * x ^ 2 from by ring]

/-- **A-1.**  Squared-standard-Gaussian MGF integral:
`∫ exp(l x²) dN(0,1) = (1 − 2l)^(−1/2)` for `l < 1/2`, i.e. `(√(1−2l))⁻¹`.
This is the closed form of the chi-squared-1 moment generating function. -/
theorem mgf_exp_sq_stdGaussian {l : ℝ} (hl : l < 1 / 2) :
    ∫ x, Real.exp (l * x ^ 2) ∂(gaussianReal 0 1) = (Real.sqrt (1 - 2 * l))⁻¹ := by
  have hb : (0 : ℝ) < 1 / 2 - l := by linarith
  rw [integral_gaussianReal_eq_integral_smul (by norm_num : (1 : ℝ≥0) ≠ 0)]
  have hpt : ∀ x : ℝ, gaussianPDFReal 0 1 x • Real.exp (l * x ^ 2)
      = (Real.sqrt (2 * π))⁻¹ * Real.exp (-(1 / 2 - l) * x ^ 2) := by
    intro x; rw [smul_eq_mul]; exact gaussianPDF_mul_exp_sq l x
  simp_rw [hpt]
  rw [integral_const_mul, integral_gaussian (1 / 2 - l)]
  -- (√(2π))⁻¹ · √(π/(1/2 − l)) = (√(1 − 2l))⁻¹
  have hπ : (0 : ℝ) ≤ (2 * π)⁻¹ := by positivity
  have h12 : (1 : ℝ) - 2 * l ≠ 0 := (ne_of_gt (by linarith : (0 : ℝ) < 1 - 2 * l))
  rw [← Real.sqrt_inv, ← Real.sqrt_mul hπ,
    show (2 * π)⁻¹ * (π / (1 / 2 - l)) = (1 - 2 * l)⁻¹ from by
      have hpine : π ≠ 0 := Real.pi_ne_zero
      have hbne : (1 / 2 - l) ≠ 0 := hb.ne'
      field_simp,
    Real.sqrt_inv]

/-- **A-1.**  Non-integrability of `x ↦ exp(l x²)` under `N(0,1)` for `l ≥ 1/2`. Pins the
chi-squared-1 MGF to `0` outside its convergence half-line. -/
theorem not_integrable_exp_sq_stdGaussian {l : ℝ} (hl : 1 / 2 ≤ l) :
    ¬ Integrable (fun x => Real.exp (l * x ^ 2)) (gaussianReal 0 1) := by
  rw [gaussianReal_of_var_ne_zero 0 (by norm_num : (1 : ℝ≥0) ≠ 0),
    integrable_withDensity_iff_integrable_smul' (measurable_gaussianPDF 0 1)
      (ae_of_all _ fun _ ↦ gaussianPDF_lt_top)]
  simp only [toReal_gaussianPDF, smul_eq_mul]
  rw [show (fun x => gaussianPDFReal 0 1 x * Real.exp (l * x ^ 2))
        = fun x => (Real.sqrt (2 * π))⁻¹ * Real.exp (-(1 / 2 - l) * x ^ 2) from
      funext (fun x => gaussianPDF_mul_exp_sq l x)]
  intro hInt
  set c : ℝ := (Real.sqrt (2 * π))⁻¹ with hcdef
  have hcpos : 0 < c := by rw [hcdef]; positivity
  have hle : ∀ x : ℝ, ‖(fun _ : ℝ => c) x‖
      ≤ ‖(fun x : ℝ => c * Real.exp (-(1 / 2 - l) * x ^ 2)) x‖ := by
    intro x
    have hexp1 : (1 : ℝ) ≤ Real.exp (-(1 / 2 - l) * x ^ 2) :=
      Real.one_le_exp (mul_nonneg (by linarith : (0 : ℝ) ≤ -(1 / 2 - l)) (sq_nonneg x))
    have hge : c ≤ c * Real.exp (-(1 / 2 - l) * x ^ 2) := by
      nlinarith [hcpos, hexp1]
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos hcpos,
      abs_of_pos (by positivity : (0 : ℝ) < c * Real.exp (-(1 / 2 - l) * x ^ 2))]
    exact hge
  have hconst : Integrable (fun _ : ℝ => c) volume :=
    hInt.mono aestronglyMeasurable_const (ae_of_all _ hle)
  rw [integrable_const_iff] at hconst
  rcases hconst with hc0 | hfin
  · exact hcpos.ne' hc0
  · haveI := hfin
    exact absurd Real.volume_univ (measure_ne_top volume Set.univ)

/-- **A-1.**  Integrability of `x ↦ exp(l x²)` under the standard Gaussian, for `l < 1/2`. -/
theorem integrable_exp_sq_stdGaussian {l : ℝ} (hl : l < 1 / 2) :
    Integrable (fun x => Real.exp (l * x ^ 2)) (gaussianReal 0 1) := by
  have hb : (0 : ℝ) < 1 / 2 - l := by linarith
  rw [gaussianReal_of_var_ne_zero 0 (by norm_num : (1 : ℝ≥0) ≠ 0),
    integrable_withDensity_iff_integrable_smul' (measurable_gaussianPDF 0 1)
      (ae_of_all _ fun _ ↦ gaussianPDF_lt_top)]
  simp only [toReal_gaussianPDF, smul_eq_mul]
  have hpt : (fun x => gaussianPDFReal 0 1 x * Real.exp (l * x ^ 2))
      = fun x => (Real.sqrt (2 * π))⁻¹ * Real.exp (-(1 / 2 - l) * x ^ 2) := by
    funext x; exact gaussianPDF_mul_exp_sq l x
  rw [hpt]
  exact (integrable_exp_neg_mul_sq hb).const_mul _

end Unit2B.FormalR1.Kappa.ProofsR5

/-! ## A-2 — the Gamma MGF and its integrability halves

Ported from StatLean's `GammaMoments.lean` and from the private lemmas of
`ChiSquared.lean`. The licence notice at the head of this file covers the provenance; the
material change is the `integrableExpSet_eq_of_mgf` re-check, which is consumed in A-3
rather than here.
-/

namespace Unit2B.FormalR1.Kappa.ProofsR5

open Set

/-- The Gamma pdf as an `ℝ≥0∞`-valued function is measurable. -/
private lemma measurable_gammaPDF' (a r : ℝ) : Measurable (gammaPDF a r) :=
  (measurable_gammaPDFReal a r).ennreal_ofReal

/-- Reduce a `gammaMeasure` integral to a Lebesgue integral over `Ioi 0`. -/
private lemma integral_gammaMeasure_eq_setIntegral {a r : ℝ} (ha : 0 < a) (hr : 0 < r)
    (g : ℝ → ℝ) :
    ∫ x, g x ∂(gammaMeasure a r) = ∫ x in Ioi (0 : ℝ), gammaPDFReal a r x * g x := by
  rw [gammaMeasure, integral_withDensity_eq_integral_toReal_smul (measurable_gammaPDF' a r)
        (ae_of_all _ (fun _ => ENNReal.ofReal_lt_top))]
  have h1 : ∀ x : ℝ, (gammaPDF a r x).toReal • g x = gammaPDFReal a r x * g x := fun x => by
    rw [show gammaPDF a r x = ENNReal.ofReal (gammaPDFReal a r x) from rfl,
        ENNReal.toReal_ofReal (gammaPDFReal_nonneg ha hr x), smul_eq_mul]
  simp_rw [h1]
  refine (setIntegral_eq_integral_of_ae_compl_eq_zero ?_).symm
  rw [MeasureTheory.ae_iff]
  refine measure_mono_null (t := {(0 : ℝ)}) ?_ Real.volume_singleton
  intro x hx
  simp only [mem_setOf_eq, Classical.not_imp, mem_Ioi, not_lt] at hx
  obtain ⟨hx1, hx2⟩ := hx
  rcases lt_or_eq_of_le hx1 with h | h
  · exact absurd (by rw [gammaPDFReal, if_neg (not_le.mpr h), zero_mul]) hx2
  · exact h

/-- **A-2.**  MGF of `Gamma(a, r)`: `∫ exp(t x) dGamma(a,r) = (r/(r-t))^a` for `t < r`. -/
theorem mgf_gammaMeasure' {a r : ℝ} (ha : 0 < a) (hr : 0 < r) {t : ℝ} (ht : t < r) :
    ∫ x, Real.exp (t * x) ∂(gammaMeasure a r) = (r / (r - t)) ^ a := by
  have hEqOn : EqOn (fun x => gammaPDFReal a r x * Real.exp (t * x))
      (fun x => r ^ a / Real.Gamma a * (x ^ (a - 1) * Real.exp (-((r - t) * x))))
      (Ioi (0 : ℝ)) := by
    intro x hx
    rw [mem_Ioi] at hx
    simp only
    rw [gammaPDFReal, if_pos hx.le, mul_assoc, ← Real.exp_add,
        show -(r * x) + t * x = -((r - t) * x) by ring]
    ring
  rw [integral_gammaMeasure_eq_setIntegral ha hr, setIntegral_congr_fun measurableSet_Ioi hEqOn,
      integral_const_mul,
      Real.integral_rpow_mul_exp_neg_mul_Ioi ha (by linarith : (0 : ℝ) < r - t)]
  rw [show r ^ a / Real.Gamma a * ((1 / (r - t)) ^ a * Real.Gamma a)
        = r ^ a * (1 / (r - t)) ^ a from by
      rw [mul_comm ((1 / (r - t)) ^ a) (Real.Gamma a), ← mul_assoc,
          div_mul_cancel₀ _ (Real.Gamma_pos_of_pos ha).ne'],
    ← Real.mul_rpow hr.le (div_nonneg zero_le_one (by linarith : (0 : ℝ) ≤ r - t)), mul_one_div]

/-- **A-2.**  Integrability of `x ↦ exp(l·x)` under `Gamma(a,r)` for `l < r`. -/
private lemma integrable_exp_mul_gammaMeasure' {a r : ℝ} (ha : 0 < a) (hr : 0 < r)
    {l : ℝ} (hl : l < r) :
    Integrable (fun x => Real.exp (l * x)) (gammaMeasure a r) := by
  have hmeasG : Measurable (gammaPDF a r) := measurable_gammaPDF' a r
  rw [gammaMeasure, integrable_withDensity_iff hmeasG
        (ae_of_all _ (fun _ => ENNReal.ofReal_lt_top))]
  rw [show (fun x => Real.exp (l * x) * (gammaPDF a r x).toReal)
        = (fun x => Real.exp (l * x) * gammaPDFReal a r x) from by
      funext x
      rw [show gammaPDF a r x = ENNReal.ofReal (gammaPDFReal a r x) from rfl,
        ENNReal.toReal_ofReal (gammaPDFReal_nonneg ha hr x)]]
  have hmodel : IntegrableOn (fun x => x ^ (a - 1) * Real.exp (-((r - l) * x)))
      (Ioi (0 : ℝ)) volume := by
    have h := integrableOn_rpow_mul_exp_neg_mul_rpow (p := 1) (s := a - 1) (b := r - l)
      (by linarith) (by norm_num : (0 : ℝ) < 1) (by linarith)
    refine h.congr_fun (fun x hx => ?_) measurableSet_Ioi
    dsimp only
    rw [Real.rpow_one, neg_mul]
  have hIoi : IntegrableOn (fun x => Real.exp (l * x) * gammaPDFReal a r x) (Ioi 0) volume := by
    refine IntegrableOn.congr_fun (hmodel.const_mul (r ^ a / Real.Gamma a))
      (fun x hx => ?_) measurableSet_Ioi
    rw [mem_Ioi] at hx
    rw [gammaPDFReal, if_pos hx.le,
      show -((r - l) * x) = l * x + -(r * x) from by ring, Real.exp_add]
    ring
  rw [← integrableOn_univ, ← Iic_union_Ioi (a := (0 : ℝ)), integrableOn_union]
  refine ⟨?_, hIoi⟩
  refine integrableOn_zero.congr ?_
  rw [Filter.EventuallyEq, ae_restrict_iff' measurableSet_Iic, MeasureTheory.ae_iff]
  refine measure_mono_null (t := {(0 : ℝ)}) ?_ Real.volume_singleton
  intro x hx
  simp only [mem_setOf_eq, Classical.not_imp, mem_Iic] at hx
  obtain ⟨hx1, hx2⟩ := hx
  rcases lt_or_eq_of_le hx1 with h | h
  · exact absurd (show Real.exp (l * x) * gammaPDFReal a r x = 0 by
      rw [gammaPDFReal, if_neg (not_le.mpr h), mul_zero]).symm hx2
  · exact h

/-- **A-2.**  Non-integrability of `x ↦ exp(l·x)` under `Gamma(a,r)` for `r ≤ l`. -/
private lemma not_integrable_exp_mul_gammaMeasure' {a r : ℝ} (ha : 0 < a) (hr : 0 < r)
    {l : ℝ} (hl : r ≤ l) :
    ¬ Integrable (fun x => Real.exp (l * x)) (gammaMeasure a r) := by
  have hmeasG : Measurable (gammaPDF a r) := measurable_gammaPDF' a r
  rw [gammaMeasure, integrable_withDensity_iff hmeasG
        (ae_of_all _ (fun _ => ENNReal.ofReal_lt_top))]
  rw [show (fun x => Real.exp (l * x) * (gammaPDF a r x).toReal)
        = (fun x => Real.exp (l * x) * gammaPDFReal a r x) from by
      funext x
      rw [show gammaPDF a r x = ENNReal.ofReal (gammaPDFReal a r x) from rfl,
        ENNReal.toReal_ofReal (gammaPDFReal_nonneg ha hr x)]]
  intro hInt
  set C : ℝ := r ^ a / Real.Gamma a with hC
  have hCpos : 0 < C := div_pos (Real.rpow_pos_of_pos hr a) (Real.Gamma_pos_of_pos ha)
  have hmodel_ni : ¬ IntegrableOn (fun x => x ^ (a - 1)) (Ioi (1 : ℝ)) volume := by
    rw [integrableOn_Ioi_rpow_iff (by norm_num : (0 : ℝ) < 1)]
    intro hlt; linarith
  apply hmodel_ni
  have hI1 : IntegrableOn (fun x => Real.exp (l * x) * gammaPDFReal a r x) (Ioi (1 : ℝ)) volume :=
    hInt.integrableOn
  have hcmp : IntegrableOn (fun x => C * x ^ (a - 1)) (Ioi (1 : ℝ)) volume := by
    refine Integrable.mono hI1 (by fun_prop) ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with x hx
    rw [mem_Ioi] at hx
    have hx0 : (0 : ℝ) < x := lt_trans one_pos hx
    have hxa : (0 : ℝ) ≤ x ^ (a - 1) := Real.rpow_nonneg hx0.le _
    have hexpge : (1 : ℝ) ≤ Real.exp ((l - r) * x) :=
      Real.one_le_exp (mul_nonneg (by linarith) hx0.le)
    rw [gammaPDFReal, if_pos hx0.le, Real.norm_eq_abs, Real.norm_eq_abs,
      abs_of_nonneg (by positivity : (0 : ℝ) ≤ C * x ^ (a - 1)),
      abs_of_nonneg (by positivity :
        (0 : ℝ) ≤ Real.exp (l * x) * (C * x ^ (a - 1) * Real.exp (-(r * x)))),
      show Real.exp (l * x) * (C * x ^ (a - 1) * Real.exp (-(r * x)))
          = C * x ^ (a - 1) * Real.exp ((l - r) * x) from by
        rw [show (l - r) * x = l * x + -(r * x) from by ring, Real.exp_add]; ring]
    nlinarith [hexpge, mul_nonneg hCpos.le hxa]
  have hscaled := hcmp.const_mul (C⁻¹)
  refine IntegrableOn.congr_fun hscaled (fun x _ => ?_) measurableSet_Ioi
  rw [← mul_assoc, inv_mul_cancel₀ hCpos.ne', one_mul]

end Unit2B.FormalR1.Kappa.ProofsR5

/-! ## A-3 — the sum-of-squares law, by MGF matching and characteristic-function uniqueness

Ported from StatLean's `ChiSquared.lean`. The order's T0 step 4 requires that StatLean's
`chiSquared` not enter a claim-bearing statement, so `gammaMeasure (k/2) (1/2)` is written
throughout instead. The order's T0 step 3, re-checking `integrableExpSet_eq_of_mgf` against
the vendored Mathlib, is discharged by the presence of that lemma at
`Probability/Moments/ComplexMGF.lean:241` in primed form, and by the build below.
-/

namespace Unit2B.FormalR1.Kappa.ProofsR5

open Set

/-- The chi-squared law at `k` degrees of freedom, written as the gamma law so that no
external definition enters a claim-bearing statement. -/
private noncomputable def chiSqLaw (k : ℕ) : Measure ℝ := gammaMeasure ((k : ℝ) / 2) (1 / 2)

/-- `chiSqLaw k` is a probability measure for `k ≥ 1`, which also supplies the
`IsFiniteMeasure` instance that `Measure.ext_of_charFun` needs. -/
instance isProbabilityMeasure_chiSq (k : ℕ) [NeZero k] :
    IsProbabilityMeasure (chiSqLaw k) := by
  have hk : (0 : ℝ) < k := by exact_mod_cast Nat.pos_of_ne_zero (NeZero.ne k)
  unfold chiSqLaw
  exact isProbabilityMeasure_gammaMeasure (by linarith) (by norm_num)

/-- `rpow` bridge: `((√(1−2l))⁻¹)ⁿ = ((1/2)/((1/2)−l))^(n/2)`, both `(1−2l)^(−n/2)`. -/
private lemma chiSq_mgf_rpow_bridge (n : ℕ) {l : ℝ} (hl : l < 1 / 2) :
    (Real.sqrt (1 - 2 * l))⁻¹ ^ n = ((1 / 2) / ((1 / 2) - l)) ^ ((n : ℝ) / 2) := by
  have hb : (0 : ℝ) < 1 - 2 * l := by linarith
  have hRHSbase : (1 / 2 : ℝ) / ((1 / 2) - l) = (1 - 2 * l)⁻¹ := by
    rw [show (1 / 2 : ℝ) - l = (1 - 2 * l) / 2 from by ring]; field_simp
  have hL : (Real.sqrt (1 - 2 * l))⁻¹ ^ n = (1 - 2 * l) ^ (-(n : ℝ) / 2) := by
    rw [Real.sqrt_eq_rpow, ← Real.rpow_neg hb.le, ← Real.rpow_natCast _ n, ← Real.rpow_mul hb.le]
    congr 1; ring
  have hR : ((1 / 2 : ℝ) / ((1 / 2) - l)) ^ ((n : ℝ) / 2) = (1 - 2 * l) ^ (-(n : ℝ) / 2) := by
    rw [hRHSbase, show (1 - 2 * l)⁻¹ = (1 - 2 * l) ^ (-1 : ℝ) from by
        rw [Real.rpow_neg hb.le, Real.rpow_one], ← Real.rpow_mul hb.le]
    congr 1; ring
  rw [hL, hR]

/-- `x ↦ exp(l·x)` is integrable under `χ²ₙ` exactly on the MGF half-line `l < 1/2`. -/
private lemma integrable_exp_mul_chiSq {n : ℕ} (hn : 0 < n) {l : ℝ} (hl : l < 1 / 2) :
    Integrable (fun x => Real.exp (l * x)) (chiSqLaw n) := by
  have hnr : (0 : ℝ) < n := by exact_mod_cast hn
  unfold chiSqLaw
  exact integrable_exp_mul_gammaMeasure' (by linarith) (by norm_num) hl

/-- `x ↦ exp(l·x)` is non-integrable under `χ²ₙ` for `1/2 ≤ l`. -/
private lemma not_integrable_exp_mul_chiSq {n : ℕ} (hn : 0 < n) {l : ℝ} (hl : 1 / 2 ≤ l) :
    ¬ Integrable (fun x => Real.exp (l * x)) (chiSqLaw n) := by
  have hnr : (0 : ℝ) < n := by exact_mod_cast hn
  unfold chiSqLaw
  exact not_integrable_exp_mul_gammaMeasure' (by linarith) (by norm_num) hl

/-- MGF of `χ²ₖ` for `t < 1/2`. -/
theorem mgf_chiSq {k : ℕ} (hk : 0 < k) {t : ℝ} (ht : t < 1 / 2) :
    ∫ x, Real.exp (t * x) ∂(chiSqLaw k) = ((1 / 2) / ((1 / 2) - t)) ^ ((k : ℝ) / 2) := by
  have hkr : (0 : ℝ) < k := by exact_mod_cast hk
  unfold chiSqLaw
  exact mgf_gammaMeasure' (by linarith) (by norm_num) ht

/-- **A-3.**  Sum-of-squares law: for i.i.d. `Zᵢ ∼ N(0,1)` under an arbitrary probability
measure, `∑ᵢ Zᵢ²` has the chi-squared law with `n` degrees of freedom.  Proved by matching
moment generating functions and then appealing to characteristic-function uniqueness. -/
theorem map_sum_sq_eq_chiSq {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω}
    [IsProbabilityMeasure μ] {n : ℕ} (hn : 0 < n) (Z : Fin n → Ω → ℝ)
    (hmeas : ∀ i, Measurable (Z i))
    (hlaw : ∀ i, Measure.map (Z i) μ = gaussianReal 0 1)
    (hindep : iIndepFun Z μ) :
    Measure.map (fun ω => ∑ i, (Z i ω) ^ 2) μ = chiSqLaw n := by
  classical
  haveI : NeZero n := ⟨hn.ne'⟩
  set S : Ω → ℝ := fun ω => ∑ i, (Z i ω) ^ 2 with hS
  set W : Fin n → Ω → ℝ := fun i ω => (Z i ω) ^ 2 with hW
  have hWmeas : ∀ i, Measurable (W i) := fun i => (hmeas i).pow_const 2
  have hWindep : iIndepFun W μ := by
    have := hindep.comp (fun _ => fun x : ℝ => x ^ 2) (fun _ => by fun_prop)
    simpa [Function.comp_def, hW] using this
  have hSsum : S = ∑ i, W i := by
    funext ω; simp only [hS, hW, Finset.sum_apply]
  have hSmeas : Measurable S := by
    rw [hS]; exact Finset.measurable_sum _ (fun i _ => (hmeas i).pow_const 2)
  have hWval : ∀ (i : Fin n) (l : ℝ),
      mgf (W i) μ l = ∫ y, Real.exp (l * y ^ 2) ∂(gaussianReal 0 1) := by
    intro i l
    have hg_asm : AEStronglyMeasurable (fun y => Real.exp (l * y ^ 2))
        (Measure.map (Z i) μ) :=
      (Real.continuous_exp.comp (continuous_const.mul (continuous_pow 2))).aestronglyMeasurable
    simp only [mgf]
    rw [show (fun ω => Real.exp (l * W i ω))
          = (fun ω => (fun y => Real.exp (l * y ^ 2)) (Z i ω)) from by funext ω; simp [hW]]
    rw [← integral_map (hmeas i).aemeasurable hg_asm, hlaw i]
  have hidval : ∀ l : ℝ, mgf id (chiSqLaw n) l = ∫ x, Real.exp (l * x) ∂(chiSqLaw n) := by
    intro l; simp only [mgf, id_eq]
  have hA : mgf S μ = mgf id (chiSqLaw n) := by
    funext l
    rw [show mgf S μ l = ∏ i, mgf (W i) μ l from by
      rw [hSsum]; exact hWindep.mgf_sum hWmeas Finset.univ]
    rcases lt_or_ge l (1 / 2) with hl | hl
    · have hWi : ∀ i, mgf (W i) μ l = (Real.sqrt (1 - 2 * l))⁻¹ := fun i => by
        rw [hWval i l, mgf_exp_sq_stdGaussian hl]
      rw [Finset.prod_congr rfl (fun i _ => hWi i), Finset.prod_const, Finset.card_univ,
        Fintype.card_fin, hidval l, mgf_chiSq hn hl]
      exact chiSq_mgf_rpow_bridge n hl
    · have hWi0 : ∀ i, mgf (W i) μ l = 0 := fun i => by
        rw [hWval i l, integral_undef (not_integrable_exp_sq_stdGaussian hl)]
      obtain ⟨i0, _⟩ := (Finset.univ_nonempty (α := Fin n))
      rw [Finset.prod_eq_zero (Finset.mem_univ i0) (hWi0 i0), hidval l,
        integral_undef (not_integrable_exp_mul_chiSq hn hl)]
  have hB : (0 : ℝ) ∈ interior (integrableExpSet S μ) := by
    rw [integrableExpSet_eq_of_mgf hA, mem_interior]
    refine ⟨Iio (1 / 2), fun l hl => ?_, isOpen_Iio, mem_Iio.mpr (by norm_num)⟩
    rw [mem_Iio] at hl
    change Integrable (fun x => Real.exp (l * id x)) (chiSqLaw n)
    simpa only [id_eq] using integrable_exp_mul_chiSq hn hl
  haveI : IsProbabilityMeasure (Measure.map S μ) := inferInstance
  refine MeasureTheory.Measure.ext_of_charFun (funext fun t => ?_)
  rw [← complexMGF_mul_I hSmeas.aemeasurable t, ← complexMGF_id_mul_I t]
  refine eqOn_complexMGF_of_mgf hA ?_
  rw [Set.mem_setOf_eq, show ((t : ℂ) * Complex.I).re = 0 from by simp]
  exact hB

end Unit2B.FormalR1.Kappa.ProofsR5

/-! ## A-4 — the marginal of a product of standard Gaussians, then the order's D-3 statement

The marginal is proved directly from `Measure.pi_pi` on a cylinder set, which avoids
characteristic functions entirely: the pushforward of a coordinate is the measure of the
cylinder obtained by pulling a set back along it, and `pi_pi` evaluates that cylinder as the
product of the factors, of which only one is non-trivial.
-/

namespace Unit2B.FormalR1.Kappa.ProofsR5

open Set

private lemma univ_gaussianReal_one : (gaussianReal 0 1) (Set.univ : Set ℝ) = 1 :=
  (inferInstance : IsProbabilityMeasure (gaussianReal 0 1)).measure_univ

/-- **A-4, first half.**  The marginal of a product of i.i.d. standard Gaussians is again a
standard Gaussian. -/
theorem gaussianReal_pi_marginal (k : ℕ) [NeZero k] (i : Fin k) :
    (Measure.pi fun _ : Fin k => gaussianReal 0 1).map (fun x => x i) = gaussianReal 0 1 := by
  have hmeas : Measurable (fun x : Fin k → ℝ => x i) := measurable_pi_apply i
  refine Measure.ext fun A hA => ?_
  have hpull : (fun x : Fin k → ℝ => x i) ⁻¹' A
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
    have h1 : (gaussianReal 0 1) (Set.univ : Set ℝ) = 1 := univ_gaussianReal_one
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

/-- **A-4, second half.**  The order's **D-3 statement, verbatim**: for every positive
integer `k`, the sum of the squares of `k` independent standard normal random variables has
law `gammaMeasure (k/2) (1/2)`.  Obtained from A-3 by instantiating the product measure, its
independent standard-normal coordinates, and the marginal just proved. -/
theorem chiSq_eq_gamma (k : ℕ) (hk : 0 < k) :
    (Measure.pi fun _ : Fin k => gaussianReal 0 1).map (fun x => ∑ i, x i ^ 2)
      = gammaMeasure ((k : ℝ) / 2) (1 / 2) := by
  haveI : NeZero k := ⟨hk.ne'⟩
  have h := map_sum_sq_eq_chiSq hk (μ := Measure.pi fun _ : Fin k => gaussianReal 0 1)
    (Z := fun i (x : Fin k → ℝ) => x i) (fun i => measurable_pi_apply i)
    (fun i => gaussianReal_pi_marginal k i)
    (iIndepFun_pi (μ := fun _ : Fin k => gaussianReal 0 1)
      (X := fun _ => (id : ℝ → ℝ)) (fun _ => measurable_id.aemeasurable))
  unfold chiSqLaw at h
  exact h

end Unit2B.FormalR1.Kappa.ProofsR5
