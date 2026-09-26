/-
Copyright (c) 2026 Momen Ahmed. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Momen Ahmed
-/
import ExactManyToOne.Reduction

/-!
# Sections 9.8, 9.11 and 10.4: the corner-substitution enclosure and its certificates

Every certificate in the paper rests on one monotonicity fact. The factor
`Φ(s · x)` is increasing in the scale `s` where its argument `x` is positive and
decreasing where it is negative, so on a box of scales the minimum and the
maximum sit at the corner that matches the sign (equation (37)). This file proves
that step, the product enclosure built from it, the integrated lower bound, and the
logical steps that turn a lower bound on the distribution function into an upper
bound on the constant (the direction Section 9.11 warns about). It then proves the
branch-and-bound tiling argument of equation (38), the partition bound of
equation (41) with its tail term, and the three substitutions of Section 10.4.

`Φ` is any monotone function. Nothing here uses the normal distribution beyond
monotonicity and nonnegativity.
-/

open MeasureTheory Set
open scoped ENNReal

namespace ExactManyToOne

/-- The scale that minimises `Φ (s * x)` over `s ∈ [sm, sp]`. -/
noncomputable def lowCorner (sm sp x : ℝ) : ℝ := if 0 ≤ x then sm else sp

/-- The scale that maximises `Φ (s * x)` over `s ∈ [sm, sp]`. -/
noncomputable def highCorner (sm sp x : ℝ) : ℝ := if 0 ≤ x then sp else sm

section corner

variable {Φ : ℝ → ℝ} (hΦ : Monotone Φ)
include hΦ

/-- **Corner substitution, lower side (equation (37)).** -/
lemma corner_le {s sm sp : ℝ} (hs : s ∈ Icc sm sp) (x : ℝ) :
    Φ (lowCorner sm sp x * x) ≤ Φ (s * x) := by
  apply hΦ
  unfold lowCorner
  split_ifs with hx
  · exact mul_le_mul_of_nonneg_right hs.1 hx
  · exact mul_le_mul_of_nonpos_right hs.2 (le_of_lt (not_le.1 hx))

/-- **Corner substitution, upper side (equation (37)).** -/
lemma le_corner {s sm sp : ℝ} (hs : s ∈ Icc sm sp) (x : ℝ) :
    Φ (s * x) ≤ Φ (highCorner sm sp x * x) := by
  apply hΦ
  unfold highCorner
  split_ifs with hx
  · exact mul_le_mul_of_nonneg_right hs.2 hx
  · exact mul_le_mul_of_nonpos_right hs.1 (le_of_lt (not_le.1 hx))

/-- **The product enclosure `G⁻ ≤ ∏ Φ(sᵢ xᵢ) ≤ G⁺` (equation (37)).** It holds at every
scale vector in the box and every value of the arguments. -/
theorem prod_enclosure {ι : Type*} [Fintype ι] (hΦ0 : ∀ x, 0 ≤ Φ x)
    (s sm sp x : ι → ℝ) (hs : ∀ i, s i ∈ Icc (sm i) (sp i)) :
    ∏ i, Φ (lowCorner (sm i) (sp i) (x i) * x i) ≤ ∏ i, Φ (s i * x i) ∧
      ∏ i, Φ (s i * x i) ≤ ∏ i, Φ (highCorner (sm i) (sp i) (x i) * x i) :=
  ⟨Finset.prod_le_prod₀ (fun _ _ => hΦ0 _) (fun i _ => corner_le hΦ (hs i) (x i)),
   Finset.prod_le_prod₀ (fun _ _ => hΦ0 _) (fun i _ => le_corner hΦ (hs i) (x i))⟩

/-- **The integrated lower envelope bounds the distribution function from below at
every ratio vector in the box.** `F(c, s) = ∫ ∏ Φ(sᵢ xᵢ(z, c)) dμ(z)` is the
conditional distribution function of Section 9.11 with `μ` the law of the
benchmark's standardised mean. -/
theorem lintegral_lower_envelope {ι : Type*} [Fintype ι] (hΦ0 : ∀ x, 0 ≤ Φ x)
    (μ : Measure ℝ) (s sm sp : ι → ℝ) (hs : ∀ i, s i ∈ Icc (sm i) (sp i))
    (x : ι → ℝ → ℝ) :
    ∫⁻ z, ENNReal.ofReal (∏ i, Φ (lowCorner (sm i) (sp i) (x i z) * x i z)) ∂μ
      ≤ ∫⁻ z, ENNReal.ofReal (∏ i, Φ (s i * x i z)) ∂μ :=
  lintegral_mono fun z =>
    ENNReal.ofReal_le_ofReal (prod_enclosure hΦ hΦ0 s sm sp (fun i => x i z) hs).1

end corner

/-- **From a box of ratios to a box of scales.** The factor in (37) carries
`sᵢ = √(nᵢ/λᵢ)`, which is decreasing in the ratio, so `λ ∈ [λ⁻, λ⁺]` puts
`s ∈ [√(n/λ⁺), √(n/λ⁻)]`. -/
lemma scale_mem_Icc {n lam lm lp : ℝ} (hn : 0 ≤ n) (hlm : 0 < lm) (hlam : lam ∈ Icc lm lp) :
    Real.sqrt (n / lam) ∈ Icc (Real.sqrt (n / lp)) (Real.sqrt (n / lm)) := by
  have hl : 0 < lam := lt_of_lt_of_le hlm hlam.1
  constructor
  · exact Real.sqrt_le_sqrt (div_le_div_of_nonneg_left hn hl hlam.2)
  · exact Real.sqrt_le_sqrt (div_le_div_of_nonneg_left hn hlm hlam.1)

/-- **A lower bound on acceptance is an upper bound on rejection.** If the
acceptance probability at `c` is at least `1 - α`, rejection at `c` has
probability at most `α`. This is the step that makes a lower bound on the
minimum of the distribution function the quantity a guarantee needs. -/
theorem reject_le_of_accept_ge {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    [IsProbabilityMeasure P] {A : Set Ω} (hA : MeasurableSet A) {α : ℝ≥0∞} (hα : α ≤ 1)
    (h : 1 - α ≤ P A) : P Aᶜ ≤ α := by
  rw [prob_compl_eq_one_sub hA]
  calc 1 - P A ≤ 1 - (1 - α) := tsub_le_tsub_left h 1
    _ = α := ENNReal.sub_sub_cancel ENNReal.one_ne_top hα

/-- **The certified constant dominates every constant in the box (equation (35)).**
Let `F c λ` be the distribution function at ratio `λ` and `cpop λ` the smallest
constant at which it reaches `1 - α`. If a lower bound `Flo` holds at every ratio in
the box and reaches `1 - α` at `crig`, then `crig` is at least `cpop λ` for every
`λ` in the box. An upper bound on the maximum of `F` would bound the smallest
constant instead, which is the wrong direction. -/
theorem cpop_le_crig {Λ : Type*} (F : ℝ → Λ → ℝ) (Flo : ℝ → ℝ) (box : Set Λ) (α crig : ℝ)
    (hlo : ∀ lam ∈ box, Flo crig ≤ F crig lam) (hcrig : 1 - α ≤ Flo crig)
    (hbdd : ∀ lam, BddBelow {c | 1 - α ≤ F c lam}) :
    ∀ lam ∈ box, sInf {c | 1 - α ≤ F c lam} ≤ crig :=
  fun lam hlam => csInf_le (hbdd lam) (hcrig.trans (hlo lam hlam))

/-- **Branch and bound (equation (38)).** The terminal boxes cover the region and
each box's bound dominates the constant on that box, so the maximum over all
terminal boxes, pruned ones included, dominates the constant everywhere in the
region, and a witness in the region is dominated by it. The certified gap
`c_U - c_L` therefore brackets the supremum. -/
theorem tiling_bound {α ι : Type*} (region : Set α) (boxes : Finset ι) (hne : boxes.Nonempty)
    (box : ι → Set α) (hcover : region ⊆ ⋃ j ∈ boxes, box j) (ub : ι → ℝ) (f : α → ℝ)
    (hub : ∀ j ∈ boxes, ∀ x ∈ box j, f x ≤ ub j) :
    ∀ x ∈ region, f x ≤ boxes.sup' hne ub := by
  intro x hx
  obtain ⟨j, hj, hxj⟩ := mem_iUnion₂.1 (hcover hx)
  exact (hub j hj x hxj).trans (Finset.le_sup' ub hj)

/-- The certificate brackets the supremum: with `c_U` dominating the constant on the
region and a witness in the region attaining `c_L`, every value the supremum could
take lies in `[c_L, c_U]`, so a gap under the tolerance pins it to that tolerance. -/
theorem certificate_brackets {α : Type*} (region : Set α) (f : α → ℝ) (cU cL tol : ℝ)
    (hU : ∀ x ∈ region, f x ≤ cU) (w : α) (hw : w ∈ region) (hL : f w = cL)
    (hgap : cU - cL ≤ tol) :
    sSup (f '' region) ∈ Icc cL cU ∧ cU - sSup (f '' region) ≤ tol := by
  have hbdd : BddAbove (f '' region) := ⟨cU, by rintro _ ⟨x, hx, rfl⟩; exact hU x hx⟩
  have hne : (f '' region).Nonempty := ⟨f w, w, hw, rfl⟩
  have h1 : cL ≤ sSup (f '' region) := hL ▸ le_csSup hbdd ⟨w, hw, rfl⟩
  have h2 : sSup (f '' region) ≤ cU := csSup_le hne (by rintro _ ⟨x, hx, rfl⟩; exact hU x hx)
  exact ⟨⟨h1, h2⟩, by linarith⟩

/-- **The partition bound of equation (41).** Split the state space into finitely
many cells and a tail. On each cell the integrand is at most that cell's corner
value, and on the tail it is at most one, being a probability. Then the integral is
at most the mass-weighted sum of the corner values plus the tail's mass. -/
theorem partition_bound {α ι : Type*} [MeasurableSpace α] (μ : Measure α)
    (cells : Finset ι) (A : ι → Set α) (tail : Set α)
    (hA : ∀ q ∈ cells, MeasurableSet (A q)) (htail : MeasurableSet tail)
    (hcover : (univ : Set α) ⊆ (⋃ q ∈ cells, A q) ∪ tail)
    (f : α → ℝ≥0∞) (R : ι → ℝ≥0∞) (hR : ∀ q ∈ cells, ∀ x ∈ A q, f x ≤ R q)
    (hone : ∀ x ∈ tail, f x ≤ 1) :
    ∫⁻ x, f x ∂μ ≤ ∑ q ∈ cells, R q * μ (A q) + μ tail := by
  classical
  -- pointwise, `f` is dominated by the step function that takes each cell's corner
  -- value on that cell and one on the tail
  have hpt : ∀ x, f x ≤ (∑ q ∈ cells, (A q).indicator (fun _ => R q) x)
      + tail.indicator (fun _ => (1 : ℝ≥0∞)) x := by
    intro x
    rcases hcover (mem_univ x) with hx | hx
    · obtain ⟨q, hq, hxq⟩ := mem_iUnion₂.1 hx
      calc f x ≤ R q := hR q hq x hxq
        _ = (A q).indicator (fun _ => R q) x := (indicator_of_mem hxq (fun _ => R q)).symm
        _ ≤ ∑ q ∈ cells, (A q).indicator (fun _ => R q) x :=
            Finset.single_le_sum (f := fun q => (A q).indicator (fun _ => R q) x)
              (fun _ _ => bot_le) hq
        _ ≤ _ := le_self_add
    · calc f x ≤ 1 := hone x hx
        _ = tail.indicator (fun _ => (1 : ℝ≥0∞)) x :=
            (indicator_of_mem hx (fun _ => (1 : ℝ≥0∞))).symm
        _ ≤ _ := le_add_self
  calc ∫⁻ x, f x ∂μ
      ≤ ∫⁻ x, (∑ q ∈ cells, (A q).indicator (fun _ => R q) x)
          + tail.indicator (fun _ => (1 : ℝ≥0∞)) x ∂μ := lintegral_mono hpt
    _ = ∑ q ∈ cells, ∫⁻ x, (A q).indicator (fun _ => R q) x ∂μ
          + ∫⁻ x, tail.indicator (fun _ => (1 : ℝ≥0∞)) x ∂μ := by
        rw [lintegral_add_right _ (measurable_const.indicator htail),
          lintegral_finsetSum _ fun q hq => measurable_const.indicator (hA q hq)]
    _ = ∑ q ∈ cells, R q * μ (A q) + μ tail := by
        rw [lintegral_indicator_const htail, one_mul]
        congr 1
        exact Finset.sum_congr rfl fun q hq => lintegral_indicator_const (hA q hq) _

/-- **The tail term `ε_tail = 4 p_min`.** On the unit square with the uniform law
(the probability-scale coordinates of the two chi-square factors), the part outside
`(p, 1 - p)²` carries at most `4p` of mass. -/
theorem tail_mass_le (p : ℝ) (hp : 0 ≤ p) :
    (volume : Measure (ℝ × ℝ)) ((Icc (0 : ℝ) 1 ×ˢ Icc (0 : ℝ) 1) \ (Ioo p (1 - p) ×ˢ Ioo p (1 - p)))
      ≤ ENNReal.ofReal (4 * p) := by
  -- the escape set lies in two strips per coordinate, each of width `p`
  have hsub : (Icc (0 : ℝ) 1 ×ˢ Icc (0 : ℝ) 1) \ (Ioo p (1 - p) ×ˢ Ioo p (1 - p))
      ⊆ ((Icc 0 p ∪ Icc (1 - p) 1) ×ˢ Icc (0 : ℝ) 1)
        ∪ (Icc (0 : ℝ) 1 ×ˢ (Icc 0 p ∪ Icc (1 - p) 1)) := by
    rintro ⟨u, v⟩ ⟨⟨hu, hv⟩, hnot⟩
    simp only [mem_prod, mem_Ioo, not_and_or, not_lt] at hnot
    simp only [mem_union, mem_prod, mem_Icc] at hu hv ⊢
    rcases hnot with (h | h) | (h | h)
    · exact Or.inl ⟨Or.inl ⟨hu.1, h⟩, hv⟩
    · exact Or.inl ⟨Or.inr ⟨h, hu.2⟩, hv⟩
    · exact Or.inr ⟨hu, Or.inl ⟨hv.1, h⟩⟩
    · exact Or.inr ⟨hu, Or.inr ⟨h, hv.2⟩⟩
  have hstrip : (volume : Measure ℝ) (Icc 0 p ∪ Icc (1 - p) 1) ≤ ENNReal.ofReal (2 * p) := by
    calc (volume : Measure ℝ) (Icc 0 p ∪ Icc (1 - p) 1)
        ≤ volume (Icc 0 p) + volume (Icc (1 - p) 1) := measure_union_le _ _
      _ = ENNReal.ofReal p + ENNReal.ofReal p := by
          rw [Real.volume_Icc, Real.volume_Icc]; congr 1 <;> ring_nf
      _ = ENNReal.ofReal (2 * p) := by
          rw [← ENNReal.ofReal_add hp hp]; ring_nf
  have hunit : (volume : Measure ℝ) (Icc 0 1) = 1 := by simp [Real.volume_Icc]
  calc (volume : Measure (ℝ × ℝ)) _
      ≤ volume (((Icc 0 p ∪ Icc (1 - p) 1) ×ˢ Icc (0 : ℝ) 1)
          ∪ (Icc (0 : ℝ) 1 ×ˢ (Icc 0 p ∪ Icc (1 - p) 1))) := measure_mono hsub
    _ ≤ volume ((Icc 0 p ∪ Icc (1 - p) 1) ×ˢ Icc (0 : ℝ) 1)
          + volume (Icc (0 : ℝ) 1 ×ˢ (Icc 0 p ∪ Icc (1 - p) 1)) := measure_union_le _ _
    _ = volume (Icc 0 p ∪ Icc (1 - p) 1) + volume (Icc 0 p ∪ Icc (1 - p) 1) := by
          rw [Measure.volume_eq_prod, Measure.prod_prod, Measure.prod_prod, hunit, mul_one, one_mul]
    _ ≤ ENNReal.ofReal (2 * p) + ENNReal.ofReal (2 * p) := add_le_add hstrip hstrip
    _ = ENNReal.ofReal (4 * p) := by
          rw [← ENNReal.ofReal_add (by linarith) (by linarith)]; ring_nf

/-- **The three substitutions of Section 10.4.** At one candidate and a fixed
chi-square state, the conditional rejection probability is `1 - Φ(c · d(λ)/σ(λ))`,
with `d` the estimated standard error and `σ` the true one, both increasing in the
ratio. On a box `[λ⁻, λ⁺]`, evaluating `d` at `λ⁻`, `σ` at `λ⁺` and the threshold at
the smallest value it can take anywhere in the box each moves the rejection
probability up, so the substituted value bounds it at every ratio in the box. -/
theorem three_substitutions {Φ : ℝ → ℝ} (hΦ : Monotone Φ) (d σ : ℝ → ℝ)
    (hd : Monotone d) (hσ : Monotone σ) (hd0 : ∀ l, 0 ≤ d l) (hσ0 : ∀ l, 0 < σ l)
    (thr : ℝ → ℝ) (cmin lm lp : ℝ) (hcmin : 0 ≤ cmin) (hthr : ∀ l ∈ Icc lm lp, cmin ≤ thr l)
    (l : ℝ) (hl : l ∈ Icc lm lp) :
    1 - Φ (thr l * d l / σ l) ≤ 1 - Φ (cmin * d lm / σ lp) := by
  have hdl : d lm ≤ d l := hd hl.1
  have hσl : σ l ≤ σ lp := hσ hl.2
  have h1 : cmin * d lm ≤ thr l * d l :=
    mul_le_mul (hthr l hl) hdl (hd0 _) (hcmin.trans (hthr l hl))
  have h2 : cmin * d lm / σ lp ≤ cmin * d lm / σ l :=
    div_le_div_of_nonneg_left (mul_nonneg hcmin (hd0 _)) (hσ0 l) hσl
  have h3 : cmin * d lm / σ l ≤ thr l * d l / σ l := div_le_div_of_nonneg_right h1 (hσ0 l).le
  linarith [hΦ (h2.trans h3)]

/-- **Decoupling only enlarges the bound.** Letting the ratio take independent values
in the places it appears gives a supremum over a larger set, so the decoupled bound
holds at every ratio on the diagonal. -/
theorem decoupling_le {α β : Type*} (f : α → α → α → β) [Preorder β] (B : β) (box : Set α)
    (h : ∀ a ∈ box, ∀ b ∈ box, ∀ c ∈ box, f a b c ≤ B) : ∀ a ∈ box, f a a a ≤ B :=
  fun a ha => h a ha a ha a ha

end ExactManyToOne
