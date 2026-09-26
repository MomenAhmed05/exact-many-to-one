/-
Copyright (c) 2026 Momen Ahmed. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Momen Ahmed
-/
import ExactManyToOne.Reduction

/-!
# Section 9.2–9.3: the admission region and the envelope guarantee

Two results. The envelope rule's familywise error rate is at most `α` plus the
probability that the true ratios escape the admitted region (equation (32)), and
that escape probability is an exact one-dimensional integral (equation (29)),
which is again the factorisation of Section 3 with the benchmark's chi-square
factor in the role of `B`.
-/

open MeasureTheory ProbabilityTheory Set
open scoped ENNReal

namespace ExactManyToOne

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}

/-- **The envelope inequality.** If, on the event `G` that the region covers the
truth, rejection implies rejection at the true ratios' own constant, then the
error rate is at most that constant's level plus the escape probability. -/
theorem envelope_bound (reject G Apop : Set Ω) (α γ : ℝ≥0∞)
    (hsub : reject ∩ G ⊆ Apop) (hpop : P Apop ≤ α) (hmiss : P Gᶜ ≤ γ) :
    P reject ≤ α + γ := by
  calc P reject ≤ P ((reject ∩ G) ∪ Gᶜ) := by
        refine measure_mono fun ω hω => ?_
        by_cases hG : ω ∈ G
        · exact Or.inl ⟨hω, hG⟩
        · exact Or.inr hG
    _ ≤ P (reject ∩ G) + P Gᶜ := measure_union_le _ _
    _ ≤ P Apop + P Gᶜ := add_le_add (measure_mono hsub) le_rfl
    _ ≤ α + γ := add_le_add hpop hmiss

/-- **Theorem 2 (equation (32)).** The envelope rule rejects when some statistic
exceeds `cenv ω`, the largest known-ratio constant over the admitted region. On
the covering event the true ratios' constant `cpop` is at most `cenv ω`, and
`cpop` has level `α` by definition. So the familywise error rate is at most
`α` plus the escape probability, and at most `α + γ` once the coverage is
certified at `1 - γ`. -/
theorem envelope_guarantee {ι : Type*} (T : ι → Ω → ℝ) (cenv : Ω → ℝ) (cpop : ℝ)
    (G : Set Ω) (α γ : ℝ≥0∞)
    (hdom : ∀ ω ∈ G, cpop ≤ cenv ω)
    (hlevel : P {ω | ∃ i, cpop < T i ω} ≤ α)
    (hcov : P Gᶜ ≤ γ) :
    P {ω | ∃ i, cenv ω < T i ω} ≤ α + γ := by
  refine envelope_bound _ G _ α γ ?_ hlevel hcov
  rintro ω ⟨⟨i, hi⟩, hG⟩
  exact ⟨i, lt_of_le_of_lt (hdom ω hG) hi⟩

/-- **Box membership in the chi-square factors.** With `λ̂ = λ Y / W` (the
estimated ratio is the true one scaled by the candidate's and the benchmark's
chi-square factors), the true ratio lies in the box of equation (27) exactly
when the candidate's factor lies in an interval proportional to the benchmark's.
This is what turns coverage into the integral (29). -/
lemma mem_box_iff {lam W Y wlo whi ylo yhi : ℝ} (hlam : 0 < lam) (hW : 0 < W)
    (hwlo : 0 < wlo) (hwhi : 0 < whi) (hylo : 0 < ylo) (hyhi : 0 < yhi) :
    lam ∈ Icc (lam * Y / W * wlo / yhi) (lam * Y / W * whi / ylo)
      ↔ Y ∈ Icc (W * (ylo / whi)) (W * (yhi / wlo)) := by
  have e1 : lam * Y / W * wlo / yhi = lam * (Y * wlo / (W * yhi)) := by ring
  have e2 : lam * Y / W * whi / ylo = lam * (Y * whi / (W * ylo)) := by ring
  have l1 : lam * (Y * wlo / (W * yhi)) ≤ lam ↔ Y ≤ W * (yhi / wlo) := by
    rw [mul_le_iff_le_one_right hlam, div_le_one (mul_pos hW hyhi), mul_div_assoc',
      le_div_iff₀ hwlo]
  have l2 : lam ≤ lam * (Y * whi / (W * ylo)) ↔ W * (ylo / whi) ≤ Y := by
    rw [le_mul_iff_one_le_right hlam, one_le_div (mul_pos hW hylo), mul_div_assoc',
      div_le_iff₀ hwhi]
  rw [mem_Icc, mem_Icc, e1, e2, l1, l2, and_comm]

/-- **The coverage integral (equation (29)).** The candidates' chi-square factors
`Y i` are mutually independent and independent of the benchmark's `W`. Then the
probability that every candidate factor lies in its `W`-proportional interval is
a single integral over the law of `W` of a product of `K` interval
probabilities. -/
theorem coverage_eq {K : ℕ} (W : Ω → ℝ) (Y : Fin K → Ω → ℝ)
    [IsProbabilityMeasure P]
    (hW : Measurable W) (hY : ∀ i, Measurable (Y i))
    (hWY : IndepFun W (fun ω i => Y i ω) P) (hYY : iIndepFun Y P)
    (alo ahi : Fin K → ℝ) :
    P {ω | ∀ i, Y i ω ∈ Icc (W ω * alo i) (W ω * ahi i)}
      = ∫⁻ w, ∏ i, (P.map (Y i)) (Icc (w * alo i) (w * ahi i)) ∂(P.map W) := by
  have h := prob_forall_eq_lintegral_prod W Y hW hY hWY hYY
    (fun i => {p : ℝ × ℝ | p.1 ∈ Icc (p.2 * alo i) (p.2 * ahi i)})
    (fun i => by
      simp only [mem_Icc, ofPred_and]
      exact (measurableSet_le (by fun_prop) (by fun_prop)).inter
        (measurableSet_le (by fun_prop) (by fun_prop)))
  exact h

/-- Coverage can only fall as the box widens its levels: nested regions have
nested coverage events. Together with the monotonicity of the chi-square
quantiles this is why the tightest feasible scale in (30) is well defined. -/
theorem coverage_mono (G G' : Set Ω) (h : G ⊆ G') : P G ≤ P G' := measure_mono h

end ExactManyToOne
