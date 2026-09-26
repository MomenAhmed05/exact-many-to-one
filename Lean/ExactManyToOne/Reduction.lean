/-
Copyright (c) 2026 Momen Ahmed. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Momen Ahmed
-/
import Mathlib

/-!
# Section 3: the reduction

The paper's first result. Both couplings between the `K` comparisons run through
the benchmark, so conditioning on the benchmark's pair of sufficient statistics
makes the comparisons independent, and the joint probability becomes an outer
integral over that pair of a product of `K` inner probabilities (equation (16)).

`prob_forall_eq_lintegral_prod` is the measure-theoretic content and holds for any
laws at all: it needs independence and nothing else. `exact_law` specialises it
to the statistic of equation (9). The outer integral is over `ℝ × ℝ` whatever `K`
is, which is the dimension result of Section 3.4.
-/

open MeasureTheory ProbabilityTheory Set
open scoped ENNReal

namespace ExactManyToOne

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]

/-- **Factorisation under conditioning.** If `B` is independent of the vector
`(Y i)ᵢ` and the `Y i` are mutually independent, the probability that every pair
`(Y i, B)` lands in its set is the `B`-average of the product of the `K` inner
probabilities. -/
theorem prob_forall_eq_lintegral_prod
    {ι : Type*} [Fintype ι] {β : Type*} [MeasurableSpace β] {γ : ι → Type*}
    [∀ i, MeasurableSpace (γ i)]
    (B : Ω → β) (Y : ∀ i, Ω → γ i) (hB : Measurable B) (hY : ∀ i, Measurable (Y i))
    (hBY : IndepFun B (fun ω i => Y i ω) P) (hYY : iIndepFun Y P)
    (S : ∀ i, Set (γ i × β)) (hS : ∀ i, MeasurableSet (S i)) :
    P {ω | ∀ i, (Y i ω, B ω) ∈ S i}
      = ∫⁻ b, ∏ i, (P.map (Y i)) {y | (y, b) ∈ S i} ∂(P.map B) := by
  have hYv : Measurable (fun ω i => Y i ω) := measurable_pi_iff.2 hY
  have h1 : P.map (fun ω => (B ω, fun i => Y i ω))
      = (P.map B).prod (P.map fun ω i => Y i ω) :=
    (indepFun_iff_map_prod_eq_prod_map_map hB.aemeasurable hYv.aemeasurable).1 hBY
  have h2 : P.map (fun ω i => Y i ω) = Measure.pi (fun i => P.map (Y i)) :=
    (iIndepFun_iff_map_fun_eq_pi_map (fun i => (hY i).aemeasurable)).1 hYY
  set T : Set (β × ∀ i, γ i) := {p | ∀ i, (p.2 i, p.1) ∈ S i} with hTdef
  have hT : MeasurableSet T := by
    have : T = ⋂ i, (fun p : β × ∀ i, γ i => (p.2 i, p.1)) ⁻¹' S i := by
      ext p; simp [T]
    rw [this]
    exact MeasurableSet.iInter fun i =>
      (((measurable_pi_apply i).comp measurable_snd).prodMk measurable_fst) (hS i)
  have hpre : {ω | ∀ i, (Y i ω, B ω) ∈ S i} = (fun ω => (B ω, fun i => Y i ω)) ⁻¹' T := rfl
  rw [hpre, ← Measure.map_apply (hB.prodMk hYv) hT, h1, h2, Measure.prod_apply hT]
  refine lintegral_congr fun b => ?_
  have hslice : Prod.mk b ⁻¹' T = Set.univ.pi (fun i => {y | (y, b) ∈ S i}) := by
    ext y; simp [T]
  rw [hslice, Measure.pi_pi]

/-- The rearrangement of equation (14): for a positive variance term the event
`t ≤ c` is the event that the candidate's centred mean lies below a threshold
in which the benchmark appears only through `(u, v)`. -/
lemma div_sqrt_le_iff {z u c d : ℝ} (hd : 0 < d) :
    (z - u) / Real.sqrt d ≤ c ↔ z ≤ u + c * Real.sqrt d := by
  have hs : 0 < Real.sqrt d := Real.sqrt_pos.2 hd
  rw [div_le_iff₀ hs]
  constructor <;> intro h <;> linarith

/-- Under the global null the numerator is `Z i - U` (equation (13)). -/
lemma numerator_under_null (Xi X0 μ : ℝ) : Xi - X0 = (Xi - μ) - (X0 - μ) := by ring

/-- The statistic of equation (9), written in the centred variables of (12)–(13):
`Z i` is candidate `i`'s centred mean, `S i` its variance estimate, `U` the
benchmark's centred mean and `V = s₀²/n₀`. -/
noncomputable def tstat {K : ℕ} (n : Fin K → ℝ) (Z S : Fin K → Ω → ℝ) (U V : Ω → ℝ)
    (i : Fin K) (ω : Ω) : ℝ :=
  (Z i ω - U ω) / Real.sqrt (S i ω / n i + V ω)

/-- The inner function of equation (15), before the Gaussian is evaluated: the
probability, at a fixed benchmark pair `(u, v)`, that candidate `i` is accepted
at its constant. -/
noncomputable def inner {K : ℕ} (P : Measure Ω) (n : Fin K → ℝ) (Z S : Fin K → Ω → ℝ)
    (c : Fin K → ℝ) (i : Fin K) (b : ℝ × ℝ) : ℝ≥0∞ :=
  (P.map fun ω => (Z i ω, S i ω)) {y | y.1 ≤ b.1 + c i * Real.sqrt (y.2 / n i + b.2)}

lemma measurableSet_accept {K : ℕ} (n : Fin K → ℝ) (c : Fin K → ℝ) (i : Fin K) :
    MeasurableSet {p : (ℝ × ℝ) × (ℝ × ℝ) |
      p.1.1 ≤ p.2.1 + c i * Real.sqrt (p.1.2 / n i + p.2.2)} := by
  apply measurableSet_le
  · fun_prop
  · exact measurable_fst.comp measurable_snd |>.add
      (measurable_const.mul (Real.continuous_sqrt.measurable.comp (by fun_prop)))

/-- **Theorem 1 (the exact law, equation (16)).** The pairs `(Z i, S i)` are
mutually independent and independent of the benchmark pair `(U, V)`. Then the
joint acceptance probability is the average, over the law of `(U, V)` alone, of
the product of the `K` inner functions. The constants may differ by candidate. -/
theorem exact_law {K : ℕ} (n : Fin K → ℝ) (Z S : Fin K → Ω → ℝ) (U V : Ω → ℝ)
    (hZ : ∀ i, Measurable (Z i)) (hS : ∀ i, Measurable (S i))
    (hU : Measurable U) (hV : Measurable V)
    (hBY : IndepFun (fun ω => (U ω, V ω)) (fun ω i => (Z i ω, S i ω)) P)
    (hYY : iIndepFun (fun i ω => (Z i ω, S i ω)) P)
    (hpos : ∀ i ω, 0 < S i ω / n i + V ω) (c : Fin K → ℝ) :
    P {ω | ∀ i, tstat n Z S U V i ω ≤ c i}
      = ∫⁻ b, ∏ i, inner P n Z S c i b ∂(P.map fun ω => (U ω, V ω)) := by
  have hset : {ω | ∀ i, tstat n Z S U V i ω ≤ c i}
      = {ω | ∀ i, ((Z i ω, S i ω), (U ω, V ω)) ∈
          {p : (ℝ × ℝ) × (ℝ × ℝ) | p.1.1 ≤ p.2.1 + c i * Real.sqrt (p.1.2 / n i + p.2.2)}} := by
    ext ω
    simp only [Set.mem_ofPred_eq, tstat]
    exact forall_congr' fun i => div_sqrt_le_iff (hpos i ω)
  rw [hset, prob_forall_eq_lintegral_prod (fun ω => (U ω, V ω)) (fun i ω => (Z i ω, S i ω))
    (hU.prodMk hV) (fun i => (hZ i).prodMk (hS i)) hBY hYY _
    (fun i => measurableSet_accept n c i)]
  rfl

/-- **The inner function as a one-dimensional integral (equation (15)).** Within a
candidate the mean is independent of the variance estimate (the Gaussian property
of Section 3.1), so the inner function averages the mean's distribution function
over the law of the variance estimate alone. -/
theorem inner_eq_lintegral {K : ℕ} (n : Fin K → ℝ) (Z S : Fin K → Ω → ℝ)
    (hZ : ∀ i, Measurable (Z i)) (hS : ∀ i, Measurable (S i))
    (hind : ∀ i, IndepFun (Z i) (S i) P) (c : Fin K → ℝ) (i : Fin K) (b : ℝ × ℝ) :
    inner P n Z S c i b
      = ∫⁻ s, (P.map (Z i)) (Iic (b.1 + c i * Real.sqrt (s / n i + b.2))) ∂(P.map (S i)) := by
  have hprod : P.map (fun ω => (Z i ω, S i ω)) = (P.map (Z i)).prod (P.map (S i)) :=
    (indepFun_iff_map_prod_eq_prod_map_map (hZ i).aemeasurable (hS i).aemeasurable).1 (hind i)
  have hmeas : MeasurableSet {y : ℝ × ℝ | y.1 ≤ b.1 + c i * Real.sqrt (y.2 / n i + b.2)} := by
    apply measurableSet_le
    · fun_prop
    · exact measurable_const.add
        (measurable_const.mul (Real.continuous_sqrt.measurable.comp (by fun_prop)))
  rw [inner, hprod, Measure.prod_apply_symm hmeas]
  rfl

/-- With a centred Gaussian mean the integrand of `inner_eq_lintegral` is the normal
distribution function, as in equation (15). -/
theorem inner_gaussian {K : ℕ} (n : Fin K → ℝ) (Z S : Fin K → Ω → ℝ)
    (hZ : ∀ i, Measurable (Z i)) (hS : ∀ i, Measurable (S i))
    (hind : ∀ i, IndepFun (Z i) (S i) P) (w : Fin K → NNReal)
    (hlaw : ∀ i, P.map (Z i) = gaussianReal 0 (w i)) (c : Fin K → ℝ) (i : Fin K) (b : ℝ × ℝ) :
    inner P n Z S c i b
      = ∫⁻ s, ENNReal.ofReal
          (cdf (gaussianReal 0 (w i)) (b.1 + c i * Real.sqrt (s / n i + b.2))) ∂(P.map (S i)) := by
  rw [inner_eq_lintegral n Z S hZ hS hind c i b, hlaw i]
  simp_rw [ofReal_cdf]

/-- **The outer law factorises (Section 3.3).** The benchmark's mean and variance
estimate are independent, so the outer average is an iterated integral against a
product law: a normal in one variable and a scaled chi-square in the other. -/
theorem exact_law_iterated {K : ℕ} (n : Fin K → ℝ) (Z S : Fin K → Ω → ℝ) (U V : Ω → ℝ)
    (hZ : ∀ i, Measurable (Z i)) (hS : ∀ i, Measurable (S i))
    (hU : Measurable U) (hV : Measurable V)
    (hBY : IndepFun (fun ω => (U ω, V ω)) (fun ω i => (Z i ω, S i ω)) P)
    (hYY : iIndepFun (fun i ω => (Z i ω, S i ω)) P) (hUV : IndepFun U V P)
    (hpos : ∀ i ω, 0 < S i ω / n i + V ω) (c : Fin K → ℝ) :
    P {ω | ∀ i, tstat n Z S U V i ω ≤ c i}
      = ∫⁻ u, ∫⁻ v, ∏ i, inner P n Z S c i (u, v) ∂(P.map V) ∂(P.map U) := by
  rw [exact_law n Z S U V hZ hS hU hV hBY hYY hpos c,
    (indepFun_iff_map_prod_eq_prod_map_map hU.aemeasurable hV.aemeasurable).1 hUV]
  refine lintegral_prod _ ?_
  refine (Finset.measurable_prod _ fun i _ => ?_).aemeasurable
  exact measurable_measure_prodMk_right (measurableSet_accept n c i)

omit [IsProbabilityMeasure P] in
/-- **The joint law is monotone in the constants.** This is what makes the critical
constant well defined as the point where the law reaches `1 - α`. -/
theorem exact_law_mono {K : ℕ} (n : Fin K → ℝ) (Z S : Fin K → Ω → ℝ) (U V : Ω → ℝ)
    {c c' : Fin K → ℝ} (h : c ≤ c') :
    P {ω | ∀ i, tstat n Z S U V i ω ≤ c i} ≤ P {ω | ∀ i, tstat n Z S U V i ω ≤ c' i} :=
  measure_mono fun _ hω i => (hω i).trans (h i)

end ExactManyToOne
