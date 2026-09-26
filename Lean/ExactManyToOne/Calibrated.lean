/-
Copyright (c) 2026 Momen Ahmed. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Momen Ahmed
-/
import ExactManyToOne.Reduction

/-!
# Section 10: the calibrated rule, the union step, both sides and the joint projection

* The calibrated rule is valid as soon as the constraint (39) holds at the level it
  uses; nothing depends on that level being the largest feasible one (Section 10.2).
* The union step (42) certifies the per-candidate rule at any number of candidates
  from the one-candidate certificate, `FWER ≤ K · sup FWER⁽¹⁾`. Both the product form
  inside the expectation and Boole's inequality are proved.
* At one candidate the two-sided error rate is exactly twice the one-sided one
  (equation (43)), because the null law is symmetric and the two tails are disjoint,
  and this holds with a data-dependent threshold as long as the threshold reads only
  the variance estimates. Hence `FWER₂ₛ ≤ 2K · sup FWER⁽¹⁾₁ₛ`.
* The joint rule is certified by projecting the deployed table (equation (44)).
-/

open MeasureTheory ProbabilityTheory Set
open scoped ENNReal

namespace ExactManyToOne

variable {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}

/-- **Validity of the calibrated rule (equations (39)–(40)).** If the unconditional
error rate at level `a` is at most the ceiling at every ratio in the declared domain,
the rule run at `a` is valid there. The level need only be feasible, not maximal. -/
theorem calibrated_valid {Λ : Type*} (FWER : Λ → ℝ → ℝ) (D : Set Λ) (α a : ℝ)
    (hfeas : ∀ lam ∈ D, FWER lam a ≤ α) : ∀ lam ∈ D, FWER lam a ≤ α := hfeas

/-- **The Weierstrass product inequality**: for factors in `[0, 1]`,
`1 - ∏ gᵢ ≤ ∑ (1 - gᵢ)`. This is Boole's inequality at a fixed benchmark state. -/
theorem one_sub_prod_le_sum {ι : Type*} (s : Finset ι) (g : ι → ℝ)
    (h0 : ∀ i ∈ s, 0 ≤ g i) (h1 : ∀ i ∈ s, g i ≤ 1) :
    1 - ∏ i ∈ s, g i ≤ ∑ i ∈ s, (1 - g i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [Finset.prod_insert ha, Finset.sum_insert ha]
    have hp0 : 0 ≤ ∏ i ∈ s, g i := Finset.prod_nonneg fun i hi => h0 i (Finset.mem_insert_of_mem hi)
    have hp1 : ∏ i ∈ s, g i ≤ 1 :=
      Finset.prod_le_one₀ (fun i hi => h0 i (Finset.mem_insert_of_mem hi))
        (fun i hi => h1 i (Finset.mem_insert_of_mem hi))
    have ih' := ih (fun i hi => h0 i (Finset.mem_insert_of_mem hi))
      (fun i hi => h1 i (Finset.mem_insert_of_mem hi))
    have ha0 := h0 a (Finset.mem_insert_self a s)
    have ha1 := h1 a (Finset.mem_insert_self a s)
    nlinarith

/-- **No rejection, conditionally independent (the product inside (42)).** For a
rule that rejects candidate `i` when its pair with the benchmark lands in `R i`, the
probability that no candidate is rejected is the benchmark average of the product
of the `K` conditional non-rejection probabilities. -/
theorem no_rejection_eq {ι : Type*} [Fintype ι] {β : Type*} [MeasurableSpace β]
    [IsProbabilityMeasure P]
    {γ : ι → Type*} [∀ i, MeasurableSpace (γ i)]
    (B : Ω → β) (Y : ∀ i, Ω → γ i) (hB : Measurable B) (hY : ∀ i, Measurable (Y i))
    (hBY : IndepFun B (fun ω i => Y i ω) P) (hYY : iIndepFun Y P)
    (R : ∀ i, Set (γ i × β)) (hR : ∀ i, MeasurableSet (R i)) :
    P {ω | ∀ i, (Y i ω, B ω) ∈ (R i)ᶜ}
      = ∫⁻ b, ∏ i, (P.map (Y i)) {y | (y, b) ∈ (R i)ᶜ} ∂(P.map B) :=
  prob_forall_eq_lintegral_prod B Y hB hY hBY hYY _ fun i => (hR i).compl

/-- **Boole's inequality for the familywise event.** -/
theorem fwer_le_sum {ι : Type*} [Fintype ι] (Rej : ι → Set Ω) :
    P (⋃ i, Rej i) ≤ ∑ i, P (Rej i) := measure_iUnion_fintype_le P Rej

/-- **Theorem (union step, equation (42)).** Each candidate's marginal rejection
probability is the one-candidate error rate at that candidate's own ratio, and the
one-candidate error rate is certified at most `b` over the declared domain. Then the
familywise error rate at any ratio vector in the domain is at most `K · b`. -/
theorem union_step {K : ℕ} {Λ : Type*} (Rej : Fin K → Set Ω) (FWER1 : Λ → ℝ≥0∞)
    (lam : Fin K → Λ) (D : Set Λ) (hD : ∀ i, lam i ∈ D)
    (hmarg : ∀ i, P (Rej i) = FWER1 (lam i)) (b : ℝ≥0∞) (hcert : ∀ l ∈ D, FWER1 l ≤ b) :
    P (⋃ i, Rej i) ≤ K * b := by
  calc P (⋃ i, Rej i) ≤ ∑ i, P (Rej i) := fwer_le_sum Rej
    _ ≤ ∑ _i : Fin K, b := Finset.sum_le_sum fun i _ => (hmarg i).trans_le (hcert _ (hD i))
    _ = K * b := by simp

/-- **Two-sided identity at one candidate (equation (43)).** If the pair
(statistic, what the threshold reads) has the same law as (−statistic, the same),
and the threshold is nonnegative, then the two-sided rejection probability is
exactly twice the one-sided one: the tails are disjoint and mirror images. -/
theorem two_sided_eq_two_mul {L : Type*} [MeasurableSpace L] (T : Ω → ℝ) (Lv : Ω → L)
    (hT : Measurable T) (hL : Measurable Lv) (c : L → ℝ) (hc : Measurable c)
    (hc0 : ∀ l, 0 ≤ c l)
    (hsym : P.map (fun ω => (-T ω, Lv ω)) = P.map (fun ω => (T ω, Lv ω))) :
    P {ω | c (Lv ω) < |T ω|} = 2 * P {ω | c (Lv ω) < T ω} := by
  have hA : MeasurableSet {p : ℝ × L | c p.2 < p.1} :=
    measurableSet_lt (hc.comp measurable_snd) measurable_fst
  have hsplit : {ω | c (Lv ω) < |T ω|} = {ω | c (Lv ω) < T ω} ∪ {ω | c (Lv ω) < -T ω} := by
    ext ω; simp only [mem_ofPred_eq, mem_union, lt_abs]
  have hdisj : Disjoint {ω | c (Lv ω) < T ω} {ω | c (Lv ω) < -T ω} := by
    rw [Set.disjoint_left]
    intro ω h1 h2
    simp only [mem_ofPred_eq] at h1 h2
    linarith [hc0 (Lv ω)]
  have hmirror : P {ω | c (Lv ω) < -T ω} = P {ω | c (Lv ω) < T ω} := by
    have e1 : P {ω | c (Lv ω) < -T ω}
        = P.map (fun ω => (-T ω, Lv ω)) {p : ℝ × L | c p.2 < p.1} := by
      rw [Measure.map_apply (f := fun ω => (-T ω, Lv ω)) (hT.neg.prodMk hL) hA]; rfl
    have e2 : P {ω | c (Lv ω) < T ω}
        = P.map (fun ω => (T ω, Lv ω)) {p : ℝ × L | c p.2 < p.1} := by
      rw [Measure.map_apply (f := fun ω => (T ω, Lv ω)) (hT.prodMk hL) hA]; rfl
    rw [e1, e2, hsym]
  have hm2 : MeasurableSet {ω | c (Lv ω) < -T ω} :=
    measurableSet_lt (hc.comp hL) hT.neg
  rw [hsplit, measure_union hdisj hm2, hmirror, two_mul]

/-- **The sign flip that makes the null law symmetric.** Flipping the signs of the
candidate's and the benchmark's centred means flips the statistic of (9) and leaves
the variance estimates alone. -/
lemma tstat_neg (z u s v n : ℝ) :
    (-z - -u) / Real.sqrt (s / n + v) = -((z - u) / Real.sqrt (s / n + v)) := by ring

/-- **The null law of one comparison is symmetric.** If the centred means `(Z, U)` are
independent of the variance estimates `(S, V)` and their joint law is invariant under
the sign flip, then (statistic, any function of the variance estimates) has the same
law as (−statistic, the same function). Centred Gaussians independent of each other
satisfy the hypothesis (`gaussian_pair_symm`). -/
theorem tstat_symm {L : Type*} [MeasurableSpace L] (n : ℝ) (Z U S V : Ω → ℝ)
    [IsProbabilityMeasure P]
    (hZ : Measurable Z) (hU : Measurable U) (hS : Measurable S) (hV : Measurable V)
    (f : ℝ × ℝ → L) (hf : Measurable f)
    (hind : IndepFun (fun ω => (Z ω, U ω)) (fun ω => (S ω, V ω)) P)
    (hflip : P.map (fun ω => (-Z ω, -U ω)) = P.map (fun ω => (Z ω, U ω))) :
    P.map (fun ω => (-((Z ω - U ω) / Real.sqrt (S ω / n + V ω)), f (S ω, V ω)))
      = P.map (fun ω => ((Z ω - U ω) / Real.sqrt (S ω / n + V ω), f (S ω, V ω))) := by
  -- both sides are images of the law of ((Z, U), (S, V)) under one measurable map,
  -- composed on the left with the flip or not
  set g : (ℝ × ℝ) × (ℝ × ℝ) → ℝ × L :=
    fun p => ((p.1.1 - p.1.2) / Real.sqrt (p.2.1 / n + p.2.2), f p.2) with hg
  have hgm : Measurable g := by
    refine Measurable.prodMk ?_ (hf.comp measurable_snd)
    exact (by fun_prop : Measurable fun p : (ℝ × ℝ) × (ℝ × ℝ) => p.1.1 - p.1.2).div
      (Real.continuous_sqrt.measurable.comp (by fun_prop))
  have hZU : Measurable fun ω => (Z ω, U ω) := hZ.prodMk hU
  have hSV : Measurable fun ω => (S ω, V ω) := hS.prodMk hV
  have hjoint : P.map (fun ω => ((Z ω, U ω), (S ω, V ω)))
      = (P.map fun ω => (Z ω, U ω)).prod (P.map fun ω => (S ω, V ω)) :=
    (indepFun_iff_map_prod_eq_prod_map_map hZU.aemeasurable hSV.aemeasurable).1 hind
  have hflipm : Measurable fun q : ℝ × ℝ => (-q.1, -q.2) := by fun_prop
  -- the flipped pair is still independent of the variance estimates, with the same law
  have hjoint' : P.map (fun ω => ((-Z ω, -U ω), (S ω, V ω)))
      = (P.map fun ω => (Z ω, U ω)).prod (P.map fun ω => (S ω, V ω)) := by
    have : (fun ω => ((-Z ω, -U ω), (S ω, V ω)))
        = (Prod.map (fun q : ℝ × ℝ => (-q.1, -q.2)) id) ∘ (fun ω => ((Z ω, U ω), (S ω, V ω))) := rfl
    rw [this, ← Measure.map_map (hflipm.prodMap (measurable_id (α := ℝ × ℝ)))
      (hZU.prodMk hSV), hjoint,
      ← Measure.map_prod_map _ _ hflipm (measurable_id (α := ℝ × ℝ)), Measure.map_id]
    rw [Measure.map_map hflipm hZU]
    exact congrArg (fun m : Measure (ℝ × ℝ) => m.prod (P.map fun ω => (S ω, V ω))) hflip
  have lhs : (fun ω => (-((Z ω - U ω) / Real.sqrt (S ω / n + V ω)), f (S ω, V ω)))
      = g ∘ (fun ω => ((-Z ω, -U ω), (S ω, V ω))) := by
    funext ω; simp only [Function.comp, hg, tstat_neg]
  have rhs : (fun ω => ((Z ω - U ω) / Real.sqrt (S ω / n + V ω), f (S ω, V ω)))
      = g ∘ (fun ω => ((Z ω, U ω), (S ω, V ω))) := rfl
  have hflipω : Measurable fun ω => ((-Z ω, -U ω), (S ω, V ω)) :=
    (hZ.neg.prodMk hU.neg).prodMk hSV
  rw [lhs, rhs, ← Measure.map_map hgm hflipω,
    ← Measure.map_map hgm (hZU.prodMk hSV), hjoint, hjoint']

/-- Independent centred Gaussians are invariant under the joint sign flip. -/
theorem gaussian_pair_symm (Z U : Ω → ℝ) (hZ : Measurable Z) (hU : Measurable U)
    [IsProbabilityMeasure P]
    (hind : IndepFun Z U P) (w w0 : NNReal)
    (hZlaw : P.map Z = gaussianReal 0 w) (hUlaw : P.map U = gaussianReal 0 w0) :
    P.map (fun ω => (-Z ω, -U ω)) = P.map (fun ω => (Z ω, U ω)) := by
  have hprod : P.map (fun ω => (Z ω, U ω)) = (P.map Z).prod (P.map U) :=
    (indepFun_iff_map_prod_eq_prod_map_map hZ.aemeasurable hU.aemeasurable).1 hind
  have hneg : ∀ v : NNReal, (gaussianReal 0 v).map (fun x => -x) = gaussianReal 0 v := by
    intro v; rw [gaussianReal_map_neg, neg_zero]
  have : (fun ω => (-Z ω, -U ω))
      = (Prod.map (fun x : ℝ => -x) (fun x : ℝ => -x)) ∘ (fun ω => (Z ω, U ω)) := rfl
  rw [this, ← Measure.map_map (measurable_neg.prodMap measurable_neg) (hZ.prodMk hU), hprod,
    ← Measure.map_prod_map _ _ measurable_neg measurable_neg, hZlaw, hUlaw, hneg, hneg]

/-- **Theorem (two-sided certificate, equation (43)).** With the per-candidate
two-sided rule, each candidate's two-sided rejection probability is twice its
one-sided one, so a one-candidate one-sided certificate `b` gives
`FWER₂ₛ ≤ 2K · b`. The identity is used one candidate at a time, never as a
`K`-candidate statement, which fails above one candidate. -/
theorem two_sided_union {K : ℕ} {Λ : Type*} (Rej2 Rej1 : Fin K → Set Ω) (FWER1 : Λ → ℝ≥0∞)
    (lam : Fin K → Λ) (D : Set Λ) (hD : ∀ i, lam i ∈ D)
    (htwice : ∀ i, P (Rej2 i) = 2 * P (Rej1 i))
    (hmarg : ∀ i, P (Rej1 i) = FWER1 (lam i)) (b : ℝ≥0∞) (hcert : ∀ l ∈ D, FWER1 l ≤ b) :
    P (⋃ i, Rej2 i) ≤ 2 * K * b := by
  calc P (⋃ i, Rej2 i) ≤ ∑ i, P (Rej2 i) := fwer_le_sum Rej2
    _ ≤ ∑ _i : Fin K, 2 * b := Finset.sum_le_sum fun i _ => by
        rw [htwice i, hmarg i]; gcongr; exact hcert _ (hD i)
    _ = 2 * K * b := by simp; ring

/-- **The projection of the deployed joint table (equation (44)).** `ψ i t` is the
smallest threshold the table assigns anywhere with candidate `i`'s ratio at `t`.
It is a minimum over the deployed table itself, so it is at most the threshold
actually applied. -/
noncomputable def projection {K : ℕ} {B : Type*} [Fintype B] [Nonempty B] [DecidableEq B]
    (cj : (Fin K → B) → ℝ) (i : Fin K) (t : B) : ℝ :=
  (Finset.univ.filter fun l : Fin K → B => l i = t).inf'
    ⟨fun _ => t, by simp⟩ cj

lemma projection_le {K : ℕ} {B : Type*} [Fintype B] [Nonempty B] [DecidableEq B]
    (cj : (Fin K → B) → ℝ) (l : Fin K → B) (i : Fin K) : projection cj i (l i) ≤ cj l :=
  Finset.inf'_le _ (by simp)

/-- **Theorem (joint certificate).** Rejection of candidate `i` by the joint rule
implies `Tᵢ > ψᵢ(λ̂ᵢ)`, so the familywise event lies in the union of `K`
one-candidate events, each bounded by a one-candidate certificate. -/
theorem joint_projection_bound {K : ℕ} {B : Type*} [Fintype B] [Nonempty B] [DecidableEq B]
    (T : Fin K → Ω → ℝ) (lhat : Ω → Fin K → B) (cj : (Fin K → B) → ℝ)
    (bnd : Fin K → ℝ≥0∞)
    (hcert : ∀ i, P {ω | projection cj i (lhat ω i) < T i ω} ≤ bnd i) :
    P {ω | ∃ i, cj (lhat ω) < T i ω} ≤ ∑ i, bnd i := by
  calc P {ω | ∃ i, cj (lhat ω) < T i ω}
      ≤ P (⋃ i, {ω | projection cj i (lhat ω i) < T i ω}) := by
        refine measure_mono fun ω ⟨i, hi⟩ => mem_iUnion.2 ⟨i, ?_⟩
        exact lt_of_le_of_lt (projection_le cj (lhat ω) i) hi
    _ ≤ ∑ i, P {ω | projection cj i (lhat ω i) < T i ω} := fwer_le_sum _
    _ ≤ ∑ i, bnd i := Finset.sum_le_sum fun i _ => hcert i

end ExactManyToOne
