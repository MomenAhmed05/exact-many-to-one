import W4lean.Unit2BSizeR11
open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.SizeConsR11

/-!
# 2B size bound, part 2a: `r1 ≥ 1/2` and the existence of `cCons`

No numerical premise, no engine, no bridge.  Every project LEMMA is used fully qualified
(`Unit2B.FormalR1.Proofs.*`, `Unit2B.FormalR1.Nonvac.*`, `Unit2B.FormalR1.SizeR11.*`).

* `two_mu1_ge_aF`: `aF θ ≤ 2 μ₁(θ)` for EVERY θ (Rayleigh quotient at `v = (-1, 1)`),
  via `2 μ₁ - aF = 2 b₀ b₁ + √discG` and `√discG ≥ |2 b₀ b₁|` from `discG_eq`.
* `r1_ge_half`: `1/2 ≤ r1 θ` on `Box`;  `r2_pos`: `0 < r2 θ` on `Box`.
* `cCons := sInf {c | 0 < c ∧ Phi P0 N0 X10 X20 c (1/2) 0 ≤ alpha}` on the canonical triple.
  The set is nonempty by continuity from above (Γ(35/2, 1/2) gives `{X₁ ≤ 0}` measure zero);
  the infimum is attained by continuity from below, because the event is a STRICT inequality;
  and `cCons > 0` because `Phi` at `c = 0` is `1/2 > alpha`.
-/

open Unit2B.FormalR1.Nonvac

/-! ## The eigenvalue inequality -/

theorem two_mu1_sub_aF (θ : Param) :
    2 * mu1 θ - aF θ = 2 * (θ.2.1 * θ.2.2) + Real.sqrt (discG θ) := by
  unfold mu1 aF trG g00 g11
  ring

theorem sqrt_discG_ge (θ : Param) : -(2 * (θ.2.1 * θ.2.2)) ≤ Real.sqrt (discG θ) := by
  have hD : discG θ = (g00 θ - g11 θ) ^ 2 + 4 * g01 θ ^ 2 :=
    Unit2B.FormalR1.Proofs.discG_eq θ
  have hg : g01 θ = θ.2.1 * θ.2.2 := rfl
  have hlow : (2 * (θ.2.1 * θ.2.2)) ^ 2 ≤ discG θ := by
    rw [hD, hg]
    nlinarith [sq_nonneg (g00 θ - g11 θ)]
  have habs := abs_le.mp (Real.abs_le_sqrt hlow)
  linarith [habs.1]

/-- `aF θ ≤ 2 μ₁(θ)` for every `θ`: the Rayleigh quotient of `Gmat θ` at `v = (-1, 1)` is
`aF θ / 2`, and `μ₁` is the top eigenvalue. -/
theorem two_mu1_ge_aF (θ : Param) : aF θ ≤ 2 * mu1 θ := by
  have h1 := Unit2B.FormalR1.SizeConsR11.two_mu1_sub_aF θ
  have h2 := Unit2B.FormalR1.SizeConsR11.sqrt_discG_ge θ
  linarith

theorem r1_ge_half_of {θ : Param} (ha : 0 < aF θ) : 1 / 2 ≤ r1 θ := by
  have h := Unit2B.FormalR1.SizeConsR11.two_mu1_ge_aF θ
  unfold r1
  rw [le_div_iff₀ ha]
  linarith

/-- `r1 θ ≥ 1/2` on the Box. -/
theorem r1_ge_half {θ : Param} (hθ : θ ∈ Box) : 1 / 2 ≤ r1 θ :=
  Unit2B.FormalR1.SizeConsR11.r1_ge_half_of (Unit2B.FormalR1.Proofs.box_pos hθ).1

/-- `r2 θ > 0` on the Box. -/
theorem r2_pos {θ : Param} (hθ : θ ∈ Box) : 0 < r2 θ := by
  unfold r2
  exact div_pos (Unit2B.FormalR1.Proofs.box_pos hθ).2.2.1 (Unit2B.FormalR1.Proofs.box_pos hθ).1

/-! ## The conservative rate `c ↦ Phi c (1/2) 0` -/

/-- The event of `Phi P N X1 X2 c (1/2) 0`, verbatim. -/
def consEvent {Ω : Type*} (N X1 X2 : Ω → ℝ) (c : ℝ) : Set Ω :=
  {ω | c * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω}

/-- The valid set of the conservative rule, verbatim the set inside `cCons`. -/
abbrev consSet {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω) (N X1 X2 : Ω → ℝ) (ε : ℝ) :
    Set ℝ :=
  {c | 0 < c ∧ Phi P N X1 X2 c (1 / 2) 0 ≤ ε}

def consSeq {Ω : Type*} (N X1 X2 : Ω → ℝ) (n : ℕ) : Set Ω :=
  Unit2B.FormalR1.SizeConsR11.consEvent N X1 X2 ((n : ℝ) + 1)

def consUp {Ω : Type*} (N X1 X2 : Ω → ℝ) (c0 : ℝ) (n : ℕ) : Set Ω :=
  Unit2B.FormalR1.SizeConsR11.consEvent N X1 X2 (c0 + 1 / ((n : ℝ) + 1))

theorem consEvent_measurable {Ω : Type*} [MeasurableSpace Ω] {N X1 X2 : Ω → ℝ}
    (hNm : Measurable N) (hX1m : Measurable X1) (hX2m : Measurable X2) (c : ℝ) :
    MeasurableSet (Unit2B.FormalR1.SizeConsR11.consEvent N X1 X2 c) := by
  change MeasurableSet {ω | c * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω}
  exact measurableSet_lt (by fun_prop) hNm

/-- Γ(k, r) puts no mass on `(-∞, 0]`. -/
theorem gamma_Iic_zero (k r : ℝ) : gammaMeasure k r (Iic 0) = 0 := by
  have hae : (Iic (0 : ℝ)) =ᵐ[volume] Iio 0 := Iio_ae_eq_Iic.symm
  have h := setLIntegral_congr (f := gammaPDF k r) hae
  rw [gammaMeasure, withDensity_apply _ measurableSet_Iic, h]
  exact lintegral_gammaPDF_of_nonpos le_rfl

theorem null_nonpos_of_gammaLaw {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {X : Ω → ℝ} (hXm : Measurable X) {k r : ℝ} (hX : P.map X = gammaMeasure k r) :
    P {ω | X ω ≤ 0} = 0 := by
  have h : P.map X (Iic 0) = P (X ⁻¹' Iic 0) := Measure.map_apply hXm measurableSet_Iic
  rw [hX, Unit2B.FormalR1.SizeConsR11.gamma_Iic_zero] at h
  exact h.symm

/-- EXISTENCE (nonemptiness): for every `ε > 0` some `c > 0` has `Phi c (1/2) 0 ≤ ε`.
Continuity from above along `c = n + 1`; the limit event lies in `{X₁ ≤ 0}`, which is null. -/
theorem cons_tail {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
    {N X1 X2 : Ω → ℝ} (hNm : Measurable N) (hX1m : Measurable X1) (hX2m : Measurable X2)
    (hX1 : P.map X1 = gammaMeasure (35 / 2) (1 / 2)) {ε : ℝ} (hε : 0 < ε) :
    ∃ c, 0 < c ∧ Phi P N X1 X2 c (1 / 2) 0 ≤ ε := by
  have hanti : Antitone (Unit2B.FormalR1.SizeConsR11.consSeq N X1 X2) := by
    intro n n' h
    change Unit2B.FormalR1.SizeConsR11.consSeq N X1 X2 n'
      ⊆ Unit2B.FormalR1.SizeConsR11.consSeq N X1 X2 n
    intro ω hω
    have hω' : ((n' : ℝ) + 1) * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω := hω
    change ((n : ℝ) + 1) * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω
    have hs := Real.sqrt_nonneg ((1 / 2 * X1 ω + 0 * X2 ω) / nu)
    have hc : ((n : ℝ) + 1) ≤ (n' : ℝ) + 1 := by
      have h' : (n : ℝ) ≤ (n' : ℝ) := Nat.cast_le.mpr h
      linarith
    exact lt_of_le_of_lt (mul_le_mul_of_nonneg_right hc hs) hω'
  have hsub : (⋂ n, Unit2B.FormalR1.SizeConsR11.consSeq N X1 X2 n) ⊆ {ω | X1 ω ≤ 0} := by
    intro ω hω
    have hall : ∀ n : ℕ,
        ((n : ℝ) + 1) * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω :=
      fun n => mem_iInter.mp hω n
    have hs := Real.sqrt_nonneg ((1 / 2 * X1 ω + 0 * X2 ω) / nu)
    have ht0 : Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) = 0 := by
      by_contra hne
      have hpos : 0 < Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) :=
        lt_of_le_of_ne hs (Ne.symm hne)
      obtain ⟨n, hn⟩ := exists_nat_gt (N ω / Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu))
      have h1 := hall n
      rw [div_lt_iff₀ hpos] at hn
      nlinarith
    have hnu : (0 : ℝ) < nu := by norm_num [nu]
    have h0 := Real.sqrt_eq_zero'.mp ht0
    have h2 := (div_le_iff₀ hnu).mp h0
    change X1 ω ≤ 0
    linarith
  have hnull : P (⋂ n, Unit2B.FormalR1.SizeConsR11.consSeq N X1 X2 n) = 0 :=
    measure_mono_null hsub (Unit2B.FormalR1.SizeConsR11.null_nonpos_of_gammaLaw hX1m hX1)
  have hinf : P (⋂ n, Unit2B.FormalR1.SizeConsR11.consSeq N X1 X2 n)
      = ⨅ n, P (Unit2B.FormalR1.SizeConsR11.consSeq N X1 X2 n) :=
    hanti.measure_iInter
      (fun n => (Unit2B.FormalR1.SizeConsR11.consEvent_measurable hNm hX1m hX2m
        ((n : ℝ) + 1)).nullMeasurableSet)
      ⟨0, measure_ne_top P _⟩
  have hlt : ⨅ n, P (Unit2B.FormalR1.SizeConsR11.consSeq N X1 X2 n) < ENNReal.ofReal ε := by
    rw [← hinf, hnull]
    exact ENNReal.ofReal_pos.mpr hε
  obtain ⟨n, hn⟩ := iInf_lt_iff.mp hlt
  refine ⟨(n : ℝ) + 1, by positivity, ?_⟩
  exact ENNReal.toReal_le_of_le_ofReal hε.le hn.le

/-- ATTAINMENT: the infimum of a nonempty valid set is itself valid.  Every `c > sInf` is
valid (upward closure by `Phi_antitone`), and the strict-inequality event at the infimum is
the increasing union of the events at `sInf + 1/(n+1)`, so continuity from below applies. -/
theorem cons_attain {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
    {N X1 X2 : Ω → ℝ} (hX1ae : ∀ᵐ ω ∂P, 0 ≤ X1 ω) (hX2ae : ∀ᵐ ω ∂P, 0 ≤ X2 ω)
    {ε : ℝ} (hε : 0 ≤ ε)
    (hne : (Unit2B.FormalR1.SizeConsR11.consSet P N X1 X2 ε).Nonempty) :
    Phi P N X1 X2 (sInf (Unit2B.FormalR1.SizeConsR11.consSet P N X1 X2 ε)) (1 / 2) 0 ≤ ε := by
  obtain ⟨c0, hc0⟩ : ∃ c0, sInf (Unit2B.FormalR1.SizeConsR11.consSet P N X1 X2 ε) = c0 :=
    ⟨_, rfl⟩
  rw [hc0]
  have habove : ∀ n : ℕ, Phi P N X1 X2 (c0 + 1 / ((n : ℝ) + 1)) (1 / 2) 0 ≤ ε := by
    intro n
    have hpos : (0 : ℝ) < 1 / ((n : ℝ) + 1) := by positivity
    have hlt : sInf (Unit2B.FormalR1.SizeConsR11.consSet P N X1 X2 ε)
        < c0 + 1 / ((n : ℝ) + 1) := by
      rw [hc0]
      linarith
    obtain ⟨s, hs, hsc⟩ := exists_lt_of_csInf_lt hne hlt
    exact (Unit2B.FormalR1.Proofs.Phi_antitone hX1ae hX2ae hs.1.le hsc.le le_rfl le_rfl).trans
      hs.2
  have hunion : Unit2B.FormalR1.SizeConsR11.consEvent N X1 X2 c0
      = ⋃ n : ℕ, Unit2B.FormalR1.SizeConsR11.consUp N X1 X2 c0 n := by
    ext ω
    rw [mem_iUnion]
    have hs := Real.sqrt_nonneg ((1 / 2 * X1 ω + 0 * X2 ω) / nu)
    constructor
    · intro h
      have h' : c0 * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω := h
      have ht1 : 0 < Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) + 1 := by linarith
      have hgap : 0 < (N ω - c0 * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu))
          / (Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) + 1) :=
        div_pos (by linarith) ht1
      obtain ⟨n, hn⟩ := exists_nat_one_div_lt hgap
      refine ⟨n, ?_⟩
      have hpos : (0 : ℝ) < 1 / ((n : ℝ) + 1) := by positivity
      rw [lt_div_iff₀ ht1] at hn
      change (c0 + 1 / ((n : ℝ) + 1)) * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω
      nlinarith
    · rintro ⟨n, hn⟩
      have hn' : (c0 + 1 / ((n : ℝ) + 1)) * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu)
          < N ω := hn
      have hpos : (0 : ℝ) < 1 / ((n : ℝ) + 1) := by positivity
      change c0 * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω
      nlinarith [mul_nonneg hpos.le hs]
  have hmono : Monotone (Unit2B.FormalR1.SizeConsR11.consUp N X1 X2 c0) := by
    intro n n' h
    change Unit2B.FormalR1.SizeConsR11.consUp N X1 X2 c0 n
      ⊆ Unit2B.FormalR1.SizeConsR11.consUp N X1 X2 c0 n'
    intro ω hω
    have hω' : (c0 + 1 / ((n : ℝ) + 1)) * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu)
        < N ω := hω
    change (c0 + 1 / ((n' : ℝ) + 1)) * Real.sqrt ((1 / 2 * X1 ω + 0 * X2 ω) / nu) < N ω
    have hs := Real.sqrt_nonneg ((1 / 2 * X1 ω + 0 * X2 ω) / nu)
    have hle : 1 / ((n' : ℝ) + 1) ≤ 1 / ((n : ℝ) + 1) :=
      one_div_le_one_div_of_le (by positivity)
        (by have h' : (n : ℝ) ≤ (n' : ℝ) := Nat.cast_le.mpr h; linarith)
    nlinarith [mul_le_mul_of_nonneg_right hle hs]
  have hmeas : P (Unit2B.FormalR1.SizeConsR11.consEvent N X1 X2 c0) ≤ ENNReal.ofReal ε := by
    rw [hunion, hmono.measure_iUnion]
    refine iSup_le (fun n => ?_)
    exact (ENNReal.le_ofReal_iff_toReal_le (measure_ne_top P _) hε).mpr (habove n)
  change (P (Unit2B.FormalR1.SizeConsR11.consEvent N X1 X2 c0)).toReal ≤ ε
  exact ENNReal.toReal_le_of_le_ofReal hε hmeas

/-- The infimum is strictly positive: at `c = 0` the rate is `P(N > 0) = 1/2 > ε`. -/
theorem cons_sInf_pos {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
    {N X1 X2 : Ω → ℝ} (hNm : Measurable N) (hN : P.map N = gaussianReal 0 1)
    (hX1ae : ∀ᵐ ω ∂P, 0 ≤ X1 ω) (hX2ae : ∀ᵐ ω ∂P, 0 ≤ X2 ω)
    {ε : ℝ} (hε0 : 0 ≤ ε) (hε : ε < 1 / 2)
    (hne : (Unit2B.FormalR1.SizeConsR11.consSet P N X1 X2 ε).Nonempty) :
    0 < sInf (Unit2B.FormalR1.SizeConsR11.consSet P N X1 X2 ε) := by
  have hnn : 0 ≤ sInf (Unit2B.FormalR1.SizeConsR11.consSet P N X1 X2 ε) :=
    le_csInf hne (fun c hc => hc.1.le)
  rcases hnn.lt_or_eq with h | h
  · exact h
  · exfalso
    have hat := Unit2B.FormalR1.SizeConsR11.cons_attain hX1ae hX2ae hε0 hne
    rw [← h] at hat
    have h00 : Phi P N X1 X2 0 (1 / 2) 0 = Phi P N X1 X2 0 0 0 := by
      simp only [Phi, zero_mul]
    have h0 : Phi P N X1 X2 0 (1 / 2) 0 = 1 / 2 := by
      rw [h00]
      exact Unit2B.FormalR1.Nonvac.gauss_half hNm hN 0
    linarith

/-! ## `cCons` on the canonical triple -/

/-- The conservative critical value.  Stated as an infimum; no decimal appears. -/
def cCons : ℝ := sInf {c | 0 < c ∧ Phi P0 N0 X10 X20 c (1 / 2) 0 ≤ alpha}

/-- EXISTENCE: the set defining `cCons` is nonempty. -/
theorem cCons_set_nonempty :
    {c | 0 < c ∧ Phi P0 N0 X10 X20 c (1 / 2) 0 ≤ alpha}.Nonempty := by
  obtain ⟨hNm, hX1m, hX2m, -, hX1, -, -⟩ := Unit2B.FormalR1.SizeR11.canonical_model
  obtain ⟨c, hc, hle⟩ := Unit2B.FormalR1.SizeConsR11.cons_tail hNm hX1m hX2m hX1
    (show (0 : ℝ) < alpha by norm_num [alpha])
  exact ⟨c, hc, hle⟩

/-- ATTAINMENT: `cCons` is positive and is itself a valid critical value of the conservative
rate. -/
theorem cCons_spec : 0 < cCons ∧ Phi P0 N0 X10 X20 cCons (1 / 2) 0 ≤ alpha := by
  obtain ⟨hNm, hX1m, hX2m, hN, hX1, hX2, -⟩ := Unit2B.FormalR1.SizeR11.canonical_model
  have hX1ae : ∀ᵐ ω ∂P0, 0 ≤ X10 ω :=
    Unit2B.FormalR1.Proofs.ae_nonneg_of_gammaLaw hX1m (by norm_num) (by norm_num) hX1
  have hX2ae : ∀ᵐ ω ∂P0, 0 ≤ X20 ω :=
    Unit2B.FormalR1.Proofs.ae_nonneg_of_gammaLaw hX2m (by norm_num) (by norm_num) hX2
  have hne : (Unit2B.FormalR1.SizeConsR11.consSet P0 N0 X10 X20 alpha).Nonempty :=
    Unit2B.FormalR1.SizeConsR11.cCons_set_nonempty
  exact ⟨Unit2B.FormalR1.SizeConsR11.cons_sInf_pos hNm hN hX1ae hX2ae
      (by norm_num [alpha]) (by norm_num [alpha]) hne,
    Unit2B.FormalR1.SizeConsR11.cons_attain hX1ae hX2ae (by norm_num [alpha]) hne⟩

#print axioms Unit2B.FormalR1.SizeConsR11.two_mu1_ge_aF
#print axioms Unit2B.FormalR1.SizeConsR11.r1_ge_half_of
#print axioms Unit2B.FormalR1.SizeConsR11.r1_ge_half
#print axioms Unit2B.FormalR1.SizeConsR11.r2_pos
#print axioms Unit2B.FormalR1.SizeConsR11.gamma_Iic_zero
#print axioms Unit2B.FormalR1.SizeConsR11.null_nonpos_of_gammaLaw
#print axioms Unit2B.FormalR1.SizeConsR11.cons_tail
#print axioms Unit2B.FormalR1.SizeConsR11.cons_attain
#print axioms Unit2B.FormalR1.SizeConsR11.cons_sInf_pos
#print axioms Unit2B.FormalR1.SizeConsR11.cCons
#print axioms Unit2B.FormalR1.SizeConsR11.cCons_set_nonempty
#print axioms Unit2B.FormalR1.SizeConsR11.cCons_spec

end Unit2B.FormalR1.SizeConsR11
