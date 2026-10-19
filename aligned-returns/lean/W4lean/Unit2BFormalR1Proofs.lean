import W4lean.Unit2BFormalR1
open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.Proofs

open Unit2B.FormalR1 in
/-- L-A. `Phi` is non-increasing in `c` on `[0, ∞)` and non-increasing in each of `ρ₁, ρ₂`,
given `X₁, X₂ ≥ 0` almost surely. The kernel is the judge. -/
theorem Phi_antitone {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsFiniteMeasure P]
    {N X1 X2 : Ω → ℝ} (hX1 : ∀ᵐ ω ∂P, 0 ≤ X1 ω) (hX2 : ∀ᵐ ω ∂P, 0 ≤ X2 ω)
    {c c' ρ1 ρ1' ρ2 ρ2' : ℝ} (hc : 0 ≤ c) (hcc' : c ≤ c') (h1 : ρ1 ≤ ρ1') (h2 : ρ2 ≤ ρ2') :
    Phi P N X1 X2 c' ρ1' ρ2' ≤ Phi P N X1 X2 c ρ1 ρ2 := by
  have hnu : 0 < nu := by norm_num [nu]
  have hg : ∀ᵐ ω ∂P, 0 ≤ X1 ω ∧ 0 ≤ X2 ω := hX1.and hX2
  have hgP : P {ω | ¬ (0 ≤ X1 ω ∧ 0 ≤ X2 ω)} = 0 := ae_iff.mp hg
  -- pointwise: on the conull good set, the unprimed threshold is at most the primed one
  have hmain : {ω | c' * Real.sqrt ((ρ1' * X1 ω + ρ2' * X2 ω) / nu) < N ω} ∩
      {ω | 0 ≤ X1 ω ∧ 0 ≤ X2 ω}
      ⊆ {ω | c * Real.sqrt ((ρ1 * X1 ω + ρ2 * X2 ω) / nu) < N ω} := by
    intro ω hω
    have hA : ρ1 * X1 ω + ρ2 * X2 ω ≤ ρ1' * X1 ω + ρ2' * X2 ω :=
      add_le_add (mul_le_mul_of_nonneg_right h1 hω.2.1) (mul_le_mul_of_nonneg_right h2 hω.2.2)
    have hA' : (ρ1 * X1 ω + ρ2 * X2 ω) / nu ≤ (ρ1' * X1 ω + ρ2' * X2 ω) / nu :=
      div_le_div_of_nonneg_right hA (le_of_lt hnu)
    have hsq : Real.sqrt ((ρ1 * X1 ω + ρ2 * X2 ω) / nu)
        ≤ Real.sqrt ((ρ1' * X1 ω + ρ2' * X2 ω) / nu) := Real.sqrt_le_sqrt hA'
    have hc' : 0 ≤ c' := by linarith
    have hkey : c * Real.sqrt ((ρ1 * X1 ω + ρ2 * X2 ω) / nu)
        ≤ c' * Real.sqrt ((ρ1' * X1 ω + ρ2' * X2 ω) / nu) := by
      refine le_trans (mul_le_mul_of_nonneg_right hcc' (Real.sqrt_nonneg _)) ?_
      simpa [mul_comm] using mul_le_mul_of_nonneg_right hsq hc'
    exact lt_of_le_of_lt hkey hω.1
  -- the a.e. hypothesis costs nothing: intersecting the primed event with the good set
  -- leaves its measure unchanged
  have hconull : P ({ω | c' * Real.sqrt ((ρ1' * X1 ω + ρ2' * X2 ω) / nu) < N ω} ∩
        {ω | 0 ≤ X1 ω ∧ 0 ≤ X2 ω})
      = P {ω | c' * Real.sqrt ((ρ1' * X1 ω + ρ2' * X2 ω) / nu) < N ω} :=
    measure_inter_conull hgP
  have hPt : P {ω | c' * Real.sqrt ((ρ1' * X1 ω + ρ2' * X2 ω) / nu) < N ω} ≤
      P {ω | c * Real.sqrt ((ρ1 * X1 ω + ρ2 * X2 ω) / nu) < N ω} :=
    le_trans hconull.symm.le (measure_mono hmain)
  unfold Phi
  exact ENNReal.toReal_mono (measure_ne_top P _) hPt

-- L-A′
theorem ae_nonneg_of_gammaLaw {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {X : Ω → ℝ} (hXm : Measurable X) {k r : ℝ} (hk : 0 < k) (hr : 0 < r)
    (hX : P.map X = gammaMeasure k r) : ∀ᵐ ω ∂P, 0 ≤ X ω := by
  have hG : ∀ᵐ x ∂(gammaMeasure k r), 0 ≤ x := by
    rw [ae_iff]
    have hset : {x : ℝ | ¬ 0 ≤ x} = Iio 0 := by
      ext x; simp [not_le]
    rw [hset, gammaMeasure, withDensity_apply _ measurableSet_Iio]
    exact lintegral_gammaPDF_of_nonpos le_rfl
  rw [← hX] at hG
  exact ae_of_ae_map hXm.aemeasurable hG

-- L-B
theorem charQ_eq_det (m : ℝ) (θ : Param) :
    charQ m θ = (Gmat θ - m • (1 : Matrix (Fin 2) (Fin 2) ℝ)).det := by
  rw [Matrix.det_fin_two]
  simp [Gmat, charQ, trG, detG, Matrix.one_apply] <;> ring

theorem discG_eq (θ : Param) : discG θ = (g00 θ - g11 θ) ^ 2 + 4 * g01 θ ^ 2 := by
  unfold discG trG detG
  ring

theorem discG_nonneg (θ : Param) : 0 ≤ discG θ := by
  rw [discG_eq]
  positivity

theorem mu_sum_prod (θ : Param) : mu1 θ + mu2 θ = trG θ ∧ mu1 θ * mu2 θ = detG θ := by
  have hs : Real.sqrt (discG θ) ^ 2 = discG θ := Real.sq_sqrt (discG_nonneg θ)
  refine ⟨by unfold mu1 mu2; ring, ?_⟩
  unfold mu1 mu2
  have e : (trG θ + Real.sqrt (discG θ)) / 2 * ((trG θ - Real.sqrt (discG θ)) / 2)
      = (trG θ ^ 2 - Real.sqrt (discG θ) ^ 2) / 4 := by ring
  rw [e, hs]
  unfold discG
  ring

theorem charQ_factor (m : ℝ) (θ : Param) : charQ m θ = (m - mu1 θ) * (m - mu2 θ) := by
  obtain ⟨hs, hp⟩ := mu_sum_prod θ
  unfold charQ
  have e : (m - mu1 θ) * (m - mu2 θ) = m ^ 2 - (mu1 θ + mu2 θ) * m + mu1 θ * mu2 θ := by ring
  rw [e, hs, hp]

theorem box_pos {θ : Param} (hθ : θ ∈ Box) :
    0 < aF θ ∧ 0 < detG θ ∧ 0 < mu2 θ ∧ mu2 θ ≤ mu1 θ := by
  simp only [Box, Set.mem_setOf_eq] at hθ
  obtain ⟨hl, -, -, -, -, -⟩ := hθ
  have hl0 : 0 < θ.1 := lt_of_lt_of_le (by norm_num) hl
  have ha : 0 < aF θ := by
    unfold aF
    nlinarith [sq_nonneg (θ.2.2 - θ.2.1)]
  have hdet_eq : detG θ = θ.2.1 ^ 2 * θ.1 + θ.2.2 ^ 2 + θ.1 := by
    unfold detG g00 g11 g01
    ring
  have hdet : 0 < detG θ := by
    rw [hdet_eq]
    nlinarith [mul_nonneg (sq_nonneg θ.2.1) hl0.le, sq_nonneg θ.2.2]
  have htr : 0 < trG θ := by
    unfold trG g00 g11
    nlinarith [sq_nonneg θ.2.1, sq_nonneg θ.2.2]
  have hlt : Real.sqrt (discG θ) < trG θ :=
    (Real.sqrt_lt' htr).mpr (by unfold discG; linarith)
  have hs := Real.sqrt_nonneg (discG θ)
  refine ⟨ha, hdet, ?_, ?_⟩
  · unfold mu2
    linarith
  · unfold mu1 mu2
    linarith

theorem r1_lower {θ : Param} (ha : 0 < aF θ) {ρ : ℝ}
    (h : ρ * aF θ ≤ trG θ / 2 ∨ charQ (ρ * aF θ) θ ≤ 0) : ρ ≤ r1 θ := by
  unfold r1
  rw [le_div_iff₀ ha]
  have hs : 0 ≤ Real.sqrt (discG θ) := Real.sqrt_nonneg _
  rcases h with h | h
  · have hhalf : trG θ / 2 ≤ mu1 θ := by
      unfold mu1
      linarith
    linarith
  · rw [charQ_factor] at h
    have h12 : mu2 θ ≤ mu1 θ := by
      have hs := Real.sqrt_nonneg (discG θ)
      unfold mu1 mu2
      linarith
    by_contra hm
    push_neg at hm
    have e1 : 0 < ρ * aF θ - mu1 θ := by linarith
    have e2 : 0 < ρ * aF θ - mu2 θ := by linarith
    linarith [mul_pos e1 e2]

theorem r2_lower {θ : Param} (ha : 0 < aF θ) {ρ : ℝ}
    (h1 : ρ * aF θ ≤ trG θ / 2) (h2 : 0 ≤ charQ (ρ * aF θ) θ) : ρ ≤ r2 θ := by
  unfold r2
  rw [le_div_iff₀ ha]
  rw [charQ_factor] at h2
  have hs : 0 ≤ Real.sqrt (discG θ) := Real.sqrt_nonneg _
  have hm1 : mu1 θ = (trG θ + Real.sqrt (discG θ)) / 2 := rfl
  have hm2 : mu2 θ = (trG θ - Real.sqrt (discG θ)) / 2 := rfl
  have hhalf : trG θ / 2 ≤ mu1 θ := by linarith
  by_contra hm
  push_neg at hm
  have h1f : ρ * aF θ - mu1 θ ≤ 0 := by linarith
  rcases lt_or_eq_of_le h1f with hlt0 | heq0
  · have e2 : 0 < ρ * aF θ - mu2 θ := by linarith
    linarith [mul_neg_of_neg_of_pos hlt0 e2]
  · have hsq0 : Real.sqrt (discG θ) = 0 := by linarith
    linarith

-- L-C
theorem cstar_bracket {Θ : Type*} (R : ℝ → Θ → ℝ) (S : Set Θ) (α L U : ℝ)
    (hanti : ∀ θ ∈ S, ∀ c c', 0 < c → c ≤ c' → R c' θ ≤ R c θ)
    (hL : 0 < L) (hLU : L ≤ U) (hU : ∀ θ ∈ S, R U θ ≤ α)
    {θ0 : Θ} (hθ0 : θ0 ∈ S) (hW : α < R L θ0) :
    U ∈ validSet R S α ∧ (∀ c ∈ validSet R S α, ∀ c', c ≤ c' → c' ∈ validSet R S α) ∧
    L ≤ cstarBox R S α ∧ cstarBox R S α ≤ U := by
  have hbdd : BddBelow (validSet R S α) := ⟨0, fun c hc => hc.1.le⟩
  have hUV : U ∈ validSet R S α := ⟨lt_of_lt_of_le hL hLU, hU⟩
  unfold cstarBox
  refine ⟨hUV, ?_, ?_, csInf_le hbdd hUV⟩
  · intro c hc c' hcc'
    exact ⟨lt_of_lt_of_le hc.1 hcc',
      fun θ hθ => (hanti θ hθ c c' hc.1 hcc').trans (hc.2 θ hθ)⟩
  · refine le_csInf ⟨U, hUV⟩ (fun c hc => ?_)
    by_contra hcL
    push_neg at hcL
    have hbad := (hanti θ0 hθ0 c L hc.1 hcL.le).trans (hc.2 θ0 hθ0)
    linarith

-- L-D
theorem cstarBox_eq_sSup {Θ : Type*} (R : ℝ → Θ → ℝ) (S : Set Θ) (α : ℝ) (hS : S.Nonempty)
    (hanti : ∀ θ ∈ S, ∀ c c', 0 < c → c ≤ c' → R c' θ ≤ R c θ)
    (hV : (validSet R S α).Nonempty) :
    cstarBox R S α = sSup (cstarPt R α '' S) := by
  have hsub : ∀ θ ∈ S, validSet R S α ⊆ {c | 0 < c ∧ R c θ ≤ α} :=
    fun θ hθ c hc => ⟨hc.1, hc.2 θ hθ⟩
  have hbddV : BddBelow (validSet R S α) := ⟨0, fun c hc => hc.1.le⟩
  have hbdd : BddAbove (cstarPt R α '' S) := by
    obtain ⟨v, hv⟩ := hV
    refine ⟨v, ?_⟩
    rintro _ ⟨θ, hθ, rfl⟩
    exact csInf_le (⟨0, fun x hx => hx.1.le⟩) (hsub θ hθ hv)
  have hne : (cstarPt R α '' S).Nonempty := hS.image _
  unfold cstarBox
  apply le_antisymm
  · by_contra hlt
    push_neg at hlt
    obtain ⟨c, hsc, hcI⟩ := exists_between hlt
    have key : ∀ θ ∈ S, 0 < c ∧ R c θ ≤ α := by
      intro θ hθ
      have ha : sInf {c | 0 < c ∧ R c θ ≤ α} < c :=
        lt_of_le_of_lt (le_csSup hbdd (mem_image_of_mem _ hθ)) hsc
      obtain ⟨d, hd, hdc⟩ := exists_lt_of_csInf_lt (hV.mono (hsub θ hθ)) ha
      exact ⟨hd.1.trans hdc, (hanti θ hθ d c hd.1 hdc.le).trans hd.2⟩
    obtain ⟨θ0, hθ0⟩ := hS
    have hcV : c ∈ validSet R S α := ⟨(key θ0 hθ0).1, fun θ hθ => (key θ hθ).2⟩
    have hinf := csInf_le hbddV hcV
    linarith
  · refine csSup_le hne ?_
    rintro _ ⟨θ, hθ, rfl⟩
    exact le_csInf hV (fun c hc => csInf_le (⟨0, fun x hx => hx.1.le⟩) (hsub θ hθ hc))

-- L-E
theorem cover_sound {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsFiniteMeasure P]
    {N X1 X2 : Ω → ℝ} (hX1 : ∀ᵐ ω ∂P, 0 ≤ X1 ω) (hX2 : ∀ᵐ ω ∂P, 0 ≤ X2 ω)
    (cells : List Cell) {c α : ℝ} (hc : 0 ≤ c)
    (hCover : ∀ θ ∈ Box, ∃ k ∈ cells, k.mem θ)
    (hLow1 : ∀ k ∈ cells, ∀ θ ∈ Box, k.mem θ →
      (k.ρ1 * aF θ ≤ trG θ / 2 ∨ charQ (k.ρ1 * aF θ) θ ≤ 0))
    (hLow2 : ∀ k ∈ cells, ∀ θ ∈ Box, k.mem θ →
      (k.ρ2 * aF θ ≤ trG θ / 2 ∧ 0 ≤ charQ (k.ρ2 * aF θ) θ))
    (hNum : ∀ k ∈ cells, Phi P N X1 X2 c k.ρ1 k.ρ2 ≤ α) :
    ∀ θ ∈ Box, Phi P N X1 X2 c (r1 θ) (r2 θ) ≤ α := by
  intro θ hθ
  obtain ⟨k, hk, hmem⟩ := hCover θ hθ
  have ha := (box_pos hθ).1
  have hρ1 : k.ρ1 ≤ r1 θ := r1_lower ha (hLow1 k hk θ hθ hmem)
  have hρ2 : k.ρ2 ≤ r2 θ :=
    r2_lower ha (hLow2 k hk θ hθ hmem).1 (hLow2 k hk θ hθ hmem).2
  exact (Phi_antitone hX1 hX2 hc le_rfl hρ1 hρ2).trans (hNum k hk)

-- COMP-THM (F-6)
theorem comp_thm_2B {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
    {N X1 X2 : Ω → ℝ} (hNm : Measurable N) (hX1m : Measurable X1) (hX2m : Measurable X2)
    (hN : P.map N = gaussianReal 0 1)
    (hX1 : P.map X1 = gammaMeasure (35 / 2) (1 / 2))
    (hX2 : P.map X2 = gammaMeasure (35 / 2) (1 / 2))
    (hInd : iIndepFun ![N, X1, X2] P)
    (R : ℝ → Param → ℝ)
    (hKappa : ∀ c, 0 < c → ∀ θ ∈ Box, R c θ = Phi P N X1 X2 c (r1 θ) (r2 θ))
    (cells : List Cell)
    (hCover : ∀ θ ∈ Box, ∃ k ∈ cells, k.mem θ)
    (hLow1 : ∀ k ∈ cells, ∀ θ ∈ Box, k.mem θ →
      (k.ρ1 * aF θ ≤ trG θ / 2 ∨ charQ (k.ρ1 * aF θ) θ ≤ 0))
    (hLow2 : ∀ k ∈ cells, ∀ θ ∈ Box, k.mem θ →
      (k.ρ2 * aF θ ≤ trG θ / 2 ∧ 0 ≤ charQ (k.ρ2 * aF θ) θ))
    (hNum : ∀ k ∈ cells, Phi P N X1 X2 Uc k.ρ1 k.ρ2 ≤ alpha)
    (hWit : alpha < Phi P N X1 X2 Lc (r1 corner) (r2 corner)) :
    Uc ∈ validSet R Box alpha ∧ Lc ≤ cstarBox R Box alpha ∧ cstarBox R Box alpha ≤ Uc ∧
    Uc - Lc ≤ tol ∧ cstarBox R Box alpha = sSup (cstarPt R alpha '' Box) := by
  have hX1ae : ∀ᵐ ω ∂P, 0 ≤ X1 ω :=
    ae_nonneg_of_gammaLaw hX1m (by norm_num) (by norm_num) hX1
  have hX2ae : ∀ᵐ ω ∂P, 0 ≤ X2 ω :=
    ae_nonneg_of_gammaLaw hX2m (by norm_num) (by norm_num) hX2
  have hanti : ∀ θ ∈ Box, ∀ c c', 0 < c → c ≤ c' → R c' θ ≤ R c θ := by
    intro θ hθ c c' hc hcc'
    rw [hKappa c' (lt_of_lt_of_le hc hcc') θ hθ, hKappa c hc θ hθ]
    exact Phi_antitone hX1ae hX2ae hc.le hcc' le_rfl le_rfl
  have hUc : (0 : ℝ) < Uc := by norm_num [Uc]
  have hLc : (0 : ℝ) < Lc := by norm_num [Lc]
  have hU : ∀ θ ∈ Box, R Uc θ ≤ alpha := by
    intro θ hθ
    rw [hKappa Uc hUc θ hθ]
    exact cover_sound hX1ae hX2ae cells hUc.le hCover hLow1 hLow2 hNum θ hθ
  have hcorner : corner ∈ Box := by
    simp only [Box, corner, Set.mem_setOf_eq]
    norm_num
  have hW : alpha < R Lc corner := by
    rw [hKappa Lc hLc corner hcorner]
    exact hWit
  obtain ⟨hUV, -, hLow, hUp⟩ :=
    cstar_bracket R Box alpha Lc Uc hanti hLc (by norm_num [Lc, Uc]) hU hcorner hW
  exact ⟨hUV, hLow, hUp, by norm_num [Uc, Lc, tol],
    cstarBox_eq_sSup R Box alpha ⟨corner, hcorner⟩ hanti ⟨Uc, hUV⟩⟩

end Unit2B.FormalR1.Proofs
