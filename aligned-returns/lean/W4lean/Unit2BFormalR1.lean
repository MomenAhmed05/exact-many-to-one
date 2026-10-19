import Mathlib
open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1

def nu : ℝ := 35
def alpha : ℝ := 1 / 40
def Lc : ℝ := 271731 / 100000
def Uc : ℝ := 271740 / 100000
def tol : ℝ := 1 / 10000

abbrev Param := ℝ × ℝ × ℝ   -- θ = (λ, b₀, b₁)
def Box : Set Param :=
  {θ | 1 / 25 ≤ θ.1 ∧ θ.1 ≤ 25 ∧ -(3 / 2) ≤ θ.2.1 ∧ θ.2.1 ≤ 3 / 2 ∧
        -(3 / 2) ≤ θ.2.2 ∧ θ.2.2 ≤ 3 / 2}
def corner : Param := (1 / 25, -(3 / 2), 3 / 2)

def aF (θ : Param) : ℝ := (θ.2.2 - θ.2.1) ^ 2 + 1 + θ.1
def g00 (θ : Param) : ℝ := θ.2.1 ^ 2 + 1
def g11 (θ : Param) : ℝ := θ.2.2 ^ 2 + θ.1
def g01 (θ : Param) : ℝ := θ.2.1 * θ.2.2
def Gmat (θ : Param) : Matrix (Fin 2) (Fin 2) ℝ := !![g00 θ, g01 θ; g01 θ, g11 θ]
def trG (θ : Param) : ℝ := g00 θ + g11 θ
def detG (θ : Param) : ℝ := g00 θ * g11 θ - g01 θ ^ 2
def discG (θ : Param) : ℝ := trG θ ^ 2 - 4 * detG θ
def mu1 (θ : Param) : ℝ := (trG θ + Real.sqrt (discG θ)) / 2
def mu2 (θ : Param) : ℝ := (trG θ - Real.sqrt (discG θ)) / 2
def r1 (θ : Param) : ℝ := mu1 θ / aF θ
def r2 (θ : Param) : ℝ := mu2 θ / aF θ
def charQ (m : ℝ) (θ : Param) : ℝ := m ^ 2 - trG θ * m + detG θ

def Phi {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω) (N X1 X2 : Ω → ℝ)
    (c ρ1 ρ2 : ℝ) : ℝ :=
  (P {ω | c * Real.sqrt ((ρ1 * X1 ω + ρ2 * X2 ω) / nu) < N ω}).toReal

def validSet {Θ : Type*} (R : ℝ → Θ → ℝ) (S : Set Θ) (α : ℝ) : Set ℝ :=
  {c | 0 < c ∧ ∀ θ ∈ S, R c θ ≤ α}
def cstarBox {Θ : Type*} (R : ℝ → Θ → ℝ) (S : Set Θ) (α : ℝ) : ℝ := sInf (validSet R S α)
def cstarPt {Θ : Type*} (R : ℝ → Θ → ℝ) (α : ℝ) (θ : Θ) : ℝ := sInf {c | 0 < c ∧ R c θ ≤ α}

structure Cell where
  lo : Param
  hi : Param
  ρ1 : ℝ
  ρ2 : ℝ
def Cell.mem (k : Cell) (θ : Param) : Prop :=
  k.lo.1 ≤ θ.1 ∧ θ.1 ≤ k.hi.1 ∧ k.lo.2.1 ≤ θ.2.1 ∧ θ.2.1 ≤ k.hi.2.1 ∧
    k.lo.2.2 ≤ θ.2.2 ∧ θ.2.2 ≤ k.hi.2.2

-- L-A
theorem Phi_antitone {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsFiniteMeasure P]
    {N X1 X2 : Ω → ℝ} (hX1 : ∀ᵐ ω ∂P, 0 ≤ X1 ω) (hX2 : ∀ᵐ ω ∂P, 0 ≤ X2 ω)
    {c c' ρ1 ρ1' ρ2 ρ2' : ℝ} (hc : 0 ≤ c) (hcc' : c ≤ c') (h1 : ρ1 ≤ ρ1') (h2 : ρ2 ≤ ρ2') :
    Phi P N X1 X2 c' ρ1' ρ2' ≤ Phi P N X1 X2 c ρ1 ρ2 := by sorry
-- L-A′
theorem ae_nonneg_of_gammaLaw {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {X : Ω → ℝ} (hXm : Measurable X) {k r : ℝ} (hk : 0 < k) (hr : 0 < r)
    (hX : P.map X = gammaMeasure k r) : ∀ᵐ ω ∂P, 0 ≤ X ω := by sorry
-- L-B
theorem charQ_eq_det (m : ℝ) (θ : Param) :
    charQ m θ = (Gmat θ - m • (1 : Matrix (Fin 2) (Fin 2) ℝ)).det := by sorry
theorem discG_eq (θ : Param) : discG θ = (g00 θ - g11 θ) ^ 2 + 4 * g01 θ ^ 2 := by sorry
theorem charQ_factor (m : ℝ) (θ : Param) : charQ m θ = (m - mu1 θ) * (m - mu2 θ) := by sorry
theorem mu_sum_prod (θ : Param) : mu1 θ + mu2 θ = trG θ ∧ mu1 θ * mu2 θ = detG θ := by sorry
theorem box_pos {θ : Param} (hθ : θ ∈ Box) :
    0 < aF θ ∧ 0 < detG θ ∧ 0 < mu2 θ ∧ mu2 θ ≤ mu1 θ := by sorry
theorem r1_lower {θ : Param} (ha : 0 < aF θ) {ρ : ℝ}
    (h : ρ * aF θ ≤ trG θ / 2 ∨ charQ (ρ * aF θ) θ ≤ 0) : ρ ≤ r1 θ := by sorry
theorem r2_lower {θ : Param} (ha : 0 < aF θ) {ρ : ℝ}
    (h1 : ρ * aF θ ≤ trG θ / 2) (h2 : 0 ≤ charQ (ρ * aF θ) θ) : ρ ≤ r2 θ := by sorry
-- L-C
theorem cstar_bracket {Θ : Type*} (R : ℝ → Θ → ℝ) (S : Set Θ) (α L U : ℝ)
    (hanti : ∀ θ ∈ S, ∀ c c', 0 < c → c ≤ c' → R c' θ ≤ R c θ)
    (hL : 0 < L) (hLU : L ≤ U) (hU : ∀ θ ∈ S, R U θ ≤ α)
    {θ0 : Θ} (hθ0 : θ0 ∈ S) (hW : α < R L θ0) :
    U ∈ validSet R S α ∧ (∀ c ∈ validSet R S α, ∀ c', c ≤ c' → c' ∈ validSet R S α) ∧
    L ≤ cstarBox R S α ∧ cstarBox R S α ≤ U := by sorry
-- L-D
theorem cstarBox_eq_sSup {Θ : Type*} (R : ℝ → Θ → ℝ) (S : Set Θ) (α : ℝ) (hS : S.Nonempty)
    (hanti : ∀ θ ∈ S, ∀ c c', 0 < c → c ≤ c' → R c' θ ≤ R c θ)
    (hV : (validSet R S α).Nonempty) :
    cstarBox R S α = sSup (cstarPt R α '' S) := by sorry
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
    ∀ θ ∈ Box, Phi P N X1 X2 c (r1 θ) (r2 θ) ≤ α := by sorry
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
    Uc - Lc ≤ tol ∧ cstarBox R Box alpha = sSup (cstarPt R alpha '' Box) := by sorry

end Unit2B.FormalR1
