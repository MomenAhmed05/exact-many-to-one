import W4lean.Unit2BFormalR1
open MeasureTheory ProbabilityTheory Set
noncomputable section
namespace Unit2B.FormalR1.Kappa

-- Raw model, K = 1, normalised sigma_f = sigma_0 = 1 (w2_factor.py:13-36, brute :144-162).
-- Equal-mean null: both series share the mean m.  θ = (λ, b₀, b₁), λ = σ₁²/σ₀².
def nObs : ℕ := 36
abbrev PanelIdx : Type := Fin 3 × Fin nObs
abbrev Panel : Type := PanelIdx → ℝ
-- z (0,t) = f_t, z (1,t) = ε₀,t, z (2,t) = ε₁,t / √λ; all iid N(0,1)
def panelLaw : Measure Panel := Measure.pi (fun _ : PanelIdx => gaussianReal 0 1)

def ret0 (m : ℝ) (θ : Param) (z : Panel) (t : Fin nObs) : ℝ :=
  m + θ.2.1 * z (0, t) + z (1, t)
def ret1 (m : ℝ) (θ : Param) (z : Panel) (t : Fin nObs) : ℝ :=
  m + θ.2.2 * z (0, t) + Real.sqrt θ.1 * z (2, t)
def sMean (x : Fin nObs → ℝ) : ℝ := (∑ t, x t) / (nObs : ℝ)
def sVar (x : Fin nObs → ℝ) : ℝ := (∑ t, (x t - sMean x) ^ 2) / ((nObs : ℝ) - 1)
def Tstat (m : ℝ) (θ : Param) (z : Panel) : ℝ :=
  (sMean (ret1 m θ z) - sMean (ret0 m θ z)) /
    Real.sqrt (sVar (ret1 m θ z) / (nObs : ℝ) + sVar (ret0 m θ z) / (nObs : ℝ))
-- the raw rejection probability P(T > c), unconditional
def RawRate (m c : ℝ) (θ : Param) : ℝ := (panelLaw {z | c < Tstat m θ z}).toReal

-- pivots of the reduction (ENCL-01 §2.1; D0 = aF θ)
def Nstd (m : ℝ) (θ : Param) (z : Panel) : ℝ :=
  Real.sqrt (nObs : ℝ) * (sMean (ret1 m θ z) - sMean (ret0 m θ z)) / Real.sqrt (aF θ)
def Sres (m : ℝ) (θ : Param) (z : Panel) : ℝ := nu * (sVar (ret0 m θ z) + sVar (ret1 m θ z))
def Mmat (θ : Param) : Matrix (Fin 2) (Fin 3) ℝ := !![θ.2.1, 1, 0; θ.2.2, 0, Real.sqrt θ.1]
def zc (z : Panel) (t : Fin nObs) : Fin 3 → ℝ := fun k => z (k, t) - sMean (fun s => z (k, s))

-- K-0: bridge constants
theorem nu_eq : nu = (nObs : ℝ) - 1 := by sorry
-- K-1: Gram matrix of the model is §6's Gmat
theorem gram_eq {θ : Param} (hθ : θ ∈ Box) : Mmat θ * (Mmat θ).transpose = Gmat θ := by sorry
-- K-2: S = Σ_t |M x_t|², pathwise, on centred coordinates
theorem Sres_eq_gram_form (m : ℝ) (θ : Param) (z : Panel) :
    Sres m θ z = ∑ t, ∑ i, (Matrix.mulVec (Mmat θ) (zc z t) i) ^ 2 := by sorry
-- K-3: S > 0 almost surely
theorem Sres_pos_ae (m : ℝ) {θ : Param} (hθ : θ ∈ Box) :
    ∀ᵐ z ∂panelLaw, 0 < Sres m θ z := by sorry
-- K-4: the event, pathwise (standardise by D0 = aF)
theorem kappa_event {m c : ℝ} {θ : Param} (hθ : θ ∈ Box) (hc : 0 < c) {z : Panel}
    (hS : 0 < Sres m θ z) :
    c < Tstat m θ z ↔ c * Real.sqrt (Sres m θ z / (nu * aF θ)) < Nstd m θ z := by sorry
-- K-5: numerator law, N(0, D0/n) standardised
theorem kappa_numerator_law (m : ℝ) {θ : Param} (hθ : θ ∈ Box) :
    panelLaw.map (Nstd m θ) = gaussianReal 0 1 := by sorry
-- K-6: numerator independent of S (means vs centred residuals)
theorem kappa_indep (m : ℝ) {θ : Param} (hθ : θ ∈ Box) :
    IndepFun (Nstd m θ) (Sres m θ) panelLaw := by sorry
-- K-7: S =_d μ₁ X₁ + μ₂ X₂, X_j iid χ²_ν central
theorem kappa_spectral {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
    {X1 X2 : Ω → ℝ} (hX1m : Measurable X1) (hX2m : Measurable X2)
    (hX1 : P.map X1 = gammaMeasure (35 / 2) (1 / 2))
    (hX2 : P.map X2 = gammaMeasure (35 / 2) (1 / 2))
    (hInd : IndepFun X1 X2 P) (m : ℝ) {θ : Param} (hθ : θ ∈ Box) :
    panelLaw.map (Sres m θ) = P.map (fun ω => mu1 θ * X1 ω + mu2 θ * X2 ω) := by sorry
-- S-κ (P6): exactly the shape of hKappa in comp_thm_2B, with R := RawRate m
theorem kappa_reduction {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω} [IsProbabilityMeasure P]
    {N X1 X2 : Ω → ℝ} (hNm : Measurable N) (hX1m : Measurable X1) (hX2m : Measurable X2)
    (hN : P.map N = gaussianReal 0 1)
    (hX1 : P.map X1 = gammaMeasure (35 / 2) (1 / 2))
    (hX2 : P.map X2 = gammaMeasure (35 / 2) (1 / 2))
    (hInd : iIndepFun ![N, X1, X2] P) (m : ℝ) :
    ∀ c, 0 < c → ∀ θ ∈ Box, RawRate m c θ = Phi P N X1 X2 c (r1 θ) (r2 θ) := by sorry

end Unit2B.FormalR1.Kappa
