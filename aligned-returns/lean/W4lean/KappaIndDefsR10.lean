/-
Model definitions, transcribed verbatim from the written statement of the model (§1).

Mathlib rev 5ed2965256430c3649e86755f9576b54eca72435, toolchain v4.34.0.

NAMESPACES (§1 "NAMESPACES"): the statements are given in isolation but in
the sources they are preceded by

    open MeasureTheory ProbabilityTheory Set

which is required here, or nothing resolves.  In particular IndepFun and
iIndepFun live in ProbabilityTheory at this rev.

NOTE on the model: ret0 and ret1 share the SINGLE parameter m, so m cancels
identically in sMean (ret1 m th z) - sMean (ret0 m th z).  Verified symbolically
(indep/math_check.py item 2): the difference of means is
    (b1-b0) * mean_t z(0,t) + sqrt(lam) * mean_t z(2,t) - mean_t z(1,t)
with no m, and the coefficient vector (b1-b0, -1, sqrt(lam)) = u^T M for
u = e1 - e0 has squared norm exactly aF.  Hence Nstd is standard normal for
every m, which is K-5.
-/
import Mathlib.MeasureTheory.Constructions.Pi
import Mathlib.Probability.Distributions.Gamma
import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.Probability.Independence.Basic
import Mathlib.LinearAlgebra.Matrix.Notation

open MeasureTheory ProbabilityTheory Set
open scoped Topology

noncomputable section

namespace KappaInd

def nu : ℝ := 35
def alpha : ℝ := 1 / 40
def Lc : ℝ := 271731 / 100000
def Uc : ℝ := 271740 / 100000
def tol : ℝ := 1 / 10000

abbrev Param := ℝ × ℝ × ℝ

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

/-! ### The model proper -/

def nObs : ℕ := 36

abbrev PanelIdx : Type := Fin 3 × Fin nObs

abbrev Panel : Type := PanelIdx → ℝ

/- z (0,t) = f_t, z (1,t) = ε₀,t, z (2,t) = ε₁,t / √λ; all iid N(0,1) -/
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

/-- The raw rejection probability P(T > c), unconditional. -/
def RawRate (m c : ℝ) (θ : Param) : ℝ := (panelLaw {z | c < Tstat m θ z}).toReal

/- Pivots of the reduction (ENCL-01 §2.1; D0 = aF θ). -/
def Nstd (m : ℝ) (θ : Param) (z : Panel) : ℝ :=
  Real.sqrt (nObs : ℝ) * (sMean (ret1 m θ z) - sMean (ret0 m θ z)) / Real.sqrt (aF θ)

def Sres (m : ℝ) (θ : Param) (z : Panel) : ℝ := nu * (sVar (ret0 m θ z) + sVar (ret1 m θ z))

def Mmat (θ : Param) : Matrix (Fin 2) (Fin 3) ℝ := !![θ.2.1, 1, 0; θ.2.2, 0, Real.sqrt θ.1]

def zc (z : Panel) (t : Fin nObs) : Fin 3 → ℝ := fun k => z (k, t) - sMean (fun s => z (k, s))

end KappaInd