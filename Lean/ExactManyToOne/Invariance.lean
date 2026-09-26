/-
Copyright (c) 2026 Momen Ahmed. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Momen Ahmed
-/
import Mathlib

/-!
# Section 10.8: the affine group, the maximal invariant, and what Hunt–Stein adds

The testing problem is unchanged by `x ↦ a + b x` with `b > 0` applied to every
group at once (equation (47)). This file

* builds that group and proves it solvable (translations form an abelian normal
  subgroup with abelian quotient), which is the hypothesis under which the
  Hunt–Stein theorem is quoted;
* proves that the statistics and the estimated ratios are invariant, and that they
  form a **maximal** invariant: two sufficient statistics that share them are related
  by an explicit group element. The paper proves this at one candidate and records
  the extension above one candidate as unwritten. The proof here is for any number
  of candidates `K`, and `maximal_invariant_one` is the paper's statement at `K = 1`;
* shows that a statistic that is invariant has the same law at every parameter in an
  orbit when the model is equivariant, and that the parameter's orbits are the level
  sets of the standardised effects and the true ratios, so the oracle's power (the
  loss's reference point) depends on the parameter only through them;
* states the Hunt–Stein conclusion as a named hypothesis. It is **not** proved here:
  it is quoted from the literature, as the paper does.
-/

open MeasureTheory

namespace ExactManyToOne

/-- The affine group of the line: `x ↦ a + b x` with `b > 0`. -/
@[ext] structure Aff where
  a : ℝ
  b : ℝ
  hb : 0 < b

namespace Aff

instance : Mul Aff := ⟨fun g h => ⟨g.a + g.b * h.a, g.b * h.b, mul_pos g.hb h.hb⟩⟩
instance : One Aff := ⟨⟨0, 1, one_pos⟩⟩
noncomputable instance : Inv Aff := ⟨fun g => ⟨-g.a / g.b, 1 / g.b, one_div_pos.2 g.hb⟩⟩

@[simp] lemma mul_a (g h : Aff) : (g * h).a = g.a + g.b * h.a := rfl
@[simp] lemma mul_b (g h : Aff) : (g * h).b = g.b * h.b := rfl
@[simp] lemma one_a : (1 : Aff).a = 0 := rfl
@[simp] lemma one_b : (1 : Aff).b = 1 := rfl
@[simp] lemma inv_a (g : Aff) : g⁻¹.a = -g.a / g.b := rfl
@[simp] lemma inv_b (g : Aff) : g⁻¹.b = 1 / g.b := rfl

noncomputable instance : Group Aff where
  mul_assoc f g h := by ext <;> simp <;> ring
  one_mul g := by ext <;> simp
  mul_one g := by ext <;> simp
  inv_mul_cancel g := by
    have := g.hb.ne'
    ext
    · simp only [mul_a, inv_a, inv_b, one_a]
      field_simp
      ring
    · simp only [mul_b, inv_b, one_b]
      field_simp

/-- Translations `x ↦ a + x`. -/
noncomputable def translation : Multiplicative ℝ →* Aff where
  toFun a := ⟨Multiplicative.toAdd a, 1, one_pos⟩
  map_one' := by ext <;> simp
  map_mul' x y := by ext <;> simp

/-- The scale part `(a, b) ↦ b`, into the units of `ℝ`. -/
noncomputable def scale : Aff →* ℝˣ where
  toFun g := Units.mk0 g.b g.hb.ne'
  map_one' := by ext; simp
  map_mul' g h := by ext; simp

/-- **The affine group is solvable.** The kernel of the scale map is the image of the
translations, and both are abelian. Solvable groups are amenable, which is the
hypothesis of the Hunt–Stein theorem. Amenability itself is not in Mathlib and is
quoted, not proved. -/
instance : Group.IsSolvable Aff := by
  refine Group.isSolvable_of_ker_le_range translation scale ?_
  intro g hg
  have hb : g.b = 1 := by
    have := congrArg Units.val (MonoidHom.mem_ker.1 hg)
    simpa [scale] using this
  exact ⟨Multiplicative.ofAdd g.a, by ext <;> simp [translation, hb]⟩

end Aff

/-- A point of the sufficient statistic for one benchmark and `K` candidates: the
benchmark's mean and variance estimate, and each candidate's. The same type serves
for the parameter `(μ₀, σ₀², μᵢ, σᵢ²)`, on which the group acts in the same way. -/
@[ext] structure Suff (K : ℕ) where
  x0 : ℝ
  q0 : ℝ
  x : Fin K → ℝ
  q : Fin K → ℝ

namespace Suff

variable {K : ℕ}

/-- The group acts on means by `a + b x` and on variances by `b² q`. -/
instance : SMul Aff (Suff K) :=
  ⟨fun g s => ⟨g.a + g.b * s.x0, g.b ^ 2 * s.q0, fun i => g.a + g.b * s.x i,
    fun i => g.b ^ 2 * s.q i⟩⟩

@[simp] lemma smul_x0 (g : Aff) (s : Suff K) : (g • s).x0 = g.a + g.b * s.x0 := rfl
@[simp] lemma smul_q0 (g : Aff) (s : Suff K) : (g • s).q0 = g.b ^ 2 * s.q0 := rfl
@[simp] lemma smul_x (g : Aff) (s : Suff K) (i) : (g • s).x i = g.a + g.b * s.x i := rfl
@[simp] lemma smul_q (g : Aff) (s : Suff K) (i) : (g • s).q i = g.b ^ 2 * s.q i := rfl

noncomputable instance : MulAction Aff (Suff K) where
  one_smul s := by ext <;> simp
  mul_smul g h s := by ext <;> simp <;> ring

/-- The statistic of equation (9) for candidate `i`, with sample sizes `n` and `n₀`.
Read on the parameter instead of the data it is the standardised effect `δᵢ`. -/
noncomputable def t (n : Fin K → ℝ) (n0 : ℝ) (s : Suff K) (i : Fin K) : ℝ :=
  (s.x i - s.x0) / Real.sqrt (s.q i / n i + s.q0 / n0)

/-- The estimated ratio `λ̂ᵢ = sᵢ²/s₀²`. Read on the parameter it is the true `λᵢ`. -/
noncomputable def ratio (s : Suff K) (i : Fin K) : ℝ := s.q i / s.q0

/-- **Invariance of the statistic.** -/
theorem t_smul (n : Fin K → ℝ) (n0 : ℝ) (g : Aff) (s : Suff K) (i : Fin K) :
    t n n0 (g • s) i = t n n0 s i := by
  have hb := g.hb
  have hden : g.b ^ 2 * s.q i / n i + g.b ^ 2 * s.q0 / n0
      = g.b ^ 2 * (s.q i / n i + s.q0 / n0) := by
    ring
  simp only [t, smul_x, smul_x0, smul_q, smul_q0]
  rw [hden, Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hb.le,
    show g.a + g.b * s.x i - (g.a + g.b * s.x0) = g.b * (s.x i - s.x0) by ring,
    mul_div_mul_left _ _ hb.ne']

/-- **Invariance of the estimated ratio.** -/
theorem ratio_smul (g : Aff) (s : Suff K) (i : Fin K) : ratio (g • s) i = ratio s i := by
  simp only [ratio, smul_q, smul_q0]
  exact mul_div_mul_left _ _ (pow_ne_zero 2 g.hb.ne')

/-- Positive sufficient statistics: all variance estimates strictly positive. -/
def Pos (s : Suff K) : Prop := 0 < s.q0 ∧ ∀ i, 0 < s.q i

/-- **The maximal invariant, at any number of candidates.** Two positive sufficient
statistics that share every statistic and every estimated ratio are related by the
explicit group element `b = √(q₀'/q₀)`, `a = x₀' − b x₀`. -/
theorem maximal_invariant (n : Fin K → ℝ) (n0 : ℝ) (hn : ∀ i, 0 < n i) (hn0 : 0 < n0)
    (s s' : Suff K) (hs : Pos s) (hs' : Pos s')
    (ht : ∀ i, t n n0 s i = t n n0 s' i) (hr : ∀ i, ratio s i = ratio s' i) :
    ∃ g : Aff, g • s = s' := by
  obtain ⟨hq0, hq⟩ := hs
  obtain ⟨hq0', hq'⟩ := hs'
  set b := Real.sqrt (s'.q0 / s.q0) with hbdef
  have hbpos : 0 < b := Real.sqrt_pos.2 (div_pos hq0' hq0)
  have hb2 : b ^ 2 = s'.q0 / s.q0 := Real.sq_sqrt (div_pos hq0' hq0).le
  refine ⟨⟨s'.x0 - b * s.x0, b, hbpos⟩, ?_⟩
  -- variances: b² q = q' follows from the ratios
  have hqi : ∀ i, b ^ 2 * s.q i = s'.q i := by
    intro i
    have h := hr i
    simp only [ratio] at h
    rw [hb2]
    field_simp at h ⊢
    linarith [h]
  have hq0eq : b ^ 2 * s.q0 = s'.q0 := by rw [hb2]; field_simp
  -- denominators scale by b
  have hden : ∀ i, Real.sqrt (s'.q i / n i + s'.q0 / n0)
      = b * Real.sqrt (s.q i / n i + s.q0 / n0) := by
    intro i
    rw [← hqi i, ← hq0eq,
      show b ^ 2 * s.q i / n i + b ^ 2 * s.q0 / n0 = b ^ 2 * (s.q i / n i + s.q0 / n0) by ring,
      Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hbpos.le]
  ext
  · simp
  · simpa using hq0eq
  · rename_i i
    have hpos : 0 < Real.sqrt (s.q i / n i + s.q0 / n0) :=
      Real.sqrt_pos.2 (add_pos (div_pos (hq i) (hn i)) (div_pos hq0 hn0))
    have h := ht i
    simp only [t] at h
    rw [hden i, eq_div_iff (mul_pos hbpos hpos).ne'] at h
    have hkey : (s.x i - s.x0) / Real.sqrt (s.q i / n i + s.q0 / n0)
        * (b * Real.sqrt (s.q i / n i + s.q0 / n0)) = b * (s.x i - s.x0) := by
      rw [mul_comm b, ← mul_assoc, div_mul_cancel₀ _ hpos.ne', mul_comm]
    simp only [smul_x]
    linarith [h, hkey]
  · rename_i i
    simpa using hqi i

/-- **The paper's statement at one candidate (Section 10.8).** The pair
`(t, λ̂)` is a maximal invariant for the one-candidate problem. -/
theorem maximal_invariant_one (n0 n1 : ℝ) (hn0 : 0 < n0) (hn1 : 0 < n1) (s s' : Suff 1)
    (hs : Pos s) (hs' : Pos s')
    (ht : t (fun _ => n1) n0 s 0 = t (fun _ => n1) n0 s' 0) (hr : ratio s 0 = ratio s' 0) :
    ∃ g : Aff, g • s = s' :=
  maximal_invariant (fun _ => n1) n0 (fun _ => hn1) hn0 s s' hs hs'
    (fun i => by fin_cases i; exact ht) (fun i => by fin_cases i; exact hr)

/-- **An invariant function factors through the maximal invariant.** Any function
of the sufficient statistic that the group leaves unchanged takes the same value at
any two positive points with the same `(t, λ̂)`. This is why the search of
Section 10.7, which ranges over functions of `(t, λ̂)`, covers every invariant test. -/
theorem invariant_factors (n : Fin K → ℝ) (n0 : ℝ) (hn : ∀ i, 0 < n i) (hn0 : 0 < n0)
    {β : Type*} (φ : Suff K → β) (hφ : ∀ (g : Aff) s, φ (g • s) = φ s)
    (s s' : Suff K) (hs : Pos s) (hs' : Pos s')
    (ht : ∀ i, t n n0 s i = t n n0 s' i) (hr : ∀ i, ratio s i = ratio s' i) :
    φ s = φ s' := by
  obtain ⟨g, rfl⟩ := maximal_invariant n n0 hn hn0 s s' hs hs' ht hr
  exact (hφ g s).symm

/-- **The law of an invariant statistic is constant on orbits.** If the model is
equivariant (the law of `g • X` under `θ` is the law of `X` under `g • θ`) and the
statistic is invariant, its law is the same at `θ` and at `g • θ`. With the
parameter's orbits the level sets of `(δ, λ)` (the same `maximal_invariant` read on the
parameter), the oracle's power is a function of `(δ, λ)` only, which is what makes
worst-case regret against it an invariant loss. -/
theorem law_invariant {X Θ M : Type*} [MeasurableSpace X] [MeasurableSpace M]
    [SMul Aff X] [SMul Aff Θ] (Pθ : Θ → Measure X)
    (hequiv : ∀ (g : Aff) θ, (Pθ θ).map (g • ·) = Pθ (g • θ))
    (hmeas : ∀ g : Aff, Measurable (g • · : X → X))
    (T : X → M) (hT : Measurable T) (hinv : ∀ (g : Aff) x, T (g • x) = T x)
    (g : Aff) (θ : Θ) : (Pθ (g • θ)).map T = (Pθ θ).map T := by
  rw [← hequiv, Measure.map_map hT (hmeas g)]
  congr 1
  funext x
  exact hinv g x

/-- **The Hunt–Stein conclusion, as a named hypothesis.** For a group-invariant loss
and an amenable group, every test is matched or beaten in worst-case risk by an
invariant one. This is quoted from the literature (Lehmann and Romano 2005,
Section 8.5, with Bondar and Milnes 1981 for the passage from almost-invariance to
invariance), exactly as in the paper. It is not proved in this development. -/
def HuntSteinConclusion {Test : Type*} (risk : Test → ℝ) (invariant : Test → Prop) : Prop :=
  ∀ φ, ∃ φ', invariant φ' ∧ risk φ' ≤ risk φ

/-- **The floor extends from invariant tests to all tests.** If no invariant test
beats `floor` in worst-case regret, and the Hunt–Stein conclusion holds, then no
test at all beats it. This is the step Section 10.8 takes, and the one place where
the development uses a result it does not prove. -/
theorem floor_all_tests {Test : Type*} (risk : Test → ℝ) (invariant : Test → Prop)
    (hHS : HuntSteinConclusion risk invariant) (floor : ℝ)
    (hfloor : ∀ φ, invariant φ → floor ≤ risk φ) : ∀ φ, floor ≤ risk φ := by
  intro φ
  obtain ⟨φ', hinv, hle⟩ := hHS φ
  exact (hfloor φ' hinv).trans hle

/-- **The linear program (48) is linear.** The rejection probability of a randomised
rule on the discretised `(t, λ̂)` plane is `∑ q_bk m_bk`, which is linear in `q`, so a
convex combination of feasible rules is feasible and the program's optimum is
global. -/
theorem lp_convex {ι : Type*} (s : Finset ι) (m : ι → ℝ) (q q' : ι → ℝ) (θ α : ℝ)
    (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) (hq : ∑ k ∈ s, q k * m k ≤ α) (hq' : ∑ k ∈ s, q' k * m k ≤ α) :
    ∑ k ∈ s, (θ * q k + (1 - θ) * q' k) * m k ≤ α := by
  have : ∑ k ∈ s, (θ * q k + (1 - θ) * q' k) * m k
      = θ * ∑ k ∈ s, q k * m k + (1 - θ) * ∑ k ∈ s, q' k * m k := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun k _ => by ring
  rw [this]
  nlinarith

end Suff

end ExactManyToOne
