import W4lean.Unit2BKappaProofsR4
import W4lean.OrthonormalCoordsR1
-- Stage P′, sub-lemma P-3a `coordLaw`.
-- Stated for an arbitrary orthonormal basis `b`, so it does not depend on D-6.
-- Route of R2 §2: the banked `panelLaw_map_toLp`, the banked `stdGaussian_repr`, `map_map`.
-- Signatures confirmed by source inspection at Mathlib 5ed2965256430c3649e86755f9576b54eca72435:
--   `Measure.map_map` (MeasureTheory/Measure/Map.lean:208), `Measure.map_id` (:199),
--   `WithLp.measurable_toLp` (Analysis/Normed/Lp/MeasurableSpace.lean:33, `@[fun_prop]`),
--   `WithLp.measurable_ofLp` (:30, `@[fun_prop]`), `WithLp.ofLp_toLp` (Normed/Lp/WithLp.lean:97,
--   `rfl`), `Fin.addCommMonoid` (Algebra/Group/Fin/Basic.lean:36, needs `[NeZero n]`).
-- Banked in this project:
--   `panelLaw_map_toLp` (Unit2BKappaProofsR4.lean:39, namespace `ProofsR4`),
--   `stdGaussian_repr` (OrthonormalCoordsR1.lean:41, namespace `ProofsR4`),
--   `panelLaw`, `Panel`, `PanelIdx`, `nObs` (Unit2BKappaStatementR1.lean:8-12).
-- `stdGaussian_repr` carries the auto-included section variables
-- `[DecidableEq ι] [Zero ι] [AddCommMonoid ι]` (reported by the linter in the
-- audit build of AxAll.lean).  At `ι = Fin 3 × Fin nObs` the last two
-- need `NeZero nObs`, which instance search cannot see through the `def`, so it is supplied.

open MeasureTheory ProbabilityTheory
noncomputable section
namespace Unit2B.FormalR1.Kappa.ProofsR5

/-- **P-3a `coordLaw`.**  For any orthonormal basis `b` of the Euclidean carrier, the
coordinates of the panel in `b`, read back as a plain function `PanelIdx → ℝ`, are i.i.d.
standard normal. -/
theorem coordLaw (b : OrthonormalBasis PanelIdx ℝ (EuclideanSpace ℝ PanelIdx)) :
    panelLaw.map (fun z : Panel => WithLp.ofLp (b.repr (WithLp.toLp 2 z)))
      = Measure.pi fun _ : PanelIdx => gaussianReal 0 1 := by
  haveI : NeZero nObs := ⟨by decide⟩
  have hto : Measurable (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx) := by fun_prop
  have hrepr : Measurable (b.repr : EuclideanSpace ℝ PanelIdx → EuclideanSpace ℝ PanelIdx) := by
    fun_prop
  have hof : Measurable (WithLp.ofLp : EuclideanSpace ℝ PanelIdx → Panel) := by fun_prop
  have h1 : panelLaw.map (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx)
      = stdGaussian (EuclideanSpace ℝ PanelIdx) :=
    Unit2B.FormalR1.Kappa.ProofsR4.panelLaw_map_toLp
  have h2 : (stdGaussian (EuclideanSpace ℝ PanelIdx)).map
        (b.repr : EuclideanSpace ℝ PanelIdx → EuclideanSpace ℝ PanelIdx)
      = (Measure.pi fun _ : PanelIdx => gaussianReal 0 1).map
        (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx) :=
    Unit2B.FormalR1.Kappa.ProofsR4.stdGaussian_repr (ι := PanelIdx) b
  have h3 : ((Measure.pi fun _ : PanelIdx => gaussianReal 0 1).map
        (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx)).map
        (WithLp.ofLp : EuclideanSpace ℝ PanelIdx → Panel)
      = Measure.pi fun _ : PanelIdx => gaussianReal 0 1 := by
    -- `ofLp ∘ toLp 2` is `id` by structure eta (`WithLp.ofLp_toLp` is `rfl`), so `map_id`
    -- closes it by defeq; no `rw` on the composite, which would need a syntactic match.
    calc ((Measure.pi fun _ : PanelIdx => gaussianReal 0 1).map
          (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx)).map
          (WithLp.ofLp : EuclideanSpace ℝ PanelIdx → Panel)
        = (Measure.pi fun _ : PanelIdx => gaussianReal 0 1).map
          ((WithLp.ofLp : EuclideanSpace ℝ PanelIdx → Panel) ∘
            (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx)) := Measure.map_map hof hto
      _ = Measure.pi fun _ : PanelIdx => gaussianReal 0 1 := Measure.map_id
  calc panelLaw.map (fun z : Panel => WithLp.ofLp (b.repr (WithLp.toLp 2 z)))
      = (panelLaw.map ((b.repr : EuclideanSpace ℝ PanelIdx → EuclideanSpace ℝ PanelIdx) ∘
          (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx))).map
          (WithLp.ofLp : EuclideanSpace ℝ PanelIdx → Panel) :=
        (Measure.map_map hof (hrepr.comp hto)).symm
    _ = ((panelLaw.map (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx)).map
          (b.repr : EuclideanSpace ℝ PanelIdx → EuclideanSpace ℝ PanelIdx)).map
          (WithLp.ofLp : EuclideanSpace ℝ PanelIdx → Panel) := by
        rw [Measure.map_map hrepr hto]
    _ = ((stdGaussian (EuclideanSpace ℝ PanelIdx)).map
          (b.repr : EuclideanSpace ℝ PanelIdx → EuclideanSpace ℝ PanelIdx)).map
          (WithLp.ofLp : EuclideanSpace ℝ PanelIdx → Panel) := by
        rw [h1]
    _ = ((Measure.pi fun _ : PanelIdx => gaussianReal 0 1).map
          (WithLp.toLp 2 : Panel → EuclideanSpace ℝ PanelIdx)).map
          (WithLp.ofLp : EuclideanSpace ℝ PanelIdx → Panel) := by
        rw [h2]
    _ = Measure.pi fun _ : PanelIdx => gaussianReal 0 1 := h3

-- P-3b.  No pull-back lemma of this shape exists at the pinned rev: a search of
-- `Mathlib/Probability` for `IndepFun` lemmas taking a hypothesis under `μ.map _` found none.
-- Built from `indepFun_iff_measure_inter_preimage_eq_mul`
-- (Probability/Independence/Basic.lean:644) and `Measure.map_apply` (Measure/Map.lean:170).
-- The final `exact` uses `R ⁻¹' (A ∩ B) = R ⁻¹' A ∩ R ⁻¹' B` and
-- `(f ∘ R) ⁻¹' s = R ⁻¹' (f ⁻¹' s)` by defeq; both are `rfl` in Mathlib (Data/Set/Image.lean:75, :141).

/-- **P-3b, the pull-back lemma.**  Independence under the pushforward `μ.map R` transfers to
independence of the composites under `μ`. -/
theorem indepFun_comp_of_map {Ω α β γ : Type*} [MeasurableSpace Ω] [MeasurableSpace α]
    [MeasurableSpace β] [MeasurableSpace γ] {μ : Measure Ω} {R : Ω → α} {f : α → β}
    {g : α → γ} (hR : Measurable R) (hf : Measurable f) (hg : Measurable g)
    (h : IndepFun f g (μ.map R)) : IndepFun (f ∘ R) (g ∘ R) μ := by
  rw [indepFun_iff_measure_inter_preimage_eq_mul] at h ⊢
  intro s t hs ht
  have hst := h s t hs ht
  rw [Measure.map_apply hR ((hf hs).inter (hg ht)), Measure.map_apply hR (hf hs),
    Measure.map_apply hR (hg ht)] at hst
  exact hst

end Unit2B.FormalR1.Kappa.ProofsR5

#print axioms Unit2B.FormalR1.Kappa.ProofsR5.coordLaw
#print axioms Unit2B.FormalR1.Kappa.ProofsR5.indepFun_comp_of_map
