import W4lean.KappaIndDefsR10
import W4lean.Unit2BKappaStatementR1

/-
Port, step 0: are the separate Lean project's model definitions DEFINITIONALLY EQUAL, in the
kernel, to this project's?

`W4lean/KappaIndDefsR10.lean` is a copy of that project's `KappaInd/Defs.lean`, byte-identical
except for its header comment (the original's sha256 is
56fa63a17384d77b5bb8e8a77f8084c823ba232ed67bd9f949326b3d80b77f44), compiled as its own module
with its own narrow Mathlib imports.  This file imports it next to the project's
`Unit2BFormalR1` / `Unit2BKappaStatementR1` (which import all of Mathlib) and asks, for each
definition, for a proof by `rfl`.  Each is a `theorem`, so it is checked by the kernel and stays
in the environment.  Each is elaborated independently: a failure of one does not hide the others.

Namespaces.  The project splits the model: `nu` .. `Phi` live in `Unit2B.FormalR1`, and
`nObs` .. `zc` live in `Unit2B.FormalR1.Kappa`.  The separate project has all of them in `KappaInd`.
Every name below is fully qualified; nothing is opened, so no bare name can resolve to a
`sorry` statement.

The two `negctl_*` declarations at the end are negative controls: they show that the same
`rfl` check REJECTS pairs that are not definitionally equal, so the positive results are not
vacuous.
-/

open MeasureTheory ProbabilityTheory Set

namespace KappaPortR10.DefEq

universe u

-- constants
theorem defeq_nu : KappaInd.nu = Unit2B.FormalR1.nu := rfl
theorem defeq_alpha : KappaInd.alpha = Unit2B.FormalR1.alpha := rfl
theorem defeq_Lc : KappaInd.Lc = Unit2B.FormalR1.Lc := rfl
theorem defeq_Uc : KappaInd.Uc = Unit2B.FormalR1.Uc := rfl
theorem defeq_tol : KappaInd.tol = Unit2B.FormalR1.tol := rfl

-- parameter space
theorem defeq_Param : KappaInd.Param = Unit2B.FormalR1.Param := rfl
theorem defeq_Box : KappaInd.Box = Unit2B.FormalR1.Box := rfl
theorem defeq_corner : KappaInd.corner = Unit2B.FormalR1.corner := rfl

-- Gram / eigenvalue layer
theorem defeq_aF : KappaInd.aF = Unit2B.FormalR1.aF := rfl
theorem defeq_g00 : KappaInd.g00 = Unit2B.FormalR1.g00 := rfl
theorem defeq_g11 : KappaInd.g11 = Unit2B.FormalR1.g11 := rfl
theorem defeq_g01 : KappaInd.g01 = Unit2B.FormalR1.g01 := rfl
theorem defeq_Gmat : KappaInd.Gmat = Unit2B.FormalR1.Gmat := rfl
theorem defeq_trG : KappaInd.trG = Unit2B.FormalR1.trG := rfl
theorem defeq_detG : KappaInd.detG = Unit2B.FormalR1.detG := rfl
theorem defeq_discG : KappaInd.discG = Unit2B.FormalR1.discG := rfl
theorem defeq_mu1 : KappaInd.mu1 = Unit2B.FormalR1.mu1 := rfl
theorem defeq_mu2 : KappaInd.mu2 = Unit2B.FormalR1.mu2 := rfl
theorem defeq_r1 : KappaInd.r1 = Unit2B.FormalR1.r1 := rfl
theorem defeq_r2 : KappaInd.r2 = Unit2B.FormalR1.r2 := rfl
theorem defeq_charQ : KappaInd.charQ = Unit2B.FormalR1.charQ := rfl

-- Phi is universe-polymorphic and takes an implicit type and an instance
theorem defeq_Phi {Ω : Type u} [inst : MeasurableSpace Ω] :
    @KappaInd.Phi Ω inst = @Unit2B.FormalR1.Phi Ω inst := rfl

-- the model proper (project namespace Unit2B.FormalR1.Kappa)
theorem defeq_nObs : KappaInd.nObs = Unit2B.FormalR1.Kappa.nObs := rfl
theorem defeq_PanelIdx : KappaInd.PanelIdx = Unit2B.FormalR1.Kappa.PanelIdx := rfl
theorem defeq_Panel : KappaInd.Panel = Unit2B.FormalR1.Kappa.Panel := rfl
theorem defeq_panelLaw : KappaInd.panelLaw = Unit2B.FormalR1.Kappa.panelLaw := rfl
theorem defeq_ret0 : KappaInd.ret0 = Unit2B.FormalR1.Kappa.ret0 := rfl
theorem defeq_ret1 : KappaInd.ret1 = Unit2B.FormalR1.Kappa.ret1 := rfl
theorem defeq_sMean : KappaInd.sMean = Unit2B.FormalR1.Kappa.sMean := rfl
theorem defeq_sVar : KappaInd.sVar = Unit2B.FormalR1.Kappa.sVar := rfl
theorem defeq_Tstat : KappaInd.Tstat = Unit2B.FormalR1.Kappa.Tstat := rfl
theorem defeq_RawRate : KappaInd.RawRate = Unit2B.FormalR1.Kappa.RawRate := rfl
theorem defeq_Nstd : KappaInd.Nstd = Unit2B.FormalR1.Kappa.Nstd := rfl
theorem defeq_Sres : KappaInd.Sres = Unit2B.FormalR1.Kappa.Sres := rfl
theorem defeq_Mmat : KappaInd.Mmat = Unit2B.FormalR1.Kappa.Mmat := rfl
theorem defeq_zc : KappaInd.zc = Unit2B.FormalR1.Kappa.zc := rfl

-- negative controls: the same check must REJECT non-equal pairs
theorem negctl_nat : True := by
  fail_if_success
    have _h : KappaInd.nObs = 35 := rfl
  trivial

theorem negctl_cross : True := by
  fail_if_success
    have _h : KappaInd.g00 = Unit2B.FormalR1.g11 := rfl
  trivial

end KappaPortR10.DefEq

#print axioms KappaPortR10.DefEq.defeq_Phi
#print axioms KappaPortR10.DefEq.defeq_RawRate
