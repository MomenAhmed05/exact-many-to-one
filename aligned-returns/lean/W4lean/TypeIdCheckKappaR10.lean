import W4lean.Unit2BKappaStatementR1
import W4lean.Unit2BKappaProofsR10

/-
Kernel-level type-identity check for the nine S-kappa statements of
`Unit2BKappaStatementR1.lean`, against the ported proofs of `Unit2BKappaProofsR10.lean`.
Same idea as `TypeIdCheck.lean` (the thirteen FORMAL-01 theorems), extended to S-kappa.

For each statement `S` and proof `P` this file declares

    theorem typeid_S : type_of% @Unit2B.FormalR1.Kappa.S := @Unit2B.FormalR1.Kappa.ProofsR10.P

which type-checks only if `P` proves exactly the proposition `S` states, and the kernel checks
it.  Because each is a named `theorem` whose TYPE is the statement's type, `#print axioms` on it
is the axiom footprint of a proof of the project's own statement.

The `run_cmd` block additionally compares the two stored types as `Expr`s, SYNTACTICALLY (after
identifying universe parameters), which is stronger than the definitional check above.  It
reports through `logInfo` only, so a syntactic difference is visible but cannot mask the kernel
result.

Every name is fully qualified and nothing is opened, so no bare name can resolve to a `sorry`
statement.
-/

open MeasureTheory ProbabilityTheory Set

namespace KappaPortR10.TypeId

universe u

-- the headline one
theorem typeid_kappa_reduction : type_of% @Unit2B.FormalR1.Kappa.kappa_reduction.{u}
    := @Unit2B.FormalR1.Kappa.ProofsR10.kappa_reduction.{u}

theorem typeid_nu_eq : type_of% @Unit2B.FormalR1.Kappa.nu_eq
    := @Unit2B.FormalR1.Kappa.ProofsR10.nu_eq

theorem typeid_gram_eq : type_of% @Unit2B.FormalR1.Kappa.gram_eq
    := @Unit2B.FormalR1.Kappa.ProofsR10.gram_eq

theorem typeid_Sres_eq_gram_form : type_of% @Unit2B.FormalR1.Kappa.Sres_eq_gram_form
    := @Unit2B.FormalR1.Kappa.ProofsR10.Sres_eq_gram_form_verbatim

theorem typeid_Sres_eq_gram_form_step4 : type_of% @Unit2B.FormalR1.Kappa.Sres_eq_gram_form
    := @Unit2B.FormalR1.Kappa.ProofsR10.Sres_eq_gram_form

theorem typeid_Sres_pos_ae : type_of% @Unit2B.FormalR1.Kappa.Sres_pos_ae
    := @Unit2B.FormalR1.Kappa.ProofsR10.Sres_pos_ae

theorem typeid_kappa_event : type_of% @Unit2B.FormalR1.Kappa.kappa_event
    := @Unit2B.FormalR1.Kappa.ProofsR10.kappa_event

theorem typeid_kappa_numerator_law : type_of% @Unit2B.FormalR1.Kappa.kappa_numerator_law
    := @Unit2B.FormalR1.Kappa.ProofsR10.kappa_numerator_law

theorem typeid_kappa_indep : type_of% @Unit2B.FormalR1.Kappa.kappa_indep
    := @Unit2B.FormalR1.Kappa.ProofsR10.kappa_indep

theorem typeid_kappa_spectral : type_of% @Unit2B.FormalR1.Kappa.kappa_spectral.{u}
    := @Unit2B.FormalR1.Kappa.ProofsR10.kappa_spectral.{u}

end KappaPortR10.TypeId

open Lean Elab Command in
run_cmd do
  let pairs : List (Name × Name) := [
    (`Unit2B.FormalR1.Kappa.kappa_reduction, `Unit2B.FormalR1.Kappa.ProofsR10.kappa_reduction),
    (`Unit2B.FormalR1.Kappa.nu_eq, `Unit2B.FormalR1.Kappa.ProofsR10.nu_eq),
    (`Unit2B.FormalR1.Kappa.gram_eq, `Unit2B.FormalR1.Kappa.ProofsR10.gram_eq),
    (`Unit2B.FormalR1.Kappa.Sres_eq_gram_form,
      `Unit2B.FormalR1.Kappa.ProofsR10.Sres_eq_gram_form_verbatim),
    (`Unit2B.FormalR1.Kappa.Sres_eq_gram_form, `Unit2B.FormalR1.Kappa.ProofsR10.Sres_eq_gram_form),
    (`Unit2B.FormalR1.Kappa.Sres_pos_ae, `Unit2B.FormalR1.Kappa.ProofsR10.Sres_pos_ae),
    (`Unit2B.FormalR1.Kappa.kappa_event, `Unit2B.FormalR1.Kappa.ProofsR10.kappa_event),
    (`Unit2B.FormalR1.Kappa.kappa_numerator_law,
      `Unit2B.FormalR1.Kappa.ProofsR10.kappa_numerator_law),
    (`Unit2B.FormalR1.Kappa.kappa_indep, `Unit2B.FormalR1.Kappa.ProofsR10.kappa_indep),
    (`Unit2B.FormalR1.Kappa.kappa_spectral, `Unit2B.FormalR1.Kappa.ProofsR10.kappa_spectral)]
  for (a, b) in pairs do
    let ia ← getConstInfo a
    let ib ← getConstInfo b
    let tb := ib.type.instantiateLevelParams ib.levelParams (ia.levelParams.map Level.param)
    logInfo s!"SYNTAX-ID {a} vs {b}: univ-params {ia.levelParams.length}/{ib.levelParams.length}, eqv={ia.type == tb}, equal={ia.type.equal tb}"

-- axiom footprint of proofs of the PROJECT's statement types
#print axioms KappaPortR10.TypeId.typeid_kappa_reduction
#print axioms KappaPortR10.TypeId.typeid_nu_eq
#print axioms KappaPortR10.TypeId.typeid_gram_eq
#print axioms KappaPortR10.TypeId.typeid_Sres_eq_gram_form
#print axioms KappaPortR10.TypeId.typeid_Sres_eq_gram_form_step4
#print axioms KappaPortR10.TypeId.typeid_Sres_pos_ae
#print axioms KappaPortR10.TypeId.typeid_kappa_event
#print axioms KappaPortR10.TypeId.typeid_kappa_numerator_law
#print axioms KappaPortR10.TypeId.typeid_kappa_indep
#print axioms KappaPortR10.TypeId.typeid_kappa_spectral

-- and of the ported declarations themselves, by their project names
#print axioms Unit2B.FormalR1.Kappa.ProofsR10.kappa_reduction
#print axioms Unit2B.FormalR1.Kappa.ProofsR10.Sres_pos_ae
#print axioms Unit2B.FormalR1.Kappa.ProofsR10.kappa_event
#print axioms Unit2B.FormalR1.Kappa.ProofsR10.kappa_numerator_law
#print axioms Unit2B.FormalR1.Kappa.ProofsR10.kappa_indep
#print axioms Unit2B.FormalR1.Kappa.ProofsR10.kappa_spectral

-- control: the statement file's own declaration is, by design, still `sorry`
#print axioms Unit2B.FormalR1.Kappa.kappa_reduction
