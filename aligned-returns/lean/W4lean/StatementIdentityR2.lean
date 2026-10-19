import W4lean.Unit2BFormalR1R2

/-
Kernel-level statement identity between the statement file of record, `Unit2BFormalR1.lean`
(namespace `Unit2B.FormalR1`), and its self-contained successor `Unit2BFormalR1R2.lean`
(namespace `Unit2B.FormalR1.R2`).

Three checks per theorem, all by fully qualified name; nothing is opened.

  (1) EQUALITY OF TYPES, closed by `rfl`:
        example : type_of% @Unit2B.FormalR1.X = type_of% @Unit2B.FormalR1.R2.X := rfl
      The six universe-polymorphic theorems are instantiated at one shared universe `u` on both
      sides, so the two sides are compared at the same universe rather than at two independent
      universe parameters.
  (2) INHABITATION, the form already used by `TypeIdCheck.lean`:
        example : type_of% @Unit2B.FormalR1.X := @Unit2B.FormalR1.R2.X
  (3) SYNTACTIC `Expr` equality of the two stored types (after identifying universe
      parameters), reported by `logInfo`, which is stronger than definitional equality.
      `==` is `Expr.eqv` (alpha-equivalence: binder names ignored); `Expr.equal` additionally
      compares binder names, so each theorem's differing binder-name pairs are listed too.
-/

open MeasureTheory ProbabilityTheory Set

namespace Unit2B.FormalR1.StatementIdentityR2

universe u

-- (1) equality of types, by `rfl`
example : type_of% @Unit2B.FormalR1.Phi_antitone.{u}
    = type_of% @Unit2B.FormalR1.R2.Phi_antitone.{u} := rfl
example : type_of% @Unit2B.FormalR1.ae_nonneg_of_gammaLaw.{u}
    = type_of% @Unit2B.FormalR1.R2.ae_nonneg_of_gammaLaw.{u} := rfl
example : type_of% @Unit2B.FormalR1.charQ_eq_det
    = type_of% @Unit2B.FormalR1.R2.charQ_eq_det := rfl
example : type_of% @Unit2B.FormalR1.discG_eq
    = type_of% @Unit2B.FormalR1.R2.discG_eq := rfl
example : type_of% @Unit2B.FormalR1.charQ_factor
    = type_of% @Unit2B.FormalR1.R2.charQ_factor := rfl
example : type_of% @Unit2B.FormalR1.mu_sum_prod
    = type_of% @Unit2B.FormalR1.R2.mu_sum_prod := rfl
example : type_of% @Unit2B.FormalR1.box_pos
    = type_of% @Unit2B.FormalR1.R2.box_pos := rfl
example : type_of% @Unit2B.FormalR1.r1_lower
    = type_of% @Unit2B.FormalR1.R2.r1_lower := rfl
example : type_of% @Unit2B.FormalR1.r2_lower
    = type_of% @Unit2B.FormalR1.R2.r2_lower := rfl
example : type_of% @Unit2B.FormalR1.cstar_bracket.{u}
    = type_of% @Unit2B.FormalR1.R2.cstar_bracket.{u} := rfl
example : type_of% @Unit2B.FormalR1.cstarBox_eq_sSup.{u}
    = type_of% @Unit2B.FormalR1.R2.cstarBox_eq_sSup.{u} := rfl
example : type_of% @Unit2B.FormalR1.cover_sound.{u}
    = type_of% @Unit2B.FormalR1.R2.cover_sound.{u} := rfl
example : type_of% @Unit2B.FormalR1.comp_thm_2B.{u}
    = type_of% @Unit2B.FormalR1.R2.comp_thm_2B.{u} := rfl

-- (2) inhabitation: the successor's declaration proves the statement file's type
example : type_of% @Unit2B.FormalR1.Phi_antitone.{u}
    := @Unit2B.FormalR1.R2.Phi_antitone.{u}
example : type_of% @Unit2B.FormalR1.ae_nonneg_of_gammaLaw.{u}
    := @Unit2B.FormalR1.R2.ae_nonneg_of_gammaLaw.{u}
example : type_of% @Unit2B.FormalR1.charQ_eq_det
    := @Unit2B.FormalR1.R2.charQ_eq_det
example : type_of% @Unit2B.FormalR1.discG_eq
    := @Unit2B.FormalR1.R2.discG_eq
example : type_of% @Unit2B.FormalR1.charQ_factor
    := @Unit2B.FormalR1.R2.charQ_factor
example : type_of% @Unit2B.FormalR1.mu_sum_prod
    := @Unit2B.FormalR1.R2.mu_sum_prod
example : type_of% @Unit2B.FormalR1.box_pos
    := @Unit2B.FormalR1.R2.box_pos
example : type_of% @Unit2B.FormalR1.r1_lower
    := @Unit2B.FormalR1.R2.r1_lower
example : type_of% @Unit2B.FormalR1.r2_lower
    := @Unit2B.FormalR1.R2.r2_lower
example : type_of% @Unit2B.FormalR1.cstar_bracket.{u}
    := @Unit2B.FormalR1.R2.cstar_bracket.{u}
example : type_of% @Unit2B.FormalR1.cstarBox_eq_sSup.{u}
    := @Unit2B.FormalR1.R2.cstarBox_eq_sSup.{u}
example : type_of% @Unit2B.FormalR1.cover_sound.{u}
    := @Unit2B.FormalR1.R2.cover_sound.{u}
example : type_of% @Unit2B.FormalR1.comp_thm_2B.{u}
    := @Unit2B.FormalR1.R2.comp_thm_2B.{u}

-- every binder (name, binder-info) of an `Expr`, in pre-order; used by (3) to show WHICH
-- binders differ when `Expr.equal` (which compares binder names) is false but `==` is true
open Lean in
partial def binders : Expr → Array (Name × String)
  | .forallE n t b i => #[(n, bi i)] ++ binders t ++ binders b
  | .lam n t b i => #[(n, bi i)] ++ binders t ++ binders b
  | .letE n t v b _ => #[(n, "let")] ++ binders t ++ binders v ++ binders b
  | .app f a => binders f ++ binders a
  | .mdata _ e => binders e
  | .proj _ _ e => binders e
  | _ => #[]
where
  bi : BinderInfo → String
    | .default => "default"
    | .implicit => "implicit"
    | .strictImplicit => "strictImplicit"
    | .instImplicit => "instImplicit"

end Unit2B.FormalR1.StatementIdentityR2

-- (3) syntactic identity of the stored types
open Lean Elab Command in
run_cmd do
  let names : List Name := [`Phi_antitone, `ae_nonneg_of_gammaLaw, `charQ_eq_det, `discG_eq,
    `charQ_factor, `mu_sum_prod, `box_pos, `r1_lower, `r2_lower, `cstar_bracket,
    `cstarBox_eq_sSup, `cover_sound, `comp_thm_2B]
  for n in names do
    let a := `Unit2B.FormalR1 ++ n
    let b := `Unit2B.FormalR1.R2 ++ n
    let ia ← getConstInfo a
    let ib ← getConstInfo b
    let tb := ib.type.instantiateLevelParams ib.levelParams (ia.levelParams.map Level.param)
    logInfo s!"SYNTAX-ID {a} vs {b}: univ-params {ia.levelParams.length}/{ib.levelParams.length}, eqv={ia.type == tb}, equal={ia.type.equal tb}"
    let ba := Unit2B.FormalR1.StatementIdentityR2.binders ia.type
    let bb := Unit2B.FormalR1.StatementIdentityR2.binders tb
    let diffs := (ba.zip bb).filter (fun (x, y) => x != y)
    logInfo s!"BINDERS {n}: {ba.size}/{bb.size} binders, differing pairs = {diffs.size}: {diffs.toList}"
