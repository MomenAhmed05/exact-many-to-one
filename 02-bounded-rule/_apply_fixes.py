"""Apply the Stage 1-5 verifier findings to the reports and the replication kit.

Every replacement is asserted: a missing anchor is reported, never silently skipped.
Also deposits the K=1 naive reference that the Stage-4 verifier flagged as prose-only.
"""

import json
import pathlib

from scipy.stats import t as student_t

HERE = pathlib.Path(__file__).resolve().parent


def patch(fn, pairs):
 p = HERE / fn
 s = p.read_text(encoding="utf-8")
 applied, missed = [], []
 for old, new in pairs:
  if old in s:
   s = s.replace(old, new, 1)
   applied.append(old[:70])
  else:
   missed.append(old[:70])
 p.write_text(s, encoding="utf-8")
 print(f"--- {fn}: applied {len(applied)} / {len(applied) + len(missed)}")
 for m in missed:
  print("    MISS:", m)


GAMMA_FIX = (
 "spends **no `gamma` at all** — its defining property is per-draw and per-lambda.",
 "has a `gamma`-free **sup step**: the admitted region is still built from the same `gamma` budget "
 "as the envelope, so its unconditional bound is `alpha + gamma` and its defining property — "
 "conditional rejection `<= alpha` — holds per-draw and per-`lambda` *within the region*.",
)
patch(
 "STAGE1-CERTIFICATE.md",
 [
  GAMMA_FIX,
  (
   "no `gamma` term; no candidate set; no covering radius",
   "the **sup step** carries no `gamma`: no candidate set, no covering radius, no Lipschitz slack "
   "(the region itself is the same `gamma`-budgeted one, so the unconditional bound is `alpha + gamma`)",
  ),
 ],
)

patch(
 "STAGE2-MAP.md",
 [
  (
   "MCSE at `alpha = 0.025` is 0.0035-0.0077, so a doubling of the FWER is a >5 sigma event.",
   "MCSE at `alpha = 0.025` is 0.0035-0.0077; at the cell carrying the boundary claim (n=36, `phi`=0.1) "
   "the doubling is 5.1 sigma, while at the largest MCSE in the table it is 3.2 sigma — the claim "
   "rests on the former, not on every cell.",
  ),
  (
   "It is never within 5% at any phi",
   "Within the tested grid (n in {36,60}, `phi` in {0,0.1,...,0.4}) no cell puts it within 5%",
  ),
  (
   "(`bootstrap_boundary_phi_max: null` in the deposit), so it does not rescue the guarantee either —",
   "(`bootstrap_boundary_phi_max: null`; that field means *no qualifying cell in this grid* — and it is "
   'also the value the writer emits on exception — so the statement is about the grid, not "any phi"; '
   "the bootstrap also carries an undisclosed extra Monte-Carlo layer, `R_boot = 200` per draw). "
   "It does not rescue the guarantee:",
  ),
  (
   "The corrections move the excess by 1-3% of its size:",
   "The corrections move the excess by 0-3% of its size (Satterthwaite route; the sample-size route "
   "reaches 7% at n=60, phi=0.4):",
  ),
  (
   "max 0.071 — that is Monte-Carlo noise, not\n  a quadrature error;",
   "max 0.071 against a Monte-Carlo SE of ~0.019\n  (600 reps), i.e. a few SE — consistent with noise "
   "rather than a quadrature error;",
  ),
  (
   "inflated the headline to 0.385; the corrected definition above is the one in the deposit.",
   "inflated the headline to 0.385; the corrected definition above is the one in the deposit.\n\n"
   "**One cell needs a caveat.** A df3's `nu80` sits at the grid endpoint (4.000), where the oracle's\n"
   'power is 0.7561 — not 0.80 — so 0.0439 of that cell\'s 0.0546 "loss" is oracle shortfall rather\n'
   "than envelope cost (its interior gap is 0.0107). The worst (0.185, B df3) and the minimum (0.017,\n"
   "A df5) are interior and sound, so the A-family range is quoted from interior cells only.",
  ),
 ],
)

patch(
 "STAGE3-INSTRUMENTS.md",
 [
  (
   "**The control-variate calibration works.** Every VR residual is an order of magnitude inside its",
   "**The control-variate calibration works, with one qualification.** At `B >= 1500` every VR residual",
  ),
  (
   "   own seed spread at `B >= 1500`, and the residual *shrinks* with `B` (S5: -0.0017 -> -0.0007 ->",
   "   clears the pre-registered criterion and the residual *shrinks* with `B` (S5: -0.0017 -> -0.0007 ->",
  ),
  (
   "**The v1 comparison is grid-sensitive.** v2's reference implementation of v1 gives 0.703 on this\n"
   "   split, not the 0.752 reported earlier, because the earlier fit searched a finer spread grid\n"
   "   (including `spread <= 25`). Both numbers are honest; the difference is the search grid, and quoting\n"
   "   only the better one would have flattered v1.",
   "**The v1 comparison cannot be attributed to a documented cause.** v2's re-implementation of v1 gives\n"
   "   0.703 on this split against the earlier 0.752. `selector.py`'s spread grid *includes* 25.0 (so it is\n"
   "   coarser, not finer) and the re-implementation used `min_cand_df` where v1 used `df`. The 0.752 figure\n"
   "   is **not recoverable from `selector2.json`** (it lives in `deposits/selector/rule.json`), so the honest\n"
   "   statement is: both numbers are reported, the difference is unexplained, and neither is allowed to\n"
   "   carry weight for the instrument.",
  ),
 ],
)

patch(
 "STAGE4-FRONTIER.md",
 [
  (
   "At B df3 the HPD shape is *narrower in width*",
   "At B df3 the HPD shape has a **larger** width at its own grid-best point",
  ),
  (
   "(9.3662 vs 8.8949 uniform is the *width* column — the HPD's best-width point is the same 8.6831 as\n"
   "   eqtail) yet it *raises* the critical value by +0.39%, the largest effect in the table. A pure width\n"
   "   objective would have picked a shape that costs power. The floor statement must therefore be made in\n"
   "   `c`, not in width — which is what the table does.",
   "(9.3662) than the equal-tail shape (8.6831) and than the shipped uniform split (8.8949), yet it *raises*\n"
   "   the critical value by +0.39% — the largest effect in the table. Width and `c` disagree in sign there,\n"
   "   so a pure width objective would have picked a shape that costs power; the floor statement must be made\n"
   "   in `c`, which is what the table does. (The A and B widths are identical by construction — region\n"
   "   geometry depends on n, K and `gamma`, not on `lambda_true` — so the two families are not independent\n"
   "   evidence about shape.)",
  ),
  (
   "in these cells (0.0185-0.0239 at 40\n   draws per cell)",
   "in these cells (0.0103-0.0238 at 40\n   draws per cell)",
  ),
  (
   "| rho | exact `c` (K=1) | naive `c` (`~ t_70`) | required factor | my `w` |",
   "| rho | exact `c` (K=1) | naive `c` (`~ t_70`) | required factor | my `w` |\n"
   "| --- | --- | --- | --- | --- |",
  ),
 ],
)

patch(
 "PRIOR-ART.md",
 [
  (
   "so the penalty is\nsmaller than the union-bound allocation by roughly 17× in miss probability.",
   "which *spends* the `γ` budget instead of\nover-covering it. The direction matters and an earlier draft had it backwards: at the A df3 design the\n"
   "union-bound allocation's miss probability is ~17× *smaller* than the tight region's (2.9e-5 versus\n"
   "5.0e-4) — it over-covers — and that surplus coverage is what makes its region wider and its critical\n"
   "value larger (union +7.24% in `c` versus tight +6.01% on the desk sweep). The tight scale takes the\n"
   "coverage at the target and spends the difference on width and power.",
  ),
 ],
)

patch(
 "STAGE5-PUBLICATION.md",
 [
  (
   "Status: **S5.1 complete** (prior art, full text), **S5.2 complete** (replication kit).",
   "Status: **S5.1 complete** — prior art gathered with **one source read in full** (Li 2009 dissertation);\n"
   "the JSTA version unretrieved and the two recent competitors abstract-level, all recorded in the\n"
   "retrieval log; **S5.2 complete** (replication kit).",
  ),
 ],
)

patch(
 "REPORT-PROVEN-BOUNDARY.md",
 [
  (
   "**no `gamma` term**, no candidate set, no Lipschitz constant, no quadrature-error assumption",
   "the **sup step carries no `gamma`** (no candidate set, no Lipschitz constant, no quadrature-error\n"
   "assumption); the region is the same `gamma`-budgeted one, so the unconditional bound stays `alpha + gamma`",
  ),
 ],
)

patch(
 "replicate/METHOD.md",
 [
  (
   "Report the scale and `cov(s)` per cell.",
   "Report the scale and `cov(s)` per cell — the original values are withheld with the other outputs.",
  ),
  (
   "at `K = 8` use `{lo, mid, hi}` on the first three axes and `{lo, hi}` on the rest, 3+2K points",
   "at `K = 8` the original used a reduced axis-wise set (`3 + 2K` points); the full `3^K` product or the "
   "`2^K` corners are equally valid candidate sets, and a replication should report which one it used",
  ),
  (
   "     decisive quantity: it must be `<= 0` up to Monte-Carlo error;",
   "     decisive quantity; the acceptance rule used here is `mean + 2*SE <= 0`;",
  ),
 ],
)

patch(
 "replicate/README.md",
 [
  (
   "- Do not read `../*.py` or `../deposits/**` before reporting your numbers.",
   "- Do not read `../*.py`, `../deposits/**`, `../PRIOR-ART.md`, `../STAGE*.md` or\n"
   "  `../REPORT-PROVEN-BOUNDARY.md` before reporting your numbers — several carry the expected answers.",
  ),
 ],
)

dpath = HERE / "replicate" / "inputs" / "designs.json"
kit = json.loads(dpath.read_text(encoding="utf-8"))
kit["constants"].pop("seed_base_formula", None)
kit["seed_reconstruction"] = {
 "note": (
  "seeds are given EXPLICITLY per cell, copied from the original deposits' per_seed "
  "records; no formula is provided or needed."
 )
}
dpath.write_text(json.dumps(kit, indent=2) + "\n", encoding="utf-8")
print("--- replicate/inputs/designs.json: seed formula removed")

q = float(student_t.ppf(0.975, 70))
ref = {
 "alpha": 0.025,
 "df": 70,
 "design": "K=1, n_i = n_0 = 36, lambda = 1 (the statistic reduces to the pooled two-sample t)",
 "naive_reference_c": q,
 "source": "Student t quantile t_{0.975,70}; the value the desk rule's plug-in law returns at K=1",
 "used_for": "the required-factor column of STAGE4-FRONTIER.md section 2",
 "required_factors": {
  str(r): round(c / q, 3)
  for r, c in (
   (0.25, 2.5746049685016748),
   (0.5, 3.4536662303841377),
   (0.75, 5.27329942174772),
  )
 },
 "note": "deposited after the Stage-4 verifier flagged that this reference existed only in prose",
}
(HERE / "deposits" / "cert" / "k1_reference.json").write_text(
 json.dumps(ref, indent=2) + "\n", encoding="utf-8"
)
print(
 f"--- deposits/cert/k1_reference.json: naive = {q:.6f} | factors {ref['required_factors']}"
)
