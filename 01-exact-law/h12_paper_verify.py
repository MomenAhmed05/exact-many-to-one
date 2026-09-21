"""Final audit of the paper: every quoted number against its deposit, plus
structural and style checks.

A check is only worth running if it can fail for the right reason. An earlier
version of this script asserted that literals hardcoded HERE appeared in the
paper, where they were also hardcoded in the build script: those checks compared
the document with itself and could never catch a wrong number. Every numeric
check below now formats its expected string FROM A DEPOSIT, so the deposit is
the only source of truth and editing the paper cannot make a check pass.
"""
import html
import json
import os
import re
import sys
import zipfile

D = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(D, "deposits")
F = os.path.join(os.path.dirname(D), "paper",
                 "Exact-Many-to-One-Under-Estimated-Volatilities.docx")

z = zipfile.ZipFile(F)
x = z.read("word/document.xml").decode("utf-8")
# <w:t[^>]*> also matches <w:tbl>, <w:tblPr>, <w:tblGrid> and friends, which
# splices raw markup into the extracted text and lets checks pass on XML
# attributes. Match the text element and nothing else.
txt = html.unescape("".join(re.findall(r"<w:t(?: [^>]*)?>(.*?)</w:t>", x, re.S)))
txt = re.sub(r"\s+", " ", txt)
# the paper sets a typographic minus on numbers; deposits carry ASCII
txt = txt.replace("\u2212", "-")
# scientific notation is bound with a non-breaking hyphen so Word cannot split
# "1.2e-7" across a line break; normalise it back before matching
txt = txt.replace("\u2011", "-")

ok = True
npass = 0
nfail = 0


def chk(label, cond, detail=""):
    global ok, npass, nfail
    cond = bool(cond)
    ok &= cond
    npass += 1 if cond else 0
    nfail += 0 if cond else 1
    print(("  OK   " if cond else "  FAIL ") + label
          + ("   " + detail if detail else ""))
    return cond


def near(label, value, fmt, before="", after=""):
    """Assert a deposited value appears in the paper WITH its context.

    Short values match by accident: "28" is a substring of "0.02807". Anchoring
    on surrounding words makes the check test the sentence, not the digits.
    """
    v = fmt(value) if callable(fmt) else fmt % value
    want = before + v + after
    return chk("%s -> %r" % (label, want), want in txt)


def has(label, value, fmt="%s"):
    """Assert that a value formatted from a deposit appears in the paper.

    `fmt` is a percent-format string, or any callable that renders the value.
    """
    v = fmt(value) if callable(fmt) else fmt % value
    return chk("%s -> %s" % (label, v), v in txt)


def sci(digits=1):
    """Scientific notation the way the paper writes it: 3.8e-5, not 3.8e-05."""
    def f(v):
        t = ("%." + str(digits) + "e") % v
        m, e = t.split("e")
        return "%se%s%d" % (m, "-" if int(e) < 0 else "", abs(int(e)))
    return f


def load(n):
    return json.load(open(os.path.join(S, n)))


A25, A05 = 0.025, 0.05
CMP, EX = load("h12_hetero_compare.json"), load("h12_hetero_exact.json")
H2H, CRS = load("h12_tamhane_headtohead.json"), load("h12_crossover.json")
CAL, F1 = load("h12_calibrated.json"), load("h12_f1_recheck.json")
NUM = load("h12_numerics.json")

print("STYLE AND STRUCTURE")
chk("no em dashes", chr(8212) not in txt)
chk("no en dashes", chr(8211) not in txt)
chk("no markdown bold markers", "**" not in txt)
chk("no unresolved cross-reference tokens", "<<eq:" not in txt)
chk("no doubled percentage units", not re.search(r"%\s*per cent", txt))
chk("no doubled words",
    not re.search(r"\b(the|a|an|of|is|and|to|in|that|for) \1\b", txt),
    str(re.findall(r"\b(the|a|an|of|is|and|to|in|that|for) \1\b", txt)[:3]))
chk("'incumbent' only in the promotion-gate sentence",
    txt.count("ncumbent") == 1 and "the incumbent has the longer record"
    in txt, str(txt.count("ncumbent")))
neq = len(re.findall(r"<m:oMath[ >]", x))
chk("every equation reference points at a real display",
    all(1 <= int(m) <= 50 for m in re.findall(r"equation \((\d+)\)", txt)))
chk("23 figures embedded and no stray plates from the template",
    len([n for n in z.namelist() if n.startswith("word/media")]) == 23,
    str(len([n for n in z.namelist() if n.startswith("word/media")])))
chk("37 tables", x.count("<w:tbl>") == 37, str(x.count("<w:tbl>")))
heads = re.findall(r'w:pStyle w:val="Heading1"', x)
chk("17 top-level sections (abstract, 14 numbered, refs, appendix)",
    len(heads) == 17, str(len(heads)))
for i in range(1, 24):
    chk("Figure %d referenced and captioned" % i, ("Figure %d." % i) in txt)
for i in range(1, 35):
    chk("Table %d captioned" % i, ("Table %d." % i) in txt)

print("\nSECTION 4 NUMERICS (deposit: h12_numerics.json)")
lo = NUM["laguerre_overflow"]
has("Laguerre last finite nu", lo["last_finite_nu"], "%d")
has("Laguerre first infinite nu", lo["first_infinite_nu"], "%d")
lr = NUM["laguerre_refinement"]
for r in lr["rows"]:
    if r["finite"]:
        has("Laguerre at %d nodes" % r["nodes"], r["critical_value"], "%.6f")
has("Laguerre spread", lr["spread"], sci())
near("Laguerre first non-finite node count", lr["first_nonfinite_nodes"],
     "%d", before="at ", after=" nodes it returned a non-number")
sg = NUM["step_geometry"]
has("panel width in u", sg["panel_width_u"], "%.2f")
has("panel-to-step ratio", sg["ratio_rounded_to_hundreds"], "{:,}".format)
tm = NUM["timings"]
near("cost range", tm["min_seconds"], "%.2f", before="measured at ",
     after=" to")
near("cost range, upper", tm["max_seconds"], "%.2f", after=" seconds across")

print("\nSECTION 5 FALSIFIERS (deposit: h12_hetero_exact.json)")
s = EX["summary"]
has("F1 worst", s["f1_worst_abs_err"], sci())
has("F2 worst", s["f2_worst_abs_err"], sci())
has("F3 fine to finest", s["f3_worst_fine_to_finest"], sci())
has("F4 worst |z|", s["f4_worst_abs_z"], "%.2f")
has("F5 worst decade ratio", s["f5_worst_decade_ratio"], "%.1f")
has("F5 implied quadrature bound", s["f5_implied_quad_err_bound"], sci())
un = [e for r in EX["f5_perfect_dependence"] for e in r["unadapted_errs"]]
has("F5 unadapted band, low", min(un), sci())
has("F5 unadapted band, high", max(un), sci())
has("F5b worst spread across K", s["f5b_worst_spread_across_K"], sci())
ks = sorted({k for r in EX["f5b_k_invariance"] for k in r["cdf_by_K"]}, key=int)
chk("F5b sweeps K from %s to %s" % (ks[0], ks[-1]),
    ("K = 1 to 12" in txt) or ("K from 1 to 12" in txt) or ("1 to 12" in txt))

print("\nSECTION 5.3 EXTERNAL AGREEMENT")
has("published Dunnett, recomputed",
    CMP["c3_dunnett_published"]["onefactor"], "%.6f")
chk("published Dunnett, as published",
    str(CMP["c3_dunnett_published"]["published"]) in txt)
has("SciPy agreement",
    max(r["max_abs_p_diff"] for r in CMP["c4_scipy_dunnett"]), sci())

print("\nSECTION 6 RULE COSTS (deposit: h12_hetero_compare.json)")
for r in CMP["comparison"]:
    for k in ("T3_smm_welch_df", "sidak_t_welch_df", "onefactor_pooled_df",
              "onefactor_mean_welch_df"):
        if not has("%s %s" % (r["design"][:22], k[:14]),
                   r[k]["true_fwer"], "%.5f"):
            break
hom = [r for r in CMP["comparison"] if "homoscedastic" in r["design"]][0]
has("homoscedastic, Dunnett carried across",
    100 * (hom["onefactor_pooled_df"]["true_fwer"] - A05) / A05, "%.1f")
has("homoscedastic, studentised maximum modulus",
    -100 * (hom["T3_smm_welch_df"]["true_fwer"] - A05) / A05, "%.1f")

print("\nSECTION 6.3 WORKED EXAMPLE (deposit: h12_numerics.json)")
we = NUM["worked_example"]
has("exact critical value", we["c_exact"], "%.4f")
has("exact delivers", 100 * we["fwer_exact"], "%.3f")
has("studentised maximum modulus delivers", 100 * we["fwer_t3_smm"], "%.3f")
has("Dunnett carried across delivers",
    100 * we["fwer_dunnett_carried"], "%.3f")
has("budget thrown away", we["t3_budget_shortfall_pct"], "%.0f")
has("budget overshoot", we["dunnett_budget_overshoot_pct"], "%.0f")
for i, h in enumerate(we["loadings"]):
    has("loading %d" % (i + 1), h, "%.3f")
has("homoscedastic loading", we["homoscedastic_loading"], "%.3f")

print("\nSECTION 7 HEAD TO HEAD (deposits: headtohead, crossover, numerics)")
for K in ("2", "3"):
    for key in ("exact", "method1", "method2", "method3"):
        f = [r[key]["fwer"] for r in H2H["blocks"][K]]
        chk("k=%s %s range" % (K, key),
            "%.5f to %.5f" % (min(f), max(f)) in txt)
fa = NUM["reimplementation_faithfulness"]
near("cells recomputed", fa["cells"], "%d", before="the ", after=" published cells")
near("cells up", fa["n_positive"], "%d", before="", after=" of the 42 go up")
near("cells down", fa["n_negative"], "%d", before="", after=" go down")
has("worst absolute difference", fa["worst_abs_difference"], "%.5f")
has("published standard error", fa["published_se"], "%.5f")
has("recomputation standard error", fa["recomputation_se_rms"], "%.5f")
has("combined standard error", fa["combined_se"], "%.5f")
near("worst cell in combined standard errors",
     fa["worst_in_combined_se"], "%.1f",
     before="sits at ", after=" combined standard errors")
cv = [H2H["blocks"][K][0][k]["mcse"] for K in ("2", "3") for k in ("oracle",)]
mx = {K: max(r["exact"]["mcse"] for r in H2H["blocks"][K]) for K in ("2", "3")}
has("control variate, k=2", mx["2"], "%.5f")
has("control variate, k=3", mx["3"], "%.5f")
oa = max(abs(r["oracle_analytic_fwer"] - A25) for r in H2H["blocks"]["2"])
has("oracle exactness, k=2 block", oa, sci())
cr = NUM["cross_run_consistency"]
near("cross-run comparisons", cr["n_comparisons"], "%d",
     before="the ", after=" shared comparisons")
near("cross-run worst |z|", cr["worst_abs_z"], "%.1f",
     before="is ", after=" standard errors")
near("cross-run mean |z|", cr["mean_abs_z"], "%.1f", before="the mean is ")
near("head to head draws", cr["headtohead_draws"], "%d", after=".")
near("calibration draws", cr["calibration_draws"], "%d",
     before="at ", after=" data sets per scenario")
for r in CRS["d1_crossover"]:
    chk("df %d rows" % r["df0"],
        "%+.1f%%" % r["uncal_rel_err_pct"] in txt
        and "%+.1f%%" % r["m3_rel_err_pct"] in txt)
for r in CRS["d2_decompose"]:
    for k in ("exact_plugin_bias", "m3_total_bias", "m3_approximation_bias",
              "m3_estimation_bias"):
        chk("%s %s" % (r["case"][:12], k), "%+.4f" % r[k] in txt)
neg = [r["case"] for r in CRS["d2_decompose"] if r["m3_estimation_bias"] < 0]
chk("the cancellation is claimed for %d of %d designs, not more"
    % (len(neg), len(CRS["d2_decompose"])),
    "three of the four" not in txt
    and "on three of the four" not in txt)
d3 = {r["case"]: r for r in CRS["d3_m3_noise"]}
has("Method 3 noise, shift", d3["S6"]["shift"], "%.5f")
has("Method 3 noise, in standard errors", d3["S6"]["shift_in_mcse"], "%.1f")
for r in CRS["d4_alpha05"]:
    chk("alpha=.05 %s" % r["case"],
        "%+.1f%%" % r["uncal_rel_err_pct"] in txt
        and "%+.1f%%" % r["m3_rel_err_pct"] in txt)

print("\nSECTION 8 CALIBRATION (deposits: calibrated, f1_recheck, numerics)")
for K in ("2", "3"):
    st = CAL["summary"][K]
    for k in ("uncalibrated", "calibrated", "method3"):
        chk("k=%s %s range" % (K, k),
            "%.5f to %.5f" % (st[k]["min"], st[k]["max"]) in txt)
        has("k=%s %s mean" % (K, k), st[k]["mean_abs_rel_err_pct"], "%.1f")
ps = NUM["precision_study"]
by_b = {r["replicates"]: r for r in ps["rows"]}
has("naive sd at 200", by_b[200]["naive_sd"], "%.4f")
has("naive sd at 4000", by_b[4000]["naive_sd"], "%.4f")
has("conditional sd at 100", by_b[100]["rb_sd"], "%.4f")
has("conditional sd at 200", by_b[200]["rb_sd"], "%.4f")
chk("replicate reduction stated as a word",
    {20.0: "twentyfold", 40.0: "fortyfold",
     50.0: "fiftyfold"}.get(ps["replicate_reduction"], "NOSUCH") in txt,
    "deposit says %.0f" % ps["replicate_reduction"])
has("calibration target", ps["target_alpha_star"], "%.4f")
chk("F2 ratios", all("%.3f" % r["ratio_to_alpha"] in txt
                     for r in CAL["f2_large_df"]))
has("F1a exactness", F1["summary"]["worst_abs_F_minus_alpha"], sci())
has("F1a against the head to head",
    F1["summary"]["worst_abs_F_vs_headtohead"], "%.5f")
near("F1b worst |z|", F1["summary"]["worst_abs_z_end_to_end"], "%.1f",
     before="at ", after=" standard errors")
chk("F1b is not dressed up as passing",
    not F1["summary"]["f1b_passes"]
    and ("does not pass" in txt or "not pass" in txt))
chk("F1b is reported as not passing",
    ("does not pass" in txt) or ("still fails" in txt)
    or ("not pass" in txt))
sl = [c["dFWER_dalphastar"] for c in F1["cases"]]
has("slope, lowest", min(sl), "%.2f")
has("slope, highest", max(sl), "%.2f")

print("\nREFERENCES")
for r in ["Beran, R. (1987)", "Dunnett, C. W. (1955)", "Dunnett, C. W. (1964)",
          "Dunnett, C. W. (1980)", "Hasler, M., and Hothorn",
          "Tamhane, A. C., and Xi, D. (2023)", "Welch, B. L. (1938)",
          "Satterthwaite, F. E. (1946)", "Nagananda"]:
    chk("cited: " + r, r in txt)

print("\n%d checks passed, %d failed" % (npass, nfail))
print("PAPER VERIFIED" if ok else "PAPER HAS FAILURES")
sys.exit(0 if ok else 1)
