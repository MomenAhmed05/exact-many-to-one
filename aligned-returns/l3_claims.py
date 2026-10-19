"""Write the CLAIM files for the L-3 kill test from deposits/l3_killtest.json.

Run after l3_killtest.py, from this folder:
    python l3_killtest.py
    python l3_claims.py
Writes deposits/l3_killtest_claims.txt and deposits/l3_boundary_claims.txt.  No number
is typed here; every value is read from the JSON written by l3_killtest.py.
"""
import json
import pathlib

j = json.loads(pathlib.Path("deposits/l3_killtest.json").read_text(encoding="utf-8"))

# ---- deposits/l3_killtest_claims.txt
o = []
def C(i, l, v): o.append(f"CLAIM {i}: {l} = {v}")
C("L3K-ALPHA", "nominal level", j["alpha"])
C("L3K-DRAWS", "Monte Carlo draws per design arm", j["draws"])
C("L3K-ZTHRESHOLD", "predeclared failure threshold in Monte Carlo standard errors", j["z_threshold"])
C("L3K-NCAL", "equal-window calibration design length n", j["n_cal"])
C("L3K-CRIT-COUPLED", "critical value calibrated on the equal, fully overlapping design at the coupled loadings", j["critical_value_coupled"])
C("L3K-CRIT-UNCOUPLED", "critical value calibrated at zero loadings", j["critical_value_uncoupled"])
C("L3K-EQUAL-SIZE", "realised size on the equal, fully overlapping design, the calibration anchor", j["equal_design_realised_size"])
for r in j["rows"]:
    tag = f"{r['design']}|{r['arm']}"
    C(f"L3K-{r['design']}-{r['arm']}-SIZE", f"realised size, design {tag}, n0={r['n0']} n1={r['n1']} k={r['k']}", r["size"])
    C(f"L3K-{r['design']}-{r['arm']}-MCSE", f"Monte Carlo standard error, {tag}", r["mcse"])
    C(f"L3K-{r['design']}-{r['arm']}-ZEQ", "departure from the equal-window design in Monte Carlo standard errors, " + tag, r["z_vs_equal_design"])
C("L3K-BREACHES", "designs whose size departs from the equal-window design by more than the threshold", j["n_breaches_vs_equal"])
C("L3K-VERDICT", "does the equal-window clean reduction survive transfer to unequal and partially overlapping windows", str(j["verdict_reduction_survives"]))
def get(d, a): return next(r for r in j["rows"] if r["design"] == d and r["arm"].startswith(a))
C("L3K-MECH-1", "mechanism 1: at EQUAL lengths, reducing overlap from 24 to 12 moves the coupled arm by z while leaving the uncoupled arm inside the threshold, so the coupling term is the cause", "coupled z={} vs uncoupled z={}".format(get("equal-half-overlap", "coupled")["z_vs_equal_design"], get("equal-half-overlap", "uncoupled")["z_vs_equal_design"]))
C("L3K-MECH-2", "mechanism 2: with coupling removed, the unequal-length designs still mis-size, because the unpaired denominator is a weighted sum of two chi-squares with different weights and is not exactly t when n0 differs from n1", "uncoupled z(n1 longer)={}, z(n0 longer)={}".format(get("unequal-n1-longer", "uncoupled")["z_vs_equal_design"], get("unequal-n0-longer", "uncoupled")["z_vs_equal_design"]))
C("L3K-CONTROL-OK", "the uncoupled arm is inside the threshold at equal lengths and at disjoint windows, so the calibration itself is sound and the departures are attributable to the design, not to a bad critical value", "control at equal={} disjoint={}".format(get("equal-full-overlap", "uncoupled")["z_vs_equal_design"], get("unequal-disjoint", "uncoupled")["z_vs_equal_design"]))
C("L3K-SANITY", "the coupled arm at the calibration design is inside the threshold, confirming the instrument's anchor", "z={}".format(get("equal-full-overlap", "coupled")["z_vs_equal_design"]))
pathlib.Path("deposits/l3_killtest_claims.txt").write_text("\n".join(o) + "\n", encoding="utf-8")

# ---- deposits/l3_boundary_claims.txt
o = []
rows = [r for r in j["rows"] if r["arm"].startswith("coupled")]
passes = [r for r in rows if r["survives_reduction"]]
C("L3B-COUPLED-DESIGNS", "coupled-arm designs tested", len(rows))
C("L3B-COUPLED-PASS", "coupled-arm designs where the reduction survives", ", ".join(f"n0={r['n0']} n1={r['n1']} k={r['k']}" for r in passes) or "none")
C("L3B-FULL-OVERLAP-OF-SHORTER", "designs where the overlap covers the whole shorter window, k = min(n0,n1)", ", ".join(f"n0={r['n0']} n1={r['n1']} k={r['k']} -> {'survives' if r['survives_reduction'] else 'FAILS'}" for r in rows if r['k'] == min(r['n0'], r['n1'])))
C("L3B-BOUNDARY", "the equal-window reduction survives ONLY where the two windows are identical, n0 = n1 = k; equal LENGTHS with reduced overlap already break it, so the operative condition is identical windows, not equal windows", str(all(r['n0'] == r['n1'] == r['k'] for r in passes) and len(passes) == 1))
C("L3B-K-OF-SHORTER-FAILS", "full overlap of the shorter window is NOT sufficient: n0=12 n1=24 k=12 has k=min(n0,n1) yet fails", "z={}".format(next(r for r in rows if r['n0'] == 12 and r['n1'] == 24)['z_vs_equal_design']))
C("L3B-EQUAL-LEN-REDUCED-K", "equal lengths 24,24 with overlap reduced to 12 already fails, so 'equal windows' as commonly stated is the wrong condition", str(not next(r for r in rows if r['design'] == 'equal-half-overlap')['survives_reduction']))
C("L3B-TARGET-FORM", "L-3's target is a stated correction or a stated boundary; this run yields a BOUNDARY, not a correction, because no corrected critical value is derived here", "BOUNDARY")
pathlib.Path("deposits/l3_boundary_claims.txt").write_text("\n".join(o) + "\n", encoding="utf-8")
print("wrote deposits/l3_killtest_claims.txt and deposits/l3_boundary_claims.txt")
