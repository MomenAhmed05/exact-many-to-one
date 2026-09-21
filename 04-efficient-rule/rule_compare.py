"""Compare the two certified general-K rules on the quantity a desk sees.

`certify_joint.json` carries a field `level_gain_vs_percand_x`, the ratio of the
joint rule's certified level to the per-candidate rule's.  That ratio is not a
comparison of rules.  The two levels feed DIFFERENT threshold functions: the
per-candidate level is a level on one comparison, and the joint level is a level
on the maximum over K.  At the same nominal level the joint constant is the
larger of the two, so a joint rule running at a higher level can still ask for a
higher threshold.

The comparable quantities are the deployed threshold and the realised error
rate, and both are computed here.  Nothing in this module is a certificate; it
reads the deposited certified levels and prices the rules they define.

Deposit: deposits/rule_compare.json
"""

from __future__ import annotations

import time

import numpy as np

import certify_joint as J
import fcore as F
import k1lab as L

# the estimated ratio vectors to price the thresholds at, spanning the domain
LAM_POINTS = ([1.0, 1.0], [0.25, 4.0], [4.0, 0.25], [0.04, 25.0])


def percand_threshold(n, level, lam_hat):
    """The deployed per-candidate threshold: the one-candidate constant at each
    candidate's own estimated ratio, at the per-candidate level."""
    tab = L.PluginTable(np.array([n[0], n[1]]), level, m=41, npanel=12,
                        order=8)
    return [float(np.max(tab(np.array([t])))) for t in lam_hat]


def main():
    t0 = time.time()
    dep = F.read_deposit("certify_joint.json")
    anyk = F.read_deposit("certify_anyk.json")
    pc = F.read_deposit("wb_percand.json")
    out = {"abar": dep["abar"], "designs": [],
           "note": ("the joint and per-candidate levels are levels on "
                    "different statistics, so they do not compare; the "
                    "threshold and the realised error rate do")}

    for rec in dep["designs"]:
        n = np.asarray(rec["n"], float)
        K = rec["K"]
        nbin = rec["nbins"]
        rows, index = J.JointStepTable.tabulate(n, nbin)
        tj = J.JointStepTable(n, rec["alpha0_certified_joint"], nbin,
                              rows=rows, index=index)
        tm = J.JointStepTable(n, rec["alpha0_measured_joint"], nbin,
                              rows=rows, index=index)

        # the per-candidate certificate for the same design, by design name
        pcert = next((x for x in anyk["designs"]
                      if x["design"] == rec["design"] and not x["two_sided"]),
                     None)
        pmeas = next((x for x in pc["designs"]
                      if x["design"] == rec["design"]), None)

        pts = []
        for lam in LAM_POINTS:
            lh = np.array((lam * K)[:K], float)
            cj = float(tj(lh))
            cm = float(tm(lh))
            cp = (percand_threshold(n, pcert["alpha0"], lh)
                  if pcert else [float("nan")] * K)
            pts.append({"lam_hat": lh.tolist(),
                        "c_joint_certified": cj,
                        "c_joint_measured": cm,
                        "c_percand_certified": cp,
                        "joint_certified_is_higher":
                            bool(cj > max(cp)) if pcert else None})

        rec_out = {
            "design": rec["design"], "K": K, "n": rec["n"], "nbins": nbin,
            "alpha0_joint_certified": rec["alpha0_certified_joint"],
            "alpha0_joint_measured": rec["alpha0_measured_joint"],
            "joint_certificate_costs_pct": rec["level_vs_measured_pct"],
            "points": pts,
            "fwer_joint_certified": rec["exact_sup_over_patterns"],
            "fwer_joint_certified_mcse": rec["exact_sup_mcse"],
        }
        if pcert:
            rec_out.update({
                "alpha0_percand_certified": pcert["alpha0"],
                "alpha0_percand_measured": pmeas["alpha0"] if pmeas else None,
                "percand_certificate_costs_pct":
                    (100.0 * (1.0 - pcert["alpha0"] / pmeas["alpha0"])
                     if pmeas else None),
                "fwer_percand_certified": pcert.get("exact_sup_over_patterns"),
            })
            rec_out["percand_certificate_dominates"] = bool(
                all(p["joint_certified_is_higher"] for p in pts))
        out["designs"].append(rec_out)
        print("%-9s joint cert thr %.4f  percand cert thr %.4f  "
              "percand wins %s"
              % (rec["design"], pts[0]["c_joint_certified"],
                 max(pts[0]["c_percand_certified"]),
                 rec_out.get("percand_certificate_dominates")))

    out["percand_certificate_dominates_everywhere"] = bool(
        all(d.get("percand_certificate_dominates") for d in out["designs"]
            if "percand_certificate_dominates" in d))
    out["seconds"] = float(time.time() - t0)
    F.write_deposit("rule_compare.json", out, "frontier_measured")
    print("\nper-candidate certificate dominates everywhere:",
          out["percand_certificate_dominates_everywhere"])


if __name__ == "__main__":
    main()
