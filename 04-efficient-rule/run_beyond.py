"""Reproduce the W1-W4 programme end to end, then rebuild and check the report.

    python run_beyond.py            run everything, in dependency order
    python run_beyond.py report     rebuild the report and verify only
    python run_beyond.py list       show the stages and which deposits exist

Each stage writes one deposit.  Stages are independent of each other except that
`w2_cert` imports the evaluator from `w2_factor`, and the report and verifier
come last.  Environment and input hashes are recorded so a rerun can be compared
with this one rather than merely repeated.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import subprocess
import sys
import time

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))

STAGES = [
    ("w1_invariant", "w1_invariant.json",
     "maximal invariant, Hunt-Stein route, bin-resolution bound"),
    ("w1_region", "w1_region.json",
     "monotonicity and determinism priced by linear program"),
    ("w1_lowdf", "w1_lowdf.json",
     "the low-degrees-of-freedom map over both group sizes"),
    ("w2_factor", "w2_factor.json",
     "the aligned latent-factor exact null law and its falsifiers"),
    ("w2_cert", "w2_cert.json",
     "what ignoring the factor costs, and a grid certificate over the box"),
    ("w2_unequal", "w2_unequal.json",
     "unequal histories: what the incumbent's prehistory buys"),
    ("w3_tiers", "w3_tiers.json",
     "validity tiers A, B and C for serial and factor nuisances"),
    ("w3_serial_repair", "w3_serial_repair.json",
     "why tiers B and C fail for serial dependence, and pooling"),
    ("w3_asymptotic", "w3_asymptotic.json",
     "HAC and block bootstraps measured against the exact rule"),
    ("w3_stress", "w3_stress.json",
     "the stress battery"),
    ("w2_identify", "w2_identify.json",
     "whether the latent loadings are knowable, and from what"),
    ("w4_manifest", "w4_manifest.json",
     "the promotion-gate verifier and its refusal probes"),
]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(65536), b""):
            h.update(c)
    return h.hexdigest()


def run(mod, extra=None):
    t0 = time.time()
    cmd = [sys.executable, os.path.join(HERE, mod + ".py")] + (extra or [])
    print("\n=== %s ===" % mod, flush=True)
    r = subprocess.run(cmd, cwd=HERE)
    return {"module": mod, "returncode": r.returncode,
            "seconds": round(time.time() - t0, 1)}


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    if arg == "list":
        for mod, dep, what in STAGES:
            mark = "have" if os.path.exists(os.path.join(F.DEP, dep)) else "--  "
            print("  [%s] %-18s %s" % (mark, mod, what))
        return 0

    results = []
    if arg == "all":
        for mod, _dep, _what in STAGES:
            extra = ["example"] if mod == "w4_manifest" else None
            results.append(run(mod, extra))

    results.append(run("w_report"))
    # the probes need a report to verify against, and the report needs their
    # deposit, so the first build is a bootstrap and the second is the real one
    results.append(run("probe_beyond"))
    results.append(run("coverage_beyond"))
    results.append(run("w_report"))
    results.append(run("w_milestones"))
    # the paper is rebuilt from the same deposits, then checked twice: once by
    # the frontier verifier for the structure it owns and once by this
    # update's, and then measured for coverage
    results.append(run("paper_figs_beyond"))
    results.append(run("paper_build_frontier", ["build"]))
    results.append(run("paper_verify_frontier"))
    results.append(run("paper_verify_beyond"))
    results.append(run("coverage_paper_beyond"))
    v = run("verify_beyond")
    results.append(v)

    import platform
    import numpy
    import scipy
    env = {
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }
    srcs = dict((m + ".py", sha(os.path.join(HERE, m + ".py")))
                for m, _d, _w in STAGES
                if os.path.exists(os.path.join(HERE, m + ".py")))
    for extra in ("fcore.py", "k1lab.py", "wc_frontier.py", "w_report.py",
                  "w_milestones.py", "verify_beyond.py", "probe_beyond.py",
                  "coverage_beyond.py", "paper_content_beyond.py",
                  "paper_build_frontier.py", "paper_figs_beyond.py",
                  "paper_verify_beyond.py",
                  "coverage_paper_beyond.py",
                  "manifest.schema.json", "PROMOTION-PROTOCOL.md",
                  "REPORT-BEYOND.md", "ROADMAP.md"):
        p = os.path.join(HERE, extra)
        if os.path.exists(p):
            srcs[extra] = sha(p)
    deps = dict((d, sha(os.path.join(F.DEP, d)))
                for d in ([x[1] for x in STAGES]
                          + ["probe_beyond.json",
                             "coverage_beyond.json",
                             "coverage_paper_beyond.json"])
                if os.path.exists(os.path.join(F.DEP, d)))

    out = {"stages": results, "environment": env, "source_sha256": srcs,
           "deposit_sha256": deps,
           "all_ok": all(r["returncode"] == 0 for r in results),
           "verify_returncode": v["returncode"]}
    io.open(os.path.join(HERE, "RUN-BEYOND.json"), "w",
            encoding="utf-8").write(json.dumps(out, indent=1) + "\n")

    print("\n%-20s %8s %8s" % ("stage", "rc", "seconds"))
    for r in results:
        print("%-20s %8d %8.1f" % (r["module"], r["returncode"], r["seconds"]))
    print("\nall ok: %s  (manifest written to RUN-BEYOND.json)" % out["all_ok"])
    return 0 if out["all_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
