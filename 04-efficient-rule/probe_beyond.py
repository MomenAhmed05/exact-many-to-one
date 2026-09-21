"""Corruption probes for verify_beyond.py.

A verifier that passes proves nothing until it has been shown to fail for the
right reason.  The preceding programme learned this the expensive way: before
whole-row matching was added, seven corruption probes all PASSED, because
substring containment cannot see a transposition -- both values are still
somewhere in the document.

Each probe here perturbs ONE deposit field, reruns the verifier in a fresh
process against the UNCHANGED report, and requires that the run fails.  The
deposit is restored afterwards whatever happens.

Two kinds are used deliberately:

  scalar        a single quantity moved by a small amount, which a substring
                check can catch;
  transposition two values swapped INSIDE one table row, which only a
                whole-row check can catch.

    python probe_beyond.py
Deposit: deposits/probe_beyond.json
"""

from __future__ import annotations

import copy
import io
import json
import os
import subprocess
import sys
import time

import fcore as F

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with io.open(os.path.join(F.DEP, name), encoding="utf-8") as fh:
        return json.load(fh)


def save(name, obj):
    with io.open(os.path.join(F.DEP, name), "w", encoding="utf-8") as fh:
        fh.write(json.dumps(obj, indent=1) + "\n")


def run_verifier():
    r = subprocess.run([sys.executable, os.path.join(HERE, "verify_beyond.py")],
                       cwd=HERE, capture_output=True, text=True)
    fails = [ln.strip()[5:].strip() for ln in r.stdout.split("\n")
             if ln.strip().startswith("FAIL")]
    return r.returncode, fails


# ------------------------------------------------------------------- mutations
def bump_scalar(key, by=1e-3):
    def f(d):
        d[key] = d[key] + by
    return f


def bump_nested(path, by=1):
    def f(d):
        node = d
        for k in path[:-1]:
            node = node[k]
        node[path[-1]] = node[path[-1]] + by
    return f


def transpose_row(listkey, index, a, b):
    """Swap two values inside one row.  Both values stay in the document."""
    def f(d):
        r = d[listkey][index]
        r[a], r[b] = r[b], r[a]
    return f


def transpose_map(mapkey, ka, kb):
    def f(d):
        m = d[mapkey]
        m[ka], m[kb] = m[kb], m[ka]
    return f


PROBES = [
    ("scalar: invariant hidden pp", "w1_invariant.json",
     bump_scalar("max_hidden_pp", 0.01)),
    ("scalar: region max gain", "w1_region.json",
     bump_scalar("max_gain_pp_from_dropping_monotonicity", 0.5)),
    ("scalar: stress breach count", "w3_stress.json",
     bump_scalar("n_breaches", 1)),
    ("scalar: manifest refusals", "w4_manifest.json",
     bump_scalar("n_refused", -1)),
    ("transposition: factor falsifier brute/exact", "w2_factor.json",
     transpose_row("falsifiers", 3, "brute", "exact")),
    # NOT monotone/randomised: the whole point of that result is that the two
    # are equal, so swapping them is a no-op and the probe would be vacuous.
    ("transposition: region monotone/det-floor", "w1_region.json",
     transpose_row("designs", 5, "loss_pp_monotone",
                   "loss_pp_deterministic_floor")),
    ("transposition: serial repair separate/pooled", "w3_serial_repair.json",
     lambda d: d.update({
         "separate_worst_at_mu0": d["pooled_worst_at_mu0"],
         "pooled_worst_at_mu0": d["separate_worst_at_mu0"]})),
    ("transposition: asymptotic two methods", "w3_asymptotic.json",
     transpose_map("worst_size", "Newey-West, normal critical value",
                   "Newey-West, t critical value")),
    ("nested: certificate box point count", "w2_cert.json",
     bump_nested(["box", "points"], 2)),
]


def main():
    t0 = time.time()
    base_rc, base_fails = run_verifier()
    print("baseline: returncode %d, %d failures" % (base_rc, len(base_fails)))
    if base_rc != 0:
        print("  (the baseline already fails; probes are still meaningful but "
              "the report is not clean)")
        for f in base_fails:
            print("   " + f)

    rows = []
    for tag, dep, mutate in PROBES:
        path = os.path.join(F.DEP, dep)
        if not os.path.exists(path):
            rows.append({"probe": tag, "deposit": dep, "skipped": True,
                         "caught": False})
            print("  %-46s SKIPPED (no %s)" % (tag, dep))
            continue
        orig = load(dep)
        try:
            bad = copy.deepcopy(orig)
            mutate(bad)
            save(dep, bad)
            rc, fails = run_verifier()
            new = [f for f in fails if f not in base_fails]
            caught = bool(new)
        finally:
            save(dep, orig)
        rows.append({"probe": tag, "deposit": dep, "skipped": False,
                     "caught": caught, "returncode": rc,
                     "new_failures": new[:4],
                     "n_new_failures": len(new)})
        print("  %-46s %s   %s"
              % (tag, "CAUGHT" if caught else "MISSED",
                 (new[0][:64] if new else "nothing failed")))

    live = [r for r in rows if not r["skipped"]]
    out = {
        "baseline_returncode": base_rc,
        "baseline_failures": base_fails,
        "probes": rows,
        "n_probes": len(rows),
        "n_live": len(live),
        "n_caught": int(sum(1 for r in live if r["caught"])),
        "all_caught": bool(live and all(r["caught"] for r in live)),
        "n_transposition_probes": int(sum(
            1 for r in live if r["probe"].startswith("transposition"))),
        "n_transposition_caught": int(sum(
            1 for r in live
            if r["probe"].startswith("transposition") and r["caught"])),
        "note": ("Each probe moves one deposit field and requires the "
                 "verifier to fail on a check it was not already failing.  "
                 "The transposition probes are the ones that matter: they "
                 "leave every value present in the document and are invisible "
                 "to substring containment, so they pass unless table rows "
                 "are matched whole."),
        "seconds": float(time.time() - t0),
    }
    F.write_deposit("probe_beyond.json", out, "diagnostic")
    print("\n%d of %d live probes caught (%d of %d transpositions); %.0f s"
          % (out["n_caught"], out["n_live"], out["n_transposition_caught"],
             out["n_transposition_probes"], out["seconds"]))
    return 0 if out["all_caught"] else 1


if __name__ == "__main__":
    sys.exit(main())
