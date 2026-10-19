"""Kill test (Section 4.5): coordinatewise monotonicity of the aligned-factor
null law in (u, v, w) = (beta_0^2, beta_1^2, beta_0*beta_1).

Declaration frozen BEFORE any grid run: declarations/W2A-DECLARATION.json.

The evaluator mirrors w2_factor.AlignedFactor.noreject with K = 1 and
(beta_0 - beta_1)^2 -> d^2 = u + v - 2w, which makes interior cone
points (w^2 < uv) well defined; on the cone surface this IS the native
law (checked by AN1), and only d^2 enters m (checked by AN2).

Commands:
  smoke     anchors AN1-AN5 at both levels + a subsample sweep; writes
            the anchor files, no rows
  run       anchors for one level, then the full declared grid into an
            append-only JSONL (crash-resumable; P-5 writes always)
  verdict   diffs both levels' JSONLs -> verdicts, H-12 survival,
            D0-2 shape -> deposits/w2a_kill.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import traceback
from multiprocessing import Pool

import numpy as np
from scipy.stats import norm
from scipy.special import ndtr

# Paper 1's modules (fcore, w2_factor and fcore's own h12_* imports) live in
# ../04-efficient-rule. Put that folder on the path without writing into it.
_REPO04 = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        os.pardir, "04-efficient-rule"))
if os.path.isdir(_REPO04) and _REPO04 not in sys.path:
    sys.path.insert(0, _REPO04)

import fcore as F
from h12_quad import chi2_over_df_nodes, normal_nodes
from w2_factor import AlignedFactor, ncx2_nodes

HERE = os.path.dirname(os.path.abspath(__file__))
DECL_PATH = os.path.join(HERE, "declarations", "W2A-DECLARATION.json")
W5_DEP = os.path.join(HERE, "deposits", "w5_limits.json")
# New-work deposits live beside this script; verdict must not write into
# the published repo (F.write_deposit uses fcore's module-global DEP).
F.DEP = os.path.join(HERE, "deposits")
os.makedirs(F.DEP, exist_ok=True)

N = 36
NU = 35.0
C = 2.20
AX = 61
DU = 0.0375
DV = 0.0375
DW = 0.075
W_LO = -2.25
WORKERS = 8
FINE = {"zf": (8, 10), "z0": (8, 10), "vf": (8, 10), "s0": 48, "si": 48}
COORDS = ("u", "v", "w")
REGIONS = ("interior", "layer", "combined")
CLASS_ORDER = {"monotone": 0, "monotone_with_exceptions": 1,
               "inconclusive": 2, "non_monotone": 3}

_ST = {}


# --------------------------------------------------------------- axes / layer
def u_at(i):
    return i * DU


def v_at(j):
    return j * DV


def w_at(k):
    return W_LO + k * DW


def feasible(u, v, w):
    return w * w <= u * v * (1.0 + 1e-12) + 1e-15


def is_layer(i, j, k):
    """Boundary layer per declaration: within one w-step of the cone surface,
    or on the cone's degenerate faces."""
    if i == 0 or j == 0:
        return True
    return (math.sqrt(u_at(i) * v_at(j)) - abs(w_at(k))) <= DW


# ---------------------------------------------------------------- evaluator
def _build_level(g):
    zf, wzf = normal_nodes(*[int(x) for x in g["zf"]])
    z0, wz0 = normal_nodes(*[int(x) for x in g["z0"]])
    q, wvf = chi2_over_df_nodes(NU, *[int(x) for x in g["vf"]])
    return dict(zf=zf, wzf=wzf, z0=z0, wz0=wz0, vf=q * NU, wvf=wvf,
                s0m=int(g["s0"]), sim=int(g["si"]))


def _s0_nodes_f(st, u):
    if u not in st["s0"]:
        lv = st["lv"]
        nodes, w = ncx2_nodes(NU, u * lv["vf"], lv["s0m"])
        st["s0"][u] = (1.0 * nodes / NU, w)
    return st["s0"][u]


def _s1_nodes_f(st, v, ratio):
    key = (v, ratio)
    if key not in st["s1"]:
        lv = st["lv"]
        nodes, w = ncx2_nodes(NU, v * lv["vf"] / ratio, lv["sim"])
        st["s1"][key] = (ratio * nodes / NU, w)
    return st["s1"][key]


def _noreject(lv, s0tab, s0_w, sitab, s_w, d, ratio, c, n):
    """P(T_1 <= c) under the K=1 aligned law with signed d = beta_0 - beta_1.

    Byte-for-byte the structure of AlignedFactor.noreject with K = 1
    (m depends on beta_0, beta_1 only through d; s0 noncentrality u*vf,
    s1 noncentrality v*vf/ratio).
    """
    zf, wzf, z0, wz0 = lv["zf"], lv["wzf"], lv["z0"], lv["wz0"]
    vf, wvf = lv["vf"], lv["wvf"]
    nf = float(n)
    wt = wzf[:, None, None] * wz0[None, :, None]
    m = (d * zf[:, None] + z0[None, :]) / math.sqrt(nf)
    sc = math.sqrt(ratio / nf)
    total = 0.0
    for jj in range(vf.size):
        s0 = s0tab[jj]
        si = sitab[jj]
        den = np.sqrt(si[None, :] / nf + s0[:, None] / nf)
        arg = (c * den[None, None, :, :] + m[:, :, None, None]) / sc
        # ndtr IS the standard normal CDF (norm.cdf delegates to it);
        # calling it directly avoids scipy's generic rv_continuous.cdf
        # dispatch, which hangs in its np.place slow branch on some
        # h12-level arrays (incident 23/09, task lam=1.0 i=44 cell
        # j=6 k=27).  Same values, no model change.
        acc = (ndtr(arg) * s_w[None, None, None, :]).sum(3)
        total += wvf[jj] * float((acc * wt * s0_w[None, None, :]).sum())
    return 1.0 - total


def _fwer_d(st, u, v, d, ratio):
    s0tab, s0_w = _s0_nodes_f(st, u)
    sitab, s_w = _s1_nodes_f(st, v, ratio)
    return _noreject(st["lv"], s0tab, s0_w, sitab, s_w, d, ratio,
                     st["c"], st["n"])


def _fwer_uvw(st, u, v, w, ratio):
    d = math.sqrt(max(u + v - 2.0 * w, 0.0))
    return _fwer_d(st, u, v, d, ratio)


def _local_state(level, decl):
    return dict(level=level, lv=_build_level(decl["levels"][level]),
                s0={}, s1={}, c=float(decl["design"]["c"]),
                n=int(decl["design"]["n"]))


def _init_worker(level):
    global _ST
    _ST = _local_state(level, load_decl())


def _task(args):
    ratio, i = args
    st = _ST
    lv = st["lv"]
    c = st["c"]
    n = st["n"]
    u = u_at(i)
    s0tab, s0_w = _s0_nodes_f(st, u)
    row = []
    for j in range(AX):
        v = v_at(j)
        uv = u * v
        sitab, s_w = _s1_nodes_f(st, v, ratio)
        rj = [None] * AX
        for k in range(AX):
            w = w_at(k)
            if not feasible(u, v, w):
                continue
            d = math.sqrt(max(u + v - 2.0 * w, 0.0))
            f = _noreject(lv, s0tab, s0_w, sitab, s_w, d, ratio, c, n)
            if not math.isfinite(f):
                raise ValueError("non-finite f at u=%r v=%r w=%r lam=%r"
                                 % (u, v, w, ratio))
            rj[k] = f
        row.append(rj)
    return {"level": st["level"], "ratio": ratio, "i": i, "f": row}


def load_decl():
    with open(DECL_PATH, encoding="utf-8") as fh:
        return json.load(fh)


def w5_tokens():
    with open(W5_DEP, encoding="utf-8") as fh:
        dep = json.load(fh)
    rows = dep["f3_ladders"]["rows"]

    def row(tag):
        return next(r for r in rows if r["tag"] == tag)

    a1 = row("A1 rho ladder, beta=(1,1)")["values"][0]
    a10 = row("A1 rho ladder, beta=(0,0)")["values"][0]
    a2 = row("A2 t ladder, equal direction (1,1)")["values"]
    return a1, a10, a2


# ----------------------------------------------------------------- verdict
def _components(idxs):
    S = set(idxs)
    comps = []
    while S:
        start = S.pop()
        stack = [start]
        comp = [start]
        while stack:
            x, y, z = stack.pop()
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        p = (x + dx, y + dy, z + dz)
                        if p in S:
                            S.discard(p)
                            stack.append(p)
                            comp.append(p)
        comps.append(comp)
    return comps


def verdict_diffs(diffs, tol, n_axis=AX):
    """diffs: [(idx_tuple, value)] -> declared verdict dict."""
    npos = nneg = 0
    for _, d in diffs:
        if d > tol:
            npos += 1
        elif d < -tol:
            nneg += 1
    ntot = len(diffs)
    conf = npos + nneg
    frac = (conf / ntot) if ntot else 0.0
    out = {"n_diffs": ntot, "n_pos": npos, "n_neg": nneg,
           "n_indet": ntot - conf, "confident_fraction": frac,
           "tol": tol}
    if ntot == 0:
        out["class"] = "inconclusive"
        out["reason"] = "no differences in region"
        return out
    if npos > 0 and nneg > 0 and npos == nneg:
        out["class"] = "non_monotone"
        out["reason"] = "confident sign tie (%d vs %d)" % (npos, nneg)
        return out
    direction = "+" if npos > nneg else ("-" if nneg > npos else None)
    if direction is None:
        out["class"] = "inconclusive"
        out["reason"] = "no confident differences"
        return out
    if direction == "+":
        viol = [idx for idx, d in diffs if d < -tol]
    else:
        viol = [idx for idx, d in diffs if d > tol]
    if not viol:
        if frac >= 0.5:
            out["class"] = "monotone"
            out["direction"] = direction
        else:
            out["class"] = "inconclusive"
            out["reason"] = "confident fraction %.3f < 0.5" % frac
        return out
    comps = _components(viol)
    axes = list(zip(*viol))
    ext = [max(a) - min(a) + 1 for a in axes]
    bmin = [min(a) for a in axes]
    bmax = [max(a) for a in axes]
    lim = 0.20 * n_axis
    out.update(n_viol=len(viol), n_components=len(comps),
               bbox_min=bmin, bbox_max=bmax, bbox_extents=ext,
               direction=direction)
    if len(comps) == 1 and all(e <= lim + 1e-9 for e in ext):
        out["class"] = "monotone_with_exceptions"
        out["reason"] = "single component inside 20%% bbox (extents %s)"
        out["reason"] = out["reason"] % ext
    else:
        out["class"] = "non_monotone"
        if len(comps) > 1:
            out["reason"] = ("%d separated exception components"
                             % len(comps))
        else:
            out["reason"] = "exception bbox %s wider than 20%% of %d" % (
                ext, n_axis)
    return out


def collect_diffs(rows, lams):
    """rows: {(lam, i): row} -> {(lam, coord): [(idx, value, layer)]}."""
    R = {}
    for (lam, i), row in rows.items():
        R.setdefault(lam, {})[i] = row
    out = {}
    for lam in lams:
        Rl = R.get(lam, {})
        ul, vl, wl = [], [], []
        for i in range(AX):
            row = Rl.get(i)
            if row is None:
                continue
            for j in range(AX):
                rj = row[j]
                for k in range(AX):
                    f = rj[k]
                    if f is None:
                        continue
                    idx = (i, j, k)
                    lay = is_layer(i, j, k)
                    if i + 1 < AX:
                        row2 = Rl.get(i + 1)
                        f2 = row2[j][k] if row2 is not None else None
                        if f2 is not None:
                            ul.append((idx, f2 - f,
                                       lay or is_layer(i + 1, j, k)))
                    if j + 1 < AX:
                        f2 = row[j + 1][k]
                        if f2 is not None:
                            vl.append((idx, f2 - f,
                                       lay or is_layer(i, j + 1, k)))
                    if k + 1 < AX:
                        f2 = rj[k + 1]
                        if f2 is not None:
                            wl.append((idx, f2 - f,
                                       lay or is_layer(i, j, k + 1)))
        out[(lam, "u")] = ul
        out[(lam, "v")] = vl
        out[(lam, "w")] = wl
    return out


def verdict_all(diffmap, tol, lams):
    res = {}
    for lam in lams:
        for coord in COORDS:
            entries = diffmap[(lam, coord)]
            for region in REGIONS:
                if region == "combined":
                    sel = entries
                elif region == "layer":
                    sel = [e for e in entries if e[2]]
                else:
                    sel = [e for e in entries if not e[2]]
                res[(lam, coord, region)] = verdict_diffs(
                    [(idx, val) for idx, val, _ in sel], tol)
    return res


def worst_across_lams(res, lams, coord, region):
    pick = None
    per = {}
    for lam in lams:
        v = res[(lam, coord, region)]
        per[str(lam)] = v["class"]
        if pick is None or CLASS_ORDER[v["class"]] > CLASS_ORDER[
                pick["class"]]:
            pick = dict(v)
            pick["worst_lam"] = lam
    pick["per_lam_class"] = per
    return pick


def d02_shape(classes):
    if "non_monotone" in classes:
        return "(iii) non-monotone"
    if all(c == "monotone" for c in classes):
        return "(i) monotone"
    if all(c in ("monotone", "monotone_with_exceptions")
           for c in classes):
        return "(ii) monotone-with-exceptions"
    return "UNDECIDED (inconclusive region, no non-monotone)"


# ----------------------------------------------------------------- anchors
def run_anchors(level, decl):
    t0 = time.time()
    entries = []
    try:
        st = _local_state(level, decl)

        # --- AN1: wrapper equals the native (MC-validated) evaluator
        spec = decl["anchors"]["AN1_wrapper_equals_native"]
        grid = {k: (tuple(v) if isinstance(v, list) else v)
                for k, v in decl["levels"][level].items()
                if k in ("zf", "z0", "vf", "s0", "si")}
        worst = 0.0
        detail = []
        for b0, b1 in spec["cells_beta"]:
            b0, b1 = float(b0), float(b1)
            for lam in spec["lam"]:
                lam = float(lam)
                nat = AlignedFactor(N, [b0, b1], [lam],
                                    grid=grid).fwer(C)
                wrap = _fwer_uvw(st, b0 ** 2, b1 ** 2, b0 * b1, lam)
                dev = abs(nat - wrap)
                worst = max(worst, dev)
                detail.append({"beta": [b0, b1], "lam": lam,
                               "native": nat, "wrapper": wrap,
                               "dev": float(dev)})
        entries.append({"name": "AN1_wrapper_equals_native",
                        "ok": bool(worst <= spec["tol"]),
                        "max_dev": float(worst),
                        "tol": spec["tol"], "detail": detail})

        # --- AN2: only d^2 enters (sign invariance of the interior extension)
        spec = decl["anchors"]["AN2_sign_invariance"]
        worst = 0.0
        detail = []
        for u, v, w in spec["cells_uvw"]:
            u, v, w = float(u), float(v), float(w)
            d = math.sqrt(u + v - 2.0 * w)
            for lam in spec["lam"]:
                lam = float(lam)
                fp = _fwer_d(st, u, v, d, lam)
                fm = _fwer_d(st, u, v, -d, lam)
                dev = abs(fp - fm)
                worst = max(worst, dev)
                detail.append({"uvw": [u, v, w], "lam": lam, "dev":
                               float(dev)})
        entries.append({"name": "AN2_sign_invariance",
                        "ok": bool(worst <= spec["tol"]),
                        "max_dev": float(worst),
                        "tol": spec["tol"], "detail": detail})

        # --- AN3: deposit value anchors (primary level only, per declaration)
        if level == "primary":
            a1, a10, a2 = w5_tokens()
            specd = decl["anchors"]["AN3_deposit_values"]
            stf = dict(level="fine", lv=_build_level(
                {"zf": FINE["zf"], "z0": FINE["z0"], "vf": FINE["vf"],
                 "s0": FINE["s0"], "si": FINE["si"]}),
                s0={}, s1={}, c=st["c"], n=st["n"])
            detail = []
            ok = True
            for cell in specd["cells"]:
                u, v, w = (float(x) for x in cell["uvw"])
                lam = float(cell["lam"])
                tok = a1 if u == 1.0 else a10
                drift = abs(float(cell["fine_token"]) - tok)
                fv = _fwer_uvw(stf, u, v, w, lam)
                gv = _fwer_uvw(st, u, v, w, lam)
                d_f = abs(fv - tok)
                d_g = abs(gv - tok)
                ok = ok and d_f <= cell["tol_fine"] and \
                    d_g <= cell["tol_grid"]
                detail.append({"uvw": [u, v, w], "lam": lam,
                               "deposit_token": tok,
                               "declared_token": float(
                                   cell["fine_token"]),
                               "declaration_drift": float(drift),
                               "fine_value": fv, "dev_fine": d_f,
                               "tol_fine": cell["tol_fine"],
                               "grid_value": gv, "dev_grid": d_g,
                               "tol_grid": cell["tol_grid"]})
            entries.append({"name": "AN3_deposit_values",
                            "ok": bool(ok), "detail": detail})

        # --- AN3b: deposit direction (A2 loading ladder order)
        spec = decl["anchors"]["AN3b_deposit_direction"]
        lam3 = float(spec["lam"])
        f11 = _fwer_uvw(st, 1.0, 1.0, 1.0, lam3)
        margins = {}
        ok = True
        for u, v, w in spec["cells_uvw"]:
            u, v, w = float(u), float(v), float(w)
            fv = _fwer_uvw(st, u, v, w, lam3)
            margin = f11 - fv
            ok = ok and margin > 0.0
            margins["%g,%g,%g" % (u, v, w)] = float(margin)
        entries.append({"name": "AN3b_deposit_direction",
                        "ok": bool(ok), "f_at_1": f11,
                        "margins": margins,
                        "note": ("order relation only; measured "
                                 "quadrature noise cannot flip these "
                                 "margins (see declaration)")})

        # --- AN5: synthetic verdict machinery
        an5 = an5_cases(decl)
        entries.append(an5)
    except Exception:
        entries.append({"name": "anchors_execution", "ok": False,
                        "error": traceback.format_exc()})
    rec = {"level": level, "entries": entries,
           "all_ok": all(bool(e.get("ok")) for e in entries),
           "seconds": time.time() - t0,
           "ran_at": time.strftime("%Y-%m-%d %H:%M:%S")}
    path = os.path.join(HERE, "w2a_kill_anchors_%s.json" % level)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=1, allow_nan=False)
    return rec


def an5_cases(decl):
    tol = float(decl["fd"]["tol_fd"])
    cases = []
    base = [((i, 10, 10), 0.05) for i in range(50)]
    base += [((i, 60, 60), 0.05) for i in range(50)]
    r = verdict_diffs(base, tol)
    cases.append(("all-confident-positive",
                  r["class"] == "monotone" and r.get("direction") == "+"))

    d = list(base) + [((20, 20, 20), -0.05), ((21, 20, 20), -0.05),
                      ((20, 21, 20), -0.05), ((20, 20, 21), -0.05),
                      ((21, 21, 21), -0.05)]
    r = verdict_diffs(d, tol)
    cases.append(("single-small-negative-blob",
                  r["class"] == "monotone_with_exceptions"))

    d = list(base) + [((20, 20, 20), -0.05), ((50, 50, 50), -0.05)]
    r = verdict_diffs(d, tol)
    cases.append(("two-separated-negative-clusters",
                  r["class"] == "non_monotone"))

    d = [((i, 10, 10), 0.0005 if i % 2 == 0 else -0.0005)
         for i in range(60)]
    r = verdict_diffs(d, tol)
    cases.append(("all-indeterminate", r["class"] == "inconclusive"))

    d = [((i, 10, 10), 0.05) for i in range(50)]
    d += [((i, 20, 20), -0.05) for i in range(50)]
    r = verdict_diffs(d, tol)
    cases.append(("confident-sign-tie", r["class"] == "non_monotone"
                  and "tie" in r.get("reason", "")))

    ok = all(passed for _, passed in cases)
    return {"name": "AN5_verdict_machinery", "ok": bool(ok),
            "cases": [{"case": n, "passed": bool(p)} for n, p in cases]}


# --------------------------------------------------------------------- run
def load_rows(path):
    rows = {}
    if not os.path.exists(path):
        return rows
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # trailing partial line from a crash
            rows[(float(r["ratio"]), int(r["i"]))] = r["f"]
    return rows


def cmd_run(level):
    decl = load_decl()
    print("[run %s] anchors..." % level, flush=True)
    a = run_anchors(level, decl)
    print("[run %s] anchors %s (%.1fs)"
          % (level, "ALL PASS" if a["all_ok"] else "FAIL", a["seconds"]),
          flush=True)
    if not a["all_ok"]:
        for e in a["entries"]:
            if not e.get("ok"):
                print("[run %s] FAILED: %s"
                      % (level, json.dumps(e)[:2000]), flush=True)
        return 2

    lams = [float(x) for x in decl["design"]["lam_sweep"]]
    jsonl = os.path.join(HERE, "w2a_kill_%s_rows.jsonl" % level)
    done = load_rows(jsonl)
    tasks = [(lam, i) for lam in lams for i in range(AX)
             if (lam, i) not in done]
    print("[run %s] rows: %d done, %d tasks remaining, %d workers"
          % (level, len(done), len(tasks), WORKERS), flush=True)
    if not tasks:
        print("[run %s] nothing to do" % level, flush=True)
        return 0
    t0 = time.time()
    n = len(tasks)
    pool = Pool(WORKERS, initializer=_init_worker, initargs=(level,))
    try:
        with open(jsonl, "a", encoding="utf-8") as fh:
            for cnt, res in enumerate(pool.imap_unordered(_task, tasks,
                                                          chunksize=1), 1):
                fh.write(json.dumps(res) + "\n")
                fh.flush()
                el = time.time() - t0
                eta = el / cnt * (n - cnt)
                print("[%s] %d/%d %.1f%% lam=%s i=%d  %.1fs/task  "
                      "elapsed %.2fh eta %.2fh"
                      % (level, cnt, n, 100.0 * cnt / n, res["ratio"],
                         res["i"], el / cnt, el / 3600.0, eta / 3600.0),
                      flush=True)
    finally:
        pool.close()
        pool.join()
    print("[run %s] DONE in %.2fh -> %s"
          % (level, (time.time() - t0) / 3600.0, jsonl), flush=True)
    return 0


# ----------------------------------------------------------------- verdict
def cmd_smoke():
    decl = load_decl()
    t0 = time.time()
    for level in ("primary", "h12"):
        a = run_anchors(level, decl)
        print("[smoke %s] anchors %s (%.1fs)"
              % (level, "ALL PASS" if a["all_ok"] else "FAIL",
                 a["seconds"]), flush=True)
        if not a["all_ok"]:
            for e in a["entries"]:
                if not e.get("ok"):
                    print("[smoke] FAILED entry:", json.dumps(e)[:3000],
                          flush=True)
            return 2
    st = _local_state("primary", decl)
    vals = []
    for i in range(0, AX, 15):
        for j in range(0, AX, 15):
            u, v = u_at(i), v_at(j)
            for k in range(0, AX, 15):
                w = w_at(k)
                if not feasible(u, v, w):
                    continue
                f = _fwer_uvw(st, u, v, w, 1.0)
                if not (math.isfinite(f) and 0.0 <= f <= 1.0):
                    print("[smoke] BAD VALUE", u, v, w, f, flush=True)
                    return 2
                vals.append(f)
    print("[smoke] subsample n=%d f in [%.8f, %.8f]  ALL OK  (%.1fs total)"
          % (len(vals), min(vals), max(vals), time.time() - t0), flush=True)
    return 0


def count_declared_feasible():
    cnt = 0
    for i in range(AX):
        u = u_at(i)
        for j in range(AX):
            v = v_at(j)
            for k in range(AX):
                if feasible(u, v, w_at(k)):
                    cnt += 1
    return cnt


def cmd_verdict():
    decl = load_decl()
    t0 = time.time()
    lams = [float(x) for x in decl["design"]["lam_sweep"]]
    tol = float(decl["fd"]["tol_fd"])
    expected = {(lam, i) for lam in lams for i in range(AX)}
    rows = {}
    missing = {}
    anchors = {}
    for level in ("primary", "h12"):
        path = os.path.join(HERE, "w2a_kill_%s_rows.jsonl" % level)
        rows[level] = load_rows(path)
        miss = sorted(expected - set(rows[level].keys()))
        missing[level] = [{"lam": m[0], "i": m[1]} for m in miss]
        apath = os.path.join(HERE, "w2a_kill_anchors_%s.json" % level)
        if os.path.exists(apath):
            with open(apath, encoding="utf-8") as fh:
                anchors[level] = json.load(fh)
    if any(missing.values()):
        obj = {"unit": "2A", "status": "incomplete",
               "reason": "grid rows missing; verdict not computed",
               "missing_tasks": {k: len(v) for k, v in missing.items()},
               "missing_sample": {k: v[:10] for k, v in missing.items()},
               "tol_fd": tol, "lam_sweep": lams}
        p = F.write_deposit("w2a_kill.json", obj, "measured")
        print("[verdict] INCOMPLETE rows: %s -> %s"
              % {k: len(v) for k, v in missing.items()}, flush=True)
        return 1

    declared_feasible = count_declared_feasible()
    coverage = {}
    diffs = {}
    verdicts = {}
    for level in ("primary", "h12"):
        ev = sum(1 for row in rows[level].values()
                 for rj in row for f in rj if f is not None)
        coverage[level] = {"declared_feasible_points":
                           declared_feasible,
                           "evaluated_points": ev,
                           "ratio": ev / float(declared_feasible * len(lams)),
                           "expected_total":
                           declared_feasible * len(lams),
                           "n_lams": len(lams)}
        diffs[level] = collect_diffs(rows[level], lams)
        verdicts[level] = verdict_all(diffs[level], tol, lams)

    def vjson(level):
        out = {}
        for lam in lams:
            out[str(lam)] = {}
            for coord in COORDS:
                out[str(lam)][coord] = {}
                for region in REGIONS:
                    out[str(lam)][coord][region] = verdicts[level][
                        (lam, coord, region)]
        return out

    combined = {}
    for coord in COORDS:
        combined[coord] = {}
        for region in REGIONS:
            combined[coord][region] = worst_across_lams(
                verdicts["primary"], lams, coord, region)

    # H-12 survival
    pairs = []
    mismatches = []
    for lam in lams:
        for coord in COORDS:
            for region in REGIONS:
                vp = verdicts["primary"][(lam, coord, region)]
                vh = verdicts["h12"][(lam, coord, region)]
                survived = vp["class"] == vh["class"]
                dir_note = None
                if survived and vp["class"] == "monotone":
                    survived = vp.get("direction") == vh.get("direction")
                    dir_note = [vp.get("direction"), vh.get("direction")]
                rec = {"lam": lam, "coord": coord, "region": region,
                       "primary_class": vp["class"],
                       "h12_class": vh["class"],
                       "primary_direction": vp.get("direction"),
                       "h12_direction": vh.get("direction"),
                       "directions": dir_note,
                       "survived": bool(survived)}
                pairs.append(rec)
                if not survived:
                    mismatches.append(rec)
    all_survived = not mismatches

    classes = [combined[c]["combined"]["class"] for c in COORDS]
    shape = d02_shape(classes)
    if not all_survived:
        route = ("H-12 KILL: the verdict does not survive the "
                 "refinement change; no route may be selected from 2A")
    elif shape.startswith("(i)"):
        route = "Route A candidate: monotone in each coordinate"
    elif shape.startswith("(ii)"):
        route = ("Route A candidate with exception-band caveat "
                 "(D0-2 shape (ii): certificate must exclude the "
                 "exception bbox)")
    elif shape.startswith("(iii)"):
        route = "Route B required: non-monotone in some coordinate"
    else:
        route = ("UNDECIDED: inconclusive region; extend the tolerance "
                 "study or rerun at a finer level")

    obj = {
        "unit": "2A",
        "declaration": "declarations/W2A-DECLARATION.json",
        "status": "complete",
        "design": decl["design"],
        "grid": decl["grid"],
        "levels": {k: decl["levels"][k] for k in ("primary", "h12")},
        "tol_fd": tol,
        "anchors": anchors,
        "coverage": coverage,
        "verdicts": {"primary": vjson("primary"), "h12": vjson("h12")},
        "combined_primary": combined,
        "h12_survival": {"all_survived": all_survived,
                         "n_pairs": len(pairs),
                         "n_mismatches": len(mismatches),
                         "mismatches": mismatches},
        "d02_shape": shape,
        "route_hint": route,
        "scope_notes": decl["scope_notes"],
        "seconds": time.time() - t0,
        "ran_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    p = F.write_deposit("w2a_kill.json", obj, "measured")
    print("[verdict] coverage primary=%.6f h12=%.6f"
          % (coverage["primary"]["ratio"], coverage["h12"]["ratio"]),
          flush=True)
    print("[verdict] classes combined: %s -> %s"
          % (dict(zip(COORDS, classes)), shape), flush=True)
    print("[verdict] H-12 survival: %d/%d pairs survived -> %s"
          % (len(pairs) - len(mismatches), len(pairs), route), flush=True)
    print("[verdict] deposit -> %s (%.1fs)" % (p, time.time() - t0),
          flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("smoke")
    r = sub.add_parser("run")
    r.add_argument("--level", choices=("primary", "h12"), required=True)
    sub.add_parser("verdict")
    args = ap.parse_args(argv)
    if args.cmd == "smoke":
        return cmd_smoke()
    if args.cmd == "run":
        return cmd_run(args.level)
    if args.cmd == "verdict":
        return cmd_verdict()
    return 2


if __name__ == "__main__":
    sys.exit(main())
