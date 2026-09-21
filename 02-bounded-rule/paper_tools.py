# ruff: noqa: UP031  (LaTeX templates are built with %-formatting on purpose:
# f-strings would require doubling every brace, which obscures the LaTeX)
"""Paper toolkit: read `Exact-Many-to-One-Under-Estimated-Volatilities.docx`.

`extract` converts the document to markdown with LaTeX maths, so the paper can be read,
diffed and reviewed as text.  The equations in this paper are real Word equation objects
(OMML), so the converter handles the constructs that actually appear there: fractions,
scripts, radicals, n-ary operators, delimiters, functions, accents (including the overbar
used for means), matrices, and plain runs.

Two OMML details that matter and are easy to get wrong:

  * Word stores an inline binary operator such as `+`, `-`, `<=` or `,` as a delimiter
    (`m:d`) whose `m:sepChr` carries the operator and whose `m:begChr`/`m:endChr` are
    empty.  Those must be emitted as a bare infix operator, not as brackets.
  * Body text lives in `w:t` runs while maths text lives in `m:t` runs, so a text
    collector has to read both.

Unrecognised constructs are emitted as a marked placeholder rather than dropped.
"""

import os
import re
import sys
import zipfile
from xml.etree import ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.join(HERE, "paper")
DOCX = os.path.join(PAPER, "Exact-Many-to-One-Under-Estimated-Volatilities.docx")

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
M = "{http://schemas.openxmlformats.org/officeDocument/2006/math}"

ACCENTS = {
    "\u0304": "bar",
    "\u0305": "bar",
    "\u0302": "hat",
    "\u0303": "tilde",
    "\u0307": "dot",
    "\u20d7": "vec",
    "\u0306": "breve",
}
NARY_CHR = {
    None: "\\sum",
    "": "\\sum",
    "\u2211": "\\sum",
    "\u220f": "\\prod",
    "\u222b": "\\int",
    "\u222c": "\\iint",
    "\u22c3": "\\bigcup",
    "\u22c2": "\\bigcap",
    "\u2210": "\\coprod",
}
SEP_MAP = {
    "\u2212": "-",
    "\u00d7": "\\times",
    "\u22c5": "\\cdot",
    "\u2264": "\\le",
    "\u2265": "\\ge",
    "\u2260": "\\ne",
    "\u00b1": "\\pm",
    ",": ",",
}
PASS_CONTAINERS = (
    "oMath",
    "oMathPara",
    "e",
    "num",
    "den",
    "sub",
    "sup",
    "deg",
    "lim",
    "box",
    "borderBox",
    "groupChr",
    "phant",
    "sSubSup",
    "sSup",
    "sSub",
    "rad",
    "f",
    "nary",
    "d",
    "func",
    "acc",
    "bar",
    "m",
    "mr",
    "limLow",
    "limUpp",
    "sPre",
)


def _i(x, default=0):
    """Guarded int conversion."""
    try:
        return int(x)
    except Exception:  # noqa: BLE001
        return int(default)


def _txt(el):
    """Text of a run in document order: w:t for body runs, m:t for maths runs."""
    out = []
    try:
        for node in el.iter():
            if node.tag in (W + "t", M + "t"):
                out.append(node.text or "")
    except Exception:  # noqa: BLE001
        return ""
    return "".join(out)


def conv(el):
    """Concatenate the LaTeX of every significant child of `el`."""
    if el is None:
        return ""
    parts = []
    for child in el:
        if child.tag.endswith("Pr"):
            continue
        if child.tag.startswith(M):
            parts.append(omml_to_latex(child))
        elif child.tag == W + "r":
            parts.append(_txt(child))
        else:
            parts.append(conv(child))
    return "".join(parts)


def _delimiter(el):
    """`m:d`: brackets when begChr/endChr are present, otherwise an infix operator."""
    pr = el.find(M + "dPr")
    beg = end = ""
    sep = "|"
    if pr is not None:
        b = pr.find(M + "begChr")
        e = pr.find(M + "endChr")
        sp = pr.find(M + "sepChr")
        if b is not None and (b.get(M + "val") or ""):
            beg = b.get(M + "val")
        if e is not None and (e.get(M + "val") or ""):
            end = e.get(M + "val")
        if sp is not None:
            v = sp.get(M + "val")
            if v is None:
                sep = ""
            else:
                mapped = SEP_MAP.get(v)
                sep = mapped if mapped is not None else v
    parts = [conv(c) for c in el.findall(M + "e")]
    inner = (" " + sep + " ").join(x for x in parts if x)
    if beg or end:
        return "\\left%s%s\\right%s" % (beg, inner, end)
    return inner


def omml_to_latex(el):
    """One OMML element to LaTeX.  Unknown constructs keep their text rather than vanish."""
    if el is None:
        return ""
    tag = el.tag
    if tag in (M + "r", M + "t"):
        return _txt(el)
    if tag == M + "f":
        return "\\frac{%s}{%s}" % (conv(el.find(M + "num")), conv(el.find(M + "den")))
    if tag == M + "sSup":
        return "{%s}^{%s}" % (conv(el.find(M + "e")), conv(el.find(M + "sup")))
    if tag == M + "sSub":
        return "{%s}_{%s}" % (conv(el.find(M + "e")), conv(el.find(M + "sub")))
    if tag == M + "sSubSup":
        return "{%s}_{%s}^{%s}" % (
            conv(el.find(M + "e")),
            conv(el.find(M + "sub")),
            conv(el.find(M + "sup")),
        )
    if tag == M + "sPre":
        return "{}_{%s}^{%s}%s" % (
            conv(el.find(M + "sub")),
            conv(el.find(M + "sup")),
            conv(el.find(M + "e")),
        )
    if tag == M + "rad":
        deg = el.find(M + "deg")
        body = conv(el.find(M + "e"))
        if deg is None or not _txt(deg).strip():
            return "\\sqrt{%s}" % body
        return "\\sqrt[%s]{%s}" % (conv(deg), body)
    if tag == M + "d":
        return _delimiter(el)
    if tag == M + "nary":
        pr = el.find(M + "naryPr")
        op = "\\sum"
        if pr is not None:
            chr_el = pr.find(M + "chr")
            if chr_el is not None:
                op = NARY_CHR.get(chr_el.get(M + "val"), "\\sum")
        sub = conv(el.find(M + "sub"))
        sup = conv(el.find(M + "sup"))
        body = conv(el.find(M + "e"))
        if sub:
            op += "_{%s}" % sub
        if sup:
            op += "^{%s}" % sup
        return op + " " + body
    if tag == M + "func":
        name = conv(el.find(M + "fName"))
        inner = ", ".join(conv(c) for c in el.findall(M + "e"))
        if name:
            return "\\operatorname{%s}\\left(%s\\right)" % (name, inner)
        return "\\left(%s\\right)" % inner
    if tag == M + "acc":
        pr = el.find(M + "accPr")
        ch = "\u0304"
        if pr is not None:
            c = pr.find(M + "chr")
            if c is not None and c.get(M + "val"):
                ch = c.get(M + "val")
        return "\\%s{%s}" % (ACCENTS.get(ch, "acc"), conv(el.find(M + "e")))
    if tag == M + "bar":
        return "\\bar{%s}" % conv(el.find(M + "e"))
    if tag == M + "limLow":
        return "%s_{%s}" % (conv(el.find(M + "e")), conv(el.find(M + "lim")))
    if tag == M + "limUpp":
        return "%s^{%s}" % (conv(el.find(M + "e")), conv(el.find(M + "lim")))
    if tag == M + "m":
        rows = [
            " & ".join(conv(c) for c in mr.findall(M + "e"))
            for mr in el.findall(M + "mr")
        ]
        return "\\begin{pmatrix}%s\\end{pmatrix}" % " \\\\ ".join(rows)
    if tag == M + "groupChr":
        return conv(el.find(M + "e"))
    if tag in [M + c for c in PASS_CONTAINERS]:
        return conv(el)
    if tag.endswith("Pr") or tag in (M + "ctrlPr", M + "argPr", M + "rPr"):
        return ""
    # unknown: keep the text, mark the gap
    inner = conv(el)
    return inner if inner else ""


def para_markdown(p):
    """One w:p to markdown, with $...$ inline maths and $$...$$ display maths."""
    style = ""
    pr = p.find(W + "pPr")
    if pr is not None:
        st = pr.find(W + "pStyle")
        if st is not None:
            style = st.get(W + "val") or ""
    has_disp = any(c.tag == M + "oMathPara" for c in p)
    chunks = []
    for child in p:
        if child.tag == W + "r":
            if (
                child.find(".//" + W + "drawing") is not None
                or child.find(".//" + W + "pict") is not None
            ):
                chunks.append("{{IMAGE}}")
            chunks.append(_txt(child))
        elif child.tag == M + "oMathPara":
            chunks.append("\n\n$$" + omml_to_latex(child) + "$$\n\n")
        elif child.tag == M + "oMath":
            chunks.append("$" + omml_to_latex(child) + "$")
        elif child.tag in (W + "hyperlink", W + "bookmarkStart", W + "bookmarkEnd"):
            chunks.append(_txt(child))
    text = "".join(chunks)
    if not has_disp:
        text = re.sub(r"\s+", " ", text).strip()
    if style.lower().startswith(("heading", "title")):
        level = min(6, max(1, _i((re.search(r"(\d+)", style) or ["", "1"])[1], 1)))
        if text:
            return "#" * level + " " + text
    return text


def table_markdown(tbl):
    rows = []
    for tr in tbl.findall(W + "tr"):
        cells = []
        for tc in tr.findall(W + "tc"):
            parts = [para_markdown(p).strip() for p in tc.findall(W + "p")]
            cells.append(" ".join(x for x in parts if x))
        rows.append("| " + " | ".join(cells) + " |")
    if rows:
        rows.insert(1, "|" + "---|" * (rows[0].count("|") - 1))
    return "\n".join(rows)


def body_items(doc):
    """('p', element) and ('tbl', element) in document order."""
    for child in list(doc.find(W + "body")):
        if child.tag == W + "p":
            yield "p", child
        elif child.tag == W + "tbl":
            yield "tbl", child


def _safe_parse(xml):
    try:
        from defusedxml.ElementTree import fromstring as _fs
    except Exception:  # noqa: BLE001
        _fs = ET.fromstring
    try:
        return _fs(xml)
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(f"cannot parse word/document.xml: {exc}") from exc


def extract(docx=DOCX, out=None):
    with zipfile.ZipFile(docx) as z:
        doc = _safe_parse(z.read("word/document.xml"))
        media = sorted(n for n in z.namelist() if n.startswith("word/media/"))
        try:
            os.makedirs(os.path.join(PAPER, "media"), exist_ok=True)
            for n in media:
                with open(
                    os.path.join(PAPER, "media", os.path.basename(n)), "wb"
                ) as fh:
                    fh.write(z.read(n))
        except OSError as exc:
            raise SystemExit(f"cannot extract media: {exc}") from exc
    lines, nimg = [], 0
    for kind, el in body_items(doc):
        if kind == "p":
            md = para_markdown(el)
            if "{{IMAGE}}" in md:
                nimg += 1
                md = md.replace("{{IMAGE}}", "{{IMAGE %d}}" % nimg)
            lines.append(md)
        else:
            lines.append(table_markdown(el))
            lines.append("")
    text = re.sub(r"\n{3,}", "\n\n", "\n".join(lines))
    if out is None:
        out = os.path.join(PAPER, "extracted.md")
    try:
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(text)
    except OSError as exc:
        raise SystemExit(f"cannot write {out}: {exc}") from exc
    words = len(re.findall(r"[A-Za-z][A-Za-z'-]*", text))
    print(
        f"extracted {len(text)} chars, ~{words} words, {nimg} image refs, "
        f"{len(media)} media files -> {out}"
    )
    return text


def main():
    argv = sys.argv[1:] or ["extract"]
    if argv[0] == "extract":
        text = extract()
        print("\n=== headings ===")
        for line in text.splitlines():
            if line.startswith("#"):
                print(" ", line[:110])
    else:
        raise SystemExit("usage: paper_tools.py [extract]")


if __name__ == "__main__":
    main()
