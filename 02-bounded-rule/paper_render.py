"""Render a .docx to page images so the document can be inspected visually.

Word does the layout (which also validates that every equation object is
well-formed, since Word refuses to render broken OMML), then pypdfium2 rasterises
the pages.  Both steps are local and deterministic.

Usage:
    python paper_render.py <docx> [out_dir] [pages]
      pages: comma list or range, e.g. "1,4,7" or "5-9"; omitted renders all.
"""

import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _i(x, default=0):
    try:
        return int(x)
    except Exception:  # noqa: BLE001
        return int(default)


def docx_to_pdf(docx_path, pdf_path):
    import win32com.client as win32

    word = win32.DispatchEx("Word.Application")
    word.Visible = False
    word.DisplayAlerts = 0
    try:
        doc = word.Documents.Open(
            str(Path(docx_path).resolve()), ReadOnly=False, AddToRecentFiles=False
        )
        try:
            doc.SaveAs(str(Path(pdf_path).resolve()), FileFormat=17)  # wdFormatPDF
        finally:
            doc.Close(False)
    finally:
        word.Quit()
    return pdf_path


def pdf_to_pngs(pdf_path, out_dir, pages=None, scale=3):
    import pypdfium2 as pdfium

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    doc = pdfium.PdfDocument(str(pdf_path))
    n = len(doc)
    wanted = list(range(n)) if pages is None else pages
    written = []
    for i in wanted:
        if i < 0 or i >= n:
            continue
        page = doc[i]
        img = page.render(scale=scale).to_pil()
        p = out_dir / f"page{i + 1:02d}.png"
        img.save(p)
        written.append(p)
    return n, written


def parse_pages(spec, n_total):
    if not spec:
        return None
    out = set()
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            a, _, b = part.partition("-")
            for v in range(_i(a, 1), _i(b, 1) + 1):
                out.add(v - 1)
        elif part:
            out.add(_i(part, 1) - 1)
    return sorted(x for x in out if 0 <= x < n_total) if n_total else sorted(out)


def main():
    argv = sys.argv[1:]
    if not argv:
        raise SystemExit(__doc__)
    docx = argv[0]
    out_dir = argv[1] if len(argv) > 1 else str(HERE / "paper" / "render")
    page_spec = argv[2] if len(argv) > 2 else None
    t0 = time.time()
    pdf = Path(out_dir) / (Path(docx).stem + ".pdf")
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    docx_to_pdf(docx, pdf)
    n_all, _ = pdf_to_pngs(pdf, out_dir)
    wanted = parse_pages(page_spec, n_all)
    n, files = pdf_to_pngs(pdf, out_dir, pages=wanted)
    print(f"rendered {len(files)} of {n} pages in {time.time() - t0:.1f}s -> {out_dir}")
    for f in files:
        print("  ", f.name)


if __name__ == "__main__":
    main()
