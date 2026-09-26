# -*- coding: utf-8 -*-
"""Gate for the short companion paper, Exact-Many-to-One-Proofs.docx.

Three checks, each able to fail:

1. Every number the companion quotes (anything with a decimal point, an exponent,
   a thousands separator or a percent sign) appears in the working paper, whose
   own numbers are gated against the deposits. The companion therefore cannot
   quote a number the record does not support.
2. Every Lean name the companion cites resolves to a declaration in the Lean
   development, checked by Lean itself with `#check`.
3. Theorems, lemmas, corollaries and propositions are numbered consecutively
   from 1, and every equation reference points at a displayed equation.

Run from the repository root:  python paper/verify_proofs_paper.py
"""
import os
import re
import subprocess
import sys
import tempfile

import docx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SHORT = os.path.join(HERE, 'Exact-Many-to-One-Proofs.docx')
MAIN = os.path.join(HERE, 'Exact-Many-to-One-Under-Estimated-Volatilities.docx')
LEAN = os.path.join(ROOT, 'Lean')
M = '{http://schemas.openxmlformats.org/officeDocument/2006/math}'
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'


def text_of(path):
    """All text, prose and equations, in document order."""
    d = docx.Document(path)
    out = []
    for el in d.element.body.iter():
        if el.tag in (W + 't', M + 't') and el.text:
            out.append(el.text)
        elif el.tag == W + 'p':
            out.append('\n')
    return ''.join(out)


def norm(s):
    for a, b in (('‑', '-'), ('−', '-'), ('–', '-'), ('\xa0', ' ')):
        s = s.replace(a, b)
    return s


def main():
    results = []

    def chk(name, ok, detail=''):
        results.append(ok)
        print('  %-4s %s%s' % ('OK' if ok else 'FAIL', name, ('  [%s]' % detail) if detail and not ok else ''))

    short = norm(text_of(SHORT))
    main_txt = norm(text_of(MAIN))

    # 1. numbers
    nums = set(re.findall(r'(?<![\w.])(\d[\d,]*\.\d+(?:e-?\d+)?%?|\d+e-?\d+|\d{1,3}(?:,\d{3})+|\d+(?:\.\d+)?%)(?![\w])',
                          short))
    missing = sorted(n for n in nums if n not in main_txt)
    chk('%d quoted numbers all appear in the working paper' % len(nums), not missing, ', '.join(missing))

    # 2. Lean names
    d = docx.Document(SHORT)
    names = set()
    for p in d.paragraphs:
        if p.text.startswith('Lean: '):
            names.update(n.strip() for n in p.text[len('Lean: '):].split(','))
    for t in d.tables:
        for row in t.rows[1:]:
            cells = [c.text for c in row.cells]
            if len(cells) == 3 and cells[2] in ('Reduction', 'Dependence', 'Envelope', 'Enclosure',
                                                  'Gaussian', 'Calibrated', 'Invariance'):
                names.update(n.strip() for n in cells[1].split(','))
    names = sorted(n for n in names if n)
    src = 'import ExactManyToOne\nopen ExactManyToOne\n' + ''.join('#check @%s\n' % n for n in names)
    with tempfile.NamedTemporaryFile('w', suffix='.lean', delete=False, encoding='utf-8', dir=LEAN) as fh:
        fh.write(src)
        tmp = fh.name
    try:
        r = subprocess.run(['lake', 'env', 'lean', os.path.basename(tmp)], cwd=LEAN, capture_output=True,
                           text=True, encoding='utf-8', errors='replace', shell=(os.name == 'nt'))
    finally:
        os.remove(tmp)
    errs = re.findall(r'error.*?Unknown (?:identifier|constant) `?([\w.]+)', r.stdout + r.stderr)
    chk('%d cited Lean names all resolve' % len(names), r.returncode == 0 and not errs,
        ', '.join(errs) or (r.stdout + r.stderr)[-400:])

    # 3. numbering
    for kind in ('Theorem', 'Lemma', 'Corollary', 'Proposition'):
        found = [int(n) for n in re.findall(r'(?m)^%s (\d+) \(' % kind, short)]
        chk('%s numbered 1..%d in order' % ({'Corollary': 'corollaries'}.get(kind, kind.lower() + 's'), len(found)),
            found == list(range(1, len(found) + 1)), str(found))
    displayed = set(int(n) for n in re.findall(r'(?m)^\((\d+)\)$', short))
    refs = set(int(n) for n in re.findall(r'(?:equation|\b[Ii]n|of|gives|coverage|from) \((\d+)\)', short))
    chk('%d displayed equations numbered 1..%d' % (len(displayed), len(displayed)),
        displayed == set(range(1, len(displayed) + 1)), str(sorted(displayed)))
    bad = sorted(refs - displayed)
    chk('every equation reference points at a displayed equation', not bad, str(bad))

    print()
    print('%d checks passed, %d failed' % (sum(results), len(results) - sum(results)))
    if all(results):
        print('PROOFS PAPER VERIFIED')
    else:
        print('PROOFS PAPER HAS FAILURES')
        sys.exit(1)


if __name__ == '__main__':
    main()
