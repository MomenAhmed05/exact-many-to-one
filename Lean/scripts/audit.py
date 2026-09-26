# -*- coding: utf-8 -*-
"""Check that the Lean development proves what it says it proves.

Four checks, each able to fail:

1. The certificate file is regenerated from the deposits and must be unchanged,
   so Certificates.lean cannot drift from the JSON it claims to check.
2. `lake build` succeeds with no errors and no warnings.
3. No source file contains `sorry`, `admit`, a new `axiom`, or `native_decide`
   (which would trust the compiler rather than the kernel).
4. Every headline theorem, and every Table 19 guarantee, depends on Lean's three
   standard axioms and nothing else, read from `#print axioms`.

Run from the Lean directory:  python scripts/audit.py
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
SRC = os.path.join(PROJ, 'ExactManyToOne')
STANDARD = {'propext', 'Classical.choice', 'Quot.sound'}


def run(cmd):
    return subprocess.run(cmd, cwd=PROJ, capture_output=True, text=True,
                          encoding='utf-8', errors='replace', shell=(os.name == 'nt'))


def main():
    failures = []

    # 1. certificates regenerate identically
    cert = os.path.join(SRC, 'Certificates.lean')
    before = open(cert, encoding='utf-8').read()
    r = run([sys.executable, os.path.join(HERE, 'gen_certificates.py')])
    after = open(cert, encoding='utf-8').read()
    ok = r.returncode == 0 and before == after
    print('%-4s certificates regenerate from the deposits unchanged' % ('OK' if ok else 'FAIL'))
    if not ok:
        failures.append('certificates')

    # 2. clean build
    r = run(['lake', 'build'])
    out = r.stdout + r.stderr
    nerr = len(re.findall(r'^error', out, re.M))
    nwarn = len(re.findall(r'^warning', out, re.M))
    ok = r.returncode == 0 and nerr == 0 and nwarn == 0
    print('%-4s lake build: %d errors, %d warnings' % ('OK' if ok else 'FAIL', nerr, nwarn))
    if not ok:
        failures.append('build')
        print(out[-3000:])

    # 3. no escape hatches in the sources
    bad = []
    for fn in sorted(os.listdir(SRC)):
        if not fn.endswith('.lean'):
            continue
        text = open(os.path.join(SRC, fn), encoding='utf-8').read()
        code = re.sub(r'/-.*?-/', '', text, flags=re.S)
        code = re.sub(r'--.*', '', code)
        for word in ('sorry', 'admit', 'native_decide'):
            if re.search(r'\b%s\b' % word, code):
                bad.append('%s: %s' % (fn, word))
        if re.search(r'^\s*axiom\b', code, re.M):
            bad.append('%s: axiom' % fn)
    ok = not bad
    print('%-4s no sorry, admit, axiom or native_decide in the sources%s'
          % ('OK' if ok else 'FAIL', '' if ok else ': ' + ', '.join(bad)))
    if not ok:
        failures.append('escape hatches')

    # 4. axioms of every headline
    r = run(['lake', 'env', 'lean', 'Audit.lean'])
    out = r.stdout + r.stderr
    blocks = re.findall(r"'([\w.]+)' depends on axioms: \[([^\]]*)\]", out)
    free = re.findall(r"'([\w.]+)' does not depend on any axioms", out)
    expected = len(re.findall(r'^#print axioms', open(os.path.join(PROJ, 'Audit.lean'),
                                                          encoding='utf-8').read(), re.M))
    extra = []
    for name, axs in blocks:
        used = {a.strip() for a in axs.split(',') if a.strip()}
        if not used <= STANDARD:
            extra.append('%s uses %s' % (name, ', '.join(sorted(used - STANDARD))))
    seen = len(blocks) + len(free)
    ok = r.returncode == 0 and not extra and seen == expected
    print('%-4s %d of %d audited theorems use only propext, Classical.choice, Quot.sound'
          % ('OK' if ok else 'FAIL', seen - len(extra), expected))
    for e in extra:
        print('     ' + e)
    if not ok:
        failures.append('axioms')

    print()
    if failures:
        print('AUDIT FAILED: ' + ', '.join(failures))
        sys.exit(1)
    print('LEAN AUDIT PASSED')


if __name__ == '__main__':
    main()
