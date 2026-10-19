# Figures for paper 2. Every value is read from a deposit under deposits/; nothing is typed in.
#   python make_figs.py   ->  figs/fig1_cover.png, figs/fig2_power.png, figs/fig3_mvt.png
import os, re
from fractions import Fraction as F
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
DEP = os.path.join(HERE, 'deposits')
OUT = os.path.join(HERE, 'figs'); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman', 'DejaVu Serif'], 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 220})
INK, MUTE, RED, BLUE = '#1A1A1A', '#8A8A8A', '#B3362C', '#2E4597'

def claims(name):
    out = {}
    for ln in open(os.path.join(DEP, name), encoding='utf-8'):
        m = re.match(r'CLAIM (.+?): (.*)$', ln.strip())
        if m: out[m.group(1)] = m.group(2)
    return out

# ---- Figure 1: the cover of record, first lambda layer, projected on (b0, b1)
cells = []
for ln in open(os.path.join(DEP, 'f8c_cover_witness-R1.txt'), encoding='utf-8'):
    box = ln.strip().split('#')[0]
    (l0, l1), (a0, a1), (c0, c1) = [[F(x) for x in p.strip('[]').split(',')] for p in box.split('|')]
    if l0 == F(1, 25):
        cells.append((float(a0), float(a1), float(c0), float(c1)))
fig, ax = plt.subplots(figsize=(4.2, 4.0))
for a0, a1, c0, c1 in cells:
    ax.add_patch(Rectangle((a0, c0), a1 - a0, c1 - c0, fill=False, lw=0.25, ec=INK))
ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5); ax.set_aspect('equal')
ax.set_xlabel('benchmark loading $b_0$'); ax.set_ylabel('candidate loading $b_1$')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig1_cover.png')); plt.close(fig)
print('fig1 cells', len(cells))

# ---- Figure 2: paired vs Welch power at the four points
pv = claims('paired_vs_welch_probe_claims.txt')
rows = []
for key, val in pv.items():
    m = re.match(r'PVW-(.+?) \((.+)\)', key)
    pw = re.search(r'power at paired-80% effect = ([\d.]+) / ([\d.]+) / ([\d.]+)', val)
    if m and pw: rows.append((m.group(1), m.group(2), [float(x) for x in pw.groups()]))
order = ['same sign', 'no loadings', 'opposite', 'corner']
rows.sort(key=lambda r: order.index(r[0]))
fig, ax = plt.subplots(figsize=(5.6, 2.9))
w = 0.26
for i, (name, pt, (p, wu, wc)) in enumerate(rows):
    ax.bar(i - w, p, w, color=BLUE, label='paired t' if i == 0 else None)
    ax.bar(i, wu, w, color=RED, label='Welch at 2.7174' if i == 0 else None)
    ax.bar(i + w, wc, w, color=MUTE, label='Welch at $c^*(\\theta)$ (needs $\\theta$)' if i == 0 else None)
ax.set_xticks(range(len(rows)), ['%s\n(%s)' % (n, p) for n, p, _ in rows], fontsize=7.5)
ax.set_ylabel('power'); ax.set_ylim(0, 1)
ax.legend(frameon=False, fontsize=7.5, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.13))
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig2_power.png')); plt.close(fig)

# ---- Figure 3: multivariate-t inflation against n
pa = claims('p_verify_A_claims.txt')
k2 = {}
for key, val in pa.items():
    m = re.match(r'PVA-T3-K2-n(\d+)$', key)
    if m:
        for rho, ratio in re.findall(r'([+-]\d\.\d+):[\d.]+\(([\d.]+)x\)', val):
            k2.setdefault(float(rho), []).append((int(m.group(1)), float(ratio)))
k5 = {}
for key, val in pa.items():
    m = re.match(r'PVA-T3-K5-(\w+)$', key)
    if m:
        k5[m.group(1)] = [(int(n), float(r)) for n, r in re.findall(r'n(\d+):[\d.]+\(([\d.]+)x\)', val)]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 2.8), sharey=True)
styles = {0.0: ('-', MUTE), 0.5: ('-', BLUE), 0.9: ('-', RED)}
for rho, (ls, col) in styles.items():
    pts = sorted(k2[rho]); a1.plot([p[0] for p in pts], [p[1] for p in pts], ls, color=col, marker='o', ms=3, lw=1.2, label='$\\rho$ = %.1f' % rho)
a1.set_title('two candidates (computed)', fontsize=8.5); a1.set_xscale('log'); a1.set_xlabel('dates n'); a1.set_ylabel('true rate / level')
cols5 = {'EQ0': MUTE, 'EQ5': BLUE, 'EQ9': RED, 'BM': INK, 'MIX': '#B8860B'}
for name, pts in k5.items():
    pts = sorted(pts); a2.plot([p[0] for p in pts], [p[1] for p in pts], color=cols5.get(name, INK), marker='o', ms=3, lw=1.2, label=name)
a2.set_title('five candidates (simulated)', fontsize=8.5); a2.set_xscale('log'); a2.set_xlabel('dates n')
for a in (a1, a2):
    a.axhline(1.0, ls='--', lw=0.8, color=INK); a.legend(frameon=False, fontsize=7)
    a.set_xticks([6, 8, 12, 24, 36, 60, 120]); a.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig3_mvt.png')); plt.close(fig)
print('figures written to', OUT)
