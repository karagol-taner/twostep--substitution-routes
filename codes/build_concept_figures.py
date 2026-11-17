#!/usr/bin/env python3
import collections
import csv
import itertools
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
STRUCT = os.path.join(ROOT, 'structures')
OUT = os.path.join(ROOT, 'figures')
DPI = 400

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7.5,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': DPI, 'savefig.dpi': DPI,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'axes.titlelocation': 'left', 'axes.titleweight': 'bold',
})
INK, GREY, PALE = '#1b1b1b', '#9a9a9a', '#d8d8d8'
BLUE, VERM, PLUM = '#0072B2', '#D55E00', '#CC79A7'
GREEN, SKY = '#009E73', '#56B4E9'


def die(m):
    sys.exit('FATAL build_concept_figures.py: %s' % m)


def panel(ax, letter, pad=4):
    ax.set_title(letter, loc='left', pad=pad, fontsize=8.5, fontweight='bold')


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name + '.png')
    fig.savefig(p, metadata={'Software': None})
    plt.close(fig)
    print('  %-8s %5d kB  %s' % (name, os.path.getsize(p) // 1024, p))


B = 'TCAG'
AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: AAS[i] for i, (a, b, c) in enumerate(itertools.product(B, B, B))}
AA20 = sorted(set(AAS) - {'*'})
COD = collections.defaultdict(list)
for c, a in CODON.items():
    if a != '*':
        COD[a].append(c)


def dmin(x, y):
    return min(sum(1 for i in range(3) if cx[i] != cy[i])
               for cx in COD[x] for cy in COD[y])


DIST = {(x, y): dmin(x, y) for x in sorted(AA20) for y in sorted(AA20)}
DIST_SET = set(AA20)
ORDER = list('GAVLIMFWPSTCYNQDEKRH')
if sorted(ORDER) != AA20:
    die('the residue ordering is not the twenty amino acids')


def codon_paths():
    out = os.path.join(DATA, 'codon_paths.json')
    if not os.path.isfile(out):
        die('%s is absent; run 19_dedup_profile.py first' % out)
    rec = json.load(open(out))
    for k in ('target_slots', 'direct_base_changes', 'two_via_one_and_one'):
        if k not in rec:
            die('%s carries no %s' % (out, k))
    return rec


def fig_route():
    fig = plt.figure(figsize=(7.1, 4.15))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.74, 1.26], width_ratios=[1.18, 1],
                          hspace=0.16, wspace=0.42)

    ax = fig.add_subplot(gs[0, :])
    ax.set_xlim(0, 10); ax.set_ylim(0.1, 2.6); ax.axis('off')
    panel(ax, 'a', pad=2)

    def box(x, y, w, h, label, sub, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.05',
                                    fc=fc, ec=INK, lw=0.7))
        ax.text(x + w / 2, y + h * 0.64, label, ha='center', va='center',
                fontsize=9, fontweight='bold', color=INK)
        ax.text(x + w / 2, y + h * 0.22, sub, ha='center', va='center',
                fontsize=5.8, color='#444444')

    box(0.30, 1.62, 1.45, 0.72, 'X', 'parent, scan 1', '#eef3f8')
    box(4.10, 1.62, 1.45, 0.72, 'Y', 'target', '#f7efe8')
    box(0.30, 0.32, 1.45, 0.72, 'Z', 'parent, scan 2', '#eef3f8')

    ar = dict(arrowstyle='-|>', mutation_scale=8, lw=1.1, shrinkA=2, shrinkB=2)
    ax.add_patch(FancyArrowPatch((1.85, 2.00), (4.00, 2.00), color=VERM, **ar))
    ax.text(2.92, 2.16, 'one-step effect, scan 1', ha='center', fontsize=6.4, color=VERM)

    ax.add_patch(FancyArrowPatch((1.02, 1.56), (1.02, 1.10), color=BLUE, **ar))
    ax.text(1.18, 1.32, 'first leg, scan 1', ha='left', fontsize=6.4, color=BLUE)

    ax.add_patch(FancyArrowPatch((1.85, 0.90), (4.28, 1.58), color=GREEN, **ar))
    ax.text(3.55, 0.70, 'second leg, scan 2', ha='center', fontsize=6.4, color=GREEN)
    ax.text(3.55, 0.46, 'whose parent here is Z', ha='center', fontsize=5.8, color='#4a8f7b')

    ax.plot([6.05, 6.05], [0.30, 2.40], color=PALE, lw=1.0)
    ax.text(6.25, 2.26, 'composition predicts', fontsize=6.8, fontweight='bold', color=INK)
    ax.text(6.25, 1.94, 'one-step = first leg + second leg', fontsize=7.2,
            color=INK, style='italic')
    ax.text(6.25, 1.58, 'the residual is the difference, and its null is', fontsize=6.2,
            color='#444444')
    ax.text(6.25, 1.36, 'the same two scans at a position whose', fontsize=6.2, color='#444444')
    ax.text(6.25, 1.14, 'parental residue they already share', fontsize=6.2, color='#444444')
    ax.text(6.25, 0.78, 'no single-parent library measures the second', fontsize=6.2, color=VERM)
    ax.text(6.25, 0.56, 'leg, because no variant in one carries Z as', fontsize=6.2, color=VERM)
    ax.text(6.25, 0.34, 'its parent at this position', fontsize=6.2, color=VERM)

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, 'b')
    M = np.array([[DIST[(x, y)] for y in ORDER] for x in ORDER], dtype=float)
    cmap = matplotlib.colors.ListedColormap(['#ffffff', PALE, SKY, VERM])
    im = ax.imshow(M, cmap=cmap, vmin=-0.5, vmax=3.5, interpolation='nearest')
    ax.set_xticks(range(20)); ax.set_xticklabels(ORDER, fontsize=5.0)
    ax.set_yticks(range(20)); ax.set_yticklabels(ORDER, fontsize=5.0)
    ax.set_xlabel('target residue', labelpad=1)
    ax.set_ylabel('parental residue', labelpad=1)
    ax.tick_params(length=1.4, pad=1)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.5); sp.set_color(GREY)
    cax = ax.inset_axes([0.0, -0.30, 1.0, 0.05])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[0, 1, 2, 3])
    cb.ax.set_xticklabels(['same', '1', '2', '3'], fontsize=5.8)
    cb.ax.tick_params(length=1.4, pad=1)
    cb.outline.set_linewidth(0.4)
    cb.set_label('minimum base changes between codons', fontsize=6, labelpad=2)

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, 'c')
    cp = codon_paths()
    tot = cp['target_slots']
    d = cp['direct_base_changes']
    vals = [d.get('1', 0), d.get('2', 0), d.get('3', 0)]
    xs = np.arange(3)
    pcv = [100.0 * v / tot for v in vals]
    ax.bar(xs, pcv, width=0.62, color=[PALE, SKY, VERM], edgecolor=INK, lw=0.5)
    opened = 100.0 * cp['two_via_one_and_one'] / tot
    ax.bar([1], [opened], width=0.62, color='none', edgecolor=GREEN, lw=1.2,
           hatch='////', zorder=3)
    for x, v in zip(xs, pcv):
        ax.text(x, v + 1.4, '%.0f%%' % v, ha='center', fontsize=6.4, color=INK)
    ax.annotate('%.0f%% reachable as two\nsingle base changes\nthrough Z' % opened,
                xy=(1.31, opened * 0.55), xytext=(2.02, 40),
                fontsize=6.0, color=GREEN, ha='center', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.7, color=GREEN,
                                connectionstyle='arc3,rad=-0.25'))
    ax.set_xticks(xs); ax.set_xticklabels(['1', '2', '3'])
    ax.set_xlabel('base changes for the direct substitution', labelpad=2)
    ax.set_ylabel('census target slots, per cent', labelpad=2)
    ax.set_ylim(0, 62)
    save(fig, 'fig1')


SMILES = {
    'A': 'C[C@@H](N)C(O)=O', 'R': 'N[C@@H](CCCNC(N)=N)C(O)=O',
    'N': 'N[C@@H](CC(N)=O)C(O)=O', 'D': 'N[C@@H](CC(O)=O)C(O)=O',
    'C': 'N[C@@H](CS)C(O)=O', 'Q': 'N[C@@H](CCC(N)=O)C(O)=O',
    'E': 'N[C@@H](CCC(O)=O)C(O)=O', 'G': 'NCC(O)=O',
    'H': 'N[C@@H](Cc1c[nH]cn1)C(O)=O', 'I': 'CC[C@H](C)[C@H](N)C(O)=O',
    'L': 'CC(C)C[C@H](N)C(O)=O', 'K': 'NCCCC[C@H](N)C(O)=O',
    'M': 'CSCC[C@H](N)C(O)=O', 'F': 'N[C@@H](Cc1ccccc1)C(O)=O',
    'P': 'OC(=O)[C@@H]1CCCN1', 'S': 'N[C@@H](CO)C(O)=O',
    'T': 'C[C@@H](O)[C@H](N)C(O)=O', 'W': 'N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O',
    'Y': 'N[C@@H](Cc1ccc(O)cc1)C(O)=O', 'V': 'CC(C)[C@H](N)C(O)=O',
}
FULL = {'A': 'alanine', 'R': 'arginine', 'N': 'asparagine', 'D': 'aspartate',
        'C': 'cysteine', 'Q': 'glutamine', 'E': 'glutamate', 'G': 'glycine',
        'H': 'histidine', 'I': 'isoleucine', 'L': 'leucine', 'K': 'lysine',
        'M': 'methionine', 'F': 'phenylalanine', 'P': 'proline', 'S': 'serine',
        'T': 'threonine', 'W': 'tryptophan', 'Y': 'tyrosine', 'V': 'valine'}


def best_codons(x, y):
    best = (9, None, None)
    for cx in COD[x]:
        for cy in COD[y]:
            d = sum(1 for i in range(3) if cx[i] != cy[i])
            if d < best[0]:
                best = (d, cx, cy)
    return best


def draw_aa(letter, px=460):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    m = Chem.MolFromSmiles(SMILES[letter])
    if m is None:
        die('could not parse the structure for %s' % letter)
    AllChem.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(px, px)
    o = d.drawOptions()
    o.bondLineWidth = 2
    o.clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    import io
    from PIL import Image
    im = Image.open(io.BytesIO(d.GetDrawingText())).convert('RGBA')
    bb = im.split()[-1].getbbox()
    if bb:
        im = im.crop(bb)
    return np.array(im)


def sg2(v):
    return ('%+.2f' % v).replace('-', '\u2212')


def fig_chem():
    rows = [r for r in csv.DictReader(open(os.path.join(DATA, 'residuals.csv')))]
    pick = [r for r in rows if r['phenotype'] == 'binding' and r['pair'] == 'Wuhan|BA.1'
            and r['pos_spike'] == '498' and r['res_to'] == 'S']
    if len(pick) != 1:
        die('the worked chemical route is not unique: %d' % len(pick))
    r = pick[0]
    X, Z, Y = r['res_from'], r['res_via'], r['res_to']
    one, l1, l2 = float(r['one_step']), float(r['first_leg']), float(r['second_leg'])
    resid = float(r['residual'])
    if abs(one - (l1 + l2) - resid) > 1e-6:
        die('the worked chemical route does not satisfy its own definition')

    fig = plt.figure(figsize=(7.1, 3.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.15, 1], wspace=0.20)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.0); ax.axis('off')
    panel(ax, 'a', pad=2)

    H = 1.6
    order = [X, Z, Y]
    xs0 = [0.25, 3.55, 6.95]
    for let, x0 in zip(order, xs0):
        img = draw_aa(let)
        h, w = img.shape[0], img.shape[1]
        wu = min(H * (w / float(h)), 2.9)
        ax.imshow(img, extent=(x0, x0 + wu, 0.35, 0.35 + H), zorder=2,
                  interpolation='antialiased')
        ax.text(x0 + wu / 2, 0.06, '%s, %s' % (FULL[let], let), ha='center',
                fontsize=7, fontweight='bold', color=INK)

    def leg(x0, x1, y, colour, label, val, nb, cx, cy):
        ax.add_patch(FancyArrowPatch((x0, y), (x1, y), color=colour,
                                     arrowstyle='-|>', mutation_scale=9,
                                     lw=1.2, zorder=3))
        mx = (x0 + x1) / 2
        ax.text(mx, y + 0.72, label, ha='center', fontsize=6.6,
                fontweight='bold', color=colour)
        ax.text(mx, y + 0.44, '%s log$_{10}$ $K_\\mathrm{A}$' % sg2(val),
                ha='center', fontsize=6.2, color=colour)
        ax.text(mx, y + 0.16, '%s to %s, %d base change%s'
                % (cx, cy, nb, '' if nb == 1 else 's'),
                ha='center', fontsize=5.8, color='#666666')

    dxy, cx1, cy1 = best_codons(X, Y)
    dxz, cx2, cy2 = best_codons(X, Z)
    dzy, cx3, cy3 = best_codons(Z, Y)
    leg(2.95, 3.45, 2.25, BLUE, 'first leg', l1, dxz, cx2, cy2)
    leg(6.35, 6.85, 2.25, GREEN, 'second leg', l2, dzy, cx3, cy3)

    ax.add_patch(FancyArrowPatch((1.25, 3.62), (8.25, 3.62), color=VERM,
                                 arrowstyle='-|>', mutation_scale=9, lw=1.2,
                                 connectionstyle='arc3,rad=-0.17', zorder=3))
    ax.text(4.75, 5.42, 'direct substitution', ha='center', fontsize=6.8,
            fontweight='bold', color=VERM)
    ax.text(4.75, 5.10, '%s log$_{10}$ $K_\\mathrm{A}$, %s to %s, %d base changes'
            % (sg2(one), cx1, cy1, dxy), ha='center', fontsize=6.2, color=VERM)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, 'b')
    bars = [('direct', one, VERM), ('first leg', l1, BLUE),
            ('second leg', l2, GREEN), ('sum of legs', l1 + l2, GREY)]
    ys = np.arange(len(bars))[::-1]
    ax.barh(ys, [b[1] for b in bars], color=[b[2] for b in bars],
            edgecolor=INK, lw=0.5, height=0.62)
    ax.set_yticks(ys); ax.set_yticklabels([b[0] for b in bars])
    ax.axvline(0, color=INK, lw=0.6)
    for y, b in zip(ys, bars):
        ax.text(b[1] - 0.12 if b[1] < 0 else b[1] + 0.12, y, sg2(b[1]),
                va='center', ha='right' if b[1] < 0 else 'left', fontsize=6.2)
    yb = ys[-1] - 0.62
    for v in (one, l1 + l2):
        ax.plot([v, v], [yb, ys[0] if v == one else ys[-1]], color=PLUM,
                lw=0.6, ls=':', zorder=1)
    ax.annotate('', xy=(one, yb), xytext=(l1 + l2, yb),
                arrowprops=dict(arrowstyle='<->', lw=0.9, color=PLUM))
    ax.text((one + l1 + l2) / 2, yb - 0.30,
            'residual %s' % sg2(resid), ha='center', fontsize=6.6,
            color=PLUM, fontweight='bold')
    ax.set_xlabel('effect on ACE2 binding, log$_{10}$ $K_\\mathrm{A}$', labelpad=2)
    ax.set_xlim(min(l1 + l2, one) - 1.5, max(one, 0) + 0.9)
    ax.set_ylim(ys[-1] - 1.25, ys[0] + 0.5)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    save(fig, 'fig5')


THREE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
         'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
         'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
         'TYR': 'Y', 'VAL': 'V'}


def read_pdb(path):
    ca = collections.defaultdict(dict)
    allat = collections.defaultdict(list)
    for ln in open(path):
        if not ln.startswith('ATOM'):
            continue
        if ln[16] not in ' A':
            continue
        aa = THREE.get(ln[17:20].strip())
        if aa is None:
            continue
        c, r = ln[21], int(ln[22:26])
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        allat[c].append((r, xyz))
        if ln[12:16].strip() == 'CA':
            ca[c][r] = (aa, xyz)
    return ca, allat


def place_labels(ax, items, colour, pad_frac=0.055, sep_frac=0.100):
    if not items:
        return
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    dx = pad_frac * (x1 - x0)
    minsep = sep_frac * (y1 - y0)
    pts = sorted(items, key=lambda t: t[1])
    ys = [t[1] for t in pts]
    for _ in range(500):
        moved = False
        for k in range(1, len(ys)):
            gap = ys[k] - ys[k - 1]
            if gap < minsep:
                shift = (minsep - gap) / 2.0
                ys[k - 1] -= shift
                ys[k] += shift
                moved = True
        if not moved:
            break
    for (x, y, t), yy in zip(pts, ys):
        ax.annotate(t, xy=(x, y), xytext=(x + dx, yy), fontsize=5.8,
                    color=colour, ha='left', va='center', zorder=5,
                    arrowprops=dict(arrowstyle='-', color=colour, lw=0.45,
                                    shrinkA=1.0, shrinkB=1.5))
    ax.set_ylim(min(y0, min(ys) - minsep * 0.6), max(y1, max(ys) + minsep * 0.6))
    ax.set_xlim(x0, max(x1, max(t[0] for t in pts) + dx * 3.0))


def fig_struct():
    jf = os.path.join(DATA, 'interface_distance.json')
    cart = os.path.join(OUT, 'structure_cartoon.png')
    for f in (jf, cart):
        if not os.path.isfile(f):
            die('%s is absent; run 18_interface_distance.py and '
                '20_render_structure.py first' % f)
    IFACE = json.load(open(jf))
    dist = {int(k): v for k, v in IFACE['distance_angstrom'].items()}
    BIND = set(IFACE['by_phenotype']['binding']['departing'])
    EXPR = set(IFACE['by_phenotype']['expression']['departing'])
    RNK = set(IFACE['departing_rank'])
    CUT = IFACE['contact_angstrom']

    fig = plt.figure(figsize=(7.1, 4.0))

    from PIL import Image
    img = np.asarray(Image.open(cart).convert('RGB'))
    ax = fig.add_axes([0.005, 0.145, 0.545, 0.825])
    ax.imshow(img, interpolation='antialiased')
    ax.set_axis_off()
    ax.text(-0.01, 1.01, 'a', transform=ax.transAxes, fontsize=8.5,
            fontweight='bold', ha='left', va='bottom', color=INK)
    fig.legend(handles=[
        Line2D([], [], marker='o', ls='', mfc=VERM, mec=INK, mew=0.4, ms=4.6,
               label='departs, binding'),
        Line2D([], [], marker='o', ls='', mfc=BLUE, mec=INK, mew=0.4, ms=4.6,
               label='departs, expression'),
        Line2D([], [], marker='o', ls='', mfc=PLUM, mec=INK, mew=0.4, ms=4.6,
               label='departs on both'),
        Line2D([], [], marker='o', ls='', mfc=GREEN, mec=INK, mew=0.4, ms=4.6,
               label='departs on rank'),
        Line2D([], [], marker='o', ls='', mfc='#cfcfcf', mec=INK, mew=0.4,
               ms=3.2, label='tested, no departure'),
        Line2D([], [], color='#b5b5b5', lw=2.4, label='receptor binding domain'),
        Line2D([], [], color=SKY, lw=2.4, label='ACE2'),
    ], loc='lower left', frameon=False, handletextpad=0.35, borderpad=0.0,
        labelspacing=0.26, ncol=4, columnspacing=0.8,
        bbox_to_anchor=(0.012, 0.005), fontsize=5.8)

    mag = [r for r in csv.DictReader(open(os.path.join(DATA, 'additivity_test.csv')))
           if float(r['margin']) == 0.0]
    rects = {'b': [0.645, 0.620, 0.340, 0.305], 'c': [0.645, 0.135, 0.340, 0.305]}
    for ph, dep, col, letter in [('binding', BIND, VERM, 'b'),
                                 ('expression', EXPR, BLUE, 'c')]:
        ax = fig.add_axes(rects[letter])
        ax.text(-0.185, 1.22, letter, transform=ax.transAxes, fontsize=8.5,
                fontweight='bold', ha='left', va='top', color=INK)
        sub = [r for r in mag if r['phenotype'] == ph]
        if not sub:
            die('no %s rows in the magnitude test' % ph)
        ax.axhline(0, color=GREY, lw=0.6, ls='--', zorder=1)
        for r in sub:
            p = int(r['pos_spike'])
            ax.scatter(dist[p], float(r['median_residual']),
                       s=26 if p in dep else 14,
                       c=col if p in dep else PALE,
                       edgecolors=INK, linewidths=0.4,
                       zorder=3 if p in dep else 2)
        ax.axvline(CUT, color=GREY, lw=0.6, ls=':', zorder=1)
        lab = []
        for p in sorted(dep):
            vv = [float(r['median_residual']) for r in sub if int(r['pos_spike']) == p]
            if vv:
                lab.append((dist[p], max(vv, key=abs), str(p)))
        place_labels(ax, lab, col)
        d = IFACE['by_phenotype'][ph]
        closer = d['p_departing_closer'] < d['p_departing_farther']
        ax.set_title('%s: departing positions sit %s ACE2\n'
                     'median %.1f \u00c5 against %.1f \u00c5, exact p = %.3f'
                     % (ph, 'closer to' if closer else 'farther from',
                        d['departing_median_angstrom'], d['other_median_angstrom'],
                        min(d['p_departing_closer'], d['p_departing_farther'])),
                     loc='left', fontsize=6.2, fontweight='normal', color=INK,
                     pad=3)
        ax.set_ylabel('routing residual', labelpad=2)
        if letter == 'c':
            ax.set_xlabel('minimum heavy-atom distance to ACE2, \u00c5', labelpad=2)
    save(fig, 'fig7')


if __name__ == '__main__':
    print('building the conceptual, chemical and structural figures')
    fig_route()
    fig_chem()
    fig_struct()
