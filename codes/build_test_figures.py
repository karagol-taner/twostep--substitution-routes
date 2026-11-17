#!/usr/bin/env python3
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
FIGS = os.path.join(ROOT, 'figures')
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

INK = '#1b1b1b'
GREY = '#9a9a9a'
PALE = '#d8d8d8'
BLUE = '#0072B2'
VERM = '#D55E00'
PLUM = '#CC79A7'
GREEN = '#009E73'
SKY = '#56B4E9'


def die(m):
    sys.exit('FATAL build_test_figures.py: %s' % m)


def jload(n):
    p = os.path.join(DATA, n)
    if not os.path.isfile(p):
        die('data/%s is absent' % n)
    return json.load(open(p))


def short(name):
    t = name
    for cut in (' Cold-shock protein, DNA-binding domain', ' domain-containing protein',
                ' POU-specific domain', ' Zinc finger C2H2-type domain',
                ' Homeobox domain', ' SH3 domain', ' domain'):
        t = t.replace(cut, '')
    t = t.replace('Ice-structuring protein', 'Ice-structuring')
    t = t.replace('Cold shock-like protein ', '')
    t = t.replace('Uncharacterized protein ', '')
    t = t.replace('N-acetylglucosaminidase', 'N-acetylglucosaminidase')
    return t if len(t) <= 22 else t[:21] + '\u2026'


def panel(ax, letter):
    ax.set_title(letter, loc='left', pad=4, fontsize=8.5, fontweight='bold')


def sgn(v, nd=3):
    return ('%+.*f' % (nd, v)).replace('-', '−')


def mn(v, nd=2):
    return ('%.*f' % (nd, v)).replace('-', '−')


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    p = os.path.join(FIGS, name + '.png')
    fig.savefig(p, metadata={'Software': None})
    plt.close(fig)
    b = open(p, 'rb').read()
    print('  %-8s %7.0f kB  %s' % (name, len(b) / 1024.0,
                                   p.replace(ROOT + '/', '')))
    return p


AGREE = jload('deposit_agreement.json')
NULLM = jload('null_matching.json')
ASSAY = jload('assay_match.json')
DOSE = jload('dose_response.json')
STAT = jload('statistic_choice.json')


def fig_test():
    fig = plt.figure(figsize=(7.1, 5.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.35],
                          hspace=0.42, wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, 'a')
    ex = AGREE.get('examples') or {}
    if 'passing' not in ex or 'failing' not in ex:
        die('deposit_agreement.json carries no worked examples')
    for key, col, mk in (('passing', BLUE, 'o'), ('failing', VERM, 's')):
        e = ex[key]
        x = np.array(e['one_step'], float)
        y = np.array(e['second_leg'], float)
        x = (x - x.mean()) / x.std()
        y = (y - y.mean()) / y.std()
        ax.plot(x, y, mk, ms=3.2, mfc=col, mec='white', mew=0.35,
                linestyle='none',
                label='%s, ρ = %s' % (e['target'].split(' domain')[0][:30], mn(e['rho'])))
    lim = 2.9
    ax.plot([-lim, lim], [-lim, lim], '-', color=GREY, lw=0.7, zorder=0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel('one-step effect, first deposit (z)', labelpad=2)
    ax.set_ylabel('second leg, second deposit (z)', labelpad=2)
    ax.legend(loc='upper left', frameon=False, handletextpad=0.4,
              borderpad=0.1, labelspacing=0.25)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, 'b')
    rt = np.array([r['rho_test'] for r in AGREE['per_pair']], float)
    bins = np.linspace(-1, 1, 41)
    ax.hist(rt, bins=bins, color=PALE, edgecolor=INK, linewidth=0.4)
    ax.axvline(1.0, color=VERM, lw=1.0)
    ax.axvline(float(np.median(rt)), color=BLUE, lw=0.9, ls=(0, (4, 2)))
    top = ax.get_ylim()[1]
    ax.set_ylim(0, top * 1.16)
    ax.text(1.0, top * 1.04, 'required', ha='right', va='bottom',
            color=VERM, fontsize=6.3)
    ax.text(0.02, 0.965,
            '%s of %s pairs above 0.9'
            % (format(int((rt > 0.9).sum()), ','), format(len(rt), ',')),
            transform=ax.transAxes, ha='left', va='top', color=INK,
            fontsize=6.2)
    ax.text(np.median(rt) - 0.05, top * 0.66, 'median %s' % mn(np.median(rt), 3),
            ha='right', va='center', color=BLUE, fontsize=6.3)
    ax.set_xlim(-1.02, 1.02)
    ax.set_xlabel('rank correlation at the test position', labelpad=2)
    ax.set_ylabel('deposit pairs', labelpad=2)
    ax.set_title('b', loc='left', pad=4, fontsize=8.5, fontweight='bold')

    ax = fig.add_subplot(gs[1, :])
    panel(ax, 'c')
    rows = sorted(STAT['by_target'], key=lambda r: r['median_gap_rho'])
    xs = np.arange(len(rows))
    med = np.array([r['median_rho'] for r in rows], float)
    ctl = np.array([r['median_rho'] - r['median_gap_rho'] for r in rows], float)
    cols = [VERM if r['median_gap_rho'] < 0 else BLUE for r in rows]
    ax.vlines(xs, ctl, med, color=PALE, lw=0.9, zorder=0)
    ax.scatter(xs, ctl, s=11, marker='_', color=GREY, zorder=2,
               linewidths=0.9)
    ax.scatter(xs, med, s=13, c=cols, zorder=3, linewidths=0)
    ax.axhline(1.0, color=VERM, lw=0.9)

    ax.set_xticks(xs)
    ax.set_xticklabels([short(r['target']) for r in rows], rotation=55,
                       ha='right', fontsize=6.3)
    ax.set_xlim(-0.8, len(rows) - 0.2)
    lo = min(float(np.min(ctl)), float(np.min(med)))
    ax.set_ylim(min(-0.05, lo - 0.06), 1.09)
    ax.set_ylabel('median rank correlation', labelpad=2)
    ax.text(len(rows) - 0.6, 1.015, 'required', ha='right', va='bottom',
            color=VERM, fontsize=6.3)
    ax.legend(handles=[
        Line2D([], [], marker='o', ls='none', ms=3.6, mfc=BLUE, mec=BLUE,
               label='test position, agreeing better than its own control'),
        Line2D([], [], marker='o', ls='none', ms=3.6, mfc=VERM, mec=VERM,
               label='test position, agreeing worse than its own control'),
        Line2D([], [], marker='_', ls='none', ms=6, color=GREY,
               label='control position, matched on dispersion'),
    ], loc='lower right', frameon=False, handletextpad=0.5,
        borderpad=0.1, labelspacing=0.3, bbox_to_anchor=(1.0, 0.02))
    return save(fig, 'fig3')


def fig_null():
    fig = plt.figure(figsize=(7.1, 3.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.92, 1.0, 1.22],
                          wspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, 'a')
    series = [('test\nposition', np.array([r['rho_test'] for r in
                                           NULLM['per_pair']], float), BLUE),
              ('null,\ncount', np.array(
                  [r['null_count_matched'] for r in NULLM['per_pair']],
                  float), VERM),
              ('null,\ndispersion', np.array(
                  [r['null_spread_matched'] for r in NULLM['per_pair']],
                  float), GREEN)]
    parts = ax.violinplot([s[1] for s in series], positions=[0, 1, 2],
                          widths=0.78, showextrema=False, showmedians=False)
    for body, (_, _, col) in zip(parts['bodies'], series):
        body.set_facecolor(col)
        body.set_alpha(0.35)
        body.set_edgecolor(col)
        body.set_linewidth(0.6)
    for i, (_, v, col) in enumerate(series):
        m = float(np.median(v))
        ax.plot([i - 0.28, i + 0.28], [m, m], '-', color=col, lw=1.3)
        ax.text(i, -1.08, mn(m, 3), ha='center', va='top', color=col,
                fontsize=6.2)
    ax.axhline(1.0, color=INK, lw=0.9, ls=(0, (4, 2)))
    ax.text(2.45, 1.0, 'required', ha='right', va='bottom', color=INK,
            fontsize=6.3)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([s[0] for s in series], fontsize=6.4)
    ax.set_xlim(-0.65, 2.65)
    ax.set_ylim(-1.25, 1.15)
    ax.set_ylabel('rank correlation', labelpad=2)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, 'b')
    so = np.array([r['shift_count_matched'] for r in NULLM['per_pair']], float)
    sn = np.array([r['shift_spread_matched'] for r in NULLM['per_pair']], float)
    bins = np.linspace(-1.0, 1.2, 45)
    ax.hist(so, bins=bins, color=VERM, alpha=0.45, edgecolor=VERM,
            linewidth=0.4, label='matched on count')
    ax.hist(sn, bins=bins, color=GREEN, alpha=0.45, edgecolor=GREEN,
            linewidth=0.4, label='matched on dispersion')
    ax.axvline(0.0, color=INK, lw=0.9, ls=(0, (4, 2)))
    for v, col in ((np.median(so), VERM), (np.median(sn), GREEN)):
        ax.axvline(float(v), color=col, lw=1.0)
    ax.set_xlabel('shift from the null, where it must be 0', labelpad=2)
    ax.set_ylabel('deposit pairs', labelpad=2)
    ax.legend(loc='upper right', frameon=False, handlelength=1.1,
              handletextpad=0.5, borderpad=0.1, labelspacing=0.25)

    ax = fig.add_subplot(gs[0, 2])
    panel(ax, 'c')
    keys = list(DOSE['by_differing_positions_count_matched'].keys())
    xs = np.arange(len(keys))
    cm = [DOSE['by_differing_positions_count_matched'][k]['median_shift']
          for k in keys]
    sm = [DOSE['by_differing_positions_spread_matched'][k]['median_shift']
          for k in keys]
    npair = [DOSE['by_differing_positions_count_matched'][k]['pairs']
             for k in keys]
    ax.axhline(0.0, color=INK, lw=0.7)
    ax.axvspan(-0.45, 0.45, color=PALE, alpha=0.7, zorder=0)
    ax.plot(xs, cm, '-o', color=VERM, lw=1.0, ms=3.6, label='matched on count')
    ax.plot(xs, sm, '-s', color=GREEN, lw=1.0, ms=3.4,
            label='matched on dispersion')
    ax.plot([0], [0], '*', color=INK, ms=9, zorder=5)
    ax.annotate('fixed at 0 by\nthe algebra', xy=(0.05, 0.0),
                xytext=(0.62, 0.052), fontsize=6.2, color=INK, ha='left',
                va='center',
                arrowprops=dict(arrowstyle='-', color=INK, lw=0.6,
                                shrinkA=0, shrinkB=2))
    ax.set_xticks(xs)
    ax.set_xticklabels(['%s\n%s' % (k, format(n, ',')) for k, n in
                        zip(keys, npair)], fontsize=5.8)
    ax.set_xlabel('positions differing, and deposit pairs', labelpad=2)
    ax.set_ylabel('median shift from the null', labelpad=2)
    ax.set_xlim(-0.55, len(keys) - 0.45)
    ax.legend(loc='upper right', frameon=False, handlelength=1.3,
              handletextpad=0.5, borderpad=0.1, labelspacing=0.25,
              bbox_to_anchor=(1.02, 1.03))
    return save(fig, 'fig4')


if __name__ == '__main__':
    print('building the exact-test figures at %d dpi' % DPI)
    fig_test()
    fig_null()
