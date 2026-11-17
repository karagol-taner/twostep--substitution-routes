#!/usr/bin/env python3
import csv
import json
import os
import re
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
F = lambda r, k: float(r[k])

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
SKY = '#56B4E9'


def sgn(v, nd=3):
    return ('%+.*f' % (nd, v)).replace('-', '\u2212')


def die(m):
    sys.exit('FATAL build_figures.py: %s' % m)


def need(p):
    if not os.path.isfile(p):
        die('%s is absent' % p)
    return p


def jload(n):
    return json.load(open(need(os.path.join(DATA, n))))


def cload(n):
    return list(csv.DictReader(open(need(os.path.join(DATA, n)))))


def panel(ax, letter):
    ax.set_title(letter, loc='left', pad=4, fontsize=8.5, fontweight='bold')


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    p = os.path.join(FIGS, name + '.png')
    fig.savefig(p, metadata={'Software': None})
    plt.close(fig)
    b = open(p, 'rb').read()
    print('  %-8s %7.0f kB  %s' % (name, len(b) / 1024.0, p.replace(ROOT + '/', '')))
    return p


def fig_resource():
    prof = jload('dedup_profile.json')
    fig = plt.figure(figsize=(6.4, 2.35))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.30, 1.0], wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, 'a')
    counts = prof['routes_per_cluster_sorted']
    xs = np.arange(1, len(counts) + 1)
    ax.bar(xs, counts, width=0.82, color=SKY, edgecolor=INK, lw=0.3)
    ax.set_yscale('log')
    ax.set_xlabel('protein cluster, ordered by routes')
    ax.set_ylabel('complete routes')
    ax.set_xlim(0.2, len(counts) + 0.8)
    med = prof['median_cluster_routes']
    ax.axhline(med, color=GREY, lw=0.7, ls=(0, (2, 2)))
    ax.text(len(counts) * 0.99, med * 1.25, 'median %s' % f'{med:,}',
            ha='right', va='bottom', fontsize=6.2, color=GREY)
    ax.annotate('%d clusters, %s routes;\nthe largest carries %.0f%%'
                % (prof['clusters'], f"{prof['routes']:,}",
                   100.0 * prof['largest_cluster_share']),
                xy=(1.6, counts[0]), xytext=(len(counts) * 0.34, counts[0] * 0.95),
                fontsize=6.2, color=INK, ha='left', va='top',
                arrowprops=dict(arrowstyle='-', lw=0.6, color=INK,
                                connectionstyle='arc3,rad=0.2'))

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, 'b')
    dp = np.array([int(r['differing_positions'])
                   for r in cload('background_pairs_mavedb.csv')])
    bins = np.logspace(0, np.log10(dp.max() * 1.06), 16)
    ax.hist(dp, bins=bins, color=SKY, edgecolor='white', lw=0.35)
    ax.set_xscale('log')
    ax.set_xlabel('differing positions per deposit pair')
    ax.set_ylabel('deposit pairs')
    one = int((dp == 1).sum())
    top = ax.get_ylim()[1]
    ax.set_ylim(0, top * 1.26)
    ax.text(0.97, 0.97, '%s pairs differ at one position' % f'{one:,}',
            transform=ax.transAxes, fontsize=6.2, color=INK,
            ha='right', va='top')
    return save(fig, 'fig2')


def fig_rbd(AX=None):
    mag = [r for r in cload('additivity_test.csv') if float(r['margin']) == 0.0]
    rnk = cload('rank_test.csv')
    res = cload('residuals.csv')
    own = AX is None
    if own:
        fig = plt.figure(figsize=(6.4, 4.5))
        gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.52)
        gstop = gs[0].subgridspec(1, 2, wspace=0.26)
        AX = {'a': lambda: fig.add_subplot(gstop[0, 0]),
              'b': lambda: fig.add_subplot(gstop[0, 1]),
              'c': lambda: fig.add_subplot(gs[1])}

    PAIRS = ['Wuhan|BA.1', 'Wuhan|BA.2', 'BA.1|BA.2']
    MARK = {'Wuhan|BA.1': 'o', 'Wuhan|BA.2': 's', 'BA.1|BA.2': '^'}
    OFF = {'Wuhan|BA.1': -0.24, 'Wuhan|BA.2': 0.0, 'BA.1|BA.2': 0.24}

    def slots(rows):
        pos = sorted({int(r['pos_spike']) for r in rows})
        return pos, {p: i for i, p in enumerate(pos)}

    ax = AX['a']()
    panel(ax, 'a')
    sub = [r for r in mag if r['phenotype'] == 'binding']
    pos, idx = slots(sub)
    nul = {}
    for r in sub:
        nul.setdefault(idx[int(r['pos_spike'])], []).append(F(r, 'null_median_abs'))
    xs = np.arange(len(pos))
    band = np.array([np.mean(nul[i]) for i in xs])
    ax.fill_between(xs, -band, band, color=PALE, lw=0, zorder=1)
    ax.axhline(0, color=INK, lw=0.6, zorder=2)
    for pr in PAIRS:
        rs = [r for r in sub if r['pair'] == pr]
        if not rs:
            continue
        xv = np.array([idx[int(r['pos_spike'])] + OFF[pr] for r in rs])
        md = np.array([F(r, 'median_residual') for r in rs])
        lo = np.array([F(r, 'ci_low') for r in rs])
        hi = np.array([F(r, 'ci_high') for r in rs])
        sg = np.array([F(r, 'p_bh') < 0.05 for r in rs])
        ax.errorbar(xv, md, yerr=[md - lo, hi - md], fmt='none', lw=0.7,
                    ecolor=np.where(sg, VERM, GREY), zorder=3)
        ax.scatter(xv[~sg], md[~sg], s=9, marker=MARK[pr], color=GREY, zorder=4)
        ax.scatter(xv[sg], md[sg], s=15, marker=MARK[pr], color=VERM, zorder=5)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(p) for p in pos], rotation=90, fontsize=5.6)
    ax.set_xlabel('spike position', labelpad=1)
    ax.set_ylabel('routing residual, log$_{10}$ $K_\\mathrm{A}$')
    ax.set_xlim(-0.7, len(pos) - 0.3)
    ax.legend(handles=[Line2D([], [], color=VERM, marker='o', ls='', ms=3.6,
                              label='departs'),
                       Line2D([], [], color=GREY, marker='o', ls='', ms=3.0,
                              label='does not'),
                       Line2D([], [], color=PALE, lw=5, label='null spread')],
              frameon=False, loc='upper left', handlelength=0.9,
              borderpad=0.1, labelspacing=0.25, fontsize=5.8)

    ax = AX['b']()
    panel(ax, 'b')
    pos, idx = slots(rnk)
    xs = np.arange(len(pos))
    p05 = {}
    for r in rnk:
        p05.setdefault(idx[int(r['pos_spike'])], []).append(F(r, 'null_p05_rho'))
    band = np.array([np.mean(p05[i]) for i in xs])
    ax.fill_between(xs, band, 1.0, color=PALE, lw=0, zorder=1)
    ax.axhline(1.0, color=INK, lw=0.6, ls=':', zorder=2)
    for pr in PAIRS:
        rs = [r for r in rnk if r['pair'] == pr]
        if not rs:
            continue
        xv = np.array([idx[int(r['pos_spike'])] + OFF[pr] for r in rs])
        rv = np.array([F(r, 'rho') for r in rs])
        sg = np.array([F(r, 'p_bh') < 0.05 for r in rs])
        ax.scatter(xv[~sg], rv[~sg], s=9, marker=MARK[pr], color=GREY, zorder=4)
        ax.scatter(xv[sg], rv[sg], s=20, marker=MARK[pr], color=VERM, zorder=5)
        for i in np.where(sg)[0]:
            ax.annotate(rs[i]['pos_spike'], (xv[i], rv[i]),
                        textcoords='offset points', xytext=(0, -10),
                        ha='center', fontsize=6.0, color=VERM)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(p) for p in pos], rotation=90, fontsize=5.6)
    ax.set_xlabel('spike position', labelpad=1)
    ax.set_ylabel('rank correlation across targets')
    ax.set_xlim(-0.7, len(pos) - 0.3)
    ax.set_ylim(-0.5, 1.14)
    ax.legend(handles=[Line2D([], [], color=GREY, marker=MARK[p], ls='', ms=3.2,
                              label=p.replace('|', ' vs ')) for p in PAIRS],
              frameon=False, loc='lower right', handlelength=0.9,
              borderpad=0.1, labelspacing=0.25, fontsize=5.6)

    ax = AX['c']()
    panel(ax, 'c')
    w = sorted([r for r in res if r['phenotype'] == 'binding'
                and r['pair'] == 'Wuhan|BA.1' and r['pos_spike'] == '498'],
               key=lambda r: r['res_to'])
    if not w:
        die('the worked case at spike 498 is absent from residuals.csv')
    xs = np.arange(len(w))
    one = np.array([F(r, 'one_step') for r in w])
    leg1 = np.array([F(r, 'first_leg') for r in w])
    leg2 = np.array([F(r, 'second_leg') for r in w])
    ax.plot(xs, one, 'o-', ms=3.4, lw=1.0, color=BLUE, label='one step, Q to Y')
    ax.plot(xs, leg1 + leg2, 's-', ms=3.2, lw=1.0, color=VERM,
            label='first leg + second leg, Q to R to Y')
    ax.fill_between(xs, leg1 + leg2, one, color=PALE, lw=0, zorder=0)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(xs)
    ax.set_xticklabels([r['res_to'] for r in w], fontsize=6.2)
    ax.set_xlabel('target residue Y', labelpad=1)
    ax.set_ylabel('effect, log$_{10}$ $K_\\mathrm{A}$')
    ax.set_xlim(-0.5, len(w) - 0.5)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.42 * (hi - lo))
    ax.legend(frameon=False, loc='upper center', ncol=2, handlelength=1.4,
              borderpad=0.1, labelspacing=0.3, columnspacing=1.6)
    ax.text(0.99, 0.04, 'spike 498, Wuhan against BA.1, binding',
            transform=ax.transAxes, ha='right', fontsize=6.4, color=INK)
    return save(fig, 'panel_rbd') if own else None


def fig_magnitude(AX=None, letters=('a', 'b', 'c')):
    cal = jload('calibrated_test_summary.json')
    sf = jload('census_test_scalefree_summary.json')
    by = sorted(cal['by_threshold'], key=lambda t: t['threshold'])
    own = AX is None
    if own:
        fig = plt.figure(figsize=(6.4, 2.1))
        gs = fig.add_gridspec(1, 3, wspace=0.42)
        AX = {k: (lambda i=i: fig.add_subplot(gs[0, i]))
              for i, k in enumerate(letters)}

    ax = AX[letters[0]]()
    panel(ax, letters[0])
    th = [t['threshold'] for t in by]
    ax.plot(th, [t['median_abs_residual'] for t in by], 'o-', ms=3.6, lw=1.1,
            color=VERM, label='observed residual')
    ax.plot(th, [t['median_residual_null'] for t in by], 's-', ms=3.4, lw=1.1,
            color=BLUE, label='its null')
    ax.set_yscale('log')
    ax.set_yticks([0.1, 0.2, 0.5, 1.0, 2.0])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel('fit correlation admitted')
    ax.set_ylabel('median absolute value')
    ax.legend(frameon=False, loc='center left', handlelength=1.3,
              borderpad=0.1, labelspacing=0.3, fontsize=6.0)
    r = by[2]['median_abs_residual'] / by[2]['median_residual_null']
    ax.annotate('%.0f fold apart' % r, xy=(th[2], by[2]['median_abs_residual']),
                textcoords='offset points', xytext=(2, 9), fontsize=6.2,
                color=INK)

    ax = AX[letters[1]]()
    panel(ax, letters[1])
    q = sorted(sf['quantile_sweep'], key=lambda r: r['quantile'])
    qq = [r['quantile'] for r in q]
    fr = [100 * r['significant_fraction'] for r in q]
    ax.plot(qq, fr, 'o-', ms=4.0, lw=1.2, color=VERM)
    ax.set_xlabel('centring quantile')
    ax.set_ylabel('cells departing, per cent')
    YTOP = max(fr) * 1.22
    ax.set_ylim(0, YTOP)
    for a, b2 in zip(qq, fr):
        ax.annotate('%.1f' % b2, (a, b2), textcoords='offset points',
                    xytext=(0, 5), ha='center', fontsize=5.8, color=INK)
    ax.text(0.5, 0.06, 'the answer moves with the estimate',
            transform=ax.transAxes, ha='center', fontsize=6.0, color=INK)
    ax.set_ylabel('cells departing, per cent')

    ax = AX[letters[2]]()
    panel(ax, letters[2])
    rk = [100.0 * t['rank_departing'] / t['dedup_cells'] for t in by]
    ax.plot(th, rk, 'o-', ms=4.0, lw=1.2, color=BLUE)
    ax.set_xlabel('fit correlation admitted')
    ax.set_ylabel('cells departing, per cent')
    ax.set_ylim(0, YTOP)
    for a, b2 in zip(th, rk):
        ax.annotate('%.2f' % b2, (a, b2), textcoords='offset points',
                    xytext=(0, 6), ha='center', fontsize=5.8, color=INK)
    ax.text(0.5, 0.88, 'the rank statistic stays\nbelow one per cent throughout',
            transform=ax.transAxes, ha='center', va='top', fontsize=6.0,
            color=INK)
    return save(fig, 'panel_magnitude') if own else None


def fig_domain_and_null():
    fig = plt.figure(figsize=(6.4, 6.3))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.06, 0.94, 0.86], hspace=0.60)
    top = gs[0].subgridspec(1, 2, wspace=0.26)
    bot = gs[2].subgridspec(1, 3, wspace=0.46)
    fig_rbd(AX={'a': lambda: fig.add_subplot(top[0, 0]),
                'b': lambda: fig.add_subplot(top[0, 1]),
                'c': lambda: fig.add_subplot(gs[1])})
    fig_magnitude(AX={'d': lambda: fig.add_subplot(bot[0, 0]),
                      'e': lambda: fig.add_subplot(bot[0, 1]),
                      'f': lambda: fig.add_subplot(bot[0, 2])},
                  letters=('d', 'e', 'f'))
    return save(fig, 'fig6')


if __name__ == '__main__':
    print('building figures at %d dpi' % DPI)
    fig_resource()
    fig_domain_and_null()
