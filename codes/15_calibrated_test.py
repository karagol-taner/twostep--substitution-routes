#!/usr/bin/env python3
import collections
import csv
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(
    os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE, 'data')
_sp = importlib.util.spec_from_file_location(
    'rk', os.path.join(HERE, '12_census_test_rank.py'))
rk = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rk)

DRAWS = 10000
SEED = 20261111
MIN_TARGETS = 6
MIN_FIT = 100
THRESHOLDS = [0.70, 0.75, 0.80, 0.85, 0.90]
AA20 = rk.AA20
LOG = []


def say(s):
    print(s, flush=True)
    LOG.append(s)


def system_of(name):
    n = (name or '').lower()
    for key, lab in (('homeobox', 'homeobox domain'), ('zinc finger', 'zinc finger'),
                     ('receptor binding domain', 'SARS-CoV-2 RBD'), ('np', 'influenza NP'),
                     ('mpl', 'MPL'), ('ubc', 'ubiquitin'), ('ubb', 'ubiquitin'),
                     ('spg', 'protein G'), ('ice-structuring', 'ice-structuring protein'),
                     ('cyp2c', 'CYP2C'), ('affibody', 'affibody')):
        if key in n:
            return lab
    return (name or 'unnamed').split('(')[0].strip()[:28] or 'unnamed'


def fit_map(sa, sb, pa, pb):
    xs, ys = [], []
    for p in range(1, min(len(pa), len(pb)) + 1):
        x = pa[p - 1]
        if x != pb[p - 1]:
            continue
        for y in sorted(AA20):
            if y == x:
                continue
            u, v = sa.get((p, y)), sb.get((p, y))
            if u is not None and v is not None:
                xs.append(u)
                ys.append(v)
    if len(xs) < MIN_FIT:
        return None
    X, Y = np.asarray(xs), np.asarray(ys)
    if np.std(X) < 1e-9 or np.std(Y) < 1e-9:
        return None
    sl, ic = np.polyfit(Y, X, 1)
    r = float(np.corrcoef(X, Y)[0, 1])
    return {'slope': float(sl), 'intercept': float(ic), 'pearson': r, 'n_fit': len(xs),
            'resid': (X - (sl * Y + ic))}


def bh(ps):
    n = len(ps)
    o = sorted(range(n), key=lambda i: ps[i])
    out = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(o), 1):
        k = n - rank + 1
        prev = min(prev, ps[i] * n / k)
        out[i] = prev
    return out


FIELDS = ['threshold', 'from_urn', 'via_urn', 'system', 'provenance_tier', 'pos',
          'res_from', 'res_via', 'via_is_alanine', 'n_targets', 'fit_pearson',
          'fit_slope', 'fit_n', 'rho', 'rho_null_median', 'rho_p_bh',
          'median_residual', 'median_abs_residual', 'residual_null_median',
          'intermediate_cost', 'mag_p_bh']


def main():
    seq = rk.load_sequences()
    say('deposits with a translatable target sequence: %d' % len(seq))
    tier = {}
    pth = os.path.join(DATA, 'common_scale_pairs.csv')
    if os.path.isfile(pth):
        for r in csv.DictReader(open(pth)):
            tier[(r['from_urn'], r['via_urn'])] = r['tier']
    rows = [r for r in csv.DictReader(open(os.path.join(DATA, 'route_census_mavedb.csv')))
            if int(r['complete_routes']) >= MIN_TARGETS]
    by_pair = collections.defaultdict(list)
    for r in rows:
        by_pair[(r['from_urn'], r['via_urn'])].append(r)
    say('deposit pairs in the census: %d' % len(by_pair))

    sc = rk.Scores(limit=200)
    rng = np.random.default_rng(SEED)
    fits, cells = {}, []
    aside = collections.Counter()
    done = 0
    for (ua, ub), items in sorted(by_pair.items()):
        done += 1
        if done % 500 == 0:
            say('  %d of %d pairs, %d cells so far' % (done, len(by_pair), len(cells)))
        pa, pb = seq.get(ua), seq.get(ub)
        if pa is None or pb is None or len(pa) != len(pb):
            aside['no comparable target sequence'] += len(items)
            continue
        sa, sb = sc.get(ua), sc.get(ub)
        if not sa or not sb:
            aside['deposit unreadable'] += len(items)
            continue
        fm = fit_map(sa, sb, pa, pb)
        if fm is None:
            aside['affine map not fittable'] += len(items)
            continue
        if fm['pearson'] < min(THRESHOLDS):
            aside['affine fit below the lowest threshold'] += len(items)
            continue
        fits[(ua, ub)] = fm['pearson']
        pool = rk.null_vectors(sa, sb, pa, pb, MIN_TARGETS)
        if len(pool) < rk.MIN_NULL_POSITIONS:
            aside['too few shared-parent positions'] += len(items)
            continue
        mpool = fm['resid']
        refs_rho, refs_mag = {}, {}
        got = []
        for r in items:
            pos, x, z = int(r['pos']), r['res_from'], r['res_via']
            if pos > len(pa) or pa[pos - 1] != x or pb[pos - 1] != z:
                aside['census row disagrees with the target sequence'] += 1
                continue
            va, vb, tgt, first = rk.cell_vectors(sa, sb, pos, x, z)
            if va.size < MIN_TARGETS or first is None:
                aside['too few shared targets or no first leg'] += 1
                continue
            n = va.size
            cal = fm['slope'] * vb + fm['intercept']
            resid = va - (first + cal)
            if n not in refs_rho:
                use = [(u, v) for u, v in pool if u.size >= n]
                if len(use) < rk.MIN_NULL_POSITIONS:
                    refs_rho[n] = None
                else:
                    pick = rng.integers(0, len(use), DRAWS)
                    A = np.empty((DRAWS, n)); B = np.empty((DRAWS, n))
                    for i, j in enumerate(pick):
                        u, v = use[j]
                        sel = slice(None) if u.size == n else rng.choice(u.size, n, replace=False)
                        A[i] = u[sel]; B[i] = v[sel]
                    refs_rho[n] = np.sort(rk.spearman(A, B))
                idx = rng.integers(0, mpool.size, (DRAWS, n))
                refs_mag[n] = np.sort(np.median(np.abs(mpool[idx]), axis=1))
            rref, mref = refs_rho[n], refs_mag[n]
            if rref is None:
                aside['no rank null at this number of targets'] += 1
                continue
            rho = float(rk.spearman(va[None, :], vb[None, :])[0])
            if not np.isfinite(rho):
                aside['rank correlation undefined'] += 1
                continue
            mar = float(np.median(np.abs(resid)))
            got.append({
                'from_urn': ua, 'via_urn': ub,
                'system': system_of(r['from_name']),
                'provenance_tier': tier.get((ua, ub), 'not labelled'),
                'pos': pos, 'res_from': x, 'res_via': z,
                'via_is_alanine': int(z == 'A'), 'n_targets': n,
                'fit_pearson': round(fm['pearson'], 4),
                'fit_slope': round(fm['slope'], 4), 'fit_n': fm['n_fit'],
                'rho': round(rho, 5),
                'rho_null_median': round(float(np.median(rref)), 5),
                '_rp': (float(np.searchsorted(rref, rho, side='right')) + 1.0) / (rref.size + 1.0),
                'median_residual': round(float(np.median(resid)), 5),
                'median_abs_residual': round(mar, 5),
                'residual_null_median': round(float(np.median(mref)), 5),
                'intermediate_cost': round(float(first), 5),
                '_mp': (float(mref.size - np.searchsorted(mref, mar, side='left')) + 1.0) / (mref.size + 1.0)})
        if got:
            cells.extend(got)

    say('cells with a calibrated test: %d' % len(cells))
    for k, v in sorted(aside.items()):
        say('  set aside, %-50s %d' % (k, v))
    say('')

    summ = {'draws': DRAWS, 'seed': SEED, 'min_targets': MIN_TARGETS,
            'min_fit_points': MIN_FIT, 'thresholds': THRESHOLDS,
            'calibration': 'least-squares affine map carrying the second deposit '
                           'onto the first, fitted on shared-parent positions only',
            'set_aside': dict(aside), 'by_threshold': []}
    out = []
    for t in THRESHOLDS:
        sub = [c for c in cells if c['fit_pearson'] >= t]
        if not sub:
            continue
        bypair = collections.defaultdict(list)
        for c in sub:
            bypair[(c['from_urn'], c['via_urn'])].append(c)
        for grp in bypair.values():
            for c, a in zip(grp, bh([c['_rp'] for c in grp])):
                c['rho_p_bh'] = round(a, 6)
            for c, a in zip(grp, bh([c['_mp'] for c in grp])):
                c['mag_p_bh'] = round(a, 6)
        sid, best = {}, {}
        for c in sub:
            ka = sid.setdefault(seq[c['from_urn']], len(sid))
            kb = sid.setdefault(seq[c['via_urn']], len(sid))
            k = (ka, kb, c['pos'], c['res_from'], c['res_via'])
            if k not in best or c['rho'] > best[k]['rho']:
                best[k] = c
        ded = list(best.values())
        rec = {'threshold': t, 'pairs': len(bypair), 'cells': len(sub),
               'dedup_cells': len(ded),
               'rank_departing': sum(1 for c in ded if c['rho_p_bh'] < 0.05),
               'magnitude_departing': sum(1 for c in ded if c['mag_p_bh'] < 0.05),
               'both': sum(1 for c in ded if c['rho_p_bh'] < 0.05 and c['mag_p_bh'] < 0.05),
               'median_rho': round(float(np.median([c['rho'] for c in ded])), 4),
               'median_rho_null': round(float(np.median([c['rho_null_median'] for c in ded])), 4),
               'median_abs_residual': round(float(np.median([c['median_abs_residual'] for c in ded])), 4),
               'median_residual_null': round(float(np.median([c['residual_null_median'] for c in ded])), 4),
               'median_intermediate_cost': round(float(np.median([c['intermediate_cost'] for c in ded])), 4),
               'systems': len({c['system'] for c in ded}),
               'alanine_cells': sum(1 for c in ded if c['via_is_alanine'])}
        by_sys = {}
        for c in ded:
            by_sys.setdefault(c['system'], []).append(c)
        rec['by_system'] = sorted(
            [{'system': s, 'cells': len(v),
              'rank_departing': sum(1 for c in v if c['rho_p_bh'] < 0.05),
              'magnitude_departing': sum(1 for c in v if c['mag_p_bh'] < 0.05),
              'median_rho': round(float(np.median([c['rho'] for c in v])), 4),
              'median_rho_null': round(float(np.median([c['rho_null_median'] for c in v])), 4),
              'median_abs_residual': round(float(np.median([c['median_abs_residual'] for c in v])), 4),
              'median_residual_null': round(float(np.median([c['residual_null_median'] for c in v])), 4)}
             for s, v in by_sys.items() if len(v) >= 8],
            key=lambda d: -d['cells'])
        summ['by_threshold'].append(rec)
        say('threshold %.2f  pairs %4d  cells %5d  dedup %5d  rank %4d  magnitude %4d  both %4d  '
            'rho %7.4f vs %7.4f  |res| %7.4f vs %7.4f  systems %d'
            % (t, rec['pairs'], rec['cells'], rec['dedup_cells'], rec['rank_departing'],
               rec['magnitude_departing'], rec['both'], rec['median_rho'],
               rec['median_rho_null'], rec['median_abs_residual'],
               rec['median_residual_null'], rec['systems']))
        for c in ded:
            d = {k: c[k] for k in FIELDS if k in c}
            d['threshold'] = t
            out.append(d)
    say('')
    for rec in summ['by_threshold']:
        if abs(rec['threshold'] - 0.80) < 1e-9:
            say('by system at threshold 0.80')
            say('  %-28s %6s %6s %6s %9s %9s %9s %9s'
                % ('system', 'cells', 'rank', 'magn', 'rho', 'rho null', '|res|', 'res null'))
            for s in rec['by_system']:
                say('  %-28s %6d %6d %6d %9.4f %9.4f %9.4f %9.4f'
                    % (s['system'], s['cells'], s['rank_departing'], s['magnitude_departing'],
                       s['median_rho'], s['median_rho_null'], s['median_abs_residual'],
                       s['median_residual_null']))
    with open(os.path.join(DATA, 'calibrated_test.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction='ignore')
        w.writeheader()
        w.writerows(out)
    with open(os.path.join(DATA, 'calibrated_test_summary.json'), 'w') as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'calibrated_test_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    say('')
    say('written to calibrated_test.csv and calibrated_test_summary.json')


if __name__ == '__main__':
    main()
