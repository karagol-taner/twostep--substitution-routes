#!/usr/bin/env python3
import collections
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
_sp = importlib.util.spec_from_file_location(
    'rk', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       '12_census_test_rank.py'))
rk = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rk)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(
    os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE, 'data')
DRAWS = 10000
SEED = 20261111
MIN_TARGETS = 6
MIN_NULL = 8
COMMON = {'same experiment', 'same experiment set', 'same publication'}
AA20 = rk.AA20
LOG = []


def say(s):
    print(s, flush=True)
    LOG.append(s)


def mag_pool(sa, sb, pa, pb):
    out = []
    for p in range(1, min(len(pa), len(pb)) + 1):
        x = pa[p - 1]
        if x != pb[p - 1]:
            continue
        for y in sorted(AA20):
            if y == x:
                continue
            u, v = sa.get((p, y)), sb.get((p, y))
            if u is not None and v is not None:
                out.append(v - u)
    return np.asarray(out, dtype=float)


def main():
    seq = rk.load_sequences()
    tier = {}
    for r in csv.DictReader(open(os.path.join(DATA, 'common_scale_pairs.csv'))):
        tier[(r['from_urn'], r['via_urn'])] = r['tier']
    rows = [r for r in csv.DictReader(open(os.path.join(DATA, 'route_census_mavedb.csv')))
            if int(r['complete_routes']) >= MIN_TARGETS
            and tier.get((r['from_urn'], r['via_urn'])) in COMMON]
    by = collections.defaultdict(list)
    for r in rows:
        by[(r['from_urn'], r['via_urn'])].append(r)
    say('common-scale cells: %d over %d deposit pairs' % (len(rows), len(by)))

    sc = rk.Scores(limit=240)
    rng = np.random.default_rng(SEED)
    nulls = rk.NullCache(DRAWS, rng)
    out = []
    aside = collections.Counter()
    for (ua, ub), items in sorted(by.items()):
        pa, pb = seq.get(ua), seq.get(ub)
        if pa is None or pb is None or len(pa) != len(pb):
            aside['no comparable target sequence'] += len(items)
            continue
        sa, sb = sc.get(ua), sc.get(ub)
        if not sa or not sb:
            aside['deposit unreadable'] += len(items)
            continue
        rpool = rk.null_vectors(sa, sb, pa, pb, MIN_TARGETS)
        mpool = mag_pool(sa, sb, pa, pb)
        if len(rpool) < MIN_NULL or mpool.size < 100:
            aside['null pool too small'] += len(items)
            continue
        cells = []
        for r in items:
            pos, x, z = int(r['pos']), r['res_from'], r['res_via']
            if pos > len(pa) or pa[pos - 1] != x or pb[pos - 1] != z:
                aside['census row disagrees with the target sequence'] += 1
                continue
            va, vb, tgt, first = rk.cell_vectors(sa, sb, pos, x, z)
            if va.size < MIN_TARGETS or first is None:
                aside['too few shared targets or no first leg'] += 1
                continue
            n = int(va.size)
            ref = nulls.get((ua, ub), rpool, n)
            if ref is None:
                aside['no rank null at this number of targets'] += 1
                continue
            rho = float(rk.spearman(va[None, :], vb[None, :])[0])
            if not np.isfinite(rho):
                aside['rank correlation undefined'] += 1
                continue
            p_rank = (float(np.searchsorted(ref, rho, side='right')) + 1.0) / (ref.size + 1.0)
            res = va - (first + vb)
            obs = float(np.median(np.abs(res)))
            idx = rng.integers(0, mpool.size, (DRAWS, n))
            mref = np.sort(np.median(np.abs(mpool[idx]), axis=1))
            p_mag = (float(mref.size - np.searchsorted(mref, obs, side='left')) + 1.0) / (mref.size + 1.0)
            cells.append({
                'from_urn': ua, 'via_urn': ub, 'tier': tier[(ua, ub)],
                'from_name': r['from_name'][:60], 'via_name': r['via_name'][:60],
                'pos': pos, 'res_from': x, 'res_via': z,
                'via_is_alanine': int(z == 'A'), 'n_targets': n,
                'rho': round(rho, 5),
                'rank_null_median': round(float(np.median(ref)), 5),
                'p_rank': round(p_rank, 6),
                'median_residual': round(float(np.median(res)), 5),
                'median_abs_residual': round(obs, 5),
                'mag_null_median': round(float(np.median(mref)), 5),
                'p_mag': round(p_mag, 6),
                'intermediate_cost': round(float(first), 5)})
        if not cells:
            continue
        for key, pk in (('p_rank_bh', 'p_rank'), ('p_mag_bh', 'p_mag')):
            for c, a in zip(cells, rk.bh([c[pk] for c in cells])):
                c[key] = round(a, 6)
        out.extend(cells)

    say('cells tested: %d' % len(out))
    for k, v in sorted(aside.items()):
        say('  set aside, %-50s %d' % (k, v))

    sid, best = {}, {}
    for c in out:
        ka = sid.setdefault(seq[c['from_urn']], len(sid))
        kb = sid.setdefault(seq[c['via_urn']], len(sid))
        k = (ka, kb, c['pos'], c['res_from'], c['res_via'])
        if k not in best or c['p_mag'] > best[k]['p_mag']:
            best[k] = c
    ded = list(best.values())
    say('deduplicated cells: %d' % len(ded))
    say('')

    def rep(label, rows_):
        if not rows_:
            return {'label': label, 'cells': 0}
        mid = lambda v: sorted(v)[len(v) // 2]
        d = {'label': label, 'cells': len(rows_),
             'rank_departing': sum(1 for c in rows_ if c['p_rank_bh'] < 0.05),
             'magnitude_departing': sum(1 for c in rows_ if c['p_mag_bh'] < 0.05),
             'both': sum(1 for c in rows_ if c['p_rank_bh'] < 0.05 and c['p_mag_bh'] < 0.05),
             'median_rho': round(mid([c['rho'] for c in rows_]), 4),
             'median_rank_null': round(mid([c['rank_null_median'] for c in rows_]), 4),
             'median_abs_residual': round(mid([c['median_abs_residual'] for c in rows_]), 4),
             'median_mag_null': round(mid([c['mag_null_median'] for c in rows_]), 4),
             'median_intermediate_cost': round(mid([c['intermediate_cost'] for c in rows_]), 4)}
        say('  %-34s cells %4d  rank %3d  magnitude %3d  both %3d  '
            'rho %7.4f vs null %7.4f  |res| %7.4f vs null %7.4f'
            % (label, d['cells'], d['rank_departing'], d['magnitude_departing'],
               d['both'], d['median_rho'], d['median_rank_null'],
               d['median_abs_residual'], d['median_mag_null']))
        return d

    summ = {'draws': DRAWS, 'seed': SEED, 'min_targets': MIN_TARGETS,
            'tiers_counted_as_common_scale': sorted(COMMON),
            'cells_before_dedup': len(out), 'set_aside': dict(aside),
            'note': ('the magnitude statistic is legitimate here because both '
                     'deposits of a pair come from one laboratory and one assay '
                     'and therefore share a scale; it is the statistic fixed in '
                     'the analysis plan and nothing about it is re-chosen')}
    summ['all'] = rep('all common-scale cells', ded)
    summ['alanine'] = rep('alanine as intermediate', [c for c in ded if c['via_is_alanine']])
    sysgroup = collections.defaultdict(list)
    for c in ded:
        base = c['from_name'].split('(')[0].strip()
        key = ('influenza NP' if base.startswith('NP') else
               'homeobox domain' if 'Homeobox' in base else
               'SARS-CoV-2 RBD' if 'receptor binding' in base else
               'MPL' if base.endswith('MPL') or base == 'MPL' else
               'zinc finger' if 'Zinc finger' in base else base[:28])
        sysgroup[key].append(c)
    summ['by_system'] = []
    say('')
    for k in sorted(sysgroup, key=lambda k: -len(sysgroup[k])):
        if len(sysgroup[k]) >= 8:
            summ['by_system'].append(rep(k, sysgroup[k]))
    summ['systems_with_at_least_eight_cells'] = len(summ['by_system'])
    say('')
    say('systems carrying at least eight cells: %d' % len(summ['by_system']))

    fields = ['from_urn', 'via_urn', 'tier', 'from_name', 'via_name', 'pos',
              'res_from', 'res_via', 'via_is_alanine', 'n_targets', 'rho',
              'rank_null_median', 'p_rank', 'p_rank_bh', 'median_residual',
              'median_abs_residual', 'mag_null_median', 'p_mag', 'p_mag_bh',
              'intermediate_cost']
    with open(os.path.join(DATA, 'common_scale_test.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(ded, key=lambda c: c['p_mag_bh']))
    with open(os.path.join(DATA, 'common_scale_test_summary.json'), 'w') as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'common_scale_test_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    say('')
    say('written to common_scale_test.csv and common_scale_test_summary.json')


if __name__ == '__main__':
    main()
