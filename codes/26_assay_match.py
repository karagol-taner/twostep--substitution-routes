#!/usr/bin/env python3
import collections
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')


def die(m):
    sys.exit('FATAL 26_assay_match.py: %s' % m)


def norm(title, target):
    t = (title or '').lower()
    for w in re.split(r'[^a-z0-9]+', (target or '').lower()):
        if len(w) >= 3:
            t = t.replace(w, ' ')
    t = re.sub(r'\([^)]*\)', ' ', t)
    t = re.sub(r'[^a-z]+', ' ', t)
    return ' '.join(t.split())


def summ(v):
    v = np.asarray(v, dtype=float)
    if not len(v):
        return None
    return {'n': int(len(v)), 'median': round(float(np.median(v)), 4),
            'mean': round(float(v.mean()), 4),
            'q1': round(float(np.percentile(v, 25)), 4),
            'q3': round(float(np.percentile(v, 75)), 4),
            'above_0.9': int((v > 0.9).sum())}


def main():
    ap = os.path.join(DATA, 'deposit_agreement.json')
    pp = os.path.join(DATA, 'score_set_provenance.json')
    for p in (ap, pp):
        if not os.path.exists(p):
            die('%s is absent' % os.path.relpath(p, ROOT))
    agree = json.load(open(ap))
    prov = json.load(open(pp))['provenance']

    rec = agree['per_pair']
    out_rows = []
    for r in rec:
        a, b = prov.get(r['from']), prov.get(r['via'])
        if a is None or b is None:
            die('%s or %s carries no provenance record' % (r['from'], r['via']))
        na = norm(a.get('title'), r['target'])
        nb = norm(b.get('title'), r['target'])
        if a['experiment'] == b['experiment']:
            cls, why = 'same assay', 'one experiment'
        elif na and na == nb:
            cls, why = 'same assay', 'one score set title'
        elif na and nb:
            cls, why = 'different assay', 'titles differ'
        else:
            cls, why = 'undetermined', 'a title is empty once the target is removed'
        out_rows.append(dict(r, assay=cls, assay_basis=why,
                             title_from=a.get('title'), title_via=b.get('title')))

    by = collections.defaultdict(list)
    for r in out_rows:
        by[r['assay']].append(r)

    res = {'pairs': len(out_rows),
           'identity_predicts': 1.0,
           'rule': ['two deposits of one MaveDB experiment are one assay',
                    'two deposits whose score set titles agree once the target '
                    'name is removed are one assay',
                    'anything else is counted as a different assay, and a pair '
                    'whose title is empty after the removal is undetermined'],
           'by_assay': {}}
    for k in ('same assay', 'different assay', 'undetermined'):
        rs = by.get(k, [])
        res['by_assay'][k] = {
            'pairs': len(rs),
            'targets': len(sorted({x['target'] for x in rs})),
            'rho_test': summ([x['rho_test'] for x in rs]),
            'rho_control': summ([x['rho_control'] for x in rs])}
    diff = collections.Counter(
        (x['title_from'], x['title_via']) for x in by.get('different assay', []))
    res['different_assay_examples'] = [
        {'from': a, 'via': b, 'pairs': n} for (a, b), n in diff.most_common(10)]
    res['same_assay_basis'] = dict(collections.Counter(
        x['assay_basis'] for x in by.get('same assay', [])))
    oneexp = [x['rho_test'] for x in out_rows
              if x['assay_basis'] == 'one experiment']
    res['one_experiment_pairs'] = len(oneexp)
    res['one_experiment_median_rho'] = (round(float(np.median(oneexp)), 4)
                                        if oneexp else None)
    bytar = collections.defaultdict(list)
    for r in out_rows:
        for nm in sorted({r['target'], r.get('target_via') or r['target']}):
            bytar[nm].append(r)
    tab = []
    for t, rs in bytar.items():
        rt = np.array([x['rho_test'] for x in rs])
        tab.append({'target': t, 'pairs': len(rs),
                    'median_rho_test': round(float(np.median(rt)), 4),
                    'above_0.9': int((rt > 0.9).sum()),
                    'median_rho_control':
                        round(float(np.median([x['rho_control'] for x in rs])), 4),
                    'same_assay': sum(1 for x in rs if x['assay'] == 'same assay')})
    tab.sort(key=lambda r: r['median_rho_test'])
    res['by_target'] = tab
    res['targets_failing_0.9'] = [r['target'] for r in tab
                                  if r['median_rho_test'] < 0.9]
    res['targets_passing_0.9'] = [r['target'] for r in tab
                                  if r['median_rho_test'] >= 0.9]

    res['per_pair'] = [{'from': x['from'], 'via': x['via'], 'pos': x['pos'],
                        'target': x['target'], 'assay': x['assay'],
                        'assay_basis': x['assay_basis'],
                        'rho_test': x['rho_test'],
                        'rho_control': x['rho_control']} for x in out_rows]

    with open(os.path.join(DATA, 'assay_match.json'), 'w') as fh:
        json.dump(res, fh, indent=1, sort_keys=True)

    print('pairs at one differing position, where the rank correlation must '
          'return 1: %d' % len(out_rows))
    print()
    print('%-16s %6s %8s %10s %10s %12s'
          % ('', 'pairs', 'targets', 'median', 'q1', 'above 0.9'))
    for k in ('same assay', 'different assay', 'undetermined'):
        v = res['by_assay'][k]
        if v['rho_test']:
            print('%-16s %6d %8d %10.4f %10.4f %6d of %d'
                  % (k, v['pairs'], v['targets'], v['rho_test']['median'],
                     v['rho_test']['q1'], v['rho_test']['above_0.9'], v['pairs']))
        else:
            print('%-16s %6d' % (k, v['pairs']))
    print()
    print('same assay decided by: %s' % res['same_assay_basis'])
    print()
    print('the pairs the test rejects, by what the two deposits measure:')
    for e in res['different_assay_examples'][:8]:
        print('  %3d  %-44s || %s' % (e['pairs'], (e['from'] or '')[:44],
                                      (e['via'] or '')[:44]))
    print()
    print('%-38s %6s %10s %10s' % ('target', 'pairs', 'median', 'above 0.9'))
    for r in res['by_target']:
        print('%-38s %6d %10.4f %6d' % (r['target'][:38], r['pairs'],
                                        r['median_rho_test'], r['above_0.9']))
    print()
    print('targets whose median falls below 0.9: %d of %d  (%s)'
          % (len(res['targets_failing_0.9']), len(tab),
             ', '.join(res['targets_failing_0.9'])))
    print('wrote data/assay_match.json')


if __name__ == '__main__':
    main()
