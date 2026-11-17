#!/usr/bin/env python3
import collections
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
DRAWS = 10000
BOOT = 10000
SEED = 20261111


def die(m):
    sys.exit('FATAL 06_additivity_test.py: %s' % m)


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    out = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        k = n - rank + 1
        prev = min(prev, ps[i] * n / float(k))
        out[i] = prev
    return out


def main():
    log = []

    def say(s):
        print(s)
        log.append(s)

    for f in ('residuals.csv', 'controls.csv', 'censoring_sweep.json'):
        if not os.path.isfile(os.path.join(DATA, f)):
            die('data/%s is absent; run scripts 04 and 05 first' % f)
    res = list(csv.DictReader(open(os.path.join(DATA, 'residuals.csv'))))
    ctl = list(csv.DictReader(open(os.path.join(DATA, 'controls.csv'))))
    margins = json.load(open(os.path.join(DATA, 'censoring_sweep.json')))['margins']

    rng = np.random.default_rng(SEED)
    refcache = {}
    rows, summary = [], []
    for margin in margins:
        R = [r for r in res if float(r['margin_to_clip']) > margin]
        C = [c for c in ctl if float(c['margin_to_clip']) > margin]
        cells = collections.defaultdict(list)
        for r in R:
            cells[(r['pair'], r['phenotype'], int(r['pos_rbd']))].append(r)
        pool = collections.defaultdict(list)
        for c in C:
            pool[(c['pair'], c['phenotype'])].append(float(c['discrepancy']))
        per_ph = collections.defaultdict(list)
        for key, items in sorted(cells.items()):
            pair, ph, pos = key
            d = np.array([float(x['residual']) for x in items])
            n = len(d)
            null = np.asarray(pool[(pair, ph)])
            if n == 0 or null.size < n:
                continue
            obs = float(np.median(np.abs(d)))
            rk = (margin, pair, ph, n)
            ref = refcache.get(rk)
            if ref is None:
                idx = np.argpartition(rng.random((DRAWS, null.size)),
                                      n - 1, axis=1)[:, :n]
                ref = np.median(np.abs(null[idx]), axis=1)
                refcache[rk] = ref
            p = float((np.sum(ref >= obs) + 1) / (DRAWS + 1.0))
            bs = np.median(d[rng.integers(0, n, size=(BOOT, n))], axis=1)
            lo, hi = np.percentile(bs, [2.5, 97.5])
            row = {'margin': margin, 'pair': pair, 'phenotype': ph,
                   'pos_rbd': pos, 'pos_spike': int(items[0]['pos_spike']),
                   'res_from': items[0]['res_from'], 'res_via': items[0]['res_via'],
                   'n_routes': n,
                   'median_residual': round(float(np.median(d)), 4),
                   'ci_low': round(float(lo), 4), 'ci_high': round(float(hi), 4),
                   'median_abs_residual': round(obs, 4),
                   'null_median_abs': round(float(np.median(ref)), 4),
                   'null_p95_abs': round(float(np.percentile(ref, 95)), 4),
                   'n_control_pool': int(null.size),
                   'intermediate_cost': round(float(items[0]['first_leg']), 4),
                   'p_perm': round(p, 5)}
            rows.append(row)
            per_ph[(margin, ph)].append(row)
        for (m, ph), rs in per_ph.items():
            adj = bh([r['p_perm'] for r in rs])
            for r, q in zip(rs, adj):
                r['p_bh'] = round(q, 5)
            sig = sum(1 for q in adj if q < 0.05)
            summary.append({'margin': m, 'phenotype': ph, 'cells': len(rs),
                            'significant_bh_0.05': sig,
                            'median_abs_residual_over_cells':
                                round(float(np.median([r['median_abs_residual'] for r in rs])), 4),
                            'median_null_over_cells':
                                round(float(np.median([r['null_median_abs'] for r in rs])), 4),
                            'median_intermediate_cost':
                                round(float(np.median([r['intermediate_cost'] for r in rs])), 4)})

    if not rows:
        die('no testable cells were produced')
    fields = ['margin', 'pair', 'phenotype', 'pos_rbd', 'pos_spike', 'res_from',
              'res_via', 'n_routes', 'median_residual', 'ci_low', 'ci_high',
              'median_abs_residual', 'null_median_abs', 'null_p95_abs',
              'n_control_pool', 'intermediate_cost', 'p_perm', 'p_bh']
    with open(os.path.join(DATA, 'additivity_test.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(DATA, 'additivity_summary.json'), 'w') as fh:
        json.dump({'draws': DRAWS, 'bootstrap': BOOT, 'seed': SEED,
                   'by_margin_and_phenotype': summary}, fh, indent=1, sort_keys=True)

    say('%-7s %-11s %5s %9s %22s %22s %s'
        % ('margin', 'phenotype', 'cells', 'BH<0.05', 'median |residual|',
           'median null |disc|', 'median intermediate cost'))
    for s in sorted(summary, key=lambda s: (s['margin'], s['phenotype'])):
        say('%-7.2f %-11s %5d %9d %22.4f %22.4f %24.4f'
            % (s['margin'], s['phenotype'], s['cells'], s['significant_bh_0.05'],
               s['median_abs_residual_over_cells'], s['median_null_over_cells'],
               s['median_intermediate_cost']))
    with open(os.path.join(DATA, 'additivity_test_log.txt'), 'w') as fh:
        fh.write('\n'.join(log) + '\n')
    print('\nwritten to data/additivity_test.csv and data/additivity_summary.json')


if __name__ == '__main__':
    main()
