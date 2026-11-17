#!/usr/bin/env python3
import collections
import csv
import importlib.util
import json
import math
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
MIN_TARGETS = 6
N_DRAW = 10
MAX_POOL = 400
SEED = 20261112


def die(m):
    sys.exit('FATAL 28_statistic_choice.py: %s' % m)


def _rank(v):
    v = np.asarray(v, dtype=float)
    o = v.argsort(kind='mergesort')
    r = np.empty(len(v), dtype=float)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[o[j + 1]] == v[o[i]]:
            j += 1
        r[o[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return r


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    da, db = a - a.mean(), b - b.mean()
    d = (da * da).sum() * (db * db).sum()
    return None if d <= 0 else float((da * db).sum() / math.sqrt(d))


def spearman(a, b):
    return pearson(_rank(a), _rank(b))


def kendall(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = len(a)
    c = d = ta = tb = 0
    for i in range(n - 1):
        da = a[i + 1:] - a[i]
        db = b[i + 1:] - b[i]
        s = np.sign(da) * np.sign(db)
        c += int((s > 0).sum())
        d += int((s < 0).sum())
        ta += int((da == 0).sum())
        tb += int((db == 0).sum())
    n0 = n * (n - 1) / 2.0
    den = math.sqrt(max(1e-12, (n0 - ta) * (n0 - tb)))
    return None if den <= 0 else (c - d) / den


def summ(v):
    v = np.asarray([x for x in v if x is not None], float)
    if not len(v):
        return None
    return {'n': int(len(v)), 'median': round(float(np.median(v)), 4),
            'mean': round(float(v.mean()), 4),
            'q1': round(float(np.percentile(v, 25)), 4),
            'q3': round(float(np.percentile(v, 75)), 4),
            'above_0.9': int((v > 0.9).sum())}


def main():
    spec = importlib.util.spec_from_file_location(
        'census', os.path.join(HERE, '08_mavedb_route_census.py'))
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)

    seq = {}
    for r in json.load(open(os.path.join(DATA, 'target_genes.json'))):
        ts = r.get('targetSequence') or {}
        s = (ts.get('sequence') or '').upper()
        t = ts.get('sequenceType')
        p = (C.translate(s) if (t == 'dna' and len(s) % 3 == 0)
             else (s if t == 'protein' else None))
        if p:
            seq[r['scoreSetUrn']] = p.rstrip('*')

    rows = [r for r in csv.DictReader(
        open(os.path.join(DATA, 'route_census_mavedb.csv')))
        if int(r['complete_routes']) >= 1]
    one = []
    for r in rows:
        a, b = r['from_urn'], r['via_urn']
        if a not in seq or b not in seq or len(seq[a]) != len(seq[b]):
            continue
        d = [i for i, (x, y) in enumerate(zip(seq[a], seq[b])) if x != y]
        if len(d) == 1 and d[0] == int(r['pos']) - 1:
            one.append((r, d[0]))
    if not one:
        die('no deposit pair differs at exactly one position')
    print('pairs at one differing position: %d' % len(one))
    seen_un = set()
    dedup = []
    for r, p0 in sorted(one, key=lambda t: (t[0]['from_urn'], t[0]['via_urn'],
                                            int(t[0]['pos']))):
        a, b = r['from_urn'], r['via_urn']
        k = (min(a, b), max(a, b), p0)
        if k in seen_un:
            continue
        seen_un.add(k)
        dedup.append((r, p0))
    print('  of those, independent comparisons after collapsing direction: %d'
          % len(dedup))
    one_directed = len(one)
    one = dedup

    scores = {}
    urns = sorted({u for r, _ in one for u in (r['from_urn'], r['via_urn'])})
    for i, u in enumerate(urns, 1):
        m, _ = C.load_scores(u)
        if m is not None:
            scores[u] = m
        if i % 150 == 0 or i == len(urns):
            print('  read %d of %d' % (i, len(urns)))

    rng = random.Random(SEED)
    rec = []
    for r, p0 in one:
        A, B = scores.get(r['from_urn']), scores.get(r['via_urn'])
        if A is None or B is None:
            continue
        pos = p0 + 1
        sa, sb = seq[r['from_urn']], seq[r['via_urn']]
        X, Z = sa[p0], sb[p0]
        ys = [y for y in sorted(C.AA20) if y not in (X, Z)
              and (pos, y) in A and (pos, y) in B]
        if len(ys) < MIN_TARGETS:
            continue
        a = np.array([A[(pos, y)] for y in ys])
        b = np.array([B[(pos, y)] for y in ys])
        if a.std() <= 0 or b.std() <= 0:
            continue
        k = len(ys)
        t_rho, t_r, t_tau = spearman(a, b), pearson(a, b), kendall(a, b)
        if None in (t_rho, t_r, t_tau):
            continue

        shared = [q for q in range(1, min(len(sa), len(sb)) + 1)
                  if q != pos and sa[q - 1] == sb[q - 1]]
        rng.shuffle(shared)
        pool = []
        for q in shared[:MAX_POOL]:
            yq = [y for y in sorted(C.AA20) if y != sa[q - 1]
                  and (q, y) in A and (q, y) in B]
            if len(yq) < k:
                continue
            pick = yq if len(yq) == k else rng.sample(yq, k)
            aq = np.array([A[(q, y)] for y in pick])
            bq = np.array([B[(q, y)] for y in pick])
            if aq.std() <= 0 or bq.std() <= 0:
                continue
            pool.append((q, aq, bq))
        if len(pool) < N_DRAW:
            continue
        la, lb = math.log(a.std()), math.log(b.std())
        near = sorted(pool, key=lambda p: (math.log(p[1].std()) - la) ** 2
                      + (math.log(p[2].std()) - lb) ** 2)[:N_DRAW]
        c_rho = float(np.median([spearman(x[1], x[2]) for x in near]))
        c_r = float(np.median([pearson(x[1], x[2]) for x in near]))
        c_tau = float(np.median([kendall(x[1], x[2]) for x in near]))

        rec.append({'target': r['from_name'].split('(')[0].strip(),
                    'target_via': r['via_name'].split('(')[0].strip(),
                    'from': r['from_urn'], 'via': r['via_urn'], 'pos': pos,
                    'targets': k,
                    'rho': round(t_rho, 5), 'r': round(t_r, 5),
                    'tau': round(t_tau, 5),
                    'rho_control': round(c_rho, 5),
                    'r_control': round(c_r, 5),
                    'tau_control': round(c_tau, 5)})

    if not rec:
        die('no pair yielded all three statistics')

    rho = np.array([x['rho'] for x in rec])
    rr = np.array([x['r'] for x in rec])
    tau = np.array([x['tau'] for x in rec])
    gap = rho - np.array([x['rho_control'] for x in rec])
    gapr = rr - np.array([x['r_control'] for x in rec])

    def ranksum(x, y):
        nx, ny = len(x), len(y)
        if nx < 3 or ny < 3:
            return None
        allv = np.concatenate([np.asarray(x, float), np.asarray(y, float)])
        r = _rank(allv)
        rx = r[:nx].sum()
        u = rx - nx * (nx + 1) / 2.0
        mu = nx * ny / 2.0
        _, counts = np.unique(allv, return_counts=True)
        tie = float(((counts ** 3 - counts).sum()))
        n = nx + ny
        var = nx * ny / 12.0 * ((n + 1) - tie / (n * (n - 1.0)))
        if var <= 0:
            return None
        z = (u - mu) / math.sqrt(var)
        return math.erfc(abs(z) / math.sqrt(2.0))

    def sign_test(diffs):
        pos = sum(1 for d in diffs if d > 0)
        neg = sum(1 for d in diffs if d < 0)
        n = pos + neg
        if n == 0:
            return None
        k = min(pos, neg)
        tail = sum(math.comb(n, i) for i in range(0, k + 1)) / float(2 ** n)
        return min(1.0, 2.0 * tail)

    bytar = collections.defaultdict(list)
    for x in rec:
        for nm in sorted({x['target'], x['target_via']}):
            bytar[nm].append(x)
    tab = []
    for t, rs in sorted(bytar.items()):
        g = np.array([x['rho'] - x['rho_control'] for x in rs])
        pv = sign_test(list(g))
        other = [x['rho'] - x['rho_control'] for x in rec
                 if t not in (x['target'], x['target_via'])]
        pw = ranksum(list(g), other)
        tab.append({'target': t, 'pairs': len(rs),
                    'pairs_with_a_negative_difference':
                        int(sum(1 for x in g if x < 0)),
                    'sign_test_p': None if pv is None else round(pv, 6),
                    'ranksum_p': None if pw is None else float('%.3g' % pw),
                    'median_rho': round(float(np.median([x['rho'] for x in rs])), 4),
                    'median_r': round(float(np.median([x['r'] for x in rs])), 4),
                    'median_tau': round(float(np.median([x['tau'] for x in rs])), 4),
                    'median_gap_rho': round(float(np.median(g)), 4)})
    tab.sort(key=lambda x: x['median_gap_rho'])

    ps = sorted([(r['ranksum_p'], r['target']) for r in tab
                 if r['ranksum_p'] is not None])
    m = len(ps)
    thresh, passing = 0.0, []
    for i, (pv, t) in enumerate(ps, 1):
        if pv <= 0.05 * i / m:
            thresh, passing = pv, [x[1] for x in ps[:i]]
    bh = {r['target']: (r['ranksum_p'] is not None
                        and r['ranksum_p'] <= thresh and
                        r['median_gap_rho'] < 0)
          for r in tab}
    for r in tab:
        r['departs_after_bh'] = bool(bh.get(r['target']))

    out = {'proteins_involved': len(bytar),
           'bh_threshold_p': round(thresh, 6),
           'targets_departing_after_bh':
               [r['target'] for r in tab if r['departs_after_bh']],
           'pairs_measured': len(rec),
           'directed_rows_before_collapsing': one_directed,
           'distinct_targets': len(bytar),
           'comparisons': len(rec),
           'min_targets': MIN_TARGETS, 'control_positions_per_null': N_DRAW,
           'seed': SEED,
           'all_three_predicted': 1.0,
           'spearman': summ(rho), 'pearson': summ(rr), 'kendall': summ(tau),
           'spearman_control': summ([x['rho_control'] for x in rec]),
           'pearson_control': summ([x['r_control'] for x in rec]),
           'kendall_control': summ([x['tau_control'] for x in rec]),
           'gap_spearman': summ(gap), 'gap_pearson': summ(gapr),
           'pairs_passing_rho_failing_r':
               int(((rho > 0.9) & (rr <= 0.9)).sum()),
           'pairs_passing_r_failing_rho':
               int(((rr > 0.9) & (rho <= 0.9)).sum()),
           'targets_with_a_negative_gap':
               [x['target'] for x in tab if x['median_gap_rho'] < 0],
           'by_target': tab,
           'per_pair': rec}
    with open(os.path.join(DATA, 'statistic_choice.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print()
    print('%-26s %8s %8s %8s' % ('', 'median', 'q1', 'above 0.9'))
    for lab, v in (('Spearman, test', out['spearman']),
                   ('Pearson, test', out['pearson']),
                   ('Kendall, test', out['kendall']),
                   ('Spearman, control', out['spearman_control']),
                   ('Pearson, control', out['pearson_control']),
                   ('Kendall, control', out['kendall_control'])):
        print('%-26s %8.4f %8.4f %8d' % (lab, v['median'], v['q1'], v['above_0.9']))
    print()
    print('pairs that pass on Spearman and fail on Pearson: %d of %d'
          % (out['pairs_passing_rho_failing_r'], len(rec)))
    print('pairs that pass on Pearson and fail on Spearman: %d'
          % out['pairs_passing_r_failing_rho'])
    print('paired gap, test minus its dispersion-matched control: '
          'Spearman median %+.4f, Pearson median %+.4f'
          % (out['gap_spearman']['median'], out['gap_pearson']['median']))
    print()
    print('proteins whose test positions agree less well than their own '
          'control positions (%d of %d):'
          % (len(out['targets_with_a_negative_gap']), len(tab)))
    for x in tab[:12]:
        print('  %-38s pairs %3d  rho %6.3f  gap %+6.3f'
              % (x['target'][:38], x['pairs'], x['median_rho'],
                 x['median_gap_rho']))
    print()
    print('paired difference against the rest of the census, '
          'Benjamini-Hochberg at 0.05 over %d proteins (threshold p %.2e):'
          % (len(tab), thresh))
    for r in tab:
        if r['departs_after_bh']:
            print('  %-38s pairs %3d, negative in %3d, rank-sum p %.2e, '
                  'sign-test p %s, gap %+.3f'
                  % (r['target'][:38], r['pairs'],
                     r['pairs_with_a_negative_difference'],
                     r['ranksum_p'],
                     'NA' if r['sign_test_p'] is None
                     else '%.2e' % r['sign_test_p'],
                     r['median_gap_rho']))
    print('wrote data/statistic_choice.json')


if __name__ == '__main__':
    main()
