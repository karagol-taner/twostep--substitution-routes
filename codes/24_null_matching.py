#!/usr/bin/env python3
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
MAX_POOL = 400
N_DRAW = 10
SEED = 20261112


def die(m):
    sys.exit('FATAL 24_null_matching.py: %s' % m)


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


def spearman(a, b):
    if len(a) < MIN_TARGETS:
        return None
    ra, rb = _rank(a), _rank(b)
    da, db = ra - ra.mean(), rb - rb.mean()
    v = (da * da).sum() * (db * db).sum()
    return None if v <= 0 else float((da * db).sum() / v ** 0.5)


def summ(v):
    v = np.asarray(v, dtype=float)
    return {'n': int(len(v)), 'median': round(float(np.median(v)), 4),
            'mean': round(float(v.mean()), 4),
            'q1': round(float(np.percentile(v, 25)), 4),
            'q3': round(float(np.percentile(v, 75)), 4),
            'min': round(float(v.min()), 4), 'max': round(float(v.max()), 4)}


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
    print('deposit pairs at one differing position: %d' % len(one))
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

    urns = sorted({u for r, _ in one for u in (r['from_urn'], r['via_urn'])})
    scores = {}
    for i, u in enumerate(urns, 1):
        m, how = C.load_scores(u)
        if m is not None:
            scores[u] = m
        if i % 100 == 0 or i == len(urns):
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
        ys = [y for y in sorted(C.AA20) if y not in (X, Z) and (pos, y) in A and (pos, y) in B]
        if len(ys) < MIN_TARGETS:
            continue
        a = np.array([A[(pos, y)] for y in ys])
        b = np.array([B[(pos, y)] for y in ys])
        rt = spearman(a, b)
        if rt is None or a.std() <= 0 or b.std() <= 0:
            continue
        k = len(ys)

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
            rq = spearman(aq, bq)
            if rq is None or aq.std() <= 0 or bq.std() <= 0:
                continue
            pool.append((rq, float(aq.std()), float(bq.std())))
        if len(pool) < N_DRAW:
            continue

        rng.shuffle(pool)
        null_count = float(np.median([p[0] for p in pool[:N_DRAW]]))

        la, lb = math.log(a.std()), math.log(b.std())
        near = sorted(pool, key=lambda p: (math.log(p[1]) - la) ** 2
                      + (math.log(p[2]) - lb) ** 2)[:N_DRAW]
        null_spread = float(np.median([p[0] for p in near]))
        ratio = float(np.median([p[1] for p in near])) / a.std()

        rec.append({'from': r['from_urn'], 'via': r['via_urn'],
                    'target': r['from_name'].split('(')[0].strip(),
                    'pos': pos, 'targets': k, 'pool': len(pool),
                    'rho_test': round(rt, 5),
                    'null_count_matched': round(null_count, 5),
                    'null_spread_matched': round(null_spread, 5),
                    'shift_count_matched': round(rt - null_count, 5),
                    'shift_spread_matched': round(rt - null_spread, 5),
                    'sd_ratio_after_matching': round(ratio, 4)})

    if not rec:
        die('no pair yielded a null of either kind')

    so = np.array([x['shift_count_matched'] for x in rec])
    sn = np.array([x['shift_spread_matched'] for x in rec])
    out = {'pairs_at_one_differing_position': len(one),
           'directed_rows_before_collapsing': one_directed,
           'pairs_measured': len(rec),
           'distinct_targets': len(sorted({x['target'] for x in rec})),
           'targets': sorted({x['target'] for x in rec}),
           'min_targets': MIN_TARGETS, 'control_positions_per_null': N_DRAW,
           'pool_cap': MAX_POOL, 'seed': SEED,
           'identity_predicts_shift': 0.0,
           'rho_test': summ([x['rho_test'] for x in rec]),
           'null_count_matched': summ([x['null_count_matched'] for x in rec]),
           'null_spread_matched': summ([x['null_spread_matched'] for x in rec]),
           'shift_count_matched': summ(so),
           'shift_spread_matched': summ(sn),
           'median_abs_shift_count_matched': round(float(np.median(np.abs(so))), 4),
           'median_abs_shift_spread_matched': round(float(np.median(np.abs(sn))), 4),
           'sd_ratio_after_matching': summ([x['sd_ratio_after_matching'] for x in rec]),
           'pairs_shift_positive_count_matched': int((so > 0).sum()),
           'pairs_shift_positive_spread_matched': int((sn > 0).sum()),
           'per_pair': rec}
    with open(os.path.join(DATA, 'null_matching.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print()
    print('pairs measured: %d over %d targets' % (len(rec), out['distinct_targets']))
    print('the identity predicts a shift of exactly 0 at every one of them')
    print('  rank correlation at the test position: median %.4f'
          % out['rho_test']['median'])
    print('  null matched on count alone:  median %.4f, shift %+.4f, '
          'median |shift| %.4f'
          % (out['null_count_matched']['median'],
             out['shift_count_matched']['median'],
             out['median_abs_shift_count_matched']))
    print('  null matched on dispersion:   median %.4f, shift %+.4f, '
          'median |shift| %.4f'
          % (out['null_spread_matched']['median'],
             out['shift_spread_matched']['median'],
             out['median_abs_shift_spread_matched']))
    print('  dispersion ratio achieved by the matching: median %.3f'
          % out['sd_ratio_after_matching']['median'])
    print('wrote data/null_matching.json')


if __name__ == '__main__':
    main()
