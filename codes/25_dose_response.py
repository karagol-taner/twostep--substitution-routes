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
MAX_POOL = 300
N_DRAW = 10
MAX_TEST = 12
SEED = 20261112


def die(m):
    sys.exit('FATAL 25_dose_response.py: %s' % m)


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


def pearson(u, v):
    u, v = np.asarray(u, float), np.asarray(v, float)
    du, dv = u - u.mean(), v - v.mean()
    d = (du * du).sum() * (dv * dv).sum()
    return None if d <= 0 else float((du * dv).sum() / d ** 0.5)


def partial(u, v, w):
    ruv, ruw, rvw = pearson(u, v), pearson(u, w), pearson(v, w)
    if None in (ruv, ruw, rvw):
        return None
    d = math.sqrt(max(0.0, (1 - ruw ** 2) * (1 - rvw ** 2)))
    return None if d <= 0 else round((ruv - ruw * rvw) / d, 4)


def summ(v):
    v = np.asarray(v, dtype=float)
    if not len(v):
        return None
    return {'n': int(len(v)), 'median': round(float(np.median(v)), 4),
            'mean': round(float(v.mean()), 4),
            'q1': round(float(np.percentile(v, 25)), 4),
            'q3': round(float(np.percentile(v, 75)), 4)}


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
    bypair = collections.defaultdict(list)
    for r in rows:
        bypair[(r['from_urn'], r['via_urn'])].append(r)
    pairs = []
    for (a, b), rs in bypair.items():
        if a not in seq or b not in seq or len(seq[a]) != len(seq[b]):
            continue
        nd = sum(1 for x, y in zip(seq[a], seq[b]) if x != y)
        if nd >= 1:
            pairs.append((a, b, nd, rs))
    if not pairs:
        die('no deposit pair carries a usable target sequence')
    print('deposit pairs with both sequences and a complete route: %d' % len(pairs))
    seen_un, dd = set(), []
    for a_, b_, nd_, rs_ in sorted(pairs, key=lambda t: (t[0], t[1])):
        k = (min(a_, b_), max(a_, b_))
        if k in seen_un:
            continue
        seen_un.add(k)
        dd.append((a_, b_, nd_, rs_))
    print('  of those, independent comparisons after collapsing direction: %d'
          % len(dd))
    pairs_directed = len(pairs)
    pairs = dd

    urns = sorted({u for a, b, _, _ in pairs for u in (a, b)})
    print('score sets to read: %d' % len(urns))
    scores = {}
    for i, u in enumerate(urns, 1):
        m, how = C.load_scores(u)
        if m is not None:
            scores[u] = m
        if i % 100 == 0 or i == len(urns):
            print('  read %d of %d, %d usable' % (i, len(urns), len(scores)))

    rng = random.Random(SEED)
    rec = []
    for n, (a, b, nd, rs) in enumerate(pairs, 1):
        if n % 1000 == 0:
            print('  measured %d of %d pairs' % (n, len(pairs)))
        A, B = scores.get(a), scores.get(b)
        if A is None or B is None:
            continue
        sa, sb = seq[a], seq[b]

        use = rs if len(rs) <= MAX_TEST else rng.sample(rs, MAX_TEST)
        tests = []
        for r in use:
            pos = int(r['pos'])
            if pos < 1 or pos > min(len(sa), len(sb)):
                continue
            X, Z = sa[pos - 1], sb[pos - 1]
            if X == Z:
                continue
            ys = [y for y in sorted(C.AA20) if y not in (X, Z)
                  and (pos, y) in A and (pos, y) in B]
            if len(ys) < MIN_TARGETS:
                continue
            va = np.array([A[(pos, y)] for y in ys])
            vb = np.array([B[(pos, y)] for y in ys])
            rt = spearman(va, vb)
            if rt is None or va.std() <= 0 or vb.std() <= 0:
                continue
            tests.append((rt, len(ys), float(va.std()), float(vb.std())))
        if not tests:
            continue

        shared = [q for q in range(1, min(len(sa), len(sb)) + 1)
                  if sa[q - 1] == sb[q - 1]]
        rng.shuffle(shared)
        pool = collections.defaultdict(list)
        for q in shared[:MAX_POOL]:
            yq = [y for y in sorted(C.AA20) if y != sa[q - 1]
                  and (q, y) in A and (q, y) in B]
            if len(yq) < MIN_TARGETS:
                continue
            pool['all'].append((q, yq))
        if len(pool['all']) < N_DRAW:
            continue

        t_rho, n_cnt, n_spr = [], [], []
        for rt, k, s1, s2 in tests:
            cand = []
            for q, yq in pool['all']:
                if len(yq) < k:
                    continue
                pick = yq if len(yq) == k else rng.sample(yq, k)
                aq = np.array([A[(q, y)] for y in pick])
                bq = np.array([B[(q, y)] for y in pick])
                rq = spearman(aq, bq)
                if rq is None or aq.std() <= 0 or bq.std() <= 0:
                    continue
                cand.append((rq, float(aq.std()), float(bq.std())))
            if len(cand) < N_DRAW:
                continue
            rng.shuffle(cand)
            nc = float(np.median([c[0] for c in cand[:N_DRAW]]))
            l1, l2 = math.log(s1), math.log(s2)
            near = sorted(cand, key=lambda c: (math.log(c[1]) - l1) ** 2
                          + (math.log(c[2]) - l2) ** 2)[:N_DRAW]
            ns = float(np.median([c[0] for c in near]))
            t_rho.append(rt)
            n_cnt.append(nc)
            n_spr.append(ns)
        if not t_rho:
            continue

        nm = rs[0]['from_name'].split('(')[0].strip()
        nmv = rs[0]['via_name'].split('(')[0].strip()
        rec.append({'from': a, 'via': b, 'target': nm, 'target_via': nmv,
                    'differing_positions': nd, 'test_positions': len(t_rho),
                    'rho_test': round(float(np.mean(t_rho)), 5),
                    'null_count_matched': round(float(np.mean(n_cnt)), 5),
                    'null_spread_matched': round(float(np.mean(n_spr)), 5),
                    'shift_count_matched': round(float(np.mean(t_rho) - np.mean(n_cnt)), 5),
                    'shift_spread_matched': round(float(np.mean(t_rho) - np.mean(n_spr)), 5)})

    if not rec:
        die('no deposit pair yielded both a test and a null')

    nd = np.array([x['differing_positions'] for x in rec], float)
    lg = np.log10(nd)
    npos = np.array([x['test_positions'] for x in rec], float)
    so = np.array([x['shift_count_matched'] for x in rec])
    sn = np.array([x['shift_spread_matched'] for x in rec])

    def bins(v):
        b = collections.OrderedDict()
        for lab, lo, hi in (('1', 1, 1), ('2', 2, 2), ('3', 3, 3),
                            ('4-6', 4, 6), ('7-20', 7, 20), ('>20', 21, 10 ** 9)):
            m = (nd >= lo) & (nd <= hi)
            if m.any():
                b[lab] = {'pairs': int(m.sum()),
                          'median_shift': round(float(np.median(v[m])), 4)}
        return b

    out = {'pairs_measured': len(rec),
           'directed_rows_before_collapsing': pairs_directed,
           'distinct_targets': len({n for x in rec
                                   for n in (x['target'], x['target_via'])}),
           'min_targets': MIN_TARGETS, 'control_positions_per_null': N_DRAW,
           'test_positions_cap': MAX_TEST, 'pool_cap': MAX_POOL, 'seed': SEED,
           'identity_predicts_shift_at_one_position': 0.0,
           'shift_count_matched': summ(so),
           'shift_spread_matched': summ(sn),
           'at_one_differing_position': {
               'pairs': int((nd == 1).sum()),
               'median_shift_count_matched': round(float(np.median(so[nd == 1])), 4),
               'median_shift_spread_matched': round(float(np.median(sn[nd == 1])), 4)},
           'by_differing_positions_count_matched': bins(so),
           'by_differing_positions_spread_matched': bins(sn),
           'trend_excluding_the_known_bin': {
               'pairs': int((nd >= 2).sum()),
               'count_matched': round(pearson(so[nd >= 2], lg[nd >= 2]), 4),
               'count_matched_partial':
                   partial(so[nd >= 2], lg[nd >= 2], npos[nd >= 2]),
               'spread_matched': round(pearson(sn[nd >= 2], lg[nd >= 2]), 4),
               'spread_matched_partial':
                   partial(sn[nd >= 2], lg[nd >= 2], npos[nd >= 2])},
           'trend_count_matched': {
               'shift_vs_log_positions': round(pearson(so, lg), 4),
               'partial_given_test_positions': partial(so, lg, npos)},
           'trend_spread_matched': {
               'shift_vs_log_positions': round(pearson(sn, lg), 4),
               'partial_given_test_positions': partial(sn, lg, npos)},
           'per_pair': rec}
    with open(os.path.join(DATA, 'dose_response.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print()
    print('pairs measured: %d over %d targets' % (len(rec), out['distinct_targets']))
    print('%-6s %6s %14s %14s' % ('diff', 'pairs', 'shift (count)', 'shift (spread)'))
    for k in out['by_differing_positions_count_matched']:
        print('%-6s %6d %14.4f %14.4f'
              % (k, out['by_differing_positions_count_matched'][k]['pairs'],
                 out['by_differing_positions_count_matched'][k]['median_shift'],
                 out['by_differing_positions_spread_matched'][k]['median_shift']))
    print()
    print('at one differing position the shift must be exactly 0:')
    print('  count-matched null  %+.4f' % out['at_one_differing_position']['median_shift_count_matched'])
    print('  spread-matched null %+.4f' % out['at_one_differing_position']['median_shift_spread_matched'])
    print('trend of shift against log differing positions:')
    print('  count-matched  %.4f (partial %s)'
          % (out['trend_count_matched']['shift_vs_log_positions'],
             out['trend_count_matched']['partial_given_test_positions']))
    print('  spread-matched %.4f (partial %s)'
          % (out['trend_spread_matched']['shift_vs_log_positions'],
             out['trend_spread_matched']['partial_given_test_positions']))
    print('with the bin whose value is known removed, over %d pairs:'
          % out['trend_excluding_the_known_bin']['pairs'])
    print('  count-matched  %.4f (partial %s)'
          % (out['trend_excluding_the_known_bin']['count_matched'],
             out['trend_excluding_the_known_bin']['count_matched_partial']))
    print('  spread-matched %.4f (partial %s)'
          % (out['trend_excluding_the_known_bin']['spread_matched'],
             out['trend_excluding_the_known_bin']['spread_matched_partial']))
    print('wrote data/dose_response.json')


if __name__ == '__main__':
    main()
