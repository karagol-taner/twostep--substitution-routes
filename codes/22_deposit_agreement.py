#!/usr/bin/env python3
import collections
import csv
import importlib.util
import itertools
import json
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
MIN_TARGETS = 6
MAX_CONTROL = 40
SEED = 20261112


def die(m):
    sys.exit('FATAL 22_deposit_agreement.py: %s' % m)


def spearman(a, b):
    n = len(a)
    if n < MIN_TARGETS:
        return None
    ra, rb = _rank(a), _rank(b)
    ma, mb = ra.mean(), rb.mean()
    da, db = ra - ma, rb - mb
    v = (da * da).sum() * (db * db).sum()
    if v <= 0:
        return None
    return float((da * db).sum() / v ** 0.5)


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


def expset(u):
    return u.split(':')[2].split('-')[0]


def main():
    C = importlib.util.module_from_spec(
        importlib.util.spec_from_file_location(
            'census', os.path.join(HERE, '08_mavedb_route_census.py')))
    importlib.util.spec_from_file_location(
        'census', os.path.join(HERE, '08_mavedb_route_census.py')).loader.exec_module(C)

    seq = {}
    for r in json.load(open(os.path.join(DATA, 'target_genes.json'))):
        ts = r.get('targetSequence') or {}
        s = (ts.get('sequence') or '').upper()
        t = ts.get('sequenceType')
        p = C.translate(s) if (t == 'dna' and len(s) % 3 == 0) else (s if t == 'protein' else None)
        if p:
            seq[r['scoreSetUrn']] = p.rstrip('*')

    rows = [r for r in csv.DictReader(open(os.path.join(DATA, 'route_census_mavedb.csv')))
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
    print('deposit pairs differing at exactly one position: %d' % len(one))
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
    print('score sets to read: %d' % len(urns))
    os.makedirs(C.CACHE, exist_ok=True)
    scores, unreadable = {}, []
    for i, u in enumerate(urns, 1):
        m, how = C.load_scores(u)
        if m is None:
            unreadable.append(u)
        else:
            scores[u] = m
        if i % 50 == 0 or i == len(urns):
            print('  read %d of %d, %d unreadable' % (i, len(urns), len(unreadable)))

    rng = random.Random(SEED)
    rec = []
    for r, p0 in one:
        A, Bq = scores.get(r['from_urn']), scores.get(r['via_urn'])
        if A is None or Bq is None:
            continue
        pos = p0 + 1
        X, Z = seq[r['from_urn']][p0], seq[r['via_urn']][p0]
        ys = [y for y in sorted(C.AA20) if y not in (X, Z) and (pos, y) in A and (pos, y) in Bq]
        rt = spearman([A[(pos, y)] for y in ys], [Bq[(pos, y)] for y in ys])
        if rt is None:
            continue
        k = len(ys)

        shared = [q for q in range(1, len(seq[r['from_urn']]) + 1)
                  if q != pos and seq[r['from_urn']][q - 1] == seq[r['via_urn']][q - 1]]
        rng.shuffle(shared)
        ctrl = []
        for q in shared:
            yq = [y for y in sorted(C.AA20) if y != seq[r['from_urn']][q - 1]
                  and (q, y) in A and (q, y) in Bq]
            if len(yq) < k:
                continue
            pick = yq if len(yq) == k else rng.sample(yq, k)
            rc = spearman([A[(q, y)] for y in pick], [Bq[(q, y)] for y in pick])
            if rc is not None:
                ctrl.append(rc)
            if len(ctrl) >= MAX_CONTROL:
                break
        if not ctrl:
            continue

        leg = A.get((pos, Z))
        rec.append({'from': r['from_urn'], 'via': r['via_urn'],
                    'target': r['from_name'].split('(')[0].strip(),
                    'target_via': r['via_name'].split('(')[0].strip(),
                    'pos': pos, 'targets': k,
                    'same_experiment_set': expset(r['from_urn']) == expset(r['via_urn']),
                    'rho_test': round(rt, 5),
                    'rho_control': round(float(np.median(ctrl)), 5),
                    'control_positions': len(ctrl),
                    'first_leg': None if leg is None else round(float(leg), 5)})

    if not rec:
        die('no pair yielded both a test and a control measurement')

    def vectors(r):
        A, B = scores[r['from']], scores[r['via']]
        pos = r['pos']
        X = seq[r['from']][pos - 1]
        Z = seq[r['via']][pos - 1]
        ys = [y for y in sorted(C.AA20) if y not in (X, Z)
              and (pos, y) in A and (pos, y) in B]
        return {'from': r['from'], 'via': r['via'], 'target': r['target'],
                'pos': pos, 'parent_from': X, 'parent_via': Z,
                'rho': r['rho_test'], 'targets': ys,
                'one_step': [round(float(A[(pos, y)]), 5) for y in ys],
                'second_leg': [round(float(B[(pos, y)]), 5) for y in ys]}

    kmax = max(r['targets'] for r in rec)
    mid = float(np.median([r['rho_test'] for r in rec]))
    passing = sorted([r for r in rec if r['targets'] == kmax
                      and r['rho_test'] > 0.9],
                     key=lambda r: abs(r['rho_test'] - mid))
    examples = {}
    if passing:
        examples['passing'] = vectors(passing[0])
    worst = min(set(r['target'] for r in rec),
                key=lambda t: np.median([r['rho_test'] for r in rec
                                         if r['target'] == t]))
    fam = [r for r in rec if r['target'].endswith('MPL')] or \
          [r for r in rec if r['target'] == worst]
    fmid = float(np.median([r['rho_test'] for r in fam]))
    examples['failing'] = vectors(
        sorted([r for r in fam if r['targets'] == max(x['targets'] for x in fam)],
               key=lambda r: abs(r['rho_test'] - fmid))[0])

    rt = np.array([x['rho_test'] for x in rec])
    rc = np.array([x['rho_control'] for x in rec])
    same = np.array([x['same_experiment_set'] for x in rec])
    names = sorted({n for x in rec for n in (x['target'], x['target_via'])})

    def summ(v):
        return {'n': int(len(v)), 'median': round(float(np.median(v)), 4),
                'mean': round(float(v.mean()), 4),
                'q1': round(float(np.percentile(v, 25)), 4),
                'q3': round(float(np.percentile(v, 75)), 4),
                'min': round(float(v.min()), 4), 'max': round(float(v.max()), 4),
                'above_0.9': int((v > 0.9).sum())}

    out = {'pairs_at_one_differing_position': len(one),
           'directed_rows_before_collapsing': one_directed,
           'pairs_measured': len(rec),
           'score_sets_read': len(scores), 'score_sets_unreadable': unreadable,
           'distinct_targets': len(names), 'targets': names,
           'min_targets': MIN_TARGETS, 'control_positions_max': MAX_CONTROL,
           'seed': SEED,
           'identity_predicts': 1.0,
           'control': summ(rc), 'test': summ(rt),
           'control_same_experiment_set': summ(rc[same]) if same.any() else None,
           'control_across_experiment_sets': summ(rc[~same]) if (~same).any() else None,
           'test_same_experiment_set': summ(rt[same]) if same.any() else None,
           'test_across_experiment_sets': summ(rt[~same]) if (~same).any() else None,
           'paired_control_minus_test': summ(rc - rt),
           'pairs_where_test_below_control': int((rt < rc).sum()),
           'attenuation_ceiling_from_control_median':
               round(float(np.median(rc)) ** 0.5, 4),
           'examples': examples,
           'per_pair': rec}

    legs = np.array([abs(x['first_leg']) for x in rec if x['first_leg'] is not None])
    gaps = np.array([x['rho_control'] - x['rho_test'] for x in rec
                     if x['first_leg'] is not None])
    if len(legs) >= 20:
        g = gaps - gaps.mean(); l = legs - legs.mean()
        out['gap_vs_abs_first_leg'] = round(
            float((g * l).sum() / ((g * g).sum() * (l * l).sum()) ** 0.5), 4)
        out['gap_vs_abs_first_leg_n'] = int(len(legs))

    with open(os.path.join(DATA, 'deposit_agreement.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print()
    print('pairs measured: %d over %d distinct targets' % (len(rec), len(names)))
    print('control, same substitution measured twice: median %.4f, %d of %d above 0.9'
          % (out['control']['median'], out['control']['above_0.9'], len(rc)))
    print('test, parents differ at that one position: median %.4f, %d of %d above 0.9'
          % (out['test']['median'], out['test']['above_0.9'], len(rt)))
    print('test below control in %d of %d pairs' % (out['pairs_where_test_below_control'], len(rec)))
    if out['control_same_experiment_set']:
        print('control within one experiment set: median %.4f over %d pairs'
              % (out['control_same_experiment_set']['median'],
                 out['control_same_experiment_set']['n']))
    if out['control_across_experiment_sets']:
        print('control across experiment sets:    median %.4f over %d pairs'
              % (out['control_across_experiment_sets']['median'],
                 out['control_across_experiment_sets']['n']))
    print('wrote data/deposit_agreement.json')


if __name__ == '__main__':
    main()
