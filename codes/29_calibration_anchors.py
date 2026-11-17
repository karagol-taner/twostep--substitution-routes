#!/usr/bin/env python3
import csv
import importlib.util
import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_sp = importlib.util.spec_from_file_location(
    'rk', os.path.join(HERE, '12_census_test_rank.py'))
rk = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rk)

DATA = os.path.join(
    os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE, 'data')

AA20 = rk.AA20
MIN_FIT = 100

RBD = [
    ('urn:mavedb:00000044-a-1', 'Wuhan Hu-1', 'binding'),
    ('urn:mavedb:00000044-b-2', 'Wuhan Hu-1', 'expression'),
    ('urn:mavedb:00000668-a-1', 'BA.1', 'binding'),
    ('urn:mavedb:00000669-a-1', 'BA.1', 'expression'),
    ('urn:mavedb:00000670-a-1', 'BA.2', 'binding'),
    ('urn:mavedb:00000671-a-1', 'BA.2', 'expression'),
]


def die(m):
    sys.stderr.write('FATAL 29_calibration_anchors.py: %s\n' % m)
    raise SystemExit(1)


def score_path(urn):
    cached = os.path.join(HERE, 'scorecache', urn.replace(':', '_') + '.csv')
    if os.path.isfile(cached):
        return cached
    beside = os.path.join(DATA, 'scores_' + urn.replace(':', '_') + '.csv')
    if os.path.isfile(beside):
        return beside
    die('%s has no score file; run 01_fetch_scores.py to populate codes/scorecache' % urn)


def deposited_scores(urn):
    path = score_path(urn)
    m = {}
    with open(path) as fh:
        rd = csv.DictReader(fh)
        cols = rd.fieldnames or []
        col = 'average' if 'average' in cols else ('score' if 'score' in cols else None)
        if col is None:
            die('%s carries neither an average nor a score column' % path)
        for row in rd:
            g = rk.PRO.match((row.get('hgvs_pro') or '').strip())
            if not g:
                continue
            alt = rk.THREE.get(g.group(3))
            if alt not in AA20:
                continue
            v = (row.get(col) or '').strip()
            if v in ('', 'NA'):
                continue
            try:
                m[(int(g.group(2)), alt)] = float(v)
            except ValueError:
                continue
    if not m:
        die('%s yielded no usable scores' % path)
    return m


def fit_pearson(sa, sb, pa, pb):
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
        return None, len(xs)
    X, Y = np.asarray(xs), np.asarray(ys)
    if np.std(X) < 1e-9 or np.std(Y) < 1e-9:
        return None, len(xs)
    return float(np.corrcoef(X, Y)[0, 1]), len(xs)


def main():
    seq = rk.load_sequences()
    sco = {u: deposited_scores(u) for u, _, _ in RBD}
    out = []
    for (ua, ba, ka), (ub, bb, kb) in itertools.combinations(RBD, 2):
        pa, pb = seq.get(ua), seq.get(ub)
        if pa is None or pb is None or len(pa) != len(pb):
            die('%s and %s have no comparable target sequence' % (ua, ub))
        r, nfit = fit_pearson(sco[ua], sco[ub], pa, pb)
        if r is None:
            die('%s against %s gave no fit over %d matched substitutions' % (ua, ub, nfit))
        out.append({
            'from_urn': ua, 'via_urn': ub,
            'from_background': ba, 'via_background': bb,
            'from_phenotype': ka, 'via_phenotype': kb,
            'same_phenotype': ka == kb,
            'fit_pearson': round(r, 4), 'fit_n': nfit,
        })

    known = {}
    pth = os.path.join(DATA, 'calibrated_test.csv')
    if os.path.isfile(pth):
        for row in csv.DictReader(open(pth)):
            key = tuple(sorted((row['from_urn'], row['via_urn'])))
            known[key] = float(row['fit_pearson'])
    checked = 0
    for c in out:
        key = tuple(sorted((c['from_urn'], c['via_urn'])))
        if key in known:
            if abs(known[key] - c['fit_pearson']) > 5e-4:
                die('%s against %s: this script gives %.4f, calibrated_test.csv %.4f'
                    % (c['from_urn'], c['via_urn'], c['fit_pearson'], known[key]))
            checked += 1
    if not checked:
        die('no pair could be checked against calibrated_test.csv')

    same = sorted(c['fit_pearson'] for c in out if c['same_phenotype'])
    cross = sorted(c['fit_pearson'] for c in out if not c['same_phenotype'])
    rec = {
        'pairs': out,
        'same_phenotype_fits': same,
        'cross_phenotype_fits': cross,
        'n_same_phenotype': len(same),
        'n_cross_phenotype': len(cross),
        'same_phenotype_min': min(same),
        'same_phenotype_max': max(same),
        'cross_phenotype_min': min(cross),
        'cross_phenotype_max': max(cross),
        'min_fit_points': MIN_FIT,
        'agree_with_calibrated_test': checked,
    }
    p = os.path.join(DATA, 'calibration_anchors.json')
    with open(p, 'w') as fh:
        json.dump(rec, fh, indent=1, sort_keys=True)
    print('wrote %s' % p)
    for c in out:
        print('  %-14s %-14s %-11s %-11s %-6s r=%.4f  n=%d' % (
            c['from_urn'].split(':')[-1], c['via_urn'].split(':')[-1],
            c['from_phenotype'], c['via_phenotype'],
            'same' if c['same_phenotype'] else 'cross', c['fit_pearson'], c['fit_n']))
    print('\nsame-phenotype  (%d): %s' % (len(same), same))
    print('cross-phenotype (%d): %s' % (len(cross), cross))
    print('pairs agreeing with calibrated_test.csv: %d' % checked)


if __name__ == '__main__':
    main()
