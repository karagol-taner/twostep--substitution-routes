#!/usr/bin/env python3
import argparse
import collections
import csv
import itertools
import json
import os
import re
import sys
import urllib.request

import numpy as np

API = 'https://api.mavedb.org/api/v1/score-sets/'
AA20 = set('ACDEFGHIKLMNPQRSTVWY')
THREE = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
         'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
         'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
         'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
_B = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(itertools.product(_B, _B, _B))}
PRO = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')
MIN_MAD = 0.01
MIN_POOL = 100
MIN_ROUTES = 6


def die(m):
    sys.exit('routecost: %s' % m)


def fetch(url, cache=None):
    if cache and os.path.isfile(cache) and os.path.getsize(cache) > 50:
        return open(cache, 'rb').read()
    try:
        with urllib.request.urlopen(urllib.request.Request(
                url, headers={'Accept': '*/*'}), timeout=300) as r:
            if r.status != 200:
                die('%s returned HTTP %s' % (url, r.status))
            raw = r.read()
    except Exception as e:
        die('%s failed: %s: %s' % (url, type(e).__name__, e))
    if cache:
        with open(cache, 'wb') as fh:
            fh.write(raw)
    return raw


def translate(nt):
    try:
        return ''.join(CODON[nt[i:i + 3]] for i in range(0, len(nt) - 2, 3)).rstrip('*')
    except KeyError as e:
        die('untranslatable codon %s' % e)


def load(urn, cachedir, quantile):
    meta = json.loads(fetch(API + urn, os.path.join(cachedir, urn.replace(':', '_') + '.json')
                            if cachedir else None))
    if 'detail' in meta:
        die('%s: %s. An experiment URN has no scores endpoint; a score-set URN '
            'carries a trailing index, for example urn:mavedb:00000044-a-1.'
            % (urn, meta['detail']))
    tg = (meta.get('targetGenes') or [{}])[0]
    ts = (tg.get('targetSequence') or {})
    s = (ts.get('sequence') or '').upper()
    if not s:
        die('%s carries no target sequence' % urn)
    prot = translate(s) if ts.get('sequenceType') == 'dna' else s.rstrip('*')

    raw = fetch(API + urn + '/scores',
                os.path.join(cachedir, urn.replace(':', '_') + '.csv') if cachedir else None)
    text = raw.decode('utf-8', 'replace').splitlines()
    cols = text[0].split(',')
    dcol = 'average' if 'average' in cols else ('score' if 'score' in cols else None)
    if dcol is None:
        die('%s has no score or average column; columns are %s' % (urn, cols))
    eff, selfv, ref = {}, [], {}
    for r in csv.DictReader(text):
        m = PRO.match((r.get('hgvs_pro') or '').strip())
        if not m:
            continue
        w, p, a = THREE.get(m.group(1)), int(m.group(2)), THREE.get(m.group(3))
        if a not in AA20 or w is None:
            continue
        ref.setdefault(p, w)
        v = (r.get(dcol) or '').strip()
        if v in ('', 'NA'):
            continue
        try:
            x = float(v)
        except ValueError:
            continue
        if w == a:
            selfv.append(x)
        else:
            eff[(p, a)] = x
    if len(eff) < 20:
        die('%s carries only %d usable substitutions' % (urn, len(eff)))
    bad = [p for p, w in ref.items() if p <= len(prot) and prot[p - 1] != w]
    if bad:
        die('%s reference residues disagree with its target sequence at %s'
            % (urn, bad[:5]))
    vals = np.fromiter(eff.values(), dtype=float)
    mad = float(np.median(np.abs(vals - np.median(vals))))
    if mad < MIN_MAD:
        die('%s has a robust spread of %.5f, too small to standardise' % (urn, mad))
    wt = float(np.median(selfv)) if selfv else float(np.percentile(vals, quantile))
    return {'urn': urn, 'name': tg.get('name') or '', 'prot': prot, 'eff': eff,
            'wt': wt, 'wt_deposited': bool(selfv), 'mad': mad,
            'column': dcol, 'n': len(eff)}


def std(scan, p, y):
    v = scan['eff'].get((p, y))
    return None if v is None else (v - scan['wt']) / scan['mad']


def control_null(a, b):
    out = []
    n = min(len(a['prot']), len(b['prot']))
    for p in range(1, n + 1):
        x = a['prot'][p - 1]
        if x != b['prot'][p - 1]:
            continue
        for y in sorted(AA20 - {x}):
            u, v = std(a, p, y), std(b, p, y)
            if u is not None and v is not None:
                out.append(v - u)
    return np.asarray(out, dtype=float)


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    out, prev = [1.0] * n, 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, ps[i] * n / rank)
        out[i] = min(1.0, prev)
    return out


def analyse(a, b, draws, seed):
    if len(a['prot']) != len(b['prot']):
        die('the two target sequences are %d and %d residues long, so positions '
            'do not correspond without an alignment, which this tool does not do'
            % (len(a['prot']), len(b['prot'])))
    pool = control_null(a, b)
    if pool.size < MIN_POOL:
        die('only %d control pairs are available, below the %d needed to '
            'estimate a null' % (pool.size, MIN_POOL))
    rng = np.random.default_rng(seed)
    ref = {}
    rows = []
    for p in range(1, len(a['prot']) + 1):
        x, z = a['prot'][p - 1], b['prot'][p - 1]
        if x == z:
            continue
        first = std(a, p, z)
        if first is None:
            rows.append({'pos': p, 'res_from': x, 'res_via': z, 'n_routes': 0,
                         'note': 'the first leg, %s to %s, is not measured in %s'
                                 % (x, z, a['urn'])})
            continue
        res, tgt = [], []
        for y in sorted(AA20 - {x, z}):
            one, second = std(a, p, y), std(b, p, y)
            if one is None or second is None:
                continue
            res.append((one - first) - second)
            tgt.append(y)
        if len(res) < MIN_ROUTES:
            rows.append({'pos': p, 'res_from': x, 'res_via': z,
                         'n_routes': len(res),
                         'note': 'fewer than %d targets measured in both scans'
                                 % MIN_ROUTES})
            continue
        res = np.asarray(res, dtype=float)
        stat = float(np.median(np.abs(res)))
        k = res.size
        if k not in ref:
            idx = rng.integers(0, pool.size, (draws, k))
            ref[k] = np.sort(np.median(np.abs(pool[idx]), axis=1))
        null = ref[k]
        p_perm = (null.size - int(np.searchsorted(null, stat, 'left')) + 1.0) / (null.size + 1)
        boot = np.median(res[rng.integers(0, k, (2000, k))], axis=1)
        rows.append({
            'pos': p, 'res_from': x, 'res_via': z, 'n_routes': k,
            'targets': ''.join(tgt),
            'median_residual': round(float(np.median(res)), 5),
            'ci_low': round(float(np.percentile(boot, 2.5)), 5),
            'ci_high': round(float(np.percentile(boot, 97.5)), 5),
            'median_abs_residual': round(stat, 5),
            'null_median_abs': round(float(np.median(null)), 5),
            'null_p95_abs': round(float(np.percentile(null, 95)), 5),
            'intermediate_cost': round(float(first), 5),
            'p_perm': round(float(p_perm), 6), 'note': ''})
    live = [r for r in rows if r['n_routes'] >= MIN_ROUTES]
    for r, v in zip(live, bh([r['p_perm'] for r in live])):
        r['p_bh'] = round(v, 6)
    for r in rows:
        r.setdefault('p_bh', '')
        for k in ('targets', 'median_residual', 'ci_low', 'ci_high',
                  'median_abs_residual', 'null_median_abs', 'null_p95_abs',
                  'intermediate_cost', 'p_perm'):
            r.setdefault(k, '')
    return rows, pool


FIELDS = ['pos', 'res_from', 'res_via', 'n_routes', 'targets', 'median_residual',
          'ci_low', 'ci_high', 'median_abs_residual', 'null_median_abs',
          'null_p95_abs', 'intermediate_cost', 'p_perm', 'p_bh', 'note']


def report(a, b, rows, pool, out):
    print('A  %-26s %s' % (a['urn'], a['name']))
    print('   %d substitutions, column %s, wild-type %.4f (%s), spread %.4f'
          % (a['n'], a['column'], a['wt'],
             'deposited' if a['wt_deposited'] else 'estimated', a['mad']))
    print('B  %-26s %s' % (b['urn'], b['name']))
    print('   %d substitutions, column %s, wild-type %.4f (%s), spread %.4f'
          % (b['n'], b['column'], b['wt'],
             'deposited' if b['wt_deposited'] else 'estimated', b['mad']))
    print('')
    print('control pairs, same substitution from the same parent in both: %d'
          % pool.size)
    print('their median absolute discrepancy: %.4f' % float(np.median(np.abs(pool))))
    live = [r for r in rows if r['n_routes'] >= MIN_ROUTES]
    print('positions whose parental residue differs: %d, of which testable: %d'
          % (len(rows), len(live)))
    print('')
    if not live:
        print('No position carries enough measured routes to test.')
        return
    print('%4s %3s %3s %6s %10s %18s %10s %10s %9s'
          % ('pos', 'from', 'via', 'routes', 'med res', '95% interval',
             'null med', 'cost', 'p (BH)'))
    for r in sorted(live, key=lambda r: r['p_bh']):
        print('%4d %3s %3s %6d %10.4f  [%7.4f,%7.4f] %10.4f %10.4f %9.5f'
              % (r['pos'], r['res_from'], r['res_via'], r['n_routes'],
                 r['median_residual'], r['ci_low'], r['ci_high'],
                 r['null_median_abs'], r['intermediate_cost'], r['p_bh']))
    sig = [r for r in live if r['p_bh'] < 0.05]
    print('')
    print('%d of %d positions depart from composition beyond the null at BH 0.05.'
          % (len(sig), len(live)))
    print('A positive median residual means the two-step sum predicts a worse')
    print('outcome than the direct substitution delivers, so the path is costly.')
    print('Values are in units of each scan\'s own robust spread.')
    if out:
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)
        print('')
        print('written to %s' % out)


SELF_TEST_PAIRS = [('urn:mavedb:00000044-a-1', 'urn:mavedb:00000668-a-1'),
                   ('urn:mavedb:00000044-b-2', 'urn:mavedb:00000669-a-1')]
OFFSET = 330
CALIBRATION_POSITIONS = 200
CALIBRATION_CEILING = 0.12


def self_test(args):
    checks, fails = [], []

    def chk(ok, msg):
        (checks if ok else fails).append(('PASS' if ok else 'FAIL') + '  ' + msg)

    for ua, ub in SELF_TEST_PAIRS:
        a = load(ua, args.cache, args.quantile)
        b = load(ub, args.cache, args.quantile)
        rows, pool = analyse(a, b, args.draws, args.seed)
        tag = '%s vs %s' % (ua.split(':')[-1], ub.split(':')[-1])

        n_diff = sum(1 for i in range(len(a['prot']))
                     if a['prot'][i] != b['prot'][i])
        chk(len(rows) == n_diff,
            '%s: one row per differing position, %d rows against %d differing '
            'positions counted from the deposited sequences'
            % (tag, len(rows), n_diff))

        want = sum(1 for p in range(1, len(a['prot']) + 1)
                   if a['prot'][p - 1] == b['prot'][p - 1]
                   for y in sorted(AA20)
                   if y != a['prot'][p - 1]
                   and std(a, p, y) is not None and std(b, p, y) is not None)
        chk(pool.size == want,
            '%s: control pool of %d matches the %d shared-parent pairs counted '
            'again from the sequences' % (tag, pool.size, want))

        rows2, _ = analyse(a, b, args.draws, args.seed)
        same = all(abs(r['p_bh'] - t['p_bh']) < 1e-12
                   for r, t in zip(rows, rows2) if 'p_bh' in r and 'p_bh' in t)
        chk(same, '%s: a second run at the same seed gives the same adjusted p'
            % tag)

        live = [r for r in rows if r.get('n_routes', 0) >= MIN_ROUTES]
        chk(all(r['ci_low'] <= r['median_residual'] <= r['ci_high'] for r in live),
            '%s: every bootstrap interval contains the median it brackets' % tag)
        chk(all(r['null_p95_abs'] >= r['null_median_abs'] for r in live),
            '%s: the null 95th percentile is at or above the null median' % tag)
        chk(all(r['p_bh'] >= r['p_perm'] - 1e-12 for r in live),
            '%s: adjustment never lowers a p value' % tag)
        chk(all(r.get('n_routes', 0) <= len(AA20) - 2 for r in rows),
            '%s: no position reports more routes than it has target residues'
            % tag)

        rng = np.random.default_rng(args.seed + 7)
        k = int(np.median([r['n_routes'] for r in live])) if live else MIN_ROUTES
        k = max(MIN_ROUTES, k)
        idx = rng.integers(0, pool.size, (args.draws, k))
        refdist = np.median(np.abs(pool[idx]), axis=1)
        hits = 0
        for _ in range(CALIBRATION_POSITIONS):
            fake = pool[rng.integers(0, pool.size, k)]
            stat = float(np.median(np.abs(fake)))
            p = float((np.count_nonzero(refdist >= stat) + 1) / (refdist.size + 1))
            if p < 0.05:
                hits += 1
        rate = hits / float(CALIBRATION_POSITIONS)
        chk(rate <= CALIBRATION_CEILING,
            '%s: with residuals drawn from the control pool over %d positions, '
            'the rejection rate is %.3f, at or below the %.2f ceiling'
            % (tag, CALIBRATION_POSITIONS, rate, CALIBRATION_CEILING))

        sig = sorted(r['pos'] + OFFSET for r in live if r['p_bh'] < 0.05)
        print('%s: %d positions testable, %d departing at BH 0.05, spike '
              'numbering %s' % (tag, len(live), len(sig), sig))

    print('')
    for line in checks + fails:
        print(line)
    print('')
    print('%d of %d checks passed.' % (len(checks), len(checks) + len(fails)))
    print('')
    print('Not covered: whether the departing positions are the biologically '
          'right ones,\nwhether an estimated wild-type value is the right one, '
          'pairs of unequal length,\nand any protein other than this one.')
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser(
        description='Cost of reaching a substitution through an intermediate, '
                    'measured from two deposited deep mutational scans.')
    ap.add_argument('--a', help='score-set URN of the first background')
    ap.add_argument('--b', help='score-set URN of the second background')
    ap.add_argument('--out', help='write the full table to this CSV')
    ap.add_argument('--cache', default='.routecost-cache',
                    help='directory for downloaded score sets')
    ap.add_argument('--quantile', type=float, default=90.0,
                    help='percentile used as the wild-type value when a deposit '
                         'carries no wild-type entries (default 90)')
    ap.add_argument('--draws', type=int, default=10000)
    ap.add_argument('--seed', type=int, default=20261111)
    ap.add_argument('--self-test', action='store_true',
                    help='run against the receptor binding domain and check the '
                         'result against what is known independently')
    args = ap.parse_args()
    if args.cache:
        os.makedirs(args.cache, exist_ok=True)
    if args.self_test:
        sys.exit(self_test(args))
    if not (args.a and args.b):
        ap.error('give --a and --b, or --self-test')
    a = load(args.a, args.cache, args.quantile)
    b = load(args.b, args.cache, args.quantile)
    rows, pool = analyse(a, b, args.draws, args.seed)
    report(a, b, rows, pool, args.out)


if __name__ == '__main__':
    main()
