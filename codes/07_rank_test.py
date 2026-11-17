#!/usr/bin/env python3
import collections
import csv
import itertools
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
DRAWS = 10000
SEED = 20261111
MIN_TARGETS = 6
AA20 = set('ACDEFGHIKLMNPQRSTVWY')
THREE = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
         'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
         'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
         'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
_B = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(itertools.product(_B, _B, _B))}
PRO = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')
LOG = []


def say(s):
    print(s)
    LOG.append(s)


def die(m):
    sys.exit('FATAL 07_rank_test.py: %s' % m)


def avg_rank(a):
    n = a.shape[-1]
    idx = np.arange(n, dtype=float)
    order = np.argsort(a, axis=-1, kind='stable')
    srt = np.take_along_axis(a, order, axis=-1)
    neq = srt[..., 1:] != srt[..., :-1]
    head = np.concatenate([np.ones(a.shape[:-1] + (1,), bool), neq], axis=-1)
    tail = np.concatenate([neq, np.ones(a.shape[:-1] + (1,), bool)], axis=-1)
    first = np.maximum.accumulate(np.where(head, idx, -1.0), axis=-1)
    last = np.minimum.accumulate(np.where(tail, idx, float(n))[..., ::-1],
                                 axis=-1)[..., ::-1]
    out = np.empty(a.shape, dtype=float)
    np.put_along_axis(out, order, (first + last) / 2.0, axis=-1)
    return out


def spearman(x, y):
    rx, ry = avg_rank(x), avg_rank(y)
    rx = rx - rx.mean(axis=-1, keepdims=True)
    ry = ry - ry.mean(axis=-1, keepdims=True)
    num = (rx * ry).sum(axis=-1)
    den = np.sqrt((rx * rx).sum(axis=-1) * (ry * ry).sum(axis=-1))
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0, num / den, np.nan)


def bh(ps):
    m = len(ps)
    o = sorted(range(m), key=lambda i: ps[i])
    out = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(o), 1):
        k = m - rank + 1
        prev = min(prev, ps[i] * m / k)
        out[i] = prev
    return out


def load():
    man = json.load(open(os.path.join(DATA, 'fetch_manifest.json')))
    prot = {}
    for bg, t in man['targets'].items():
        nt = open(os.path.join(ROOT, t['file'])).read().strip()
        prot[bg] = ''.join(CODON[nt[i:i + 3]] for i in range(0, len(nt), 3))
    S = {}
    for e in man['sets']:
        col = 'average' if 'average' in e['columns'] else 'score'
        d = {}
        for r in csv.DictReader(open(os.path.join(ROOT, e['file']))):
            g = PRO.match((r.get('hgvs_pro') or '').strip())
            if not g:
                die('%s carries an unparsed hgvs_pro' % e['urn'])
            alt = THREE.get(g.group(3))
            if alt not in AA20:
                continue
            v = (r.get(col) or '').strip()
            if v in ('', 'NA'):
                continue
            d[(int(g.group(2)), alt)] = float(v)
        S[(e['background'], e['phenotype'])] = d
        say('%-6s %-11s %4d substitutions from column %s'
            % (e['background'], e['phenotype'], len(d), col))
    return prot, S


def main():
    prot, S = load()
    say('')
    bp = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(DATA, 'background_positions.csv'))):
        bp[(r['bg_a'], r['bg_b'])].append(
            (int(r['pos_rbd']), int(r['pos_spike']), r['res_a'], r['res_b']))

    rng = np.random.default_rng(SEED)
    phen = sorted({p for _, p in S})
    out, calib = [], []
    for (a, b), items in sorted(bp.items()):
        for ph in phen:
            ka, kb = (a, ph), (b, ph)
            sa, sb = S[ka], S[kb]
            pool = []
            for p in range(1, len(prot[a]) + 1):
                x = prot[a][p - 1]
                if x != prot[b][p - 1]:
                    continue
                va, vb = [], []
                for y in sorted(AA20 - {x}):
                    u, v = sa.get((p, y)), sb.get((p, y))
                    if u is not None and v is not None:
                        va.append(u)
                        vb.append(v)
                if len(va) >= MIN_TARGETS:
                    pool.append((np.asarray(va), np.asarray(vb)))
            if len(pool) < 8:
                die('%s against %s on %s has only %d shared-parent positions'
                    % (a, b, ph, len(pool)))
            cells = []
            refs = {}
            for pos, sp, X, Z in sorted(items):
                va, vb, tgt = [], [], []
                for y in sorted(AA20 - {X, Z}):
                    one, second = sa.get((pos, y)), sb.get((pos, y))
                    if one is None or second is None:
                        continue
                    va.append(one)
                    vb.append(second)
                    tgt.append(y)
                if len(va) < MIN_TARGETS:
                    continue
                n = len(va)
                if n not in refs:
                    use = [(u, v) for u, v in pool if u.size >= n]
                    pick = rng.integers(0, len(use), DRAWS)
                    A = np.empty((DRAWS, n))
                    B = np.empty((DRAWS, n))
                    for i, j in enumerate(pick):
                        u, v = use[j]
                        sel = slice(None) if u.size == n else rng.choice(u.size, n, replace=False)
                        A[i] = u[sel]
                        B[i] = v[sel]
                    refs[n] = np.sort(spearman(A, B))
                ref = refs[n]
                rho = float(spearman(np.asarray(va)[None, :], np.asarray(vb)[None, :])[0])
                p = (float(np.searchsorted(ref, rho, side='right')) + 1.0) / (ref.size + 1.0)
                first = sa.get((pos, Z))
                cells.append({'pair': '%s|%s' % (a, b), 'phenotype': ph,
                              'pos_rbd': pos, 'pos_spike': sp, 'res_from': X,
                              'res_via': Z, 'n_targets': n, 'rho': round(rho, 5),
                              'null_median_rho': round(float(np.median(ref)), 5),
                              'null_p05_rho': round(float(np.percentile(ref, 5)), 5),
                              'first_leg': ('' if first is None else round(first, 5)),
                              'p_perm': round(p, 6)})
                q = float(ref[rng.integers(0, ref.size)])
                calib.append((float(np.searchsorted(ref, q, side='right')) + 1.0)
                             / (ref.size + 1.0))
            out.extend(cells)

    for ph in phen:
        sub = [c for c in out if c['phenotype'] == ph]
        for c, adj in zip(sub, bh([c['p_perm'] for c in sub])):
            c['p_bh'] = round(adj, 6)

    fields = ['pair', 'phenotype', 'pos_rbd', 'pos_spike', 'res_from', 'res_via',
              'n_targets', 'rho', 'null_median_rho', 'null_p05_rho', 'first_leg',
              'p_perm', 'p_bh']
    with open(os.path.join(DATA, 'rank_test.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(out, key=lambda c: c['p_bh']))

    rate = sum(1 for p in calib if p < 0.05) / float(len(calib))
    say('calibration: a value drawn from the null rejects at %.4f against a '
        'nominal 0.05, over %d draws' % (rate, len(calib)))
    say('')
    summ = {'draws': DRAWS, 'seed': SEED, 'min_targets': MIN_TARGETS,
            'calibration_rejection_at_nominal_05': round(rate, 4),
            'calibration_draws': len(calib), 'by_phenotype': []}
    for ph in phen:
        sub = [c for c in out if c['phenotype'] == ph]
        sig = [c for c in sub if c['p_bh'] < 0.05]
        mid = lambda v: sorted(v)[len(v) // 2]
        rec = {'phenotype': ph, 'cells': len(sub), 'departing': len(sig),
               'median_rho': round(mid([c['rho'] for c in sub]), 4),
               'median_null_rho': round(mid([c['null_median_rho'] for c in sub]), 4),
               'departing_spike_positions': sorted({c['pos_spike'] for c in sig})}
        summ['by_phenotype'].append(rec)
        say('%-11s cells %2d  departing %2d  median rho %7.4f  median null rho %7.4f'
            % (ph, rec['cells'], rec['departing'], rec['median_rho'],
               rec['median_null_rho']))
        say('%-11s departing spike positions: %s'
            % ('', ', '.join(str(x) for x in rec['departing_spike_positions']) or 'none'))
    allsig = sorted({c['pos_spike'] for c in out if c['p_bh'] < 0.05})
    summ['departing_spike_positions_either_phenotype'] = allsig
    summ['cells'] = len(out)
    summ['departing'] = sum(1 for c in out if c['p_bh'] < 0.05)
    say('')
    say('%d of %d cells depart; spike positions %s'
        % (summ['departing'], summ['cells'], ', '.join(str(x) for x in allsig)))
    with open(os.path.join(DATA, 'rank_test_summary.json'), 'w') as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'rank_test_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    print('\nwritten to data/rank_test.csv and data/rank_test_summary.json')


if __name__ == '__main__':
    main()
