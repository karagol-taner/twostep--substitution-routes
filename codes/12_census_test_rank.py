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
CACHE = os.path.join(HERE, 'scorecache')
DATA = os.path.join(
    os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE, 'data')
DRAWS = 10000
SEED = 20261111
MIN_TARGETS = 6
MIN_NULL_POSITIONS = 8
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
    print(s, flush=True)
    LOG.append(s)


def die(m):
    sys.exit('FATAL 12_census_test_rank.py: %s' % m)


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
    avg = (first + last) / 2.0
    out = np.empty(a.shape, dtype=float)
    np.put_along_axis(out, order, avg, axis=-1)
    return out


def spearman(x, y):
    rx, ry = avg_rank(x), avg_rank(y)
    rx = rx - rx.mean(axis=-1, keepdims=True)
    ry = ry - ry.mean(axis=-1, keepdims=True)
    num = (rx * ry).sum(axis=-1)
    den = np.sqrt((rx * rx).sum(axis=-1) * (ry * ry).sum(axis=-1))
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0, num / den, np.nan)


def load_sequences():
    seq = {}
    for r in json.load(open(os.path.join(DATA, 'target_genes.json'))):
        ts = r.get('targetSequence') or {}
        s = (ts.get('sequence') or '').upper()
        t = ts.get('sequenceType')
        if not s:
            continue
        if t == 'dna':
            if len(s) % 3:
                continue
            try:
                p = ''.join(CODON[s[i:i + 3]] for i in range(0, len(s) - 2, 3))
            except KeyError:
                continue
        elif t == 'protein':
            p = s
        else:
            continue
        p = p.rstrip('*')
        if '*' in p or len(p) < 20:
            continue
        seq[r['scoreSetUrn']] = p
    return seq


class Scores(object):
    def __init__(self, limit=240):
        self.limit = limit
        self.d = collections.OrderedDict()

    def get(self, urn):
        if urn in self.d:
            self.d.move_to_end(urn)
            return self.d[urn]
        path = os.path.join(CACHE, urn.replace(':', '_') + '.csv')
        m = {}
        if os.path.isfile(path):
            with open(path) as fh:
                rd = csv.DictReader(fh)
                cols = rd.fieldnames or []
                col = 'average' if 'average' in cols else ('score' if 'score' in cols else None)
                if col:
                    for row in rd:
                        g = PRO.match((row.get('hgvs_pro') or '').strip())
                        if not g:
                            continue
                        alt = THREE.get(g.group(3))
                        if alt not in AA20:
                            continue
                        v = (row.get(col) or '').strip()
                        if v in ('', 'NA'):
                            continue
                        try:
                            m[(int(g.group(2)), alt)] = float(v)
                        except ValueError:
                            continue
        self.d[urn] = m
        while len(self.d) > self.limit:
            self.d.popitem(last=False)
        return m


def null_vectors(sa, sb, pa, pb, n_min):
    out = []
    n = min(len(pa), len(pb))
    for p in range(1, n + 1):
        x = pa[p - 1]
        if x != pb[p - 1]:
            continue
        va, vb = [], []
        for y in sorted(AA20 - {x}):
            u, v = sa.get((p, y)), sb.get((p, y))
            if u is not None and v is not None:
                va.append(u)
                vb.append(v)
        if len(va) >= n_min:
            out.append((np.asarray(va), np.asarray(vb)))
    return out


def cell_vectors(sa, sb, pos, x, z):
    first = sa.get((pos, z))
    va, vb, tgt = [], [], []
    for y in sorted(AA20 - {x, z}):
        one, second = sa.get((pos, y)), sb.get((pos, y))
        if one is None or second is None:
            continue
        va.append(one)
        vb.append(second)
        tgt.append(y)
    return np.asarray(va), np.asarray(vb), tgt, first


class NullCache(object):

    def __init__(self, draws, rng):
        self.draws = draws
        self.rng = rng
        self.d = {}

    def get(self, key, pool, n):
        k = (key, n)
        if k in self.d:
            return self.d[k]
        usable = [(a, b) for a, b in pool if a.size >= n]
        if len(usable) < MIN_NULL_POSITIONS:
            self.d[k] = None
            return None
        pick = self.rng.integers(0, len(usable), self.draws)
        A = np.empty((self.draws, n))
        B = np.empty((self.draws, n))
        for i, j in enumerate(pick):
            a, b = usable[j]
            if a.size == n:
                sel = slice(None)
            else:
                sel = self.rng.choice(a.size, n, replace=False)
            A[i] = a[sel]
            B[i] = b[sel]
        ref = np.sort(spearman(A, B))
        self.d[k] = ref
        if len(self.d) > 4000:
            self.d.clear()
            self.d[k] = ref
        return ref


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


FIELDS = ['from_urn', 'via_urn', 'from_name', 'via_name', 'pos', 'res_from',
          'res_via', 'via_is_alanine', 'n_targets', 'rho', 'null_median_rho',
          'null_p05_rho', 'first_leg', 'p_perm', 'p_bh']
RBD = {'urn:mavedb:00000044-a-1', 'urn:mavedb:00000044-b-2',
       'urn:mavedb:00000668-a-1', 'urn:mavedb:00000669-a-1',
       'urn:mavedb:00000670-a-1', 'urn:mavedb:00000671-a-1'}


def summarise(label, rows):
    if not rows:
        return {'label': label, 'cells': 0}
    sig = [r for r in rows if r['p_bh'] < 0.05]
    rho = sorted(r['rho'] for r in rows)
    nul = sorted(r['null_median_rho'] for r in rows)
    mid = lambda v: v[len(v) // 2]
    d = {'label': label, 'cells': len(rows), 'significant': len(sig),
         'significant_fraction': round(len(sig) / float(len(rows)), 4),
         'median_rho': round(mid(rho), 4),
         'median_null_rho': round(mid(nul), 4),
         'cells_with_rho_below_null_median':
             sum(1 for r in rows if r['rho'] < r['null_median_rho'])}
    say('  %-30s cells %5d  departing %5d (%5.1f%%)  median rho %7.4f  '
        'median null rho %7.4f'
        % (label, d['cells'], d['significant'], 100 * d['significant_fraction'],
           d['median_rho'], d['median_null_rho']))
    return d


def main():
    seq = load_sequences()
    say('score sets with a translatable target sequence: %d' % len(seq))
    rows = list(csv.DictReader(open(os.path.join(DATA, 'route_census_mavedb.csv'))))
    live = [r for r in rows if int(r['complete_routes']) >= MIN_TARGETS]
    by_pair = collections.defaultdict(list)
    for r in live:
        by_pair[(r['from_urn'], r['via_urn'])].append(r)
    say('census rows with at least %d routes: %d, over %d score-set pairs'
        % (MIN_TARGETS, len(live), len(by_pair)))

    sc = Scores(limit=240)
    rng = np.random.default_rng(SEED)
    nulls = NullCache(DRAWS, rng)
    out = []
    aside = collections.Counter()
    calib = []
    done = 0
    for (ua, ub), items in sorted(by_pair.items()):
        done += 1
        if done % 500 == 0:
            say('  %d of %d score-set pairs, %d cells so far' % (done, len(by_pair), len(out)))
        pa, pb = seq.get(ua), seq.get(ub)
        if pa is None or pb is None or len(pa) != len(pb):
            aside['no comparable target sequence'] += len(items)
            continue
        sa, sb = sc.get(ua), sc.get(ub)
        if not sa or not sb:
            aside['deposit unreadable'] += len(items)
            continue
        pool = null_vectors(sa, sb, pa, pb, MIN_TARGETS)
        if len(pool) < MIN_NULL_POSITIONS:
            aside['too few shared-parent positions for a null'] += len(items)
            continue
        cells = []
        for r in items:
            pos, x, z = int(r['pos']), r['res_from'], r['res_via']
            if pos > len(pa) or pa[pos - 1] != x or pb[pos - 1] != z:
                aside['census row disagrees with the target sequence'] += 1
                continue
            va, vb, tgt, first = cell_vectors(sa, sb, pos, x, z)
            if va.size < MIN_TARGETS:
                aside['too few shared targets'] += 1
                continue
            ref = nulls.get((ua, ub), pool, va.size)
            if ref is None:
                aside['no null at this number of targets'] += 1
                continue
            rho = float(spearman(va[None, :], vb[None, :])[0])
            if not np.isfinite(rho):
                aside['rank correlation undefined, a vector is constant'] += 1
                continue
            p = (float(np.searchsorted(ref, rho, side='right')) + 1.0) / (ref.size + 1.0)
            cells.append({'from_urn': ua, 'via_urn': ub,
                          'from_name': r['from_name'], 'via_name': r['via_name'],
                          'pos': pos, 'res_from': x, 'res_via': z,
                          'via_is_alanine': int(z == 'A'),
                          'n_targets': int(va.size), 'rho': round(rho, 5),
                          'null_median_rho': round(float(np.median(ref)), 5),
                          'null_p05_rho': round(float(np.percentile(ref, 5)), 5),
                          'first_leg': ('' if first is None else round(first, 5)),
                          'p_perm': round(p, 6)})
            if len(calib) < 4000:
                q = float(ref[rng.integers(0, ref.size)])
                calib.append((float(np.searchsorted(ref, q, side='right')) + 1.0)
                             / (ref.size + 1.0))
        if not cells:
            continue
        adj = bh([c['p_perm'] for c in cells])
        for c, a in zip(cells, adj):
            c['p_bh'] = round(a, 6)
        out.extend(cells)

    say('cells tested: %d' % len(out))
    for k, v in sorted(aside.items()):
        say('  set aside, %-52s %d' % (k, v))

    seqid, best = {}, {}
    for c in out:
        a, b = seq.get(c['from_urn']), seq.get(c['via_urn'])
        ka = seqid.setdefault(a, len(seqid))
        kb = seqid.setdefault(b, len(seqid))
        k = (ka, kb, c['pos'], c['res_from'], c['res_via'])
        if k not in best or c['rho'] > best[k]['rho']:
            best[k] = c
    ded = list(best.values())
    say('deduplicated cells, one per sequence pair and position: %d' % len(ded))
    say('')

    rate = sum(1 for p in calib if p < 0.05) / float(len(calib)) if calib else float('nan')
    say('calibration: a value drawn from the null rejects at %.4f against a '
        'nominal 0.05, over %d draws' % (rate, len(calib)))
    say('')
    summ = {'draws': DRAWS, 'seed': SEED, 'min_targets': MIN_TARGETS,
            'min_null_positions': MIN_NULL_POSITIONS,
            'statistic': 'Spearman rank correlation between the one-step vector '
                         'and the second-leg vector over the target residues at '
                         'one position; additivity predicts 1',
            'null': 'the same statistic at shared-parent positions of the same '
                    'two deposits, matched on the number of targets',
            'invariance': 'unchanged by any increasing affine rescaling of '
                          'either deposit, so no offset or scale is estimated',
            'calibration_rejection_at_nominal_05': round(rate, 4),
            'calibration_draws': len(calib),
            'cells_before_dedup': len(out), 'set_aside': dict(aside)}
    summ['all'] = summarise('all cells', ded)
    summ['alanine'] = summarise('alanine as intermediate',
                                [c for c in ded if c['via_is_alanine']])
    summ['rbd'] = summarise('receptor binding domain',
                            [c for c in out if c['from_urn'] in RBD and c['via_urn'] in RBD])
    pairs = {(c['from_urn'], c['via_urn']) for c in ded}
    sigp = {(c['from_urn'], c['via_urn']) for c in ded if c['p_bh'] < 0.05}
    summ['sequence_pairs'] = len(pairs)
    summ['pairs_with_a_departing_cell'] = len(sigp)
    say('')
    say('deposit pairs: %d; with at least one departing position: %d' % (len(pairs), len(sigp)))

    with open(os.path.join(DATA, 'census_test_rank.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(sorted(ded, key=lambda c: c['p_bh']))
    with open(os.path.join(DATA, 'census_test_rank_summary.json'), 'w') as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'census_test_rank_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    say('')
    say('written to census_test_rank.csv and census_test_rank_summary.json')


if __name__ == '__main__':
    main()
