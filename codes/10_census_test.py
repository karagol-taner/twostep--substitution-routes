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
MIN_CONTROL_POOL = 100
MIN_ROUTES = 6
SEED = 20261111
AA20 = set('ACDEFGHIKLMNPQRSTVWY')
THREE = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
         'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
         'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
         'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
PRO = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')
_B = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(itertools.product(_B, _B, _B))}

LOG = []


def say(s):
    print(s, flush=True)
    LOG.append(s)


def die(m):
    sys.exit('FATAL 10_census_test.py: %s' % m)


def translate(nt):
    out = []
    for i in range(0, len(nt) - 2, 3):
        c = CODON.get(nt[i:i + 3])
        if c is None:
            return None
        out.append(c)
    return ''.join(out)


def load_sequences():
    path = os.path.join(DATA, 'target_genes.json')
    if not os.path.isfile(path):
        die('target_genes.json is absent')
    seq = {}
    for r in json.load(open(path)):
        ts = r.get('targetSequence') or {}
        s = (ts.get('sequence') or '').upper()
        t = ts.get('sequenceType')
        if not s:
            continue
        if t == 'dna':
            if len(s) % 3:
                continue
            p = translate(s)
        elif t == 'protein':
            p = s
        else:
            continue
        if not p:
            continue
        p = p.rstrip('*')
        if '*' in p or len(p) < 20:
            continue
        seq[r['scoreSetUrn']] = p
    return seq


class Scores(object):

    def __init__(self, limit=60):
        self.limit = limit
        self.store = collections.OrderedDict()

    def get(self, urn):
        if urn in self.store:
            self.store.move_to_end(urn)
            return self.store[urn]
        path = os.path.join(CACHE, urn.replace(':', '_') + '.csv')
        if not os.path.isfile(path):
            self.store[urn] = None
            return None
        with open(path) as fh:
            head = fh.readline().strip().split(',')
            dcol = 'average' if 'average' in head else ('score' if 'score' in head else None)
            if dcol is None:
                self.store[urn] = None
                return None
            fh.seek(0)
            m = {}
            for r in csv.DictReader(fh):
                g = PRO.match((r.get('hgvs_pro') or '').strip())
                if not g:
                    continue
                w, pos, alt = THREE.get(g.group(1)), int(g.group(2)), THREE.get(g.group(3))
                if alt not in AA20 or w not in AA20 or alt == w:
                    continue
                v = (r.get(dcol) or '').strip()
                if v in ('', 'NA', 'nan'):
                    continue
                try:
                    m[(pos, alt)] = float(v)
                except ValueError:
                    continue
        self.store[urn] = m
        while len(self.store) > self.limit:
            self.store.popitem(last=False)
        return m


def control_pool(sa, sb, pa, pb):
    out = []
    n = min(len(pa), len(pb))
    for p in range(1, n + 1):
        x = pa[p - 1]
        if x != pb[p - 1]:
            continue
        for y in sorted(AA20):
            if y == x:
                continue
            u, v = sa.get((p, y)), sb.get((p, y))
            if u is not None and v is not None:
                out.append(v - u)
    return np.asarray(out, dtype=float)


def residuals(sa, sb, pos, x, z):
    first = sa.get((pos, z))
    if first is None:
        return None, None
    rs, tgt = [], []
    for y in sorted(AA20 - {x, z}):
        one, second = sa.get((pos, y)), sb.get((pos, y))
        if one is None or second is None:
            continue
        rs.append(one - (first + second))
        tgt.append(y)
    if not rs:
        return None, None
    return np.asarray(rs, dtype=float), (first, tgt)


class Reference(object):

    def __init__(self, draws, seed):
        self.draws = draws
        self.rng = np.random.default_rng(seed)
        self.cache = {}

    def get(self, pool_id, pool, n):
        key = (pool_id, n)
        hit = self.cache.get(key)
        if hit is not None:
            return hit
        if pool.size < n:
            return None
        idx = self.rng.integers(0, pool.size, size=(self.draws, n))
        ref = np.median(np.abs(pool[idx]), axis=1)
        ref.sort()
        self.cache[key] = ref
        if len(self.cache) > 4000:
            self.cache.clear()
            self.cache[key] = ref
        return ref


def bh(pvals):
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        val = min(prev, pvals[i] * n / float(rank + 1))
        adj[i] = val
        prev = val
    return adj


def boot_ci(x, rng, draws=2000):
    if x.size < 2:
        return float(x[0]), float(x[0])
    idx = rng.integers(0, x.size, size=(draws, x.size))
    med = np.median(x[idx], axis=1)
    return float(np.percentile(med, 2.5)), float(np.percentile(med, 97.5))


def main():
    seq = load_sequences()
    say('score sets with a translatable target sequence: %d' % len(seq))

    rows = list(csv.DictReader(open(os.path.join(DATA, 'route_census_mavedb.csv'))))
    live = [r for r in rows if int(r['complete_routes']) >= MIN_ROUTES]
    say('census rows: %d; with at least %d routes: %d' % (len(rows), MIN_ROUTES, len(live)))

    by_pair = collections.defaultdict(list)
    for r in live:
        by_pair[(r['from_urn'], r['via_urn'])].append(r)
    say('score-set pairs to test: %d' % len(by_pair))

    sc = Scores(limit=150)
    ref = Reference(DRAWS, SEED)
    rng = np.random.default_rng(SEED + 1)
    sid, out = {}, []
    skipped = collections.Counter()
    done = 0

    for (ua, ub), items in sorted(by_pair.items()):
        done += 1
        if done % 200 == 0:
            say('  %d of %d score-set pairs, %d cells so far' % (done, len(by_pair), len(out)))
        if done % 1000 == 0 and out:
            with open(os.path.join(DATA, 'census_test_partial.csv'), 'w', newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=sorted(out[0].keys()))
                w.writeheader()
                w.writerows(out)
        pa, pb = seq.get(ua), seq.get(ub)
        if pa is None or pb is None or len(pa) != len(pb):
            skipped['no comparable target sequence'] += len(items)
            continue
        sa, sb = sc.get(ua), sc.get(ub)
        if not sa or not sb:
            skipped['score set unreadable'] += len(items)
            continue
        pool = control_pool(sa, sb, pa, pb)
        if pool.size < MIN_CONTROL_POOL:
            skipped['control pool under %d' % MIN_CONTROL_POOL] += len(items)
            continue
        pool_id = (ua, ub)
        seq_pair = tuple(sorted((sid.setdefault(pa, len(sid)), sid.setdefault(pb, len(sid)))))
        for r in items:
            pos, x, z = int(r['pos']), r['res_from'], r['res_via']
            if pa[pos - 1] != x or pb[pos - 1] != z:
                skipped['position disagrees with the target sequence'] += 1
                continue
            rs, extra = residuals(sa, sb, pos, x, z)
            if rs is None or rs.size < MIN_ROUTES:
                skipped['too few complete routes'] += 1
                continue
            reference = ref.get(pool_id, pool, rs.size)
            if reference is None:
                skipped['control pool smaller than the route count'] += 1
                continue
            stat = float(np.median(np.abs(rs)))
            p = float((np.searchsorted(reference, stat, side='left') * -1 + reference.size)
                      + 1) / (reference.size + 1)
            lo, hi = boot_ci(rs, rng)
            first, _ = extra
            out.append({
                'from_urn': ua, 'via_urn': ub, 'seq_pair': '%d|%d' % seq_pair,
                'from_name': r['from_name'][:70], 'via_name': r['via_name'][:70],
                'pos': pos, 'res_from': x, 'res_via': z,
                'via_is_alanine': int(z == 'A'),
                'n_routes': int(rs.size),
                'median_residual': round(float(np.median(rs)), 4),
                'ci_low': round(lo, 4), 'ci_high': round(hi, 4),
                'median_abs_residual': round(stat, 4),
                'null_median_abs': round(float(np.median(reference)), 4),
                'null_p95_abs': round(float(np.percentile(reference, 95)), 4),
                'n_control_pool': int(pool.size),
                'intermediate_cost': round(float(first), 4),
                'p_perm': round(p, 6),
            })

    say('cells tested: %d' % len(out))
    for k, v in sorted(skipped.items()):
        say('  set aside, %-46s %d' % (k, v))
    if not out:
        die('no cell survived the filters')

    groups = collections.defaultdict(list)
    for i, c in enumerate(out):
        groups[(c['from_urn'], c['via_urn'])].append(i)
    for idxs in groups.values():
        adj = bh([out[i]['p_perm'] for i in idxs])
        for i, a in zip(idxs, adj):
            out[i]['p_bh'] = round(a, 6)

    best = {}
    for c in out:
        k = (c['seq_pair'], c['pos'], c['res_from'], c['res_via'])
        if k not in best or c['n_control_pool'] > best[k]['n_control_pool']:
            best[k] = c
    dedup = sorted(best.values(), key=lambda c: c['p_bh'])
    say('deduplicated cells, one per sequence pair and position: %d' % len(dedup))

    sig = [c for c in dedup if c['p_bh'] < 0.05]
    ala = [c for c in dedup if c['via_is_alanine']]
    ala_sig = [c for c in ala if c['p_bh'] < 0.05]
    say('')
    say('significant at BH 0.05: %d of %d cells (%.1f%%)'
        % (len(sig), len(dedup), 100.0 * len(sig) / len(dedup)))
    say('alanine as the intermediate: %d cells, %d significant (%.1f%%)'
        % (len(ala), len(ala_sig), 100.0 * len(ala_sig) / len(ala) if ala else 0.0))
    say('distinct sequence pairs tested: %d' % len({c['seq_pair'] for c in dedup}))
    say('distinct sequence pairs with a significant cell: %d'
        % len({c['seq_pair'] for c in sig}))
    med_cost = float(np.median([c['intermediate_cost'] for c in dedup]))
    say('median intermediate cost over cells: %.4f' % med_cost)
    neg = sum(1 for c in dedup if c['median_residual'] < 0)
    say('cells whose median residual is negative: %d of %d' % (neg, len(dedup)))

    fields = ['from_urn', 'via_urn', 'seq_pair', 'from_name', 'via_name', 'pos',
              'res_from', 'res_via', 'via_is_alanine', 'n_routes',
              'median_residual', 'ci_low', 'ci_high', 'median_abs_residual',
              'null_median_abs', 'null_p95_abs', 'n_control_pool',
              'intermediate_cost', 'p_perm', 'p_bh']
    with open(os.path.join(DATA, 'census_test.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(dedup)
    summary = {
        'draws': DRAWS, 'seed': SEED, 'min_control_pool': MIN_CONTROL_POOL,
        'min_routes': MIN_ROUTES, 'censoring_margin': 0.0,
        'cells_tested_before_dedup': len(out),
        'cells': len(dedup), 'significant_bh_0.05': len(sig),
        'alanine_cells': len(ala), 'alanine_significant': len(ala_sig),
        'sequence_pairs': len({c['seq_pair'] for c in dedup}),
        'sequence_pairs_with_significant_cell': len({c['seq_pair'] for c in sig}),
        'median_intermediate_cost': round(med_cost, 4),
        'negative_median_residual_cells': neg,
        'set_aside': dict(skipped),
        'departures_from_the_plan': [
            'multiplicity controlled within a score-set pair, not within a '
            'phenotype, because the census spans assays with no common '
            'phenotype label',
            'censoring margin held at zero, because a published limit of '
            'detection is not available for most deposits; the receptor '
            'binding domain sweep was flat across the whole margin range',
        ],
    }
    with open(os.path.join(DATA, 'census_test_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'census_test_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    print('\nwritten to census_test.csv, census_test_summary.json')


if __name__ == '__main__':
    main()
