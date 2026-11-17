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
QUANTILES = [75, 85, 90, 95]
HEADLINE_Q = 90
DRAWS_SWEEP = 2000
DRAWS_HEADLINE = 10000
SEED = 20261111
MIN_ROUTES = 6
MIN_CONTROL_POOL = 100
MIN_MAD = 0.01
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
    sys.exit('FATAL 11_census_test_scalefree.py: %s' % m)


def load_sequences():
    out = {}
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
        out[r['scoreSetUrn']] = p
    return out


class Scores(object):

    def __init__(self, limit=240):
        self.limit = limit
        self.d = collections.OrderedDict()

    def get(self, urn):
        if urn in self.d:
            self.d.move_to_end(urn)
            return self.d[urn]
        v = self._parse(urn)
        self.d[urn] = v
        if len(self.d) > self.limit:
            self.d.popitem(last=False)
        return v

    def _parse(self, urn):
        f = os.path.join(CACHE, urn.replace(':', '_') + '.csv')
        if not os.path.isfile(f):
            return None
        try:
            cols = open(f).readline().strip().split(',')
        except Exception:
            return None
        dc = 'average' if 'average' in cols else ('score' if 'score' in cols else None)
        if dc is None:
            return None
        eff, selfv = {}, collections.defaultdict(list)
        for r in csv.DictReader(open(f)):
            m = PRO.match((r.get('hgvs_pro') or '').strip())
            if not m:
                continue
            w, p, a = THREE.get(m.group(1)), int(m.group(2)), THREE.get(m.group(3))
            if a not in AA20 or w is None:
                continue
            v = (r.get(dc) or '').strip()
            if v in ('', 'NA'):
                continue
            try:
                x = float(v)
            except ValueError:
                continue
            if w == a:
                selfv[p].append(x)
            else:
                eff[(p, a)] = x
        if len(eff) < 20:
            return None
        vals = np.fromiter(eff.values(), dtype=float)
        med = float(np.median(vals))
        mad = float(np.median(np.abs(vals - med)))
        if mad < MIN_MAD:
            return None
        wt_dep = (float(np.median([x for v in selfv.values() for x in v]))
                  if selfv else None)
        wt_est = {q: float(np.percentile(vals, q)) for q in QUANTILES}
        return {'eff': eff, 'mad': mad, 'wt_dep': wt_dep, 'wt_est': wt_est,
                'n': len(eff), 'column': dc}


def centre(sc, q):
    return sc['wt_dep'] if sc['wt_dep'] is not None else sc['wt_est'][q]


def control_pool(sa, sb, pa, pb, ca, cb):
    out = []
    n = min(len(pa), len(pb))
    ea, eb, ma, mb = sa['eff'], sb['eff'], sa['mad'], sb['mad']
    for p in range(1, n + 1):
        x = pa[p - 1]
        if x != pb[p - 1]:
            continue
        for y in sorted(AA20):
            if y == x:
                continue
            u, v = ea.get((p, y)), eb.get((p, y))
            if u is not None and v is not None:
                out.append((v - cb) / mb - (u - ca) / ma)
    return np.asarray(out, dtype=float)


def residuals(sa, sb, pos, x, z, ca, cb):
    ea, eb, ma, mb = sa['eff'], sb['eff'], sa['mad'], sb['mad']
    az = ea.get((pos, z))
    if az is None:
        return None, None, None
    first = (az - ca) / ma
    rs, tgt = [], []
    for y in sorted(AA20 - {x, z}):
        one, second = ea.get((pos, y)), eb.get((pos, y))
        if one is None or second is None:
            continue
        rs.append(((one - ca) / ma - first) - ((second - cb) / mb))
        tgt.append(y)
    if not rs:
        return None, None, None
    return np.asarray(rs, dtype=float), first, tgt


class Reference(object):

    def __init__(self, seed):
        self.rng = np.random.default_rng(seed)
        self.cache = {}

    def get(self, key, pool, n, draws):
        k = (key, n, draws)
        if k not in self.cache:
            idx = self.rng.integers(0, pool.size, (draws, n))
            self.cache[k] = np.sort(np.median(np.abs(pool[idx]), axis=1))
        return self.cache[k]

    def clear(self):
        self.cache.clear()


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    out = [1.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, ps[i] * n / rank)
        out[i] = min(1.0, prev)
    return out


def run_quantile(q, draws, by_pair, seq, sc, ref, seq_id):
    ref.clear()
    out = []
    skipped = collections.Counter()
    for (ua, ub), items in sorted(by_pair.items()):
        pa, pb = seq.get(ua), seq.get(ub)
        if pa is None or pb is None or len(pa) != len(pb):
            skipped['no comparable target sequence'] += len(items)
            continue
        sa, sb = sc.get(ua), sc.get(ub)
        if sa is None or sb is None:
            skipped['score set unusable or degenerate spread'] += len(items)
            continue
        ca, cb = centre(sa, q), centre(sb, q)
        pool = control_pool(sa, sb, pa, pb, ca, cb)
        if pool.size < MIN_CONTROL_POOL:
            skipped['control pool below %d' % MIN_CONTROL_POOL] += len(items)
            continue
        cells = []
        for r in items:
            pos = int(r['pos'])
            rs, first, tgt = residuals(sa, sb, pos, r['res_from'], r['res_via'],
                                       ca, cb)
            if rs is None or rs.size < MIN_ROUTES:
                skipped['fewer than %d routes after centring' % MIN_ROUTES] += 1
                continue
            stat = float(np.median(np.abs(rs)))
            null = ref.get((ua, ub), pool, rs.size, draws)
            p = float((np.searchsorted(null, stat, side='left').item() * -1
                       + null.size) + 1) / (null.size + 1)
            cells.append({
                'from_urn': ua, 'via_urn': ub,
                'from_name': r['from_name'], 'via_name': r['via_name'],
                'seq_pair': '%d|%d' % (seq_id(pa), seq_id(pb)),
                'pos': pos, 'res_from': r['res_from'], 'res_via': r['res_via'],
                'via_is_alanine': int(r['res_via'] == 'A'),
                'n_routes': int(rs.size),
                'median_residual': round(float(np.median(rs)), 5),
                'median_abs_residual': round(stat, 5),
                'null_median_abs': round(float(np.median(null)), 5),
                'n_control_pool': int(pool.size),
                'intermediate_cost': round(float(first), 5),
                'p_perm': round(p, 6)})
        if cells:
            ps = bh([c['p_perm'] for c in cells])
            for c, v in zip(cells, ps):
                c['p_bh'] = round(v, 6)
            out.extend(cells)
    return out, skipped


def dedup(cells):
    best = {}
    for c in cells:
        k = (c['seq_pair'], c['pos'], c['res_from'], c['res_via'])
        if k not in best or c['n_routes'] > best[k]['n_routes']:
            best[k] = c
    return list(best.values())


def describe(cells, label):
    if not cells:
        say('%-28s no cells' % label)
        return {}
    sig = sum(1 for c in cells if c['p_bh'] < 0.05)
    cost = sorted(c['intermediate_cost'] for c in cells)
    nul = sorted(c['null_median_abs'] for c in cells)
    mar = sorted(c['median_abs_residual'] for c in cells)
    neg = sum(1 for c in cells if c['median_residual'] < 0)
    md = lambda v: v[len(v) // 2]
    d = {'cells': len(cells), 'significant_bh_0.05': sig,
         'significant_fraction': round(sig / float(len(cells)), 4),
         'median_intermediate_cost': round(md(cost), 4),
         'median_null_abs': round(md(nul), 4),
         'median_abs_residual': round(md(mar), 4),
         'negative_median_residual_cells': neg,
         'negative_fraction': round(neg / float(len(cells)), 4)}
    say('%-28s cells %5d  sig %5d (%5.1f%%)  med|res| %7.3f  med null %7.3f  '
        'med cost %7.3f  negative %5.1f%%'
        % (label, d['cells'], sig, 100 * d['significant_fraction'],
           d['median_abs_residual'], d['median_null_abs'],
           d['median_intermediate_cost'], 100 * d['negative_fraction']))
    return d


RBD = {'urn:mavedb:00000044-a-1', 'urn:mavedb:00000044-b-2',
       'urn:mavedb:00000668-a-1', 'urn:mavedb:00000669-a-1',
       'urn:mavedb:00000670-a-1', 'urn:mavedb:00000671-a-1'}


def main():
    seq = load_sequences()
    say('score sets with a translatable target sequence: %d' % len(seq))

    rows = [r for r in csv.DictReader(open(os.path.join(DATA, 'route_census_mavedb.csv')))
            if int(r['complete_routes']) >= MIN_ROUTES]
    by_pair = collections.defaultdict(list)
    for r in rows:
        by_pair[(r['from_urn'], r['via_urn'])].append(r)
    say('census rows with at least %d routes: %d, over %d score-set pairs'
        % (MIN_ROUTES, len(rows), len(by_pair)))

    ids = {}
    def seq_id(p):
        return ids.setdefault(p, len(ids))

    sc = Scores(limit=240)
    ref = Reference(SEED)

    say('')
    say('wild-type value used for the six receptor binding domain score sets, '
        'as a visible check on the estimator')
    for u in sorted(RBD):
        s = sc.get(u)
        if s is None:
            say('  %-24s unusable' % u)
            continue
        say('  %-24s deposited %-9s estimated at q%d %8.4f   spread %7.4f'
            % (u, ('%.4f' % s['wt_dep']) if s['wt_dep'] is not None else 'none',
               HEADLINE_Q, s['wt_est'][HEADLINE_Q], s['mad']))

    say('')
    say('sweep of the wild-type quantile, %d draws per cell' % DRAWS_SWEEP)
    sweep, cells_by_q = [], {}
    for q in QUANTILES:
        cells, skipped = run_quantile(q, DRAWS_SWEEP, by_pair, seq, sc, ref, seq_id)
        dd = dedup(cells)
        cells_by_q[q] = dd
        d = describe(dd, '  quantile %d' % q)
        d['quantile'] = q
        d['draws'] = DRAWS_SWEEP
        d['set_aside'] = dict(skipped)
        sweep.append(d)

    say('')
    say('headline pass at quantile %d, %d draws' % (HEADLINE_Q, DRAWS_HEADLINE))
    cells, skipped = run_quantile(HEADLINE_Q, DRAWS_HEADLINE, by_pair, seq, sc,
                                  ref, seq_id)
    head = dedup(cells)
    hd = describe(head, '  headline')

    say('')
    ala = [c for c in head if c['via_is_alanine']]
    describe(ala, '  alanine as intermediate')
    rbd = [c for c in head if c['from_urn'] in RBD and c['via_urn'] in RBD]
    describe(rbd, '  receptor binding domain')

    pairs = {c['seq_pair'] for c in head}
    sig_pairs = {c['seq_pair'] for c in head if c['p_bh'] < 0.05}
    say('')
    say('distinct sequence pairs: %d; with at least one significant cell: %d'
        % (len(pairs), len(sig_pairs)))

    fields = ['from_urn', 'via_urn', 'from_name', 'via_name', 'seq_pair', 'pos',
              'res_from', 'res_via', 'via_is_alanine', 'n_routes',
              'median_residual', 'median_abs_residual', 'null_median_abs',
              'n_control_pool', 'intermediate_cost', 'p_perm', 'p_bh']
    with open(os.path.join(DATA, 'census_test_scalefree.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(head, key=lambda c: (c['p_bh'], c['seq_pair'], c['pos'])))

    summary = {
        'headline_quantile': HEADLINE_Q,
        'draws_headline': DRAWS_HEADLINE,
        'draws_sweep': DRAWS_SWEEP,
        'seed': SEED,
        'min_routes': MIN_ROUTES,
        'min_control_pool': MIN_CONTROL_POOL,
        'min_mad': MIN_MAD,
        'headline': hd,
        'alanine_subset': {
            'cells': len(ala),
            'significant_bh_0.05': sum(1 for c in ala if c['p_bh'] < 0.05)},
        'receptor_binding_domain_subset': {
            'cells': len(rbd),
            'significant_bh_0.05': sum(1 for c in rbd if c['p_bh'] < 0.05)},
        'sequence_pairs': len(pairs),
        'sequence_pairs_with_significant_cell': len(sig_pairs),
        'quantile_sweep': sweep,
        'set_aside_headline': dict(skipped),
        'departures_from_the_plan': [
            'each scan is centred on its own wild-type value and divided by its '
            'own robust spread, so residuals are in units of each scan\'s own '
            'variation and not in the units of either deposit',
            'the wild-type value is deposited for a few score sets and estimated '
            'as the %dth percentile of the effect distribution for the rest; the '
            'estimator was validated against the deposits that carry the true '
            'value and the quantile is swept here' % HEADLINE_Q,
            'permutation draws are taken with replacement, not without',
            'multiplicity is controlled within a score-set pair, not within a '
            'phenotype, because the census spans assays with no common phenotype '
            'label',
            'the censoring margin is held at zero, because a published limit of '
            'detection is not available for most deposits',
        ],
    }
    with open(os.path.join(DATA, 'census_test_scalefree_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'census_test_scalefree_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    print('\nwritten to census_test_scalefree.csv and '
          'census_test_scalefree_summary.json')


if __name__ == '__main__':
    main()
