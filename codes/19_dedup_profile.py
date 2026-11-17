#!/usr/bin/env python3
import collections
import csv
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')

B = 'TCAG'
AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: AAS[i] for i, (a, b, c) in enumerate(itertools.product(B, B, B))}
AA20 = sorted(set(AAS) - {'*'})
COD = collections.defaultdict(list)
for c, a in CODON.items():
    if a != '*':
        COD[a].append(c)
DIST = {(x, y): min(sum(1 for i in range(3) if cx[i] != cy[i])
                    for cx in COD[x] for cy in COD[y])
        for x in sorted(AA20) for y in sorted(AA20)}


def die(m):
    sys.exit('FATAL 19_dedup_profile.py: %s' % m)


def translate(nt):
    out = []
    for i in range(0, len(nt) - 2, 3):
        c = CODON.get(nt[i:i + 3])
        if c is None:
            return None
        out.append(c)
    return ''.join(out)


def main():
    tgp = os.path.join(DATA, 'target_genes.json')
    cen = os.path.join(DATA, 'route_census_mavedb.csv')
    pub = os.path.join(DATA, 'mavedb_census_dedup.json')
    for p in (tgp, cen, pub):
        if not os.path.isfile(p):
            die('%s is absent' % os.path.relpath(p, ROOT))

    seq_of = {}
    for r in json.load(open(tgp)):
        ts = r.get('targetSequence') or {}
        s = (ts.get('sequence') or '').upper()
        t = ts.get('sequenceType')
        if not s:
            continue
        if t == 'dna':
            p = translate(s) if len(s) % 3 == 0 else None
        elif t == 'protein':
            p = s
        else:
            p = None
        if p:
            seq_of[r['scoreSetUrn']] = p.rstrip('*')

    sid = {}

    def sequence_id(u):
        s = seq_of.get(u)
        if s is None:
            return None
        if s not in sid:
            sid[s] = len(sid)
        return sid[s]

    best = {}
    via_ala = set()
    seq_pairs = set()
    fan_res = {}
    raw_fans = set()
    with open(cen) as fh:
        for r in csv.DictReader(fh):
            n = int(r['complete_routes'])
            if n < 1:
                continue
            raw_fans.add((r['from_urn'], r['via_urn'], r['pos']))
            a, b = sequence_id(r['from_urn']), sequence_id(r['via_urn'])
            if a is None or b is None:
                continue
            key = (a, b, r['pos'], r['res_from'], r['res_via'])
            if n > best.get(key, -1):
                best[key] = n
            if r['res_via'] == 'A':
                via_ala.add(key)
            fan_res[key] = (r['res_from'], r['res_via'])
            seq_pairs.add((min(a, b), max(a, b)))

    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in seq_pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    clusters = {find(a) for a, _ in seq_pairs}

    total = sum(best.values())
    ala_total = sum(v for k, v in best.items() if k in via_ala)

    pubd = json.load(open(pub))
    for name, got, want in (('route fans', len(best), pubd['route_fans']),
                            ('routes', total, pubd['deduplicated_routes']),
                            ('routes via alanine', ala_total,
                             pubd['deduplicated_routes_via_alanine']),
                            ('background pairs', len(seq_pairs),
                             pubd['background_pairs_by_sequence']),
                            ('protein clusters', len(clusters),
                             pubd['protein_clusters'])):
        if got != want:
            die('the deduplication does not reproduce 09_dedup_census.py: '
                '%s %d against %d' % (name, got, want))
    print('deduplication reproduces 09_dedup_census.py on all five totals')

    per_cluster = collections.Counter()
    for key, nroutes in best.items():
        per_cluster[find(key[0])] += nroutes
    counts = sorted(per_cluster.values(), reverse=True)
    if sum(counts) != total:
        die('the per-cluster counts do not sum to the deduplicated total')

    tally = collections.Counter()
    per = collections.Counter()
    xz = collections.Counter()
    opened = 0
    for key in best:
        X, Z = fan_res[key]
        if X not in DIST_SET or Z not in DIST_SET:
            continue
        dxz = DIST[(X, Z)]
        xz[dxz] += 1
        for Y in sorted(AA20):
            if Y in (X, Z):
                continue
            dxy, dzy = DIST[(X, Y)], DIST[(Z, Y)]
            tally[dxy] += 1
            per[(dxy, dxz, dzy)] += 1
            if dxy > max(dxz, dzy):
                opened += 1
    tot = sum(tally.values())
    if tot == 0:
        die('no target slots were counted')

    code_off = [DIST[(x, y)] for x in sorted(AA20) for y in sorted(AA20) if x != y]
    code = {
        'basis': ('deduplicated route fans, one per sequence pair, position and '
                  'pair of parental residues, as in 09_dedup_census.py'),
        'route_fans': len(best),
        'route_fans_raw_census_rows': len(raw_fans),
        'target_slots': tot,
        'direct_base_changes': {str(k): v for k, v in sorted(tally.items())},
        'direct_two_or_more': sum(v for k, v in tally.items() if k >= 2),
        'direct_two_or_more_fraction': round(
            sum(v for k, v in tally.items() if k >= 2) / tot, 4),
        'first_leg_base_changes': {str(k): v for k, v in sorted(xz.items())},
        'opened_by_the_route': opened,
        'opened_fraction': round(opened / tot, 4),
        'two_via_one_and_one': per.get((2, 1, 1), 0),
        'two_via_one_and_one_fraction': round(per.get((2, 1, 1), 0) / tot, 4),
        'code_pairs_total': len(code_off),
        'code_pairs_two_or_more': sum(1 for d in code_off if d >= 2),
        'code_distribution': {str(k): v for k, v in
                              sorted(collections.Counter(code_off).items())},
    }
    with open(os.path.join(DATA, 'codon_paths.json'), 'w') as fh:
        json.dump(code, fh, indent=1, sort_keys=True)

    prof = {'clusters': len(counts), 'routes': total,
            'routes_per_cluster_sorted': counts,
            'largest_cluster_routes': counts[0],
            'largest_cluster_share': round(counts[0] / total, 4),
            'median_cluster_routes': counts[len(counts) // 2]}
    with open(os.path.join(DATA, 'dedup_profile.json'), 'w') as fh:
        json.dump(prof, fh, indent=1, sort_keys=True)

    print('%d clusters carry %d routes; the largest carries %d, %.1f per cent'
          % (len(counts), total, counts[0], 100.0 * counts[0] / total))
    print('target slots %d; direct needs two or more base changes in %d, %.1f per cent'
          % (tot, code['direct_two_or_more'],
             100.0 * code['direct_two_or_more_fraction']))
    print('two by one and one: %d, %.1f per cent' %
          (code['two_via_one_and_one'], 100.0 * code['two_via_one_and_one_fraction']))
    print('wrote data/codon_paths.json and data/dedup_profile.json')


DIST_SET = set(AA20)

if __name__ == '__main__':
    main()
