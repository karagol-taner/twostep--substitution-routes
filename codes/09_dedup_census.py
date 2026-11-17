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
_B = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(itertools.product(_B, _B, _B))}


def die(m):
    sys.exit('FATAL 09_dedup_census.py: %s' % m)


def translate(nt):
    out = []
    for i in range(0, len(nt) - 2, 3):
        c = CODON.get(nt[i:i + 3])
        if c is None:
            return None
        out.append(c)
    return ''.join(out)


def main():
    log = []

    def say(s):
        print(s)
        log.append(s)

    tgp = os.path.join(DATA, 'target_genes.json')
    cen = os.path.join(DATA, 'route_census_mavedb.csv')
    for p in (tgp, cen):
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
    say('score sets with a translatable target sequence: %d' % len(seq_of))

    sid = {}

    def sequence_id(u):
        s = seq_of.get(u)
        if s is None:
            return None
        return sid.setdefault(s, len(sid))

    rows = 0
    best = collections.defaultdict(int)
    via_ala = set()
    seq_pairs = set()
    unmapped = set()
    for r in csv.DictReader(open(cen)):
        rows += 1
        n = int(r['complete_routes'])
        if n == 0:
            continue
        a, b = sequence_id(r['from_urn']), sequence_id(r['via_urn'])
        if a is None or b is None:
            unmapped.add(r['from_urn'] if a is None else r['via_urn'])
            continue
        key = (a, b, r['pos'], r['res_from'], r['res_via'])
        if n > best[key]:
            best[key] = n
        if r['res_via'] == 'A':
            via_ala.add(key)
        seq_pairs.add((min(a, b), max(a, b)))

    total = sum(best.values())
    ala_total = sum(v for k, v in best.items() if k in via_ala)
    clusters = set()
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x

    for a, b in seq_pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    for a, b in seq_pairs:
        clusters.add(find(a))

    say('census rows read: %d' % rows)
    say('score sets in the census with no mappable sequence: %d' % len(unmapped))
    say('distinct route fans, by sequence pair and position: %d' % len(best))
    say('  of those routed through alanine as the intermediate: %d' % len(via_ala))
    say('deduplicated complete routes: %d' % total)
    say('  of those through alanine: %d' % ala_total)
    say('distinct background pairs by sequence: %d' % len(seq_pairs))
    say('distinct protein clusters yielding routes: %d' % len(clusters))

    out = {'census_rows': rows, 'route_fans': len(best),
           'route_fans_via_alanine': len(via_ala),
           'deduplicated_routes': total, 'deduplicated_routes_via_alanine': ala_total,
           'background_pairs_by_sequence': len(seq_pairs),
           'protein_clusters': len(clusters),
           'score_sets_unmappable': sorted(unmapped),
           'note': ('a route fan is one sequence pair at one position in one '
                    'direction; its count is the largest achieved by any of the '
                    'redundant score-set pairings that realise it')}
    with open(os.path.join(DATA, 'mavedb_census_dedup.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'mavedb_census_dedup_log.txt'), 'w') as fh:
        fh.write('\n'.join(log) + '\n')
    print('\nwritten to data/mavedb_census_dedup.json')


if __name__ == '__main__':
    main()
