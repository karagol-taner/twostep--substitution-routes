#!/usr/bin/env python3
import collections
import csv
import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
STRUCT = os.path.join(ROOT, 'structures')

B = 'TCAG'
AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: AAS[i] for i, (a, b, c) in enumerate(itertools.product(B, B, B))}
THREE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
         'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
         'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
         'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
FIRST_RESIDUE = 331
CONTACT = 5.0


def die(m):
    sys.exit('FATAL 18_interface_distance.py: %s' % m)


def read_pdb(path):
    ca = collections.defaultdict(dict)
    allat = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.startswith('ATOM'):
                continue
            if line[16] not in (' ', 'A'):
                continue
            if line[76:78].strip() == 'H':
                continue
            chain, res = line[21], int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            allat[chain].append((res, xyz))
            if line[12:16].strip() == 'CA':
                ca[chain][res] = (THREE.get(line[17:20].strip(), 'X'), xyz)
    return ca, allat


def solve_chain(ca, seq):
    best = None
    for c in ca:
        for off in range(-30, 31):
            hit = sum(1 for i, a in enumerate(seq)
                      if ca[c].get(i + FIRST_RESIDUE + off, ('', None))[0] == a)
            if best is None or hit > best[0]:
                best = (hit, c, off)
    hit, chain, off = best
    if hit < 0.9 * len(seq):
        die('the deposited sequence does not align into the structure: %d of %d'
            % (hit, len(seq)))
    return chain, off, hit


def exact_ranksum_p(dep, oth):
    pool = list(dep) + list(oth)
    n, k = len(pool), len(dep)
    order = sorted(range(n), key=lambda i: pool[i])
    rank = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and pool[order[j + 1]] == pool[order[i]]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for t in range(i, j + 1):
            rank[order[t]] = r
        i = j + 1
    obs = sum(rank[i] for i in range(k))
    tot = 0
    hit = 0
    for comb in itertools.combinations(range(n), k):
        tot += 1
        if sum(rank[i] for i in comb) <= obs + 1e-9:
            hit += 1
    return hit / float(tot), tot


def main():
    pdb = os.path.join(STRUCT, '7WPB.pdb')
    for p in (pdb, os.path.join(DATA, 'target_BA1.txt'),
              os.path.join(DATA, 'additivity_test.csv'),
              os.path.join(DATA, 'rank_test.csv')):
        if not os.path.isfile(p):
            die('%s is absent' % p)

    nt = open(os.path.join(DATA, 'target_BA1.txt')).read().strip()
    seq = ''.join(CODON[nt[i:i + 3]] for i in range(0, len(nt), 3))
    ca, allat = read_pdb(pdb)
    chain, off, hit = solve_chain(ca, seq)
    other = [c for c in ca if c != chain]
    if len(other) != 1:
        die('expected one partner chain, found %s' % other)
    partner = other[0]
    print('chain %s is the domain, offset %+d, %d of %d residues match; '
          'partner chain %s' % (chain, off, hit, len(seq), partner))

    mag = [r for r in csv.DictReader(open(os.path.join(DATA, 'additivity_test.csv')))
           if float(r['margin']) == 0.0]
    rnk = list(csv.DictReader(open(os.path.join(DATA, 'rank_test.csv'))))
    pos = sorted({int(r['pos_spike']) for r in mag})

    pxyz = np.array([x for _, x in allat[partner]])
    dist = {}
    for p in pos:
        at = [x for r, x in allat[chain] if r == p + off]
        if not at:
            die('position %d is absent from the structure' % p)
        a = np.array(at)
        dist[p] = float(np.sqrt(((a[:, None, :] - pxyz[None, :, :]) ** 2).sum(-1)).min())

    rec = {'structure': '7WPB', 'chain': chain, 'partner_chain': partner,
           'offset': off, 'residues_matched': hit, 'residues_total': len(seq),
           'contact_angstrom': CONTACT,
           'distance_angstrom': {str(k): round(v, 2) for k, v in sorted(dist.items())},
           'by_phenotype': {}}

    for ph in sorted({r['phenotype'] for r in mag}):
        sub = [r for r in mag if r['phenotype'] == ph]
        pp = sorted({int(r['pos_spike']) for r in sub})
        sig = sorted({int(r['pos_spike']) for r in sub if float(r['p_bh']) < 0.05})
        dep = [dist[p] for p in sig]
        oth = [dist[p] for p in pp if p not in sig]
        d = {'positions': len(pp), 'departing': sig,
             'departing_median_angstrom': round(float(np.median(dep)), 2) if dep else None,
             'other_median_angstrom': round(float(np.median(oth)), 2) if oth else None,
             'departing_within_contact': int(sum(x < CONTACT for x in dep)),
             'other_within_contact': int(sum(x < CONTACT for x in oth))}
        if dep and oth:
            p_close, nsplit = exact_ranksum_p(dep, oth)
            p_far, _ = exact_ranksum_p(oth, dep)
            d['p_departing_closer'] = round(p_close, 6)
            d['p_departing_farther'] = round(p_far, 6)
            d['splits_enumerated'] = nsplit
        rec['by_phenotype'][ph] = d
        print('  %-11s %2d departing of %2d, median %s A against %s A, '
              '%d of %d within %.0f A against %d of %d; p closer %s, p farther %s'
              % (ph, len(sig), len(pp), d['departing_median_angstrom'],
                 d['other_median_angstrom'], d['departing_within_contact'], len(dep),
                 CONTACT, d['other_within_contact'], len(oth),
                 d.get('p_departing_closer'), d.get('p_departing_farther')))

    rsig = sorted({int(r['pos_spike']) for r in rnk if float(r['p_bh']) < 0.05})
    rec['departing_rank'] = rsig
    rec['rank_distance_angstrom'] = {str(p): round(dist[p], 2) for p in rsig}
    print('  rank statistic departs at %s, at %s A'
          % (', '.join(str(p) for p in rsig),
             ', '.join('%.2f' % dist[p] for p in rsig)))

    out = os.path.join(DATA, 'interface_distance.json')
    with open(out, 'w') as fh:
        json.dump(rec, fh, indent=1, sort_keys=True)
    print('wrote %s' % out.replace(ROOT + '/', ''))


if __name__ == '__main__':
    main()
