#!/usr/bin/env python3
import collections
import csv
import itertools
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
AA20 = list('ACDEFGHIKLMNPQRSTVWY')
THREE = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
         'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
         'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
         'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
_B = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(itertools.product(_B, _B, _B))}
PRO = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')


def die(m):
    sys.exit('FATAL 05_residuals.py: %s' % m)


def main():
    log = []

    def say(s):
        print(s)
        log.append(s)

    man = json.load(open(os.path.join(DATA, 'fetch_manifest.json')))
    prot = {}
    for bg, t in man['targets'].items():
        nt = open(os.path.join(ROOT, t['file'])).read().strip()
        prot[bg] = ''.join(CODON[nt[i:i + 3]] for i in range(0, len(nt), 3))

    S = {}
    for e in man['sets']:
        dcol = 'average' if 'average' in e['columns'] else 'score'
        d = {}
        for r in csv.DictReader(open(os.path.join(ROOT, e['file']))):
            m = PRO.match((r.get('hgvs_pro') or '').strip())
            if not m:
                die('%s carries an unparsed hgvs_pro' % e['urn'])
            w, pos, alt = THREE.get(m.group(1)), int(m.group(2)), THREE.get(m.group(3))
            if alt not in AA20:
                continue
            if prot[e['background']][pos - 1] != w:
                die('%s reference residue disagrees with the target at %d' % (e['urn'], pos))
            if alt == w:
                continue
            v = (r.get(dcol) or '').strip()
            if v in ('', 'NA'):
                continue
            d[(pos, alt)] = float(v)
        S[(e['background'], e['phenotype'])] = d
        say('%-6s %-11s %4d substitutions from column %s'
            % (e['background'], e['phenotype'], len(d), dcol))

    cs = json.load(open(os.path.join(DATA, 'censoring_sweep.json')))
    clip = {}
    for k, v in cs['per_set'].items():
        bg, ph = k.split('|')
        clip[(bg, ph)] = (v['detected_clip'] if v.get('detected_clip') is not None
                          else v['minimum'])

    bgpos = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(DATA, 'background_positions.csv'))):
        bgpos[(r['bg_a'], r['bg_b'])].append(
            (int(r['pos_rbd']), int(r['pos_spike']), r['res_a'], r['res_b']))

    phen = sorted({p for _, p in S})
    res, ctl = [], []
    for (a, b), items in sorted(bgpos.items()):
        for ph in phen:
            ka, kb = (a, ph), (b, ph)
            ca, cb = clip[ka], clip[kb]
            for pos, sp, X, Z in sorted(items):
                leg1 = S[ka].get((pos, Z))
                if leg1 is None:
                    continue
                for Y in sorted(AA20):
                    if Y in (X, Z):
                        continue
                    one, leg2 = S[ka].get((pos, Y)), S[kb].get((pos, Y))
                    if one is None or leg2 is None:
                        continue
                    res.append({
                        'pair': '%s|%s' % (a, b), 'phenotype': ph,
                        'pos_rbd': pos, 'pos_spike': sp,
                        'res_from': X, 'res_via': Z, 'res_to': Y,
                        'one_step': round(one, 6), 'first_leg': round(leg1, 6),
                        'second_leg': round(leg2, 6),
                        'residual': round(one - (leg1 + leg2), 6),
                        'min_leg': round(min(one, leg1, leg2), 6),
                        'margin_to_clip': round(min(one - ca, leg1 - ca, leg2 - cb), 6),
                    })
            same = [p for p in range(1, len(prot[a]) + 1)
                    if prot[a][p - 1] == prot[b][p - 1]]
            for p in same:
                X = prot[a][p - 1]
                for Y in sorted(AA20):
                    if Y == X:
                        continue
                    va, vb = S[ka].get((p, Y)), S[kb].get((p, Y))
                    if va is None or vb is None:
                        continue
                    ctl.append({
                        'pair': '%s|%s' % (a, b), 'phenotype': ph, 'pos_rbd': p,
                        'res_from': X, 'res_to': Y,
                        'effect_a': round(va, 6), 'effect_b': round(vb, 6),
                        'discrepancy': round(vb - va, 6),
                        'margin_to_clip': round(min(va - ca, vb - cb), 6),
                    })

    if not res or not ctl:
        die('no residuals or no controls were produced')
    for name, rows in (('residuals.csv', res), ('controls.csv', ctl)):
        with open(os.path.join(DATA, name), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    say('')
    for ph in phen:
        say('%-11s residuals %4d, controls %5d'
            % (ph, sum(1 for r in res if r['phenotype'] == ph),
               sum(1 for r in ctl if r['phenotype'] == ph)))
    say('residuals written: %d; controls written: %d' % (len(res), len(ctl)))
    with open(os.path.join(DATA, 'residuals_log.txt'), 'w') as fh:
        fh.write('\n'.join(log) + '\n')
    print('\nwritten to data/residuals.csv and data/controls.csv; no test applied')


if __name__ == '__main__':
    main()
