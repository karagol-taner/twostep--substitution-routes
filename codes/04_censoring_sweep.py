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
MARGINS = [0.0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.50]
TAIL_FRACTION = 0.05
MIN_REPEATS = 20
THREE = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
         'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
         'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
         'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
_B = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(itertools.product(_B, _B, _B))}
PRO = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')


def die(m):
    sys.exit('FATAL 04_censoring_sweep.py: %s' % m)


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

    S, meta = {}, {}
    for e in man['sets']:
        cols = e['columns']
        dcol = 'average' if 'average' in cols else 'score'
        d = {}
        for r in csv.DictReader(open(os.path.join(ROOT, e['file']))):
            m = PRO.match((r.get('hgvs_pro') or '').strip())
            if not m:
                die('%s carries an unparsed hgvs_pro' % e['urn'])
            pos, alt = int(m.group(2)), THREE.get(m.group(3))
            if alt not in AA20:
                continue
            v = (r.get(dcol) or '').strip()
            if v in ('', 'NA'):
                continue
            d[(pos, alt)] = float(v)
        S[(e['background'], e['phenotype'])] = d
        meta[(e['background'], e['phenotype'])] = {'urn': e['urn'], 'delta_column': dcol}

    clip = {}
    say('lower tail of each score set, counted on unique substitutions')
    for k in sorted(S):
        vals = sorted(S[k].values())
        n = len(vals)
        tail = vals[:max(1, int(n * TAIL_FRACTION))]
        rep = collections.Counter(round(x, 3) for x in tail)
        top = rep.most_common(3)
        detected = next((v for v, c in sorted(rep.items()) if c >= MIN_REPEATS), None)
        clip[k] = {'minimum': round(vals[0], 4),
                   'detected_clip': detected,
                   'tail_repeats': [[v, c] for v, c in top],
                   'n_measured': n,
                   'n_at_detected_clip': (rep[detected] if detected is not None else 0)}
        say('  %-6s %-11s min %7.3f  tail repeats %-34s  clip %s'
            % (k[0], k[1], vals[0],
               ', '.join('%.3f x%d' % (v, c) for v, c in top),
               ('%.3f, n=%d' % (detected, rep[detected])) if detected is not None
               else 'none detected'))

    bgpos = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(DATA, 'background_positions.csv'))):
        bgpos[(r['bg_a'], r['bg_b'])].append((int(r['pos_rbd']), r['res_a'], r['res_b']))

    def threshold(k, margin):
        c = clip[k]['detected_clip']
        return (c if c is not None else clip[k]['minimum']) + margin

    phen = sorted({p for _, p in S})
    sweep = []
    say('')
    say('sweep of the exclusion margin')
    for margin in MARGINS:
        row = {'margin': margin, 'by_phenotype': {}}
        for ph in phen:
            tot = exc = 0
            for (a, b), items in sorted(bgpos.items()):
                ka, kb = (a, ph), (b, ph)
                ta, tb = threshold(ka, margin), threshold(kb, margin)
                for pos, X, Z in items:
                    leg1 = S[ka].get((pos, Z))
                    for Y in sorted(AA20):
                        if Y in (X, Z):
                            continue
                        one, leg2 = S[ka].get((pos, Y)), S[kb].get((pos, Y))
                        if None in (leg1, one, leg2):
                            continue
                        tot += 1
                        if one <= ta or leg1 <= ta or leg2 <= tb:
                            exc += 1
            row['by_phenotype'][ph] = {'routes': tot, 'excluded': exc, 'usable': tot - exc}
        sweep.append(row)
        say('  margin %4.2f  ' % margin + '  '.join(
            '%s %d of %d excluded, %d usable'
            % (ph, row['by_phenotype'][ph]['excluded'], row['by_phenotype'][ph]['routes'],
               row['by_phenotype'][ph]['usable']) for ph in phen))

    out = {
        'tail_fraction_searched': TAIL_FRACTION,
        'min_repeats_to_call_a_clip': MIN_REPEATS,
        'margins': MARGINS,
        'clip_value_source': ('detected from the deposited scores, not sourced '
                              'from a published limit of detection'),
        'per_set': {'%s|%s' % k: dict(v, **meta[k]) for k, v in clip.items()},
        'sweep': sweep,
    }
    with open(os.path.join(DATA, 'censoring_sweep.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'censoring_sweep_log.txt'), 'w') as fh:
        fh.write('\n'.join(log) + '\n')
    print('\nwritten to data/censoring_sweep.json')


if __name__ == '__main__':
    main()
