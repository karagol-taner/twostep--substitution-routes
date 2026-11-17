#!/usr/bin/env python3
import csv
import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')

BACKGROUNDS = ['Wuhan', 'BA.1', 'BA.2']
BA1_SPIKE = ['G339D', 'S371L', 'S373P', 'S375F', 'K417N', 'N440K', 'G446S',
             'S477N', 'T478K', 'E484A', 'Q493R', 'G496S', 'Q498R', 'N501Y', 'Y505H']

_B = 'TCAG'
_AAS = ('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
CODON = {a + b + c: _AAS[i] for i, (a, b, c)
         in enumerate(itertools.product(_B, _B, _B))}


def die(m):
    sys.exit('FATAL 02_background_map.py: %s' % m)


def translate(nt):
    try:
        return ''.join(CODON[nt[i:i + 3]] for i in range(0, len(nt), 3))
    except KeyError as e:
        die('untranslatable codon %s' % e)


def main():
    log = []

    def say(s):
        print(s)
        log.append(s)

    prot = {}
    for bg in BACKGROUNDS:
        p = os.path.join(DATA, 'target_%s.txt' % bg.replace('.', ''))
        if not os.path.isfile(p):
            die('%s is absent; run 01_fetch_scores.py first' % os.path.relpath(p, ROOT))
        aa = translate(open(p).read().strip())
        if '*' in aa:
            die('%s translation carries a stop codon' % bg)
        prot[bg] = aa
        say('%-6s %d residues  %s...%s' % (bg, len(aa), aa[:10], aa[-5:]))
    n = len(prot[BACKGROUNDS[0]])
    if any(len(v) != n for v in prot.values()):
        die('the three targets are not the same length, so positions are not comparable')

    want = {(int(s[1:-1]), s[0], s[-1]) for s in BA1_SPIKE}
    hits = []
    for off in range(0, 600):
        got = {(i + 1 + off, prot['Wuhan'][i], prot['BA.1'][i])
               for i in range(n) if prot['Wuhan'][i] != prot['BA.1'][i]}
        if got == want:
            hits.append(off)
    if len(hits) != 1:
        die('spike offset did not resolve uniquely: candidates %s' % hits)
    off = hits[0]
    say('spike offset resolved uniquely: RBD-local i maps to spike i+%d' % off)

    rows = []
    for a, b in itertools.combinations(BACKGROUNDS, 2):
        d = [(i + 1, prot[a][i], prot[b][i]) for i in range(n)
             if prot[a][i] != prot[b][i]]
        say('%-13s %2d differing positions: %s'
            % (a + ' vs ' + b, len(d),
               ', '.join('%s%d%s' % (x, p, y) for p, x, y in d)))
        for p, x, y in d:
            rows.append({'pair': '%s|%s' % (a, b), 'bg_a': a, 'bg_b': b,
                         'pos_rbd': p, 'pos_spike': p + off,
                         'res_a': x, 'res_b': y})

    pos = sorted({r['pos_rbd'] for r in rows})
    say('distinct differing positions: %d' % len(pos))
    say('  RBD-local: %s' % ', '.join(str(p) for p in pos))
    say('  spike    : %s' % ', '.join(str(p + off) for p in pos))
    say('non-differing positions, identical parent in all three: %d' % (n - len(pos)))

    three = [p for p in pos if len({prot[b][p - 1] for b in BACKGROUNDS}) == 3]
    say('positions with three distinct parental residues: %d  %s'
        % (len(three), ', '.join('%d (spike %d, %s)'
                                 % (p, p + off,
                                    '/'.join(prot[b][p - 1] for b in BACKGROUNDS))
                                 for p in three)))

    out = os.path.join(DATA, 'background_positions.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['pair', 'bg_a', 'bg_b', 'pos_rbd',
                                           'pos_spike', 'res_a', 'res_b'])
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(DATA, 'background_map_log.txt'), 'w') as fh:
        fh.write('\n'.join(log) + '\n')
    print('\nwritten to data/background_positions.csv (%d rows)' % len(rows))


if __name__ == '__main__':
    main()
