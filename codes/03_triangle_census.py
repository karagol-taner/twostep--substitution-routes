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
TOL = 0.005
THREE = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
         'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
         'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
         'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}
_B = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(itertools.product(_B, _B, _B))}
PRO = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')


def die(m):
    sys.exit('FATAL 03_triangle_census.py: %s' % m)


def translate(nt):
    return ''.join(CODON[nt[i:i + 3]] for i in range(0, len(nt), 3))


def load(man_entry):
    path = os.path.join(ROOT, man_entry['file'])
    if not os.path.isfile(path):
        die('%s is absent; run 01_fetch_scores.py first' % man_entry['file'])
    cols = man_entry['columns']
    delta_col = 'average' if 'average' in cols else 'score'
    abs_col = 'bind' if 'bind' in cols else ('expr' if 'expr' in cols else None)
    delta, absolute, ref = {}, {}, {}
    seen = collections.defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            m = PRO.match((r.get('hgvs_pro') or '').strip())
            if not m:
                die('%s carries an unparsed hgvs_pro %r' % (man_entry['urn'], r.get('hgvs_pro')))
            w, pos, a = THREE.get(m.group(1)), int(m.group(2)), THREE.get(m.group(3))
            if w is None or a is None:
                die('%s carries an unknown residue in %s' % (man_entry['urn'], m.group(0)))
            ref.setdefault(pos, w)
            if ref[pos] != w:
                die('%s gives two reference residues at position %d' % (man_entry['urn'], pos))
            if a not in AA20:
                continue
            v = (r.get(delta_col) or '').strip()
            key = (pos, a)
            if v in ('', 'NA'):
                seen[key].add(None)
                continue
            fv = float(v)
            if key in delta and abs(delta[key] - fv) > 1e-9:
                die('%s gives two different %s values for %s%d%s'
                    % (man_entry['urn'], delta_col, w, pos, a))
            delta[key] = fv
            if abs_col:
                av = (r.get(abs_col) or '').strip()
                if av not in ('', 'NA'):
                    absolute[key] = float(av)
    return delta, absolute, ref, delta_col, abs_col


def main():
    log = []

    def say(s):
        print(s)
        log.append(s)

    man = json.load(open(os.path.join(DATA, 'fetch_manifest.json')))
    prot = {}
    for bg, t in man['targets'].items():
        prot[bg] = translate(open(os.path.join(ROOT, t['file'])).read().strip())

    S, A, scale = {}, {}, {}
    for e in man['sets']:
        k = (e['background'], e['phenotype'])
        d, ab, ref, dcol, acol = load(e)
        S[k], A[k] = d, ab
        bad = [p for p, w in ref.items() if prot[e['background']][p - 1] != w]
        if bad:
            die('%s reference residues disagree with the translated target at %s'
                % (e['urn'], bad[:5]))
        say('%-6s %-11s %-24s %4d measured deltas from column %-8s absolute column %s'
            % (e['background'], e['phenotype'], e['urn'], len(d), dcol, acol or 'none'))
        if ab:
            offs = sorted(ab[kk] - d[kk] for kk in ab if kk in d)
            med = offs[len(offs) // 2]
            dev = max(abs(o - med) for o in offs)
            n_out = sum(1 for o in offs if abs(o - med) > TOL)
            const = n_out == 0
            scale['%s|%s' % k] = {
                'parental_absolute': round(med, 5), 'column': acol,
                'max_deviation_from_median': round(dev, 5), 'tolerance': TOL,
                'constant_offset': bool(const),
                'n_rows_outside_tolerance': n_out,
                'relation': ('absolute = parental + delta, to within the deposited rounding'
                             if const else 'offset varies by more than the tolerance'),
                'n_rows_compared': len(offs)}
            say('%-6s %-11s parental absolute %s = %.5f; %d of %d rows depart by more '
                'than %.3f, largest departure %.5f'
                % (e['background'], e['phenotype'], acol, med, n_out, len(offs), TOL, dev))

    bgpos = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(DATA, 'background_positions.csv'))):
        bgpos[(r['bg_a'], r['bg_b'])].append((int(r['pos_rbd']), int(r['pos_spike']),
                                              r['res_a'], r['res_b']))

    phen = sorted({e['phenotype'] for e in man['sets']})
    routes, controls = [], []
    for (a, b), items in sorted(bgpos.items()):
        for ph in phen:
            ka, kb = (a, ph), (b, ph)
            for pos, sp, X, Z in sorted(items):
                leg1 = S[ka].get((pos, Z))
                n_ok = n_missing = 0
                miss = []
                for Y in sorted(AA20):
                    if Y in (X, Z):
                        continue
                    one = S[ka].get((pos, Y))
                    leg2 = S[kb].get((pos, Y))
                    if leg1 is not None and one is not None and leg2 is not None:
                        n_ok += 1
                    else:
                        n_missing += 1
                        miss.append(Y)
                routes.append({'pair': '%s|%s' % (a, b), 'phenotype': ph,
                               'pos_rbd': pos, 'pos_spike': sp, 'res_from': X,
                               'res_via': Z, 'first_leg_measured': int(leg1 is not None),
                               'complete_routes': n_ok, 'incomplete_routes': n_missing,
                               'targets_missing': ''.join(miss)})
            same = [p for p in range(1, len(prot[a]) + 1)
                    if prot[a][p - 1] == prot[b][p - 1]]
            n_pair = sum(1 for p in same for Y in sorted(AA20)
                         if Y != prot[a][p - 1]
                         and S[ka].get((p, Y)) is not None
                         and S[kb].get((p, Y)) is not None)
            controls.append({'pair': '%s|%s' % (a, b), 'phenotype': ph,
                             'shared_parent_positions': len(same),
                             'control_substitution_pairs': n_pair})

    for out, rows, fields in [
            ('route_census.csv', routes,
             ['pair', 'phenotype', 'pos_rbd', 'pos_spike', 'res_from', 'res_via',
              'first_leg_measured', 'complete_routes', 'incomplete_routes', 'targets_missing']),
            ('control_census.csv', controls,
             ['pair', 'phenotype', 'shared_parent_positions', 'control_substitution_pairs'])]:
        with open(os.path.join(DATA, out), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    say('')
    for ph in phen:
        rs = [r for r in routes if r['phenotype'] == ph]
        say('%-11s positions %d, complete routes %d, incomplete %d, first leg missing at %d positions'
            % (ph, len(rs), sum(r['complete_routes'] for r in rs),
               sum(r['incomplete_routes'] for r in rs),
               sum(1 for r in rs if not r['first_leg_measured'])))
        for c in [c for c in controls if c['phenotype'] == ph]:
            say('%-11s   control pairs, %-12s %d over %d shared-parent positions'
                % ('', c['pair'], c['control_substitution_pairs'], c['shared_parent_positions']))
    say('')
    say('total complete routes, both phenotypes: %d' % sum(r['complete_routes'] for r in routes))
    say('total control pairs, both phenotypes:   %d'
        % sum(c['control_substitution_pairs'] for c in controls))

    with open(os.path.join(DATA, 'scale_check.json'), 'w') as fh:
        json.dump(scale, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'triangle_census_log.txt'), 'w') as fh:
        fh.write('\n'.join(log) + '\n')
    print('\nwritten to data/route_census.csv, data/control_census.csv, data/scale_check.json')


if __name__ == '__main__':
    main()
