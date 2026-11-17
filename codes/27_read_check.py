#!/usr/bin/env python3
import collections
import csv
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')


def die(m):
    sys.exit('FATAL 27_read_check.py: %s' % m)


def main():
    spec = importlib.util.spec_from_file_location(
        'census', os.path.join(HERE, '08_mavedb_route_census.py'))
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)

    seq = {}
    for r in json.load(open(os.path.join(DATA, 'target_genes.json'))):
        ts = r.get('targetSequence') or {}
        s = (ts.get('sequence') or '').upper()
        t = ts.get('sequenceType')
        p = (C.translate(s) if (t == 'dna' and len(s) % 3 == 0)
             else (s if t == 'protein' else None))
        if p:
            seq[r['scoreSetUrn']] = p.rstrip('*')

    urns = sorted(u for u in seq
                  if os.path.isfile(os.path.join(C.CACHE,
                                                 u.replace(':', '_') + '.csv')))
    if not urns:
        die('no deposit is cached')
    print('deposits cached and translatable: %d' % len(urns))

    ref_ok = ref_bad = ref_off = 0
    bad_by_urn = collections.Counter()
    worst = []
    cols = collections.Counter()
    both_cols = []
    drop_shapes = collections.Counter()
    rows_total = rows_kept = 0
    offsets = collections.Counter()

    for i, u in enumerate(urns, 1):
        path = os.path.join(C.CACHE, u.replace(':', '_') + '.csv')
        text = open(path, encoding='utf-8', errors='replace').read()
        head = text.split('\n', 1)[0].strip().split(',')
        dcol = ('average' if 'average' in head
                else ('score' if 'score' in head else None))
        cols[dcol or '(none)'] += 1
        if dcol is None:
            continue
        other = 'score' if dcol == 'average' and 'score' in head else None

        s = seq[u]
        ok = bad = 0
        diffs = []
        shift = collections.Counter()
        for r in csv.DictReader(text.splitlines()):
            raw = (r.get('hgvs_pro') or '').strip()
            rows_total += 1
            mm = C.PRO.match(raw)
            if not mm:
                if raw:
                    shape = (raw[:2] + '...' if len(raw) > 14 else raw)
                    drop_shapes[('unparsed', shape if not raw[2:3].isalpha()
                                 else 'p.' + ('multi' if ';' in raw or ',' in raw
                                              else 'other'))] += 1
                else:
                    drop_shapes[('empty', '')] += 1
                continue
            pos = int(mm.group(2))
            refaa = C.THREE.get(mm.group(1))
            alt = C.THREE.get(mm.group(3))
            if alt not in C.AA20:
                drop_shapes[('alt not standard', mm.group(3))] += 1
                continue
            v = (r.get(dcol) or '').strip()
            if v in ('', 'NA'):
                drop_shapes[('no value', dcol)] += 1
                continue
            rows_kept += 1
            if other is not None:
                w = (r.get(other) or '').strip()
                if w not in ('', 'NA'):
                    try:
                        if abs(float(w) - float(v)) > 1e-9:
                            diffs.append(abs(float(w) - float(v)))
                    except ValueError:
                        pass
            if refaa is None or pos < 1 or pos > len(s):
                bad += 1
                continue
            if s[pos - 1] == refaa:
                ok += 1
            else:
                bad += 1
                for d in range(-5, 6):
                    q = pos - 1 + d
                    if 0 <= q < len(s) and s[q] == refaa:
                        shift[d] += 1
        ref_ok += ok
        ref_bad += bad
        if bad:
            bad_by_urn[u] = bad
            best = shift.most_common(1)
            offsets[best[0][0] if best else None] += 1
            worst.append((u, ok, bad, best[0] if best else None))
        if diffs:
            both_cols.append((u, len(diffs), round(max(diffs), 4)))
        if i % 150 == 0 or i == len(urns):
            print('  read %d of %d' % (i, len(urns)))

    worst.sort(key=lambda t: -t[2])
    out = {'deposits_checked': len(urns),
           'variant_rows_seen': rows_total,
           'variant_rows_kept': rows_kept,
           'reference_residue_matches': ref_ok,
           'reference_residue_mismatches': ref_bad,
           'deposits_with_any_mismatch': len(bad_by_urn),
           'effect_column_used': dict(cols),
           'deposits_where_a_second_column_differs': len(both_cols),
           'second_column_examples': both_cols[:10],
           'discarded_rows': {('%s %s' % k).strip(): v
                              for k, v in drop_shapes.most_common(12)},
           'mismatch_offsets_if_any': {str(k): v for k, v in offsets.items()},
           'worst_deposits': [{'urn': u, 'matched': o, 'mismatched': b,
                               'commonest_shift': (None if sh is None
                                                   else {'offset': sh[0],
                                                         'rows': sh[1]})}
                              for u, o, b, sh in worst[:15]]}
    with open(os.path.join(DATA, 'read_check.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print()
    print('variant rows seen %s, kept %s'
          % ('{:,}'.format(rows_total), '{:,}'.format(rows_kept)))
    print('reference residue agrees with the target sequence in %s rows'
          % '{:,}'.format(ref_ok))
    print('reference residue disagrees in %s rows, across %d of %d deposits'
          % ('{:,}'.format(ref_bad), len(bad_by_urn), len(urns)))
    print('effect column taken: %s' % dict(cols))
    print('deposits where another column in the same file differs: %d'
          % len(both_cols))
    print('discarded rows: %s' % out['discarded_rows'])
    if worst:
        print()
        print('deposits with the most mismatches:')
        for w in out['worst_deposits'][:10]:
            print('  %-28s matched %6d, mismatched %6d, commonest shift %s'
                  % (w['urn'], w['matched'], w['mismatched'],
                     w['commonest_shift']))
    print('wrote data/read_check.json')


if __name__ == '__main__':
    main()
