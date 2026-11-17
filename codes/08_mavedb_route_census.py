#!/usr/bin/env python3
import collections
import csv
import itertools
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request

OUT = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(OUT, 'scorecache')
PAUSE = 1.0
DATA = os.path.join(
    os.path.dirname(OUT) if os.path.basename(OUT) == 'codes' else OUT, 'data')
API = 'https://api.mavedb.org/api/v1/'
MAX_DIVERGENCE = 0.25
MIN_LEN = 20
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
    print(s, flush=True)
    LOG.append(s)


def die(m):
    say('FATAL 08_mavedb_route_census.py: %s' % m)
    with open(os.path.join(DATA, 'mavedb_census_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    sys.exit(1)


def get(url, path, timeout=300, tries=8):
    if os.path.exists(path) and os.path.getsize(path) > 40:
        return open(path, 'rb').read()
    last = ''
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                b = r.read()
            with open(path, 'wb') as fh:
                fh.write(b)
            return b
        except urllib.error.HTTPError as e:
            last = 'HTTPError: %s' % e
            if e.code == 429:
                wait = 30 * (i + 1)
                try:
                    wait = max(wait, float(e.headers.get('Retry-After')))
                except (AttributeError, TypeError, ValueError):
                    pass
            else:
                wait = 2 * (i + 1)
            time.sleep(wait)
        except Exception as e:
            last = '%s: %s' % (type(e).__name__, e)
            time.sleep(2 * (i + 1))
    return ('ERROR ' + last).encode()


def translate(nt):
    out = []
    for i in range(0, len(nt) - 2, 3):
        c = CODON.get(nt[i:i + 3])
        if c is None:
            return None
        out.append(c)
    return ''.join(out)


def load_scores(urn):
    path = os.path.join(CACHE, urn.replace(':', '_') + '.csv')
    b = get(API + 'score-sets/%s/scores' % urn, path)
    if b.startswith(b'ERROR'):
        return None, b.decode()[:120]
    try:
        text = b.decode('utf-8', 'replace')
    except Exception as e:
        return None, 'undecodable: %s' % e
    head = text.split('\n', 1)[0].strip().split(',')
    dcol = 'average' if 'average' in head else ('score' if 'score' in head else None)
    if dcol is None:
        return None, 'no effect column in %s' % head[:6]
    m = {}
    for r in csv.DictReader(text.splitlines()):
        mm = PRO.match((r.get('hgvs_pro') or '').strip())
        if not mm:
            continue
        alt = THREE.get(mm.group(3))
        if alt not in AA20:
            continue
        v = (r.get(dcol) or '').strip()
        if v in ('', 'NA'):
            continue
        try:
            m[(int(mm.group(2)), alt)] = float(v)
        except ValueError:
            continue
    return m, dcol


def main():
    os.makedirs(CACHE, exist_ok=True)
    tgp = os.path.join(DATA, 'target_genes.json')
    stamp = os.path.join(DATA, 'mavedb_access.json')
    fresh = not (os.path.exists(tgp) and os.path.getsize(tgp) > 40)
    tg = get(API + 'target-genes', tgp, timeout=600)
    if tg.startswith(b'ERROR'):
        die('target-genes could not be fetched: %s' % tg.decode()[:200])
    if fresh:
        with open(stamp, 'w') as fh:
            json.dump({'target_genes_retrieved':
                       time.strftime('%Y-%m-%d', time.gmtime())}, fh, indent=1)
    if not os.path.exists(stamp):
        die('data/mavedb_access.json records the date the deposit snapshot was '
            'retrieved and is absent; delete data/target_genes.json to '
            'retrieve the snapshot again')
    access = json.load(open(stamp))['target_genes_retrieved']
    say('deposit snapshot retrieved: %s' % access)
    recs = json.loads(tg)
    say('target-gene records deposited: %d' % len(recs))

    byseq = collections.defaultdict(list)
    skipped = collections.Counter()
    for r in recs:
        ts = r.get('targetSequence') or {}
        seq = (ts.get('sequence') or '').upper()
        st = ts.get('sequenceType')
        if not seq:
            skipped['no sequence'] += 1
            continue
        if st == 'dna':
            if len(seq) % 3:
                skipped['dna length not a multiple of three'] += 1
                continue
            p = translate(seq)
            if p is None:
                skipped['untranslatable codon'] += 1
                continue
        elif st == 'protein':
            p = seq
        else:
            skipped['sequence type %s' % st] += 1
            continue
        p = p.rstrip('*')
        if '*' in p:
            skipped['internal stop'] += 1
            continue
        if len(p) < MIN_LEN:
            skipped['shorter than %d residues' % MIN_LEN] += 1
            continue
        byseq[p].append({'urn': r['scoreSetUrn'], 'name': (r.get('name') or '').strip(),
                         'uniprot': r.get('uniprotIdFromMappedMetadata'),
                         'hgnc': r.get('mappedHgncName')})
    for k, v in skipped.most_common():
        say('  set aside, %-40s %d' % (k, v))
    say('distinct target sequences usable: %d' % len(byseq))

    bylen = collections.defaultdict(list)
    for s in byseq:
        bylen[len(s)].append(s)
    pairs, unequal = [], 0
    for L, seqs in bylen.items():
        for s1, s2 in itertools.combinations(seqs, 2):
            diff = [(i + 1, s1[i], s2[i]) for i in range(L) if s1[i] != s2[i]]
            if not diff:
                continue
            if len(diff) / float(L) > MAX_DIVERGENCE:
                unequal += 1
                continue
            pairs.append({'len': L, 'diff': diff, 's1': s1, 's2': s2})
    say('background pairs at equal length and at most %.0f%% divergence: %d'
        % (MAX_DIVERGENCE * 100, len(pairs)))
    say('equal-length pairs set aside as too divergent: %d' % unequal)
    say('differing positions summed over kept pairs: %d'
        % sum(len(p['diff']) for p in pairs))

    urns = sorted({r['urn'] for p in pairs for s in (p['s1'], p['s2'])
                   for r in byseq[s]})
    say('distinct score sets to read: %d' % len(urns))
    done = 0
    for u in urns:
        path = os.path.join(CACHE, u.replace(':', '_') + '.csv')
        fresh = not (os.path.exists(path) and os.path.getsize(path) > 40)
        get(API + 'score-sets/%s/scores' % u, path)
        if fresh:
            time.sleep(PAUSE)
        done += 1
        if done % 100 == 0:
            say('  fetched %d of %d score sets' % (done, len(urns)))
    say('  fetched %d of %d score sets' % (done, len(urns)))
    scores, bad = {}, {}
    for i, u in enumerate(urns, 1):
        m, note = load_scores(u)
        if m is None:
            bad[u] = note
        else:
            scores[u] = m
        if i % 200 == 0 or i == len(urns):
            say('  parsed %d of %d score sets, %d unreadable' % (i, len(urns), len(bad)))
    say('score sets read: %d; unreadable: %d' % (len(scores), len(bad)))
    stalled = {u: n for u, n in bad.items() if n.startswith('ERROR ')}
    if stalled:
        die('%d of %d score sets could not be retrieved, so the census would be '
            'incomplete; rerun to resume from the cache. First failures: %s'
            % (len(stalled), len(urns),
               '; '.join('%s %s' % kv for kv in list(stalled.items())[:3])))

    rows, prows = [], []
    for p in pairs:
        for ra, rb in itertools.product(byseq[p['s1']], byseq[p['s2']]):
            ua, ub = ra['urn'], rb['urn']
            if ua not in scores or ub not in scores:
                continue
            A, Bm = scores[ua], scores[ub]
            tot = 0
            for direction, (U, V, Sa, Sb) in enumerate(
                    ((ua, ub, A, Bm), (ub, ua, Bm, A))):
                for pos, x, y in p['diff']:
                    X, Z = (x, y) if direction == 0 else (y, x)
                    leg1 = Sa.get((pos, Z))
                    n = 0
                    for Y in sorted(AA20):
                        if Y in (X, Z):
                            continue
                        if (leg1 is not None and Sa.get((pos, Y)) is not None
                                and Sb.get((pos, Y)) is not None):
                            n += 1
                    tot += n
                    rows.append({'from_urn': U, 'via_urn': V, 'from_name': (ra if direction == 0 else rb)['name'],
                                 'via_name': (rb if direction == 0 else ra)['name'],
                                 'length': p['len'], 'pos': pos, 'res_from': X, 'res_via': Z,
                                 'first_leg_measured': int(leg1 is not None),
                                 'complete_routes': n,
                                 'via_is_alanine': int(Z == 'A')})
            prows.append({'urn_a': ua, 'urn_b': ub, 'name_a': ra['name'], 'name_b': rb['name'],
                          'uniprot': ra.get('uniprot') or rb.get('uniprot') or '',
                          'length': p['len'], 'differing_positions': len(p['diff']),
                          'routes_both_directions': tot,
                          'differences': ' '.join('%s%d%s' % (x, i, y) for i, x, y in p['diff'])})

    usable = [r for r in rows if r['complete_routes'] > 0]
    tot_routes = sum(r['complete_routes'] for r in rows)
    ala = sum(r['complete_routes'] for r in rows if r['via_is_alanine'])
    prot_groups = len({r['uniprot'] or r['name_a'] for r in prows if r['routes_both_directions'] > 0})
    say('')
    say('position-and-direction cells examined: %d' % len(rows))
    say('cells yielding at least one route:     %d' % len(usable))
    say('complete routes, both directions:      %d' % tot_routes)
    say('of those, routed through alanine:      %d' % ala)
    say('pairs yielding at least one route:     %d of %d'
        % (sum(1 for p in prows if p['routes_both_directions'] > 0), len(prows)))
    say('distinct proteins yielding routes:     %d' % prot_groups)

    with open(os.path.join(DATA, 'route_census_mavedb.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(DATA, 'background_pairs_mavedb.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(prows[0].keys()))
        w.writeheader()
        w.writerows(sorted(prows, key=lambda r: -r['routes_both_directions']))
    with open(os.path.join(DATA, 'mavedb_census_summary.json'), 'w') as fh:
        json.dump({'target_gene_records': len(recs),
                   'set_aside': dict(skipped),
                   'distinct_sequences': len(byseq),
                   'background_pairs': len(pairs),
                   'too_divergent_pairs': unequal,
                   'target_genes_retrieved': access,
                   'score_sets_read': len(scores),
                   'score_sets_unreadable': bad,
                   'cells_examined': len(rows),
                   'cells_with_routes': len(usable),
                   'complete_routes': tot_routes,
                   'routes_via_alanine': ala,
                   'pairs_with_routes': sum(1 for p in prows if p['routes_both_directions'] > 0),
                   'pairs_examined': len(prows),
                   'distinct_proteins_with_routes': prot_groups,
                   'max_divergence': MAX_DIVERGENCE}, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'mavedb_census_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    say('')
    say('written to route_census_mavedb.csv, background_pairs_mavedb.csv, '
        'mavedb_census_summary.json')


if __name__ == '__main__':
    main()
