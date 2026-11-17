#!/usr/bin/env python3
import collections
import csv
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(
    os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE, 'data')
META = os.path.join(HERE, 'metacache')
API = 'https://api.mavedb.org/api/v1/score-sets/'


def die(m):
    sys.exit('FATAL 13_common_scale_subset.py: %s' % m)


def fetch(urn):
    p = os.path.join(META, urn.replace(':', '_') + '.json')
    if os.path.isfile(p) and os.path.getsize(p) > 80:
        return p
    try:
        with urllib.request.urlopen(API + urn, timeout=120) as r:
            b = r.read()
        with open(p, 'wb') as fh:
            fh.write(b)
    except Exception as e:
        return None
    return p


def pubs(d):
    out = []
    for key in ('primaryPublicationIdentifiers', 'publicationIdentifiers',
                'secondaryPublicationIdentifiers'):
        for x in (d.get(key) or []):
            if not isinstance(x, dict):
                continue
            ident = x.get('identifier') or x.get('doi') or x.get('pubmedId')
            db = x.get('dbName') or ''
            if ident:
                out.append('%s:%s' % (db, ident))
    return sorted(set(out))


def main():
    os.makedirs(META, exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(DATA, 'route_census_mavedb.csv'))))
    urns = sorted({r['from_urn'] for r in rows} | {r['via_urn'] for r in rows})
    print('score sets in the census: %d' % len(urns), flush=True)
    done = [0]
    with ThreadPoolExecutor(max_workers=8) as ex:
        for _ in ex.map(fetch, urns):
            done[0] += 1
            if done[0] % 100 == 0:
                print('  fetched %d of %d' % (done[0], len(urns)), flush=True)
    print('  fetched %d of %d' % (done[0], len(urns)), flush=True)

    prov, missing = {}, []
    for u in urns:
        p = os.path.join(META, u.replace(':', '_') + '.json')
        if not os.path.isfile(p):
            missing.append(u)
            continue
        try:
            d = json.load(open(p))
        except Exception:
            missing.append(u)
            continue
        if 'detail' in d:
            missing.append(u)
            continue
        exp = d.get('experiment') or {}
        prov[u] = {
            'experiment': exp.get('urn'),
            'experiment_set': exp.get('experimentSetUrn'),
            'publications': pubs(d) or pubs(exp),
            'title': d.get('title'),
            'short': d.get('shortDescription'),
        }
    print('provenance recorded: %d; unreadable: %d' % (len(prov), len(missing)), flush=True)
    with open(os.path.join(DATA, 'score_set_provenance.json'), 'w') as fh:
        json.dump({'provenance': prov, 'unreadable': missing}, fh, indent=1, sort_keys=True)

    cnt = collections.Counter()
    out = []
    seen = set()
    for r in rows:
        if int(r['complete_routes']) < 6:
            continue
        a, b = r['from_urn'], r['via_urn']
        k = (a, b)
        if k in seen:
            continue
        seen.add(k)
        pa, pb = prov.get(a), prov.get(b)
        if pa is None or pb is None:
            tier = 'provenance unknown'
        elif pa['experiment'] and pa['experiment'] == pb['experiment']:
            tier = 'same experiment'
        elif pa['experiment_set'] and pa['experiment_set'] == pb['experiment_set']:
            tier = 'same experiment set'
        elif pa['publications'] and set(pa['publications']) & set(pb['publications']):
            tier = 'same publication'
        else:
            tier = 'different provenance'
        cnt[tier] += 1
        out.append({'from_urn': a, 'via_urn': b, 'tier': tier,
                    'from_experiment_set': (pa or {}).get('experiment_set'),
                    'via_experiment_set': (pb or {}).get('experiment_set'),
                    'shared_publications': ';'.join(
                        sorted(set((pa or {}).get('publications') or [])
                               & set((pb or {}).get('publications') or [])))})
    with open(os.path.join(DATA, 'common_scale_pairs.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['from_urn', 'via_urn', 'tier',
                                            'from_experiment_set',
                                            'via_experiment_set',
                                            'shared_publications'])
        w.writeheader()
        w.writerows(out)
    print('', flush=True)
    print('deposit pairs in the census with at least six routes: %d' % len(out), flush=True)
    for t, n in cnt.most_common():
        print('  %-24s %6d  %5.1f%%' % (t, n, 100.0 * n / len(out)), flush=True)
    common = sum(n for t, n in cnt.items()
                 if t in ('same experiment', 'same experiment set', 'same publication'))
    print('', flush=True)
    print('common-scale pairs, any tier: %d of %d (%.1f%%)'
          % (common, len(out), 100.0 * common / len(out)), flush=True)
    print('written to score_set_provenance.json and common_scale_pairs.csv', flush=True)


if __name__ == '__main__':
    main()
