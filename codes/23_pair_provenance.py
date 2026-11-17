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

PHENO = [
    ('binding', ('binding', 'affinity', 'bind ', 'ace2')),
    ('expression', ('expression', 'abundance', 'stability', 'folding',
                    'protein level', 'degradation')),
    ('growth', ('growth', 'fitness', 'proliferation', 'survival',
                'competition', 'selection coefficient')),
    ('activity', ('activity', 'catalytic', 'enzym', 'function',
                  'transcription', 'signal', 'transport', 'channel')),
    ('resistance', ('resistance', 'resistant', 'drug', 'inhibitor',
                    'antibiotic', 'mic ', 'ic50')),
    ('splicing', ('splicing', 'splice', 'exon')),
    ('localization', ('localization', 'localisation', 'traffick', 'surface')),
]


def die(m):
    sys.exit('FATAL 23_pair_provenance.py: %s' % m)


def pheno(rec):
    t = ((rec.get('title') or '') + ' ' + (rec.get('short') or '')).lower()
    for name, words in PHENO:
        if any(w in t for w in words):
            return name
    return None


def main():
    prov_path = os.path.join(DATA, 'score_set_provenance.json')
    if not os.path.exists(prov_path):
        die('data/score_set_provenance.json is absent; run the census first')
    prov = json.load(open(prov_path))['provenance']

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

    rows = [r for r in csv.DictReader(open(os.path.join(DATA,
                                                        'route_census_mavedb.csv')))
            if int(r['complete_routes']) >= 1]
    if not rows:
        die('no census row carries a complete route')

    pairs, onepos = {}, set()
    for r in rows:
        a, b = r['from_urn'], r['via_urn']
        pairs.setdefault((a, b), 0)
        pairs[(a, b)] += 1
        if a in seq and b in seq and len(seq[a]) == len(seq[b]):
            d = [i for i, (x, y) in enumerate(zip(seq[a], seq[b])) if x != y]
            if len(d) == 1 and d[0] == int(r['pos']) - 1:
                onepos.add((a, b))

    missing = [u for p in pairs for u in p if u not in prov]
    if missing:
        die('%d deposits in the census carry no provenance record' % len(missing))

    def relation(a, b):
        pa, pb = prov[a], prov[b]
        if pa['experiment'] == pb['experiment']:
            return 'same experiment'
        if pa['experiment_set'] == pb['experiment_set']:
            return 'same experiment set, different experiment'
        if set(pa.get('publications') or []) & set(pb.get('publications') or []):
            return 'shared publication, different experiment set'
        return 'unrelated deposits'

    order = ['same experiment', 'same experiment set, different experiment',
             'shared publication, different experiment set', 'unrelated deposits']

    def tally(keys):
        rel = collections.Counter(relation(a, b) for a, b in keys)
        ph = collections.Counter()
        for a, b in keys:
            x, y = pheno(prov[a]), pheno(prov[b])
            if x is None or y is None:
                ph['one or both unclassified'] += 1
            elif x == y:
                ph['same phenotype word: ' + x] += 1
            else:
                ph['different phenotype words'] += 1
        return (collections.OrderedDict((k, rel.get(k, 0)) for k in order),
                collections.OrderedDict(sorted(ph.items())))

    allrel, allph = tally(pairs.keys())
    onerel, oneph = tally(onepos)

    clean = sorted((a, b) for (a, b) in onepos
                   if relation(a, b) == 'same experiment')
    cleanrec = [{'from': a, 'via': b,
                 'experiment': prov[a]['experiment'],
                 'rows_in_census': pairs[(a, b)],
                 'title_from': prov[a]['title'], 'title_via': prov[b]['title'],
                 'publications': sorted(set(prov[a].get('publications') or []))}
                for a, b in clean]

    byname = collections.defaultdict(lambda: {'sets': set(), 'pubs': set(),
                                              'deposits': set(),
                                              'sets_with_pub': set(),
                                              'deposits_with_pub': set()})
    for r in rows:
        a2, b2 = r['from_urn'], r['via_urn']
        if (a2, b2) not in onepos:
            continue
        nm = r['from_name'].split('(')[0].strip()
        rec = byname[nm]
        for u in (a2, b2):
            rec['sets'].add(prov[u]['experiment_set'])
            ps = prov[u].get('publications') or []
            rec['pubs'].update(ps)
            rec['deposits'].add(u)
            if ps:
                rec['deposits_with_pub'].add(u)
                rec['sets_with_pub'].add(prov[u]['experiment_set'])
    provenance_by_target = {
        k: {'experiment_sets': sorted(v['sets']),
            'publications': sorted(v['pubs']),
            'experiment_sets_with_a_publication': sorted(v['sets_with_pub']),
            'deposits': len(v['deposits']),
            'deposits_with_a_publication': len(v['deposits_with_pub'])}
        for k, v in sorted(byname.items())}

    lens, ident = [], collections.Counter()
    meta = {r['scoreSetUrn']: r for r in
            json.load(open(os.path.join(DATA, 'target_genes.json')))}
    for a2, b2 in sorted({(min(x, y), max(x, y)) for (x, y) in onepos}):
        lens.append(len(seq[a2]))
        ua = (meta[a2].get('uniprotIdFromMappedMetadata')
              or meta[a2].get('targetAccession'))
        ub = (meta[b2].get('uniprotIdFromMappedMetadata')
              or meta[b2].get('targetAccession'))
        if ua and ub:
            ident['same accession' if ua == ub else 'different accessions'] += 1
        else:
            ident['one or both unidentified'] += 1
    lens.sort()
    target_extent = {
        'pairs': len(lens),
        'median_target_length': lens[len(lens) // 2] if lens else None,
        'target_length_q1': lens[len(lens) // 4] if lens else None,
        'target_length_q3': lens[3 * len(lens) // 4] if lens else None,
        'pairs_with_target_100_or_fewer': sum(1 for x in lens if x <= 100),
        'accession_match': dict(sorted(ident.items()))}

    und_all = {(min(a2, b2), max(a2, b2)) for (a2, b2) in pairs}
    und_one = {(min(a2, b2), max(a2, b2)) for (a2, b2) in onepos}
    und_exp = {(min(a2, b2), max(a2, b2)) for (a2, b2) in clean}

    out = {'target_extent': target_extent,
           'deposit_pairs_undirected': len(und_all),
           'deposit_pairs_at_one_differing_position_undirected': len(und_one),
           'one_position_and_one_experiment_undirected': len(und_exp),
           'provenance_by_target': provenance_by_target,
           'deposit_pairs': len(pairs),
           'deposit_pairs_at_one_differing_position': len(onepos),
           'phenotype_word_classes': [k for k, _ in PHENO],
           'all_pairs_by_relation': allrel,
           'all_pairs_by_phenotype': allph,
           'one_position_pairs_by_relation': onerel,
           'one_position_pairs_by_phenotype': oneph,
           'one_position_and_one_experiment': len(clean),
           'one_position_and_one_experiment_pairs': cleanrec}

    with open(os.path.join(DATA, 'pair_provenance.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print('deposit pairs in the census: %d directed, %d undirected'
          % (len(pairs), len(und_all)))
    for k in order:
        print('  %-46s %5d  (%d at one differing position)'
              % (k, allrel[k], onerel[k]))
    print()
    print('phenotype words, all pairs:')
    for k, v in allph.items():
        print('  %-40s %5d' % (k, v))
    print()
    print('phenotype words, pairs at one differing position:')
    for k, v in oneph.items():
        print('  %-40s %5d' % (k, v))
    print()
    print('pairs that are both at one differing position and from one '
          'experiment: %d' % len(clean))
    for r in cleanrec[:12]:
        print('  %s / %s  %s' % (r['from'], r['via'], r['experiment']))
    print('wrote data/pair_provenance.json')


if __name__ == '__main__':
    main()
