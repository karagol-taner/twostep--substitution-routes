#!/usr/bin/env python3
import hashlib
import json
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
API = 'https://api.mavedb.org/api/v1/score-sets/'

SETS = [
    ('Wuhan', 'binding',    'urn:mavedb:00000044-a-1'),
    ('Wuhan', 'expression', 'urn:mavedb:00000044-b-2'),
    ('BA.1',  'binding',    'urn:mavedb:00000668-a-1'),
    ('BA.1',  'expression', 'urn:mavedb:00000669-a-1'),
    ('BA.2',  'binding',    'urn:mavedb:00000670-a-1'),
    ('BA.2',  'expression', 'urn:mavedb:00000671-a-1'),
]


def die(m):
    sys.exit('FATAL 01_fetch_scores.py: %s' % m)


def get(url, timeout=300):
    req = urllib.request.Request(url, headers={'Accept': '*/*'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            if r.status != 200:
                die('%s returned HTTP %s' % (url, r.status))
            return r.read()
    except Exception as e:
        die('%s failed: %s: %s' % (url, type(e).__name__, e))


def main():
    os.makedirs(DATA, exist_ok=True)
    man = {'api': API, 'sets': [], 'targets': {}}
    seqs = {}
    for bg, pheno, urn in SETS:
        meta = json.loads(get(API + urn))
        if 'detail' in meta:
            die('%s: %s' % (urn, meta['detail']))
        tg = (meta.get('targetGenes') or [{}])[0]
        tname = tg.get('name') or ''
        seq = ((tg.get('targetSequence') or {}).get('sequence') or '').upper()
        if not seq:
            die('%s carries no target sequence' % urn)
        if bg == 'Wuhan':
            if 'Omicron' in tname:
                die('%s is labelled %s but expected the Wuhan background' % (urn, tname))
        elif bg not in tname:
            die('%s is labelled %s but expected background %s' % (urn, tname, bg))
        want = {'binding': 'binding', 'expression': 'expression'}[pheno]
        if want not in (meta.get('shortDescription') or '').lower() + (meta.get('title') or '').lower():
            die('%s does not describe itself as a %s score set' % (urn, pheno))
        if bg in seqs and seqs[bg] != seq:
            die('two score sets for %s carry different target sequences' % bg)
        seqs[bg] = seq

        raw = get(API + urn + '/scores')
        path = os.path.join(DATA, 'scores_%s.csv' % urn.replace(':', '_'))
        with open(path, 'wb') as fh:
            fh.write(raw)
        head = raw.split(b'\n', 1)[0].decode().strip()
        nlines = raw.count(b'\n')
        man['sets'].append({
            'background': bg, 'phenotype': pheno, 'urn': urn,
            'target_name': tname, 'file': os.path.relpath(path, ROOT),
            'bytes': len(raw), 'data_rows': nlines - 1 if nlines else 0,
            'columns': head.split(','),
            'md5': hashlib.md5(raw).hexdigest(),
            'num_variants_declared': meta.get('numVariants'),
        })
        print('%-6s %-11s %-24s %7d rows  %s'
              % (bg, pheno, urn, nlines - 1, head))

    for bg, seq in sorted(seqs.items()):
        p = os.path.join(DATA, 'target_%s.txt' % bg.replace('.', ''))
        with open(p, 'w') as fh:
            fh.write(seq + '\n')
        man['targets'][bg] = {'file': os.path.relpath(p, ROOT), 'nt_length': len(seq),
                              'md5': hashlib.md5(seq.encode()).hexdigest()}
        if len(seq) % 3:
            die('%s target sequence length %d is not a multiple of three' % (bg, len(seq)))
        print('%-6s target %d nt = %d codons' % (bg, len(seq), len(seq) // 3))

    with open(os.path.join(DATA, 'fetch_manifest.json'), 'w') as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print('\nwritten to data/fetch_manifest.json')


if __name__ == '__main__':
    main()
