#!/usr/bin/env python3
import asyncio
import collections
import functools
import http.server
import json
import os
import shutil
import socketserver
import sys
import tempfile
import threading

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE
DATA = os.path.join(ROOT, 'data')
STRUCT = os.path.join(ROOT, 'structures')
OUT = os.path.join(ROOT, 'figures')
REND = os.path.join(HERE, 'structure_render')
PORT = 8734
NEAR = 18.0
MIN_RUN = 4


def die(m):
    sys.exit('FATAL 20_render_structure.py: %s' % m)


def payload():
    jf = os.path.join(DATA, 'interface_distance.json')
    pdb = os.path.join(STRUCT, '7WPB.pdb')
    for f in (jf, pdb, os.path.join(REND, 'page.html'),
              os.path.join(REND, '3Dmol-min.js')):
        if not os.path.isfile(f):
            die('%s is absent' % f)
    IF = json.load(open(jf))
    off, rbd, partner = IF['offset'], IF['chain'], IF['partner_chain']
    BIND = set(IF['by_phenotype']['binding']['departing'])
    EXPR = set(IF['by_phenotype']['expression']['departing'])
    RNK = set(IF['departing_rank'])
    allpos = sorted(int(k) for k in IF['distance_angstrom'])

    ca = collections.defaultdict(dict)
    for line in open(pdb):
        if not line.startswith('ATOM') or line[16] not in (' ', 'A'):
            continue
        if line[12:16].strip() != 'CA':
            continue
        ca[line[21]][int(line[22:26])] = (float(line[30:38]), float(line[38:46]),
                                          float(line[46:54]))
    if rbd not in ca or partner not in ca:
        die('the chains named in interface_distance.json are not in the structure')

    tr = np.array([ca[rbd][k] for k in sorted(ca[rbd])])
    pk = sorted(ca[partner])
    pall = np.array([ca[partner][k] for k in pk])
    dmin = np.sqrt(((pall[:, None, :] - tr[None, :, :]) ** 2).sum(-1)).min(1)
    hits = [int(k) for k, d in zip(pk, dmin) if d < NEAR]
    runs, cur = [], []
    for k in hits:
        if cur and k == cur[-1] + 1:
            cur.append(k)
        else:
            if len(cur) >= MIN_RUN:
                runs.extend(cur)
            cur = [k]
    if len(cur) >= MIN_RUN:
        runs.extend(cur)
    if not runs:
        die('no contiguous stretch of the receptor lies within %.0f A' % NEAR)

    groups = {k: [] for k in ('binding', 'expression', 'both', 'rank', 'none')}
    for p in allpos:
        if p + off not in ca[rbd]:
            die('position %d is absent from chain %s' % (p, rbd))
        r = p + off
        if p in BIND and p in EXPR:
            groups['both'].append(r)
        elif p in BIND:
            groups['binding'].append(r)
        elif p in EXPR:
            groups['expression'].append(r)
        elif p in RNK:
            groups['rank'].append(r)
        else:
            groups['none'].append(r)
    if sum(len(v) for v in groups.values()) != len(allpos):
        die('the colour groups do not partition the tested positions')

    return {'rbd': rbd, 'partner': partner, 'offset': off, 'near': runs,
            'groups': groups,
            'labels': {str(p + off): str(p) for p in sorted(BIND | EXPR | RNK)}}, pdb


async def shoot(tmp, out):
    from playwright.async_api import async_playwright
    H = functools.partial(http.server.SimpleHTTPRequestHandler, directory=tmp)
    socketserver.TCPServer.allow_reuse_address = True
    srv = socketserver.TCPServer(('127.0.0.1', PORT), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        async with async_playwright() as pw:
            b = await pw.chromium.launch(args=['--use-gl=swiftshader',
                                               '--enable-unsafe-swiftshader'])
            pg = await b.new_page(viewport={'width': 1000, 'height': 850},
                                  device_scale_factor=3)
            await pg.goto('http://127.0.0.1:%d/page.html' % PORT)
            await pg.wait_for_function('window.ready === true', timeout=90000)
            await pg.wait_for_timeout(900)
            raw = os.path.join(tmp, 'raw.png')
            await (await pg.query_selector('#v')).screenshot(path=raw)
            await b.close()
    finally:
        srv.shutdown()

    im = Image.open(raw).convert('RGB')
    bb = Image.eval(im.convert('L'), lambda p: 255 - p).getbbox()
    if not bb:
        die('the render came back blank')
    pad = 14
    im = im.crop((max(0, bb[0] - pad), max(0, bb[1] - pad),
                  min(im.width, bb[2] + pad), min(im.height, bb[3] + pad)))
    if im.width < 1200 or im.height < 800:
        die('the cropped render is %dx%d, too small to be the complex'
            % (im.width, im.height))
    os.makedirs(OUT, exist_ok=True)
    im.save(out)
    return im.size


def main():
    rec, pdb = payload()
    with tempfile.TemporaryDirectory() as tmp:
        for f in ('page.html', '3Dmol-min.js'):
            shutil.copy(os.path.join(REND, f), tmp)
        shutil.copy(pdb, os.path.join(tmp, 'structure.pdb'))
        json.dump(rec, open(os.path.join(tmp, 'payload.json'), 'w'))
        out = os.path.join(OUT, 'structure_cartoon.png')
        size = asyncio.run(shoot(tmp, out))
    print('chain %s against chain %s, offset %+d' % (rec['rbd'], rec['partner'],
                                                     rec['offset']))
    print('receptor residues drawn: %d' % len(rec['near']))
    for k in ('binding', 'expression', 'both', 'rank', 'none'):
        print('  %-11s %d' % (k, len(rec['groups'][k])))
    print('wrote figures/structure_cartoon.png at %dx%d' % size)


if __name__ == '__main__':
    main()
