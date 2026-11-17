#!/usr/bin/env python3
import collections
import hashlib
import csv
import importlib.util
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(
    os.path.dirname(HERE) if os.path.basename(HERE) == 'codes' else HERE, 'data')
_sp = importlib.util.spec_from_file_location(
    'rk', os.path.join(HERE, '12_census_test_rank.py'))
rk = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rk)

DRAWS = 10000
SEED = 20261111
MIN_TARGETS = 6
MIN_CELLS = 8
CAL_DRAWS = 200
ENGINEERED = re.compile(
    r'\((?:[A-Z]\d+[A-Z](?:/[A-Z]\d+[A-Z])*|pross\d+|[a-z]+\d+)\)')
LOG = []


def say(s):
    print(s, flush=True)
    LOG.append(s)


def base_name(n):
    s = re.sub(r'\s*\([^)]*\)\s*', ' ', n or '').strip()
    s = re.sub(r'\s+', ' ', s)
    return s or 'unnamed'


def bh(ps):
    m = len(ps)
    o = sorted(range(m), key=lambda i: ps[i])
    out = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(o), 1):
        k = m - rank + 1
        prev = min(prev, ps[i] * m / k)
        out[i] = prev
    return out


def main():
    seq = rk.load_sequences()
    say('deposits with a translatable target sequence: %d' % len(seq))
    rows = [r for r in csv.DictReader(open(os.path.join(DATA, 'route_census_mavedb.csv')))
            if int(r['complete_routes']) >= MIN_TARGETS]
    by_pair = collections.defaultdict(list)
    for r in rows:
        by_pair[(r['from_urn'], r['via_urn'])].append(r)
    say('deposit pairs in the census: %d' % len(by_pair))

    sc = rk.Scores(limit=200)
    rng = np.random.default_rng(SEED)
    obs = collections.defaultdict(list)
    pools = collections.defaultdict(dict)
    rawnames = {}
    aside = collections.Counter()
    done = 0
    for (ua, ub), items in sorted(by_pair.items()):
        done += 1
        if done % 1000 == 0:
            say('  %d of %d pairs' % (done, len(by_pair)))
        pa, pb = seq.get(ua), seq.get(ub)
        if pa is None or pb is None or len(pa) != len(pb):
            aside['no comparable target sequence'] += len(items)
            continue
        sa, sb = sc.get(ua), sc.get(ub)
        if not sa or not sb:
            aside['deposit unreadable'] += len(items)
            continue
        pool = rk.null_vectors(sa, sb, pa, pb, MIN_TARGETS)
        if len(pool) < 8:
            aside['too few shared-parent positions'] += len(items)
            continue
        sysname = base_name(items[0]['from_name'])
        rawnames.setdefault(sysname, set()).update(
            [items[0]['from_name'], items[0]['via_name']])
        got = False
        for r in items:
            pos, x, z = int(r['pos']), r['res_from'], r['res_via']
            if pos > len(pa) or pa[pos - 1] != x or pb[pos - 1] != z:
                aside['census row disagrees with the target sequence'] += 1
                continue
            va, vb, tgt, first = rk.cell_vectors(sa, sb, pos, x, z)
            if va.size < MIN_TARGETS:
                aside['too few shared targets'] += 1
                continue
            rho = float(rk.spearman(va[None, :], vb[None, :])[0])
            if not np.isfinite(rho):
                aside['rank correlation undefined'] += 1
                continue
            obs[sysname].append({'pair': (ua, ub), 'n': int(va.size), 'rho': rho,
                                 'pos': pos, 'via_is_alanine': r['via_is_alanine']})
            got = True
        if got:
            pools[sysname][(ua, ub)] = pool

    sig = {u: hashlib.sha1(p.encode()).hexdigest()[:16] for u, p in seq.items()}
    kept = collections.Counter()
    for sysname in list(obs):
        best = {}
        for c in sorted(obs[sysname], key=lambda c: (c['pair'], c['pos'])):
            a, b = sig[c['pair'][0]], sig[c['pair'][1]]
            key = (min(a, b), max(a, b), c['pos'])
            cur = best.get(key)
            if cur is None or c['n'] > cur['n']:
                best[key] = c
        kept['before'] += len(obs[sysname])
        obs[sysname] = sorted(best.values(), key=lambda c: (c['pair'], c['pos']))
        kept['after'] += len(obs[sysname])
        live = {c['pair'] for c in obs[sysname]}
        pools[sysname] = {k: v for k, v in pools[sysname].items() if k in live}
    say('cells before collapsing redundant realisations: %d' % kept['before'])
    say('cells after:                                   %d' % kept['after'])

    say('systems with at least one cell: %d' % len(obs))
    for k, v in sorted(aside.items()):
        say('  set aside, %-50s %d' % (k, v))
    say('')

    results, calib = [], []
    for sysname in sorted(obs):
        cells = obs[sysname]
        if len(cells) < MIN_CELLS:
            continue
        pool = pools[sysname]
        k = len(cells)
        draws = np.empty(DRAWS)
        usable = True
        per = []
        for c in cells:
            p = [(a, b) for a, b in pool.get(c['pair'], []) if a.size >= c['n']]
            if len(p) < 4:
                usable = False
                break
            per.append((p, c['n']))
        if not usable:
            continue
        NDRAW = DRAWS + CAL_DRAWS
        acc = np.zeros((NDRAW, k))
        for j, (p, n) in enumerate(per):
            bylen = collections.defaultdict(list)
            for a, b in p:
                bylen[a.size].append((a, b))
            lens = sorted(bylen)
            sizes = np.array([len(bylen[L]) for L in lens], dtype=float)
            which = rng.choice(len(lens), NDRAW, p=sizes / sizes.sum())
            A = np.empty((NDRAW, n))
            B = np.empty((NDRAW, n))
            for li, L in enumerate(lens):
                msk = which == li
                c = int(msk.sum())
                if not c:
                    continue
                blk = bylen[L]
                As = np.stack([x[0] for x in blk])
                Bs = np.stack([x[1] for x in blk])
                pick = rng.integers(0, len(blk), c)
                cols = rng.random((c, L)).argsort(axis=1)[:, :n]
                A[msk] = np.take_along_axis(As[pick], cols, axis=1)
                B[msk] = np.take_along_axis(Bs[pick], cols, axis=1)
            acc[:, j] = rk.spearman(A, B)
        allmeans = np.nanmean(acc, axis=1)
        draws = np.sort(allmeans[:DRAWS])
        mean_obs = float(np.mean([c['rho'] for c in cells]))
        p = (float(np.searchsorted(draws, mean_obs, side='right')) + 1.0) / (DRAWS + 1.0)
        for q in allmeans[DRAWS:]:
            calib.append((float(np.searchsorted(draws, float(q), side='right')) + 1.0)
                         / (DRAWS + 1.0))
        kind = ('designed background'
                if any(ENGINEERED.search(n) for n in rawnames.get(sysname, ()))
                else 'natural variant')
        results.append({'system': sysname, 'kind': kind, 'cells': k,
                        'cells_via_alanine': sum(1 for c in cells if c['via_is_alanine'] == '1'),
                        'deposit_pairs': len(pool),
                        'mean_rho': round(mean_obs, 5),
                        'null_mean_rho': round(float(np.mean(draws)), 5),
                        'null_p05': round(float(np.percentile(draws, 5)), 5),
                        'shift': round(mean_obs - float(np.mean(draws)), 5),
                        'p_perm': round(p, 6)})
        with open(os.path.join(DATA, 'combined2_partial.csv'), 'a') as fh:
            fh.write('%s,%d,%.5f,%.5f,%.6f\n'
                     % (sysname.replace(',', ';'), k, mean_obs,
                        float(np.mean(draws)), p))

    if not results:
        sys.exit('FATAL 16_combined_rank_test.py: no system carried enough cells')
    for r, a in zip(results, bh([r['p_perm'] for r in results])):
        r['p_bh'] = round(a, 6)
    results.sort(key=lambda r: r['p_bh'])

    rate = sum(1 for p in calib if p < 0.05) / float(len(calib))
    say('systems tested: %d' % len(results))
    say('calibration: a mean drawn from the null rejects at %.4f against a '
        'nominal 0.05, over %d draws' % (rate, len(calib)))
    sig = [r for r in results if r['p_bh'] < 0.05]
    say('systems departing at BH 0.05: %d of %d' % (len(sig), len(results)))
    below = sum(1 for r in results if r['shift'] < 0)
    say('systems whose mean sits below its null mean: %d of %d' % (below, len(results)))
    say('')
    import numpy as _np
    sh = _np.array([r['shift'] for r in results])
    cl = _np.array([r['cells'] for r in results], dtype=float)
    pb = _np.array([r['p_bh'] for r in results])
    say('')
    say('BY SYSTEM KIND, derived from the deposit names')
    say('  %-22s %6s %10s %14s %12s' % ('kind', 'systems', 'departing',
                                        'median shift', 'median cells'))
    bykind = {}
    for kd in sorted({r['kind'] for r in results}):
        m = _np.array([r['kind'] == kd for r in results])
        bykind[kd] = {
            'systems': int(m.sum()),
            'departing': int((pb[m] < 0.05).sum()),
            'median_shift': round(float(_np.median(sh[m])), 5),
            'median_cells': float(_np.median(cl[m])),
            'shift_vs_log_cells': (round(float(_np.corrcoef(sh[m], _np.log10(cl[m]))[0, 1]), 3)
                                   if m.sum() > 2 else None),
            'systems_below_null': int((sh[m] < 0).sum())}
        say('  %-22s %6d %10d %14.4f %12.0f'
            % (kd, bykind[kd]['systems'], bykind[kd]['departing'],
               bykind[kd]['median_shift'], bykind[kd]['median_cells']))
    say('')
    say('IS DEPARTURE STILL TIED TO SIZE INSIDE A KIND?')
    say('  %-22s %s' % ('all systems',
                        'shift against log cells %+.3f'
                        % float(_np.corrcoef(sh, _np.log10(cl))[0, 1])))
    for kd in sorted(bykind):
        v = bykind[kd]['shift_vs_log_cells']
        say('  %-22s %s' % (kd, 'shift against log cells %+.3f' % v
                            if v is not None else 'too few systems'))
    say('')
    say('%-34s %6s %6s %10s %10s %9s %9s'
        % ('system', 'cells', 'alan', 'mean rho', 'null mean', 'shift', 'p (BH)'))
    for r in results[:30]:
        say('%-34s %6d %6d %10.4f %10.4f %+9.4f %9.5f'
            % (r['system'][:34], r['cells'], r['cells_via_alanine'], r['mean_rho'],
               r['null_mean_rho'], r['shift'], r['p_bh']))

    with open(os.path.join(DATA, 'combined2_bykind.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['system', 'kind', 'cells', 'cells_via_alanine',
                                           'deposit_pairs', 'mean_rho',
                                           'null_mean_rho', 'null_p05', 'shift',
                                           'p_perm', 'p_bh'])
        w.writeheader()
        w.writerows(results)
    summ = {'draws': DRAWS, 'seed': SEED, 'min_cells_per_system': MIN_CELLS,
            'cells_before_collapse': int(kept['before']),
            'cells_tested': int(kept['after']),
            'min_targets': MIN_TARGETS,
            'statistic': 'mean rank correlation over a system test cells',
            'null': 'same count of shared-parent positions from the same deposit '
                    'pairs, matched on target count and cell count',
            'calibration_rejection_at_nominal_05': round(rate, 4),
            'calibration_draws': len(calib),
            'systems_tested': len(results), 'systems_departing': len(sig),
            'systems_below_their_null': below,
            'calibration_draws_per_system': CAL_DRAWS,
            'calibration_is_disjoint_from_the_reference': True,
            'by_kind': bykind,
            'kind_rule': 'a deposit name carrying a parenthesised point '
                         'substitution or design tag marks a designed background',
            'median_shift': round(float(np.median([r['shift'] for r in results])), 5),
            'set_aside': dict(aside)}
    with open(os.path.join(DATA, 'combined2_bykind_summary.json'), 'w') as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    with open(os.path.join(DATA, 'combined2_bykind_log.txt'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    say('')
    say('written to combined2_bykind.csv and combined2_bykind_summary.json')


if __name__ == '__main__':
    main()
