#!/usr/bin/env python3
"""CPU-only n2p2/AccelNet comparison on official models and synthetic SF fixtures.

Both executables must be built without OpenMP. Public CPU batches and prepared
CPU models use identical numerical kernels. n2p2 returns energy/forces; AccelNet
also computes virial internally. Loading, the independent initial result check,
and warmup are excluded. Five samples per process; reverse process order on
alternate rounds. Every timed process must return matching energy and forces.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'test'))
from check_n2p2_extended import compare
from check_n2p2_official import MODELS
from n2p2_extended_fixtures import make_model, geometry, write_structure, TYPES


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(command, log):
    p = subprocess.run(list(map(str, command)), capture_output=True, text=True, timeout=900)
    log.write_text(p.stdout + p.stderr)
    if p.returncode:
        raise RuntimeError(f'{command}: see {log}')
    result = {}
    for line in p.stdout.splitlines():
        fields = line.split()
        if fields and fields[0] in ('ENERGY', 'FORCE', 'VIRIAL', 'TIMING'):
            result.setdefault(fields[0], []).append([float(x) for x in fields[1:]])
    if len(result.get('TIMING', [])) != 5:
        raise RuntimeError(f'expected five timing samples: {log}')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('reference', 'candidate', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--examples', type=Path, help='n2p2/examples/nnp-predict')
    p.add_argument('--suite', choices=('official', 'synthetic', 'both'), default='official')
    p.add_argument('--scopes', nargs='+', choices=('fixed', 'full'), default=['fixed', 'full'])
    p.add_argument('--models', nargs='+', default=MODELS)
    p.add_argument('--atoms', type=int, default=512)
    p.add_argument('--types', type=int, nargs='+', choices=(2, 3, 9, *TYPES), default=(2, 3, 9, *TYPES))
    p.add_argument('--seconds', type=float, default=.2)
    p.add_argument('--rounds', type=int, default=2)
    p.add_argument('--affinity', type=int, default=6)
    p.add_argument('--direct', action='store_true', help='Force G5 direct evaluation as well; compact/weighted angular SFs are always direct')
    p.add_argument('--baseline', type=Path, help='Prior report; reject CPU regressions beyond --max-slowdown')
    p.add_argument('--max-slowdown', type=float, default=1.10)
    p.add_argument('--max-n2p2-slowdown', type=float,
                   help='Reject a public CPU/n2p2 time ratio above this limit (use fixed scope to isolate inference)')
    args = p.parse_args()
    if args.seconds <= 0 or args.rounds < 2 or args.atoms < 8:
        p.error('positive seconds, at least two rounds, and at least eight atoms required')
    if args.suite != 'synthetic' and args.examples is None:
        p.error('--examples required for official models')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    binaries = {}
    for name in ('reference', 'candidate'):
        exe = getattr(args, name).resolve(); setattr(args, name, exe)
        symbols = subprocess.check_output(['nm', '-u', str(exe)], text=True)
        omp = re.findall(r'\b(?:GOMP_|__kmpc_|__nvomp|__pgi_omp|omp_get_|omp_set_)\w*', symbols)
        if omp:
            raise ValueError(f'{name} imports OpenMP runtime: {omp[:5]}')
        binaries[name] = dict(path=str(exe), sha256=sha(exe), openmp_runtime_symbols=omp)
    os.sched_setaffinity(0, {args.affinity})
    # Warm the pinned core out of its idle frequency state before comparisons.
    warmup_end = time.monotonic() + 20.0
    while time.monotonic() < warmup_end:
        sum(range(10000))
    os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    cases = []
    if args.suite != 'synthetic':
        for name in args.models:
            directory = args.examples.resolve() / name
            # The shipped prediction examples each contain one configuration.
            text = directory.joinpath('input.data').read_text()
            if sum(line.strip() == 'begin' for line in text.splitlines()) != 1:
                raise ValueError(f'{name}: select a single-frame input.data')
            rows = [line.split() for line in text.splitlines()]
            atoms = [r for r in rows if r and r[0] == 'atom']
            lattice = [r[1:] for r in rows if r and r[0] == 'lattice']
            stem = out / name
            stem.with_suffix('.data').write_text(text)
            xsf = ['CRYSTAL', 'PRIMVEC', *[' '.join(r) for r in lattice]] if lattice else []
            xsf += ['PRIMCOORD', f'{len(atoms)} 1', *[' '.join([r[4], *r[1:4]]) for r in atoms]]
            stem.with_suffix('.xsf').write_text('\n'.join(xsf) + '\n')
            sf = Counter(tuple(line.split()[1:3]) for line in directory.joinpath('input.nn').read_text().splitlines()
                         if line.strip().startswith('symfunction_short '))
            hashes = {f.name: sha(f) for f in directory.iterdir()
                      if f.name in ('input.nn', 'input.data', 'scaling.data') or f.name.startswith('weights.')}
            cases.append(dict(name=name, model=str(directory), stem=str(stem), atoms=len(atoms),
                              elements=sorted({r[4] for r in atoms}), scopes=args.scopes,
                              descriptors={':'.join(k): v for k, v in sorted(sf.items())}, model_sha256=hashes))
    if args.suite != 'official':
        positions, cell = geometry(args.atoms)
        stem = out / 'synthetic-structure'; write_structure(stem, positions, cell)
        for kind in args.types:
            directory = out / f'type{kind}'; make_model(directory, (kind,), count=3 if kind == 9 else 6)
            cases.append(dict(name=f'type{kind}', model=str(directory), stem=str(stem), atoms=args.atoms,
                              elements=['H', 'O'], scopes=['fixed']))
    report = dict(affinity=args.affinity, core_warmup_seconds=20, g5_mode='direct' if args.direct else 'auto',
                  cpu_info=Path('/proc/cpuinfo').read_text().split('model name')[1].split('\n')[0],
                  seconds_per_sample=args.seconds, rounds=args.rounds, binaries=binaries, cases=cases, rows=[],
                  extra_accelnet_work='Virial is computed by the common kernel even when only energy/forces are requested.')
    for case in cases:
        for scope in case['scopes']:
            stem = Path(case['stem']); directory = Path(case['model'])
            commands = {
                'n2p2': [args.reference, directory, stem.with_suffix('.data'), args.seconds, scope],
                'accelnet_cpu': [args.candidate, directory, stem.with_suffix('.xsf'), args.seconds, 'cpu', scope, 'n2p2'],
                'accelnet_prepared': [args.candidate, directory, stem.with_suffix('.xsf'), args.seconds, 'host', scope, 'n2p2'],
            }
            if args.direct:
                for key in ('accelnet_cpu', 'accelnet_prepared'):
                    commands[key].append(1)
            records = {key: [] for key in commands}; errors = []
            for repeat in range(args.rounds):
                keys = list(commands) if repeat % 2 == 0 else list(reversed(commands))
                result = {}
                for key in keys:
                    log = out / f"{case['name']}-{scope}-{repeat}-{key}.log"
                    print('RUN', case['name'], scope, repeat, key, flush=True)
                    result[key] = execute(commands[key], log)
                    records[key].append(result[key]['TIMING'])
                errors.append({key: compare(result[key], result['n2p2']) for key in ('accelnet_cpu', 'accelnet_prepared')})
            timings = {key: dict(median=statistics.median(row[0] for group in groups for row in group), samples=groups)
                       for key, groups in records.items()}
            speedup = {key: timings['n2p2']['median'] / timings[key]['median'] for key in ('accelnet_cpu', 'accelnet_prepared')}
            row = dict(name=case['name'], atoms=case['atoms'], scope=scope, timings=timings, speedup=speedup,
                       errors=errors, commands={key: list(map(str, cmd)) for key, cmd in commands.items()})
            report['rows'].append(row)
            (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
            print('RESULT', case['name'], scope, {key: round(v['median']*1000, 4) for key, v in timings.items()}, speedup, flush=True)
    if args.baseline:
        baseline = json.loads(args.baseline.read_text())
        if not baseline.get('passed'):
            raise ValueError('baseline did not pass its numerical checks')
        old = {(r['name'], r['atoms'], r['scope']): r for r in baseline['rows']}
        regressions = []
        for row in report['rows']:
            if row['name'] == 'type9' and baseline.get('g5_mode', 'auto') != report['g5_mode']:
                row['baseline_note'] = 'Different G5 mode; compare n2p2 directly, not the old auto timing.'
                continue
            previous = old[(row['name'], row['atoms'], row['scope'])]
            row['baseline_ratio'] = {}
            for key in ('accelnet_cpu', 'accelnet_prepared'):
                ratio = row['timings'][key]['median'] / previous['timings'][key]['median']
                row['baseline_ratio'][key] = ratio
                if ratio > args.max_slowdown:
                    regressions.append((row['name'], row['scope'], key, ratio))
        report['baseline'] = str(args.baseline.resolve())
        report['regressions'] = regressions
        if regressions:
            (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
            raise AssertionError(f'CPU regressions: {regressions}')
    if args.max_n2p2_slowdown is not None:
        report['max_n2p2_slowdown'] = args.max_n2p2_slowdown
        report['n2p2_regressions'] = [
            (row['name'], row['scope'], 1 / row['speedup']['accelnet_cpu'])
            for row in report['rows']
            if 1 / row['speedup']['accelnet_cpu'] > args.max_n2p2_slowdown
        ]
        if report['n2p2_regressions']:
            (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
            raise AssertionError(f"n2p2 parity check: {report['n2p2_regressions']}")
    report['passed'] = True
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
