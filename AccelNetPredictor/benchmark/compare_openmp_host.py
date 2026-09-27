#!/usr/bin/env python3
"""Measure identical prepared host kernels with OpenMP OFF and ON at 1/2/4/8 threads.

Use a serial target build and a GNU -fopenmp host target build. The ordinary
cpu-shared API deliberately strips OpenMP directives and is NOT the threaded
path. This benchmark uses backend=host for both builds, checks numerical
equivalence inside the driver, and saves GNU runtime affinity diagnostics.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('serial', 'openmp', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--cores', nargs='+', type=int, default=list(range(6, 14)))
    p.add_argument('--rounds', type=int, default=2)
    p.add_argument('--seconds', type=float, default=.05)
    a = p.parse_args()
    if len(set(a.cores)) < 8 or a.rounds < 2 or a.seconds <= 0:
        p.error('need eight distinct physical cores, >=2 rounds and positive seconds')
    a.serial = a.serial.resolve(); a.openmp = a.openmp.resolve()
    for exe, threaded in [(a.serial, False), (a.openmp, True)]:
        symbols = subprocess.check_output(['nm', '-u', exe], text=True)
        found = bool(re.search(r'\b(?:GOMP_|__kmpc_|__nvomp|omp_get_|omp_set_)\w*', symbols))
        if found != threaded:
            p.error(f'OpenMP symbol audit failed: {exe}')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    topology = subprocess.check_output(['lscpu', '-p=CPU,CORE,SOCKET,NODE'], text=True)
    physical = {int(x[0]): tuple(x[1:]) for line in topology.splitlines()
                if not line.startswith('#') for x in [line.split(',')]}
    selected = [physical[c] for c in a.cores[:8]]
    if len(set(selected)) != 8 or len({x[1:] for x in selected}) != 1:
        p.error('cores must be distinct physical cores on one socket/NUMA node')
    report = dict(rounds=a.rounds, seconds=a.seconds, cores=a.cores[:8],
                  topology=topology, core_warmup_seconds=20, backend='host', candidate_only=True,
                  binaries={key: dict(path=str(exe), sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
                            for key, exe in [('off', a.serial), ('on', a.openmp)]},
                  affinity_checks=[], rows=[])

    def run(exe, nt, args, tag, diagnostics=False):
        cores = a.cores[:nt]
        env = {k: v for k, v in os.environ.items() if not k.startswith(('OMP_', 'GOMP_'))}
        controls = dict(OMP_NUM_THREADS=str(nt), OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close',
                        OMP_PLACES=','.join('{'+str(c)+'}' for c in cores),
                        OMP_MAX_ACTIVE_LEVELS='1', OMP_TARGET_OFFLOAD='DISABLED',
                        OMP_WAIT_POLICY='PASSIVE', GOMP_SPINCOUNT='300000',
                        OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        env.update(controls)
        if diagnostics: env['OMP_DISPLAY_AFFINITY'] = 'TRUE'
        cmd = ['taskset', '-c', ','.join(map(str, cores)), str(exe), *map(str, args), 'candidate-only']
        result = subprocess.run(cmd, env=env, text=True, capture_output=True, timeout=900)
        (out/(tag+'.log')).write_text(result.stdout+result.stderr)
        if result.returncode: raise RuntimeError(result.stderr[-3000:])
        samples = [float(line.split()[3]) for line in result.stdout.splitlines() if line.startswith('TIMING ')]
        errors = [float(line.split()[4]) for line in result.stdout.splitlines() if line.startswith('CASE ')]
        profiles = [[float(x) for x in line.split()[3:]] for line in result.stdout.splitlines() if line.startswith('PROFILE ')]
        if len(samples) != 5 or len(errors) != 1: raise ValueError('missing timing/accuracy output')
        return dict(samples=samples, errors=errors, profiles=profiles, command=cmd, environment=controls), result.stderr

    # Confirm worker affinities from the actual numerical executable before timing.
    for nt in (1, 2, 4, 8):
        record, stderr = run(a.openmp, nt, [16, 11, .001, 3, 1.7, 'host', 'g5-high', 'no-neighbors', 16],
                             f'affinity-{nt}', diagnostics=True)
        bindings = sorted(set(map(int, re.findall(r'affinity (\d+)\s*$', stderr, re.MULTILINE))))
        # libgomp does not print an affinity line for a serialized one-thread
        # region; taskset and the explicit one-place binding still pin it.
        if bindings != sorted(a.cores[:nt]) and not (nt == 1 and not bindings):
            raise RuntimeError(f'Expected actual worker bindings {a.cores[:nt]}, got {bindings}')
        report['affinity_checks'].append(dict(threads=nt, bindings=bindings, **record))
    os.sched_setaffinity(0, {a.cores[0]})
    until = time.monotonic()+20
    while time.monotonic() < until: sum(range(10000))
    # Explicit taskset in run() widens child affinity again for multiple threads.
    cases = [('chebyshev', 2048, 8, 64), ('g5-series', 2048, 4, 64),
             ('g5-high', 256, 16, 64), ('g5-high', 2048, 16, 64)]
    paths = [('off1', a.serial, 1)] + [(f'on{nt}', a.openmp, nt) for nt in (1, 2, 4, 8)]
    for family, n, degree, neighbors in cases:
        for mode in (1, 3):
            # Chebyshev supports forced moment as mode 2 (G5 mode 3).
            actual_mode = 2 if family == 'chebyshev' and mode == 3 else mode
            tag = f'{family}-N{n}-p{degree}-K{neighbors}-'+('direct' if mode == 1 else 'moment')
            records = {key: dict(samples=[], errors=[], profiles=[], runs=[]) for key, *_ in paths}
            for repeat in range(a.rounds):
                for key, exe, nt in (paths if repeat % 2 == 0 else reversed(paths)):
                    print('RUN', tag, repeat, key, flush=True)
                    record, _ = run(exe, nt, [n, degree, a.seconds, actual_mode, 1.7, 'host', family, 'no-neighbors', neighbors],
                                    f'{tag}-{repeat}-{key}')
                    records[key]['samples'] += record.pop('samples')
                    records[key]['errors'] += record.pop('errors')
                    records[key]['profiles'] += record.pop('profiles')
                    records[key]['runs'].append(record)
            for record in records.values():
                record['seconds'] = statistics.median(record['samples'])
                record['median_stages_seconds'] = dict(zip(
                    ('neighbors', 'prepare', 'upload', 'descriptors', 'network', 'forces', 'download', 'total'),
                    map(statistics.median, zip(*record['profiles']))))
            baseline = records['off1']['seconds']
            for record in records.values(): record['speedup_over_off1'] = baseline/record['seconds']
            report['rows'].append(dict(family=family, atoms=n, order=degree, neighbors=neighbors,
                                       mode=actual_mode, records=records))
            (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
            print('RESULT', tag, {k: round(v['seconds']*1000, 4) for k, v in records.items()}, flush=True)
    report['passed'] = True
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
