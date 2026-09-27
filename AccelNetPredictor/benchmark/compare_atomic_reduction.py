#!/usr/bin/env python3
"""Paired before/after timing of atomic removal, retaining all E/F/virial checks."""
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
    for name in ('before-serial', 'before-openmp', 'before-gpu', 'serial', 'openmp', 'gpu', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--cores', nargs='+', type=int, default=list(range(6, 14)))
    p.add_argument('--rounds', type=int, default=2)
    p.add_argument('--seconds', type=float, default=.05)
    a = p.parse_args()
    if a.rounds < 2 or len(set(a.cores)) < 8 or a.seconds <= 0:
        p.error('need >=2 rounds, eight distinct physical cores, and positive seconds')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    binaries = {name: getattr(a, name).resolve() for name in
                ('before_serial', 'before_openmp', 'before_gpu', 'serial', 'openmp', 'gpu')}
    for name in ('serial', 'before_serial'):
        if re.search(r'\b(?:GOMP_|__kmpc_|__nvomp|omp_get_|omp_set_)\w*',
                     subprocess.check_output(['nm', '-u', binaries[name]], text=True)):
            p.error('serial binaries must compile OpenMP out')
    report = dict(cores=a.cores, rounds=a.rounds, seconds=a.seconds,
                  device=os.environ.get('CUDA_VISIBLE_DEVICES'), candidate_only=True,
                  binaries={k: dict(path=str(v), sha256=hashlib.sha256(v.read_bytes()).hexdigest())
                            for k, v in binaries.items()}, rows=[])
    paths = [('before_off1', 'before_serial', 1, 'host'), ('off1', 'serial', 1, 'host'),
             ('before_on1', 'before_openmp', 1, 'host'), ('on1', 'openmp', 1, 'host'),
             ('on2', 'openmp', 2, 'host'), ('on4', 'openmp', 4, 'host'),
             ('before_on8', 'before_openmp', 8, 'host'), ('on8', 'openmp', 8, 'host'),
             ('before_gpu', 'before_gpu', 1, 'gpu'), ('gpu', 'gpu', 1, 'gpu')]
    os.sched_setaffinity(0, {a.cores[0]})
    until = time.monotonic()+20
    while time.monotonic() < until: sum(range(10000))
    for family, n, degree in [('chebyshev', 2048, 8), ('g5-series', 2048, 4), ('g5-high', 256, 16)]:
        for mode in (1, 2 if family == 'chebyshev' else 3):
            tag = f'{family}-N{n}-p{degree}-mode{mode}'
            records = {key: dict(samples=[], errors=[], profiles=[], runs=[]) for key, *_ in paths}
            for repeat in range(a.rounds):
                for key, binary, nt, backend in (paths if repeat % 2 == 0 else reversed(paths)):
                    cores = a.cores[:nt]
                    env = {k: v for k, v in os.environ.items() if not k.startswith(('OMP_', 'GOMP_'))}
                    controls = dict(OMP_NUM_THREADS=str(nt), OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close',
                                    OMP_PLACES=','.join('{'+str(c)+'}' for c in cores), OMP_MAX_ACTIVE_LEVELS='1',
                                    OMP_TARGET_OFFLOAD='MANDATORY' if backend == 'gpu' else 'DISABLED',
                                    OMP_WAIT_POLICY='PASSIVE', GOMP_SPINCOUNT='300000',
                                    OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
                    env.update(controls)
                    cmd = ['taskset', '-c', ','.join(map(str, cores)), str(binaries[binary]),
                           *map(str, [n, degree, a.seconds, mode, 1.7, backend, family, 'no-neighbors', 64, 'candidate-only'])]
                    print('RUN', tag, repeat, key, flush=True)
                    result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=900)
                    (out/f'{tag}-{repeat}-{key}.log').write_text(result.stdout+result.stderr)
                    if result.returncode: raise RuntimeError(result.stderr[-3000:])
                    samples = [float(line.split()[3]) for line in result.stdout.splitlines() if line.startswith('TIMING ')]
                    errors = [float(line.split()[4]) for line in result.stdout.splitlines() if line.startswith('CASE ')]
                    profiles = [[float(x) for x in line.split()[3:]] for line in result.stdout.splitlines() if line.startswith('PROFILE ')]
                    if len(samples) != 5 or len(errors) != 1: raise RuntimeError('incomplete benchmark output')
                    record = records[key]
                    record['samples'] += samples; record['errors'] += errors; record['profiles'] += profiles
                    record['runs'].append(dict(command=cmd, environment=controls))
            for record in records.values():
                record['seconds'] = statistics.median(record['samples'])
                record['round_medians'] = [statistics.median(record['samples'][i:i+5]) for i in range(0, 5*a.rounds, 5)]
                record['median_stages_seconds'] = dict(zip(
                    ('neighbors', 'prepare', 'upload', 'descriptors', 'network', 'forces', 'download', 'total'),
                    map(statistics.median, zip(*record['profiles']))))
            report['rows'].append(dict(family=family, atoms=n, order=degree, mode=mode, neighbors=64, records=records))
            (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
            print('RESULT', tag, {k: round(v['seconds']*1000, 4) for k, v in records.items()}, flush=True)
    report['passed'] = True
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
