#!/usr/bin/env python3
"""Compare shared direct/moment CPU and GPU inference for exact high-order G5.

CPU must be compiled without OpenMP. CSR neighbors are fixed, all inside Rc;
GPU wall time includes transfers, NN and force/virial assembly. The executable
checks every E/F/virial component against the retained CPU reference before and
during timing. No assertion assumes that high-order moments are always faster.
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
    for name in ('cpu', 'gpu', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--affinity', type=int, default=6)
    p.add_argument('--rounds', type=int, default=2)
    p.add_argument('--seconds', type=float, default=.05)
    p.add_argument('--atoms', type=int, nargs='+', default=[256, 2048])
    p.add_argument('--orders', type=int, nargs='+', default=list(range(11, 17)))
    p.add_argument('--series-atoms', type=int, nargs='+', default=[256])
    a = p.parse_args()
    if a.rounds < 2 or a.seconds <= 0 or min(a.atoms+a.series_atoms) < 2 or min(a.orders) < 11 or max(a.orders) > 16:
        p.error('need >=2 rounds, positive duration, >=2 atoms, and orders 11..16')
    a.cpu = a.cpu.resolve(); a.gpu = a.gpu.resolve()
    if re.search(r'\b(?:GOMP_|__kmpc_|__nvomp|omp_get_|omp_set_)\w*', subprocess.check_output(['nm', '-u', a.cpu], text=True)):
        p.error('CPU executable must have OpenMP compiled out')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, {a.affinity})
    os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_TARGET_OFFLOAD='MANDATORY')
    until = time.monotonic()+20
    while time.monotonic() < until: sum(range(10000))
    report = dict(affinity=a.affinity, rounds=a.rounds, seconds=a.seconds, core_warmup_seconds=20,
                  device=os.environ.get('CUDA_VISIBLE_DEVICES'), binaries={key: dict(path=str(exe), sha256=hashlib.sha256(exe.read_bytes()).hexdigest()) for key, exe in [('cpu', a.cpu), ('gpu', a.gpu)]}, rows=[])
    cases = [('g5-high', n, degree, 64) for n in a.atoms for degree in a.orders]
    cases += [('g5-high-series', n, max(a.orders), 256) for n in a.series_atoms]
    for family, n, degree, neighbors in cases:
        tag = f'{family}-N{n}-p{degree}-K{neighbors}'
        paths = [('cpu_direct', a.cpu, 'cpu-shared', 1), ('cpu_moment', a.cpu, 'cpu-shared', 3),
                 ('gpu_direct', a.gpu, 'gpu', 1), ('gpu_moment', a.gpu, 'gpu', 3)]
        records = {key: dict(samples=[], errors=[], commands=[]) for key, *_ in paths}
        for repeat in range(a.rounds):
            for key, exe, backend, mode in (paths if repeat % 2 == 0 else reversed(paths)):
                cmd = list(map(str, [exe, n, degree, a.seconds, mode, 1.7, backend, family, 'no-neighbors', neighbors]))
                print('RUN', tag, repeat, key, flush=True)
                result = subprocess.run(cmd, text=True, capture_output=True, timeout=900)
                (out/f'{tag}-{repeat}-{key}.log').write_text(result.stdout+result.stderr)
                if result.returncode: raise RuntimeError(result.stderr[-2000:])
                samples = [float(line.split()[3]) for line in result.stdout.splitlines() if line.startswith('TIMING ')]
                errors = [float(line.split()[4]) for line in result.stdout.splitlines() if line.startswith('CASE ')]
                if len(samples) != 5 or len(errors) != 1: raise ValueError('missing timing/accuracy output')
                records[key]['samples'] += samples; records[key]['errors'] += errors; records[key]['commands'].append(cmd)
        for record in records.values(): record['seconds'] = statistics.median(record['samples'])
        t = {key: r['seconds'] for key, r in records.items()}
        row = dict(family=family, atoms=n, order=degree, neighbors=neighbors, descriptors=6 if family=='g5-high' else 6*degree,
                   records=records, cpu_direct_over_moment=t['cpu_direct']/t['cpu_moment'], gpu_direct_over_moment=t['gpu_direct']/t['gpu_moment'],
                   cpu_direct_over_gpu_moment=t['cpu_direct']/t['gpu_moment'])
        report['rows'].append(row)
        (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        print('RESULT', tag, {k: round(v*1000, 4) for k, v in t.items()}, flush=True)
    report['passed'] = True
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
