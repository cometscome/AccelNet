#!/usr/bin/env python3
"""Measure the SAME OpenMP kernels on host/device versus the established CPU path.

Each child validates all E/F/virial components, warms persistent buffers, and
alternates five samples. Ratios > 1 mean the target path is slower. The two
paths use libraries from one build (same Fortran compiler and optimization).
"""
import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--benchmark', required=True, type=Path)
p.add_argument('--output', required=True, type=Path)
p.add_argument('--sizes', nargs='+', type=int, default=[64, 512, 4096])
p.add_argument('--backends', nargs='+', choices=['host', 'gpu', 'cpu-shared'], default=['host', 'gpu'])
p.add_argument('--families', nargs='+', default=['chebyshev', 'lj', 'g4', 'g5', 'behler'])
p.add_argument('--order', type=int, default=5)
p.add_argument('--seconds', type=float, default=0.08)
p.add_argument('--spacing', type=float, default=1.7)
p.add_argument('--cpu', type=int, default=6)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
a.benchmark = a.benchmark.resolve()
os.sched_setaffinity(0, {a.cpu})
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_TARGET_OFFLOAD='MANDATORY')
report = []
for n in a.sizes:
    for family in a.families:
        for mode in ([1, 2] if family == 'chebyshev' else [1]):
            for backend in a.backends:
                name = f'{family}-mode{mode}-{n}-{backend}'
                result = subprocess.run([str(a.benchmark), str(n), str(a.order), str(a.seconds),
                                         str(mode), str(a.spacing), backend, family, 'no-neighbors'],
                                        text=True, capture_output=True, env=env, timeout=600)
                (a.output / f'{name}.log').write_text(result.stdout + result.stderr)
                if result.returncode:
                    raise RuntimeError(f'{name} failed; see log')
                samples = [list(map(float, line.split()[2:])) for line in result.stdout.splitlines()
                           if line.startswith('TIMING ')]
                if len(samples) != 5:
                    raise RuntimeError(f'{name}: missing timing samples')
                medians = [statistics.median(s[i] for s in samples) for i in range(4)]
                case = next(line.split() for line in result.stdout.splitlines() if line.startswith('CASE '))
                entry = dict(natoms=n, family=family, mode=mode, backend=backend, order=a.order,
                             spacing=a.spacing, neighbors=int(case[3]), max_error=float(case[4]),
                             cpu_seconds=medians[0], target_seconds=medians[1],
                             neighbor_construction_timed=False,
                             target_over_cpu=statistics.median(s[1]/s[0] for s in samples),
                             samples=samples)
                report.append(entry)
                (a.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
                print(name, 'target/CPU', round(entry['target_over_cpu'], 3), flush=True)
