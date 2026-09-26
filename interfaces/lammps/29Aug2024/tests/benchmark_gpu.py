#!/usr/bin/env python3
"""LAMMPS end-to-end timings; run after check_gpu.py generated an input.
Includes integration, neighbor rebuilds and MPI. No initialization is timed.
"""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--lammps', type=Path, required=True)
p.add_argument('--input', type=Path, required=True, help='CPU in.test from check_gpu.py')
p.add_argument('--output', type=Path, required=True)
p.add_argument('--replicates', nargs='+', type=int, default=[1, 2, 6, 10])
p.add_argument('--samples', type=int, default=3)
p.add_argument('--steps', type=int, default=50)
p.add_argument('--cpu-ranks', nargs='+', type=int, default=[1, 4])
p.add_argument('--mpiexec', default='/opt/ompi-cuda/bin/mpiexec')
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
base = a.input.read_text().split('run 0')[0]
base = base.replace('thermo 1', f'thermo {a.steps}')
base = base.replace('neigh_modify every 1 delay 0 check no', 'neigh_modify every 10 delay 0 check yes')
base = base.replace('atom_modify sort 1 0.0', 'atom_modify sort 100 0.0')
base = '\n'.join(line for line in base.splitlines() if not line.startswith('compute ea'))
report = []
for rep in a.replicates:
    for backend, ranks in [('cpu', n) for n in a.cpu_ranks]+[('no', 1), ('yes', 1)]:
        timings = []
        pair_times = []
        for sample in range(a.samples):
            directory = a.output/f'{24*rep**3}-{backend}-{ranks}rank-{sample}'
            directory.mkdir(exist_ok=True)
            text = re.sub(r'replicate \d+ \d+ \d+', f'replicate {rep} {rep} {rep}', base)
            if backend != 'cpu':
                text = text.replace('clear', f'clear\npackage gpu 1 neigh {backend} newton on split 1', 1)
                text = text.replace('pair_style accelnet ', 'pair_style accelnet/gpu ')
            text += f'\nrun 10\nrun {a.steps}\n'
            (directory/'in.benchmark').write_text(text)
            command = [str(a.lammps.resolve()), '-in', 'in.benchmark']
            if ranks > 1:
                command = [a.mpiexec, '--bind-to', 'core', '-n', str(ranks)]+command
            else:
                command = ['taskset', '-c', '6']+command
            env = dict(os.environ, OMP_NUM_THREADS='1', OMP_TARGET_OFFLOAD='MANDATORY', OPENBLAS_NUM_THREADS='1')
            with (directory/'stdout.txt').open('w') as f:
                subprocess.run(command, cwd=directory, env=env, stdout=f, stderr=subprocess.STDOUT,
                               timeout=600, check=True)
            log = (directory/'stdout.txt').read_text()
            loops = re.findall(r'Loop time of ([\d.eE+-]+)', log)
            timings.append(float(loops[-1])/a.steps)
            pairs = re.findall(r'^Pair\s+\|\s*([\d.eE+-]+)\s*\|\s*([\d.eE+-]+)', log, re.M)
            pair_times.append(float(pairs[-1][1])/a.steps)
        item = dict(atoms=24*rep**3, backend=backend, ranks=ranks, steps=a.steps,
                    seconds_per_step=statistics.median(timings), samples=timings,
                    pair_seconds_per_step=statistics.median(pair_times))
        report.append(item)
        (a.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(item), flush=True)
