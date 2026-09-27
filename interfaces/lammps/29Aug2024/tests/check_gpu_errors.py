#!/usr/bin/env python3
"""Ensure unsupported inputs fail explicitly without a crash or MPI hang."""
import argparse
import os
from pathlib import Path
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--lammps', type=Path, required=True)
p.add_argument('--input', type=Path, required=True, help='neigh no in.test from check_gpu.py')
p.add_argument('--output', type=Path, required=True)
p.add_argument('--mpiexec', default='/opt/ompi-cuda/bin/mpiexec')
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
base = a.input.read_text().split('run 0')[0]+'run 0\n'
line = next(line for line in base.splitlines() if line.startswith('pair_style'))
paths = line.split()[-2:]
malformed = (a.output/'malformed.nn.ascii').resolve()
malformed.write_text('4\n56\n')
cases = {
    'newton': (base.replace('newton on', 'newton off'), 'requires newton on'),
    'missing': (base.replace(paths[0], str(malformed)+'.missing'), 'Cannot open network'),
    'malformed': (base.replace(paths[0], str(malformed)), 'Invalid or incomplete network'),
    'ordering': (base.replace(line, line.replace(' '.join(paths), ' '.join(reversed(paths)))), 'species ordering'),
    'stress': (base.replace('run 0', 'compute stress all stress/atom NULL\ncompute sum all reduce sum c_stress[1]\nthermo_style custom step c_sum\nrun 0'), 'per-atom stress'),
    'excluded': (base.replace('run 0', 'neigh_modify exclude type 1 2\nrun 0'), 'exclusions'),
}
for name, (text, expected) in cases.items():
    directory = a.output/name
    directory.mkdir(exist_ok=True)
    (directory/'in.test').write_text(text)
    cmd = [a.mpiexec, '--bind-to', 'none', '-n', '2', str(a.lammps.resolve()), '-nonbuf', '-in', 'in.test']
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_TARGET_OFFLOAD='MANDATORY')
    proc = subprocess.run(cmd, cwd=directory, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, timeout=60)
    (directory/'stdout.txt').write_text(proc.stdout)
    assert proc.returncode != 0 and expected in proc.stdout, (name, proc.returncode, proc.stdout)
    assert 'Segmentation fault' not in proc.stdout and 'CUDA_ERROR_ILLEGAL' not in proc.stdout, name
    print(f'{name}: explicit MPI-safe rejection', flush=True)
