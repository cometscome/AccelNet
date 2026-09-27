#!/usr/bin/env python3
"""Compare LAMMPS AccelNet CPU and GPU E/F/virial and short trajectories.
Requires a CUDA GPU, MPI launcher, and the immutable Ti/O golden networks.
Outputs logs, snapshots, and a machine-readable report in --output.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--lammps', type=Path, required=True)
p.add_argument('--golden', type=Path, required=True)
p.add_argument('--elements', nargs=2, default=['Ti', 'O'], help='Embedded global species order')
p.add_argument('--output', type=Path, required=True)
p.add_argument('--mpiexec', default='/opt/ompi-cuda/bin/mpiexec')
p.add_argument('--ranks', type=int, nargs='+', default=[1, 2, 4])
p.add_argument('--gpus', type=int, default=1)
p.add_argument('--steps', type=int, default=20)
p.add_argument('--replicate', type=int, default=1)
p.add_argument('--rebuild-every', type=int, default=1)
p.add_argument('--cases', nargs='+', default=['orthogonal', 'triclinic', 'empty'])
p.add_argument('--modes', nargs='+', default=['auto'])
p.add_argument('--g5-mode', choices=['auto','direct','moment'], default='auto')
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=True)
a.output = a.output.resolve()
a.golden = a.golden.resolve()
lines = (a.golden/'structure0001.xsf').read_text().splitlines()
i = lines.index('PRIMVEC')
box = np.array([list(map(float, line.split())) for line in lines[i+1:i+4]])
i = lines.index('PRIMCOORD')
n = int(lines[i+1].split()[0])
atoms = [line.split() for line in lines[i+2:i+2+n]]

def datafile(case):
    length = np.diag(box).copy()
    if case == 'empty':
        length[0] *= 8  # ranks with no local centers, periodic boundaries preserved
    result = f'AccelNet regression\n\n{n} atoms\n2 atom types\n\n'
    for size, axis in zip(length, 'xyz'):
        result += f'0 {size:.17g} {axis}lo {axis}hi\n'
    if case == 'triclinic':
        result += '0.7 -0.4 0.3 xy xz yz\n'
    result += '\nMasses\n\n1 47.867\n2 15.999\n\nAtoms # atomic\n\n'
    for j, at in enumerate(atoms, 1):
        xyz = np.array(list(map(float, at[1:4])))
        if case == 'triclinic':
            frac = xyz/np.diag(box)
            xyz += [0.7*frac[1]-0.4*frac[2], 0.3*frac[2], 0]
        result += f'{j} {a.elements.index(at[0])+1} '+ ' '.join(f'{x:.17g}' for x in xyz)+'\n'
    path = a.output/f'{case}.data'
    path.write_text(result)
    return path

def run(case, mode, backend, ranks):
    name = f'{case}-{mode}-{backend}-{ranks}rank'
    directory = a.output/name
    directory.mkdir(exist_ok=True)
    gpu = backend != 'cpu'
    package = f'package gpu {a.gpus} neigh {backend} newton on split 1' if gpu else ''
    processors = f'processors {ranks} 1 1' if case in ['empty', 'migration'] else ''
    drift = 'velocity all set 100.0 80.0 60.0 sum yes units box' if case == 'migration' else ''
    text = f'''clear
{package}
units metal
atom_style atomic
boundary p p p
{processors}
read_data {datafile(case)}
replicate {a.replicate} {a.replicate} {a.replicate}
pair_style accelnet{'/gpu' if gpu else ''} {mode} {a.golden}/{a.elements[0]}.nn.ascii {a.golden}/{a.elements[1]}.nn.ascii g5 {a.g5_mode if gpu else "direct"}
pair_coeff * *
neighbor 0.6 bin
neigh_modify every {a.rebuild_every} delay 0 check no
atom_modify sort 1 0.0
velocity all create 100 4928459 mom yes rot no dist gaussian
{drift}
fix integrate all nve
timestep 0.0001
compute ea all pe/atom
compute vir all pressure NULL virial
thermo 1
thermo_style custom step pe c_vir[1] c_vir[2] c_vir[3] c_vir[4] c_vir[5] c_vir[6]
thermo_modify format float %.17g lost error
run 0
write_dump all custom initial.dump id type x y z fx fy fz c_ea modify sort id format float %.17g
run {a.steps}
write_dump all custom final.dump id type x y z fx fy fz c_ea modify sort id format float %.17g
'''
    (directory/'in.test').write_text(text)
    cmd = [str(a.lammps.resolve()), '-in', 'in.test']
    if ranks > 1:
        cmd = [a.mpiexec, '--bind-to', 'none', '-n', str(ranks)]+cmd
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_TARGET_OFFLOAD='MANDATORY')
    with (directory/'stdout.txt').open('w') as f:
        proc = subprocess.run(cmd, cwd=directory, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=240)
    if proc.returncode:
        raise RuntimeError(f'{name} failed: '+(directory/'stdout.txt').read_text()[-2500:])
    stdout = (directory/'stdout.txt').read_text()
    thermo = []
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) == 8 and re.fullmatch(r'\d+', fields[0]):
            try:
                thermo.append(list(map(float, fields)))
            except ValueError:
                pass
    return directory, np.array(thermo)

def read_dump(path):
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('ITEM: ATOMS'))+1
    return np.array([list(map(float, line.split())) for line in lines[start:]])

report = []
for case in a.cases:
    for mode in a.modes:
        reference, rt = run(case, mode, 'cpu', 1)
        for ranks in a.ranks:
            for backend in ['cpu', 'no', 'yes', 'hybrid']:
                if backend == 'cpu' and ranks == 1:
                    continue
                directory, thermo = run(case, mode, backend, ranks)
                assert np.isfinite(thermo).all()
                np.testing.assert_allclose(thermo, rt, atol=2e-5, rtol=2e-9)
                diffs = {}
                for snapshot in ['initial', 'final']:
                    ref = read_dump(reference/f'{snapshot}.dump')
                    got = read_dump(directory/f'{snapshot}.dump')
                    assert np.isfinite(got).all()
                    np.testing.assert_allclose(got[:, :2], ref[:, :2], atol=0, rtol=0)
                    np.testing.assert_allclose(got[:, 2:], ref[:, 2:], atol=2e-8, rtol=2e-9)
                    diffs[snapshot] = {'max_force_error': float(np.max(np.abs(got[:, 5:8]-ref[:, 5:8]))),
                                      'max_atom_energy_error': float(np.max(np.abs(got[:, 8]-ref[:, 8])))}
                entry = dict(case=case, mode=mode, g5_mode=a.g5_mode, backend=backend, ranks=ranks,
                             max_total_energy_error=float(np.max(np.abs(thermo[:, 1]-rt[:, 1]))),
                             max_virial_pressure_error=float(np.max(np.abs(thermo[:, 2:]-rt[:, 2:]))), **diffs)
                report.append(entry)
                print(json.dumps(entry), flush=True)
                (a.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
