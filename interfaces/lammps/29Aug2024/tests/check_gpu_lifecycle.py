#!/usr/bin/env python3
"""Check workspace growth and complete LAMMPS clear/model reinitialization."""
import argparse
import os
from pathlib import Path
import subprocess
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--lammps',type=Path,required=True)
p.add_argument('--input',type=Path,required=True,help='CPU in.test from check_gpu.py')
p.add_argument('--output',type=Path,required=True)
p.add_argument('--mpiexec',default='/opt/ompi-cuda/bin/mpiexec')
a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
base=a.input.read_text().split('run 0')[0]
def dump(name):
    return f'write_dump all custom {name}.dump id type x y z fx fy fz c_ea modify sort id format float %.17g\n'
for backend in ['cpu','gpu']:
    directory=a.output/backend
    directory.mkdir(exist_ok=True)
    def setup(neigh):
        if backend=='cpu': return base
        return base.replace('clear',f'clear\npackage gpu 1 neigh {neigh} newton on split 1',1).replace('pair_style accelnet ','pair_style accelnet/gpu ')
    text=setup('yes')+'run 0\n'+dump('initial')
    text+='unfix integrate\nreplicate 2 2 2\nfix integrate all nve\nrun 2\n'+dump('grown')
    text+=setup('no')+'run 2\n'+dump('reloaded')
    (directory/'in.test').write_text(text)
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY')
    with (directory/'stdout.txt').open('w') as f:
        subprocess.run([a.mpiexec,'--bind-to','none','-n','2',str(a.lammps.resolve()),'-in','in.test'],
                       cwd=directory,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=120,check=True)
def read(path):
    lines=path.read_text().splitlines()
    start=next(i for i,line in enumerate(lines) if line.startswith('ITEM: ATOMS'))+1
    return np.array([list(map(float,line.split())) for line in lines[start:]])
for snapshot in ['initial','grown','reloaded']:
    ref=read(a.output/'cpu'/f'{snapshot}.dump')
    got=read(a.output/'gpu'/f'{snapshot}.dump')
    np.testing.assert_allclose(got,ref,atol=2e-8,rtol=2e-9)
    print(f'{snapshot}: {len(ref)} atoms, max force error {np.max(np.abs(got[:,5:8]-ref[:,5:8])):.3e}',flush=True)
