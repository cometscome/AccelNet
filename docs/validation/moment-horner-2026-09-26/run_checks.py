from pathlib import Path
import os, subprocess
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone'); out=Path('/tmp/accelnet-moment-horner')
h100='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3'; blackwell='GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9'
env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES=h100)
def run(name,cmd,extra=None):
 print('START',name,flush=True)
 with (out/(name+'.log')).open('w') as f:
  p=subprocess.run(cmd,cwd=root,env=dict(env,**(extra or {})),stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:
  print((out/(name+'.log')).read_text()[-8000:],flush=True)
  raise SystemExit(f'FAILED {name}: {p.returncode}')
 print('PASS',name,flush=True)
def compare(name,before,after,backend,uuid=h100,orders=['5','10'],sizes=['512','4096'],modes=['2'],rounds='2'):
 run(name,['python3','AccelNetPredictor/benchmark/compare_chebyshev_variants.py','--variant','before',backend,str(before),'--variant','horner',backend,str(after),'--output',str(out/name),'--orders',*orders,'--sizes',*sizes,'--modes',*modes,'--rounds',rounds,'--seconds','0.08'],{'CUDA_VISIBLE_DEVICES':uuid})
for compiler in ['gnu','nvhpc']:
 run('tests-'+compiler,['ctest','--test-dir',f'/tmp/accelnet-other-descriptors/build-{compiler}-serial','-R','predictor_batch|predictor_target_host','--output-on-failure'])
run('tests-h100',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-L','gpu|target-host','--output-on-failure'])
run('tests-blackwell',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-L','gpu','--output-on-failure'],{'CUDA_VISIBLE_DEVICES':blackwell})
for name,uuid in [('h100',h100),('blackwell',blackwell)]:
 compare(name,out/'before-gpu/accelnet-target-benchmark','/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark','gpu',uuid)
for compiler in ['gnu','nvhpc']:
 before=out/f'before-{compiler}/accelnet-target-benchmark'; after=f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/accelnet-target-benchmark'
 compare(compiler,before,after,'host')
 compare('default-'+compiler,before,after,'cpu-shared',orders=['5'],rounds='1')
for label,exe in [('before',str(out/'lmp-before')),('horner','/tmp/accelnet-lammps-gpu/build/lmp')]:
 run('lammps-'+label,['python3','interfaces/lammps/29Aug2024/tests/benchmark_gpu_modes.py','--lammps',exe,'--input','/tmp/accelnet-lammps-gpu/validation-h100/orthogonal-auto-cpu-1rank/in.test','--output',str(out/('lammps-'+label)),'--replicates','10'])
