from pathlib import Path
import os,subprocess,time
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone')
out=Path('/tmp/accelnet-cheb-common')
env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3')
def run(name,cmd,extra=None):
 print('START',name,flush=True)
 with (out/(name+'.log')).open('w') as f:
  p=subprocess.run(cmd,cwd=root,env=dict(env,**(extra or {})),stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:
  print((out/(name+'.log')).read_text()[-6000:],flush=True)
  raise SystemExit(f'FAILED {name}: {p.returncode}')
 print('PASS',name,flush=True)
for compiler in ['gnu','nvhpc']:
 run('tests-'+compiler,['ctest','--test-dir',f'/tmp/accelnet-other-descriptors/build-{compiler}-serial','-R','predictor_batch|predictor_target_host','--output-on-failure'])
run('tests-cpu',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-current','-LE','performance','--output-on-failure'])
for name,uuid in [('h100','GPU-2644154d-7268-af42-6631-59e1f3c6e7f3'),('blackwell','GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9')]:
 run('tests-'+name,['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-L','gpu','--output-on-failure'],{'CUDA_VISIBLE_DEVICES':uuid})
run('h100-variants',['python3','AccelNetPredictor/benchmark/compare_chebyshev_variants.py','--variant','before','gpu',str(out/'before-gpu/accelnet-target-benchmark'),'--variant','geometry','gpu',str(out/'geometry-gpu-bin/accelnet-target-benchmark'),'--variant','clenshaw','gpu','/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark','--output',str(out/'h100'),'--modes','1','--rounds','2','--seconds','0.08'])
run('nvhpc-variants',['python3','AccelNetPredictor/benchmark/compare_chebyshev_variants.py','--variant','before','host',str(out/'before-nvhpc/accelnet-target-benchmark'),'--variant','clenshaw','host','/tmp/accelnet-other-descriptors/build-nvhpc-serial/bin/accelnet-target-benchmark','--variant','default','cpu-shared','/tmp/accelnet-other-descriptors/build-nvhpc-serial/bin/accelnet-target-benchmark','--output',str(out/'nvhpc'),'--modes','1','--rounds','2','--seconds','0.08'])
for compiler in ['gnu','nvhpc']:
 run('moment-'+compiler,['python3','AccelNetPredictor/benchmark/compare_chebyshev_variants.py','--variant','before','host',str(out/f'before-{compiler}/accelnet-target-benchmark'),'--variant','default','cpu-shared',f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/accelnet-target-benchmark','--output',str(out/f'moment-{compiler}'),'--modes','2','--rounds','1','--seconds','0.12'])
for label,exe in [('before',str(out/'lmp-before')),('after','/tmp/accelnet-lammps-gpu/build/lmp')]:
 run('lammps-'+label,['python3','interfaces/lammps/29Aug2024/tests/benchmark_gpu_modes.py','--lammps',exe,'--input','/tmp/accelnet-lammps-gpu/validation-h100/orthogonal-auto-cpu-1rank/in.test','--output',str(out/('lammps-'+label)),'--replicates','10'])
