from pathlib import Path
import os,subprocess
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-moment-horner')
env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3')
def run(name,cmd,extra=None):
 print('START',name,flush=True)
 with (out/(name+'.log')).open('w') as f:
  p=subprocess.run(cmd,cwd=root,env=dict(env,**(extra or {})),stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:
  print((out/(name+'.log')).read_text()[-8000:],flush=True);raise SystemExit(p.returncode)
 print('PASS',name,flush=True)
run('tests-checked',['ctest','--test-dir','/tmp/accelnet-gpu-stage23/build-gnu-target','-R','predictor_batch|predictor_target_host','--output-on-failure'])
for compiler in ['gnu','nvhpc']:
 run('mixed-degrees-'+compiler,[f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/test_batch_target','--host','--moment'])
for name,uuid in [('h100','GPU-2644154d-7268-af42-6631-59e1f3c6e7f3'),('blackwell','GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9')]:
 run('mixed-degrees-'+name,['/tmp/accelnet-gpu-research/build-nvhpc/bin/test_batch_target','--moment'],{'CUDA_VISIBLE_DEVICES':uuid})
run('memcheck',['/usr/local/cuda-13.1/bin/compute-sanitizer','--tool','memcheck','--leak-check','no','--error-exitcode','1','/tmp/accelnet-gpu-research/build-nvhpc/bin/test_batch_target','--quick'])
run('tests-cpu',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-current','-LE','performance','--output-on-failure'])
run('tests-cpu-performance',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-current','-L','performance','--output-on-failure'])
