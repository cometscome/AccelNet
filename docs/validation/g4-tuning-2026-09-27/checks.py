from pathlib import Path
import os,subprocess,json,hashlib,re,sys
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-g4-tuning/checks');out.mkdir(exist_ok=True)
env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3')
def run(name,cmd,extra=None):
 if name.startswith('build-') and name != 'build-lammps':
  help_text=subprocess.check_output(['cmake','--build',cmd[2],'--target','help'],text=True)
  candidates=re.findall(r'^\.\.\. ([A-Za-z0-9_-]+)$',help_text,re.M)
  targets=[t for t in candidates if t.startswith(('test_','write_','dump_','compare_','original_descriptor')) or t in ['accelnet-target-benchmark','accelnet-cpu-regression-benchmark','accelnet-predict','accelnet-descriptor','accelnet-setup-descriptor','accelnet-model-converter-fortran']]
  if name != 'build-cpu':
   required={'test_batch','test_batch_target','test_target_c','test_openmp_target','write_target_fixtures','accelnet-target-benchmark','accelnet-model-converter-fortran'}
   targets=[t for t in targets if t in required]
  assert targets
  cmd+=['--target',*targets]
  (out/(name+'-command.json')).write_text(json.dumps(cmd,indent=2)+'\n')
 print('START',name,flush=True)
 with (out/(name+'.log')).open('w') as f:p=subprocess.run(cmd,cwd=root,env=dict(env,**(extra or {})),stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise SystemExit('FAILED '+name+'\n'+(out/(name+'.log')).read_text()[-5000:])
 print('PASS',name,flush=True)
for name,build in [('gnu','/tmp/accelnet-other-descriptors/build-gnu-serial'),('nvhpc','/tmp/accelnet-other-descriptors/build-nvhpc-serial'),('gpu','/tmp/accelnet-gpu-research/build-nvhpc')]:
 if name == 'gnu' and '--resume' in sys.argv:continue
 run('build-'+name,['cmake','--build',build,'-j','8'])
 if name != 'gpu':run('tests-'+name,['ctest','--test-dir',build,'-R','predictor_batch|predictor_target_host','--output-on-failure'])
run('tests-h100',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-L','gpu|target-host','--output-on-failure'])
run('tests-blackwell',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-L','gpu','--output-on-failure'],{'CUDA_VISIBLE_DEVICES':'GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9'})
run('build-checked',['cmake','--build','/tmp/accelnet-gpu-stage23/build-gnu-target','-j','8'])
run('tests-checked',['ctest','--test-dir','/tmp/accelnet-gpu-stage23/build-gnu-target','-R','predictor_batch|predictor_target_host','--output-on-failure'])
for family in ['g4-distinct','g4-series','behler']:
 run('memcheck-'+family,['/usr/local/cuda-13.1/bin/compute-sanitizer','--tool','memcheck','--leak-check','no','--error-exitcode','99','/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark','8','8','0.001','1','1.7','gpu',family,'no-neighbors'])
run('build-lammps',['cmake','--build','/tmp/accelnet-lammps-gpu/build','-j','8'])
run('tests-lammps',['python3','interfaces/lammps/29Aug2024/tests/check_gpu_descriptors.py','--lammps','/tmp/accelnet-lammps-gpu/build/lmp','--fixtures','/tmp/accelnet-gpu-research/build-nvhpc/AccelNetPredictor/target-fixtures','--converter','/tmp/accelnet-gpu-research/build-current/bin/accelnet-model-converter-fortran','--output',str(out/'lammps')])
run('build-cpu',['cmake','--build','/tmp/accelnet-gpu-research/build-current','-j','8'])
run('tests-cpu',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-current','-LE','performance','--output-on-failure'])
run('tests-cpu-performance',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-current','-L','performance','--output-on-failure'])
