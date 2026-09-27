from pathlib import Path
import subprocess,os,json,re
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-g5-moment')
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3');commands=[]
def run(name,cmd,extra=None):
 commands.append(dict(name=name,command=cmd,environment=extra or {}));(out/'check-commands.json').write_text(json.dumps(commands,indent=2)+'\n')
 with (out/(name+'.log')).open('w') as f:p=subprocess.run(cmd,cwd=root,env=dict(env,**(extra or {})),stdout=f,stderr=subprocess.STDOUT)
 print(name,p.returncode,flush=True)
 if p.returncode:raise RuntimeError(name+'\n'+(out/(name+'.log')).read_text()[-4000:])
for name,build in [('gnu','/tmp/accelnet-other-descriptors/build-gnu-serial'),('nvhpc','/tmp/accelnet-other-descriptors/build-nvhpc-serial'),('checked','/tmp/accelnet-gpu-stage23/build-gnu-target')]:
 run('verified-tests-'+name,['ctest','--test-dir',build,'-R','predictor_batch|predictor_target_host','--output-on-failure','-V'])
run('verified-tests-h100',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-L','gpu|target-host','--output-on-failure','-V'])
run('verified-tests-blackwell',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-R','predictor_target_(g5_moments|equivalence|other_descriptors|c_api)$','--output-on-failure','-V'],{'CUDA_VISIBLE_DEVICES':'GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9'})
run('verified-build-lammps',['cmake','--build','/tmp/accelnet-lammps-gpu/build','-j','6'])
run('verified-lammps-moment',['python3','interfaces/lammps/29Aug2024/tests/check_gpu_descriptors.py','--lammps','/tmp/accelnet-lammps-gpu/build/lmp','--fixtures','/tmp/accelnet-gpu-research/build-nvhpc/AccelNetPredictor/target-fixtures','--converter','/tmp/accelnet-gpu-research/build-current/bin/accelnet-model-converter-fortran','--g5-mode','moment','--output',str(out/'lammps-moment')])
run('verified-memcheck',['/usr/local/cuda-13.1/bin/compute-sanitizer','--tool','memcheck','--leak-check','no','--error-exitcode','99','/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark','8','4','0.001','3','1.7','gpu','g5-scaling','no-neighbors','64'])
build='/tmp/accelnet-gpu-research/build-current'
help_text=subprocess.check_output(['cmake','--build',build,'--target','help'],text=True)
targets=[t for t in re.findall(r'^\.\.\. ([A-Za-z0-9_-]+)$',help_text,re.M) if t.startswith(('test_','write_','dump_','compare_','original_descriptor')) or t in ['accelnet-predict','accelnet-descriptor','accelnet-setup-descriptor','accelnet-model-converter-fortran']]
run('verified-build-ordinary',['cmake','--build',build,'-j','6','--target',*targets])
run('verified-tests-ordinary',['ctest','--test-dir',build,'-LE','performance','--output-on-failure'])
for compiler in ['gnu','nvhpc']:
 build=f'/tmp/accelnet-other-descriptors/build-{compiler}-serial'
 run('verified-gate-'+compiler,['ctest','--test-dir',build,'-R','^predictor_common_cpu_performance$','--output-on-failure'])
 (out/('gate-'+compiler+'.json')).write_text(Path(build,'AccelNetPredictor/common-cpu-performance.json').read_text())
