from pathlib import Path
import os,subprocess,json,hashlib
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-angular-contraction')
h100='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3';blackwell='GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9'
env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES=h100)
def run(name,cmd,extra=None):
 print('START',name,flush=True)
 with (out/(name+'.log')).open('w') as f:p=subprocess.run(cmd,cwd=root,env=dict(env,**(extra or {})),stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise SystemExit('FAILED '+name+'\n'+(out/(name+'.log')).read_text()[-5000:])
 print('PASS',name,flush=True)
for name,build in [('gnu','/tmp/accelnet-other-descriptors/build-gnu-serial'),('nvhpc','/tmp/accelnet-other-descriptors/build-nvhpc-serial')]:
 run('tests-final-'+name,['ctest','--test-dir',build,'-R','predictor_batch|predictor_target_host','--output-on-failure'])
run('tests-final-h100',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-L','gpu|target-host','--output-on-failure'])
run('tests-final-blackwell',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-nvhpc','-L','gpu','--output-on-failure'],{'CUDA_VISIBLE_DEVICES':blackwell})
run('build-checked',['cmake','--build','/tmp/accelnet-gpu-stage23/build-gnu-target','-j','8','--target','test_batch','test_batch_target'])
run('tests-checked',['ctest','--test-dir','/tmp/accelnet-gpu-stage23/build-gnu-target','-R','predictor_batch|predictor_target_host','--output-on-failure'])
run('memcheck',['/usr/local/cuda-13.1/bin/compute-sanitizer','--tool','memcheck','--leak-check','no','--error-exitcode','99','/tmp/accelnet-gpu-research/build-nvhpc/bin/test_batch_target','--descriptors'])
run('build-lammps',['cmake','--build','/tmp/accelnet-lammps-gpu/build','-j','8'])
run('tests-lammps',['python3','interfaces/lammps/29Aug2024/tests/check_gpu_descriptors.py','--lammps','/tmp/accelnet-lammps-gpu/build/lmp','--fixtures','/tmp/accelnet-gpu-research/build-nvhpc/AccelNetPredictor/target-fixtures','--converter','/tmp/accelnet-gpu-research/build-current/bin/accelnet-model-converter-fortran','--output',str(out/'lammps')])

run('build-cpu',['cmake','--build','/tmp/accelnet-gpu-research/build-current','-j','8','--target','test_batch','accelnet-cpu-regression-benchmark'])
run('tests-cpu',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-current','-LE','performance','--output-on-failure'])
run('tests-cpu-performance',['ctest','--test-dir','/tmp/accelnet-gpu-research/build-current','-L','performance','--output-on-failure'])
for compiler,backend in [('gnu','host'),('nvhpc','host'),('gpu','gpu')]:
 before=out/f'baseline-{compiler}/accelnet-target-benchmark'
 after=('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark' if compiler=='gpu' else f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/accelnet-target-benchmark')
 run('chebyshev-regression-'+compiler,['python3','AccelNetPredictor/benchmark/compare_chebyshev_variants.py','--variant','before',backend,str(before),'--variant','after',backend,after,'--sizes','512','4096','--orders','5','--modes','1','2','--rounds','1','--seconds','0.08','--output',str(out/('chebyshev-'+compiler))])
paths={c:Path(f'/tmp/accelnet-other-descriptors/build-{c}-serial/bin/accelnet-target-benchmark') for c in ['gnu','nvhpc']}
paths['gpu']=Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark')
paths['lammps']=Path('/tmp/accelnet-lammps-gpu/build/lmp')
(out/'final-binaries.json').write_text(json.dumps({c:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for c,p in paths.items()},indent=2)+'\n')
