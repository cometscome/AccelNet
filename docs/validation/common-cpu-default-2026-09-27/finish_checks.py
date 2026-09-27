from pathlib import Path
import subprocess, os, json
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-common-default')
commands=[]
def run(name,cmd):
 commands.append(dict(name=name,command=cmd));(out/'check-commands.json').write_text(json.dumps(commands,indent=2)+'\n')
 with (out/(name+'.log')).open('w') as f:
  p=subprocess.run(cmd,cwd=root,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
 print(name,p.returncode,flush=True)
 if p.returncode:print((out/(name+'.log')).read_text()[-5000:],flush=True);raise SystemExit(p.returncode)
for compiler in ['gnu','nvhpc']:
 build=f'/tmp/accelnet-other-descriptors/build-{compiler}-serial'
 run('configure-gate-'+compiler,['cmake','-S',str(root),'-B',build,'-DACCELNET_TEST_COMMON_CPU_PERFORMANCE=ON'])
 run('build-gate-'+compiler,['cmake','--build',build,'-j','6','--target','accelnet-cpu-regression-benchmark'])
 run('gate-'+compiler,['ctest','--test-dir',build,'-R','^predictor_common_cpu_performance$','--output-on-failure'])
 (out/('gate-'+compiler+'.json')).write_text(Path(build,'AccelNetPredictor/common-cpu-performance.json').read_text())
# Rebuild the ordinary non-target configuration: tests must run on current code.
build='/tmp/accelnet-gpu-research/build-current'
help_text=subprocess.check_output(['cmake','--build',build,'--target','help'],text=True)
import re
targets=[t for t in re.findall(r'^\.\.\. ([A-Za-z0-9_-]+)$',help_text,re.M) if t.startswith(('test_','write_','dump_','compare_','original_descriptor')) or t in ['accelnet-predict','accelnet-descriptor','accelnet-setup-descriptor','accelnet-model-converter-fortran']]
run('build-ordinary',['cmake','--build',build,'-j','6','--target',*targets])
run('tests-ordinary',['ctest','--test-dir',build,'-LE','performance','--output-on-failure'])
