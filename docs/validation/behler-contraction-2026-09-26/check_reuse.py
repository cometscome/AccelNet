from pathlib import Path
import json,subprocess,time,os
out=Path('/tmp/accelnet-angular-contraction');root=Path('/home/nagai/AccelNetGPU/AccelNet-clone')
while True:
 try:done=len(json.loads((out/'final-timings.json').read_text()))==96
 except (FileNotFoundError,json.JSONDecodeError):done=False
 if done:break
 time.sleep(5)
for name,cmd in [
 ('configure-reuse',['cmake','-S',str(out/'source-reuse'),'-B',str(out/'build-reuse-gnu'),'-DCMAKE_BUILD_TYPE=Release','-DACCELNET_BUILD_OPENMP_TARGET=ON','-DACCELNET_TARGET_SERIAL=ON']),
 ('build-reuse',['cmake','--build',str(out/'build-reuse-gnu'),'-j','8','--target','accelnet-target-benchmark','test_batch_target']),
 ('test-reuse',['taskset','-c','6',str(out/'build-reuse-gnu/bin/test_batch_target'),'--host','--descriptors']),
 ('repeat-chebyshev',['python3',str(root/'AccelNetPredictor/benchmark/compare_chebyshev_variants.py'),'--variant','before','host',str(out/'baseline-gnu/accelnet-target-benchmark'),'--variant','current','host','/tmp/accelnet-other-descriptors/build-gnu-serial/bin/accelnet-target-benchmark','--variant','reuse','host',str(out/'build-reuse-gnu/bin/accelnet-target-benchmark'),'--sizes','4096','--modes','1','--rounds','2','--seconds','0.08','--output',str(out/'chebyshev-repeat-gnu')])]:
 print('START',name,flush=True)
 with (out/(name+'.log')).open('w') as f:p=subprocess.run(cmd,cwd=root,stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise SystemExit('FAIL '+name)
 print('PASS',name,flush=True)
