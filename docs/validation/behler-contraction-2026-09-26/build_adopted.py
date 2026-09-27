from pathlib import Path
import subprocess
out=Path('/tmp/accelnet-angular-contraction');root=Path('/home/nagai/AccelNetGPU/AccelNet-clone')
for name,build in [('gnu','/tmp/accelnet-other-descriptors/build-gnu-serial'),('nvhpc','/tmp/accelnet-other-descriptors/build-nvhpc-serial'),('gpu','/tmp/accelnet-gpu-research/build-nvhpc')]:
 cmd=['cmake','--build',build,'-j','8','--target','accelnet-target-benchmark','test_batch_target','test_batch','test_target_c','write_target_fixtures']
 with (out/f'build-final-{name}.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise SystemExit('FAIL '+name)
 print('BUILT',name,flush=True)
