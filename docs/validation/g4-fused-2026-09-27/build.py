from pathlib import Path
import subprocess,sys,shutil
out=Path('/tmp/accelnet-g4-fused')
for n in sys.argv[1:] or ['gnu','nvhpc','gpu']:
 b='/tmp/accelnet-gpu-research/build-nvhpc' if n=='gpu' else f'/tmp/accelnet-other-descriptors/build-{n}-serial'
 with (out/f'build-final-{n}.log').open('w') as f:
  subprocess.run(['cmake','--build',b,'-j','8','--target','accelnet-target-benchmark','test_batch_target','test_batch','test_target_c','write_target_fixtures'],stdout=f,stderr=subprocess.STDOUT,check=True)
 print('BUILT',n,flush=True)
