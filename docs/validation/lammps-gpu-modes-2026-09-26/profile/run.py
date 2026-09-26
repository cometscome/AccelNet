from pathlib import Path
import os,subprocess
base=Path('/tmp/accelnet-modes-profile')
env=dict(os.environ,CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3',OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',OPENBLAS_NUM_THREADS='1')
for atoms,mode in [(24000,'direct'),(5184,'moment'),(5184,'direct'),(192,'moment'),(192,'direct')]:
 d=base/f'{atoms}-{mode}';d.mkdir(exist_ok=True)
 command=['nsys','profile','--trace=cuda','--sample=none','--cpuctxsw=none','--force-overwrite=true','--output=trace','taskset','-c','6','/home/nagai/AccelNetGPU/lammps-accelnet-gpu/lmp','-in',f'/tmp/accelnet-lammps-modes-h100/{atoms}-{mode}-0/in.benchmark']
 with (d/'stdout.txt').open('w') as f: subprocess.run(command,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=300)
 print('done',atoms,mode,flush=True)
