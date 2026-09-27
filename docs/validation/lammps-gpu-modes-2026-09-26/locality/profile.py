from pathlib import Path
import subprocess,os
base=Path('/tmp/accelnet-moment-analysis')
env=dict(os.environ,CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3',OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',OPENBLAS_NUM_THREADS='1')
for atoms,variant in [(24000,'tile1'),(24000,'tile32'),(24000,'tile128'),(24000,'tile512'),(5184,'tile1'),(192,'tile1')]:
 d=base/'profile'/f'{atoms}-{variant}';d.mkdir(parents=True,exist_ok=True)
 cmd=['nsys','profile','--trace=cuda','--sample=none','--cpuctxsw=none','--force-overwrite=true','--output=trace','taskset','-c','6',str(base/variant/'lmp'),'-in',f'/tmp/accelnet-lammps-modes-h100/{atoms}-moment-0/in.benchmark']
 with (d/'stdout.txt').open('w') as f:subprocess.run(cmd,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=240)
 print('profiled',atoms,variant,flush=True)
