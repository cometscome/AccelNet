from pathlib import Path
import os,subprocess,statistics,json,hashlib
out=Path('/tmp/accelnet-angular-contraction');os.sched_setaffinity(0,{6});rows=[]
for label,uuid in [('gnu',''),('nvhpc',''),('h100','GPU-2644154d-7268-af42-6631-59e1f3c6e7f3'),('blackwell','GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9')]:
 host=label in ['gnu','nvhpc']; backend='host' if host else 'gpu';compiler=label if host else 'gpu'
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES=uuid)
 before=out/f'baseline-{compiler}/accelnet-target-benchmark'
 after=Path(f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/accelnet-target-benchmark') if host else Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark')
 for n in [512,4096]:
  for family in ['lj','g4','g5','behler','g4-series','g5-series']:
   # Dense bases at 4096 expose a different cost balance from the four-output models.
   for variant,exe in [('before',before),('after',after)]:
    cmd=[str(exe),str(n),'8','0.06','1','1.7',backend,family,'no-neighbors']
    p=subprocess.run(cmd,env=env,text=True,capture_output=True);name=f'final-{label}-{n}-{family}-{variant}'
    (out/(name+'.log')).write_text(p.stdout+p.stderr)
    if p.returncode:raise SystemExit(name+'\n'+p.stdout+p.stderr)
    samples=[list(map(float,l.split()[2:])) for l in p.stdout.splitlines() if l.startswith('TIMING ')]
    phases=[list(map(float,l.split()[3:])) for l in p.stdout.splitlines() if l.startswith('PROFILE ')]
    case=next(l.split() for l in p.stdout.splitlines() if l.startswith('CASE '))
    row=dict(backend=label,n=n,family=family,variant=variant,command=cmd,sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),max_error=float(case[4]),edges=int(case[3]),samples=samples,phases=phases,ms=1000*statistics.median(s[1] for s in samples),ratio=statistics.median(s[1]/s[0] for s in samples))
    rows.append(row);(out/'final-timings.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,round(row['ms'],3),round(row['ratio'],3),flush=True)
