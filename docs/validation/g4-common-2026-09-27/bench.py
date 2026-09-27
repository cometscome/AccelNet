from pathlib import Path
import os,subprocess,json,statistics,sys
out=Path('/tmp/accelnet-g4-improve');os.sched_setaffinity(0,{6})
label=sys.argv[1];tag=sys.argv[2];rows=[]
backend='host' if label in ('gnu','nvhpc') else 'gpu'
exe=Path(f'/tmp/accelnet-other-descriptors/build-{label}-serial/bin/accelnet-target-benchmark') if backend=='host' else Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark')
uuid={'h100':'GPU-2644154d-7268-af42-6631-59e1f3c6e7f3','blackwell':'GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9'}.get(label,'')
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES=uuid)
for n in ([512] if tag == 'regression' else [512,4096]):
 for family in sys.argv[3:] or ['g4','g4-series','g5','behler']:
  for variant,binary in [('before',out/'before'/('gpu' if backend=='gpu' else label)),('after',exe)]:
   if family=='g4-distinct' and variant=='before': continue
   cmd=[str(binary),str(n),'8','0.10','1','1.7',backend,family,'no-neighbors']
   p=subprocess.run(cmd,text=True,capture_output=True,env=env);log=f'{tag}-{label}-{n}-{family}-{variant}.log';(out/log).write_text(p.stdout+p.stderr)
   if p.returncode:raise SystemExit(p.stderr)
   samples=[list(map(float,l.split()[2:])) for l in p.stdout.splitlines() if l.startswith('TIMING')]
   phases=[list(map(float,l.split()[3:])) for l in p.stdout.splitlines() if l.startswith('PROFILE')]
   row=dict(backend=label,n=n,family=family,variant=variant,command=cmd,max_error=float(next(l.split()[4] for l in p.stdout.splitlines() if l.startswith('CASE '))),samples=samples,phases=phases,ms=1000*statistics.median(x[1] for x in samples),reference_ms=1000*statistics.median(x[0] for x in samples),ratio=statistics.median(x[1]/x[0] for x in samples))
   rows.append(row);print(n,family,variant,round(row['ms'],3),round(row['ratio'],3),flush=True)
   (out/f'{tag}-{label}.json').write_text(json.dumps(rows,indent=2)+'\n')
