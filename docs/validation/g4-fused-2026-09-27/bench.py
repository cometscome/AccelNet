from pathlib import Path
import os,subprocess,json,statistics,hashlib,sys
out=Path('/tmp/accelnet-g4-fused');os.sched_setaffinity(0,{6})
label=sys.argv[1];tag=sys.argv[2];rows=[];host=label in ('gnu','nvhpc');backend='host' if host else 'gpu'
exe=Path(f'/tmp/accelnet-other-descriptors/build-{label}-serial/bin/accelnet-target-benchmark') if host else Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark')
uuid={'h100':'GPU-2644154d-7268-af42-6631-59e1f3c6e7f3','blackwell':'GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9'}.get(label,'')
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES=uuid)
variants=[('before',out/'before'/('gpu' if not host else label)),('after',exe)]
for n in ([512] if tag=='regression' else [512,4096]):
 for family in sys.argv[3:] or ['g4','g4-distinct','g4-series']:
  for rep in range(1 if tag=='regression' else 2):
   for variant,binary in variants[::1 if rep%2==0 else -1]:
    cmd=[str(binary),str(n),'8','0.10','1','1.7',backend,family,'no-neighbors']
    p=subprocess.run(cmd,text=True,capture_output=True,env=env,timeout=600)
    name=f'{tag}-{label}-{n}-{family}-{rep}-{variant}';(out/(name+'.log')).write_text(p.stdout+p.stderr)
    if p.returncode:raise SystemExit(name+'\n'+p.stdout+p.stderr)
    samples=[list(map(float,l.split()[2:])) for l in p.stdout.splitlines() if l.startswith('TIMING')]
    phases=[list(map(float,l.split()[3:])) for l in p.stdout.splitlines() if l.startswith('PROFILE')]
    case=next(l.split() for l in p.stdout.splitlines() if l.startswith('CASE '))
    row=dict(backend=label,n=n,family=family,round=rep,variant=variant,command=cmd,sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),max_error=float(case[4]),edges=int(case[3]),samples=samples,phases=phases,ms=1000*statistics.median(x[1] for x in samples),reference_ms=1000*statistics.median(x[0] for x in samples),ratio=statistics.median(x[1]/x[0] for x in samples))
    rows.append(row);print(name,round(row['ms'],3),round(row['ratio'],3),flush=True)
    (out/f'{tag}-{label}.json').write_text(json.dumps(rows,indent=2)+'\n')
