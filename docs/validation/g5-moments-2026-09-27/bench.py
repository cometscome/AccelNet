from pathlib import Path
import os,subprocess,json,statistics,hashlib,re
out=Path('/tmp/accelnet-g5-moment'); os.sched_setaffinity(0,{6})
rows=[];audit=[]
backends=[('gnu','/tmp/accelnet-other-descriptors/build-gnu-serial/bin/accelnet-target-benchmark','cpu-shared',''),('nvhpc','/tmp/accelnet-other-descriptors/build-nvhpc-serial/bin/accelnet-target-benchmark','cpu-shared',''),('h100','/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark','gpu','GPU-2644154d-7268-af42-6631-59e1f3c6e7f3'),('blackwell','/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark','gpu','GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9')]
for name,exe,backend,gpu in backends:
 if not gpu:
  build=Path(exe).parents[1]
  flags=[(build/f'AccelNetPredictor/CMakeFiles/{t}.dir/flags.make').read_text() for t in ['AccelNetPredictor','AccelNetTarget','accelnet-target-benchmark']]
  symbols=subprocess.check_output(['nm','-u',exe],text=True)
  assert not re.search(r'\b(?:GOMP_|__kmpc_|__nvomp_|omp_)[A-Za-z0-9_]*',symbols)
  assert not any(re.search(r'(?:^|\s)(?:-fopenmp|-mp)(?:=\S+)?(?:\s|$)',x) for x in flags)
 else: flags=[]
 audit.append(dict(backend=name,path=exe,sha256=hashlib.sha256(Path(exe).read_bytes()).hexdigest(),flags=flags,gpu=gpu,cpu_affinity=[6]))
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY')
 if gpu:env['CUDA_VISIBLE_DEVICES']=gpu
 cases=[(n,nb,mode) for n in [64,512] for nb in [16,32,64,128] for mode in [1,3]]
 for rnd in range(2):
  for n,nb,mode in (cases if rnd==0 else cases[::-1]):
   cmd=[exe,str(n),'4','0.05',str(mode),'1.7',backend,'g5-scaling','no-neighbors',str(nb)]
   p=subprocess.run(cmd,env=env,text=True,capture_output=True,timeout=600)
   tag=f'{name}-n{n}-nb{nb}-mode{mode}-round{rnd}'
   (out/(tag+'.log')).write_text(p.stdout+p.stderr)
   if p.returncode:raise RuntimeError(tag+' failed '+p.stderr[-1000:])
   samples=[list(map(float,l.split()[2:])) for l in p.stdout.splitlines() if l.startswith('TIMING ')]
   case=next(l.split() for l in p.stdout.splitlines() if l.startswith('CASE '))
   phases=[list(map(float,l.split()[3:])) for l in p.stdout.splitlines() if l.startswith('PROFILE ')]
   assert len(samples)==5
   row=dict(backend=name,natoms=n,neighbors=nb,mode=mode,round=rnd,command=cmd,samples=samples,old_ms=1000*statistics.median(t[0] for t in samples),new_ms=1000*statistics.median(t[1] for t in samples),error=float(case[4]),phases_ms=[1000*statistics.median(x[i] for x in phases) for i in range(8)])
   rows.append(row);(out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
   print(tag,'old',round(row['old_ms'],4),'new',round(row['new_ms'],4),flush=True)
(out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
