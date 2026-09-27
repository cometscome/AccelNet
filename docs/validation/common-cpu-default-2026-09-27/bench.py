from pathlib import Path
import os, subprocess, json, statistics, hashlib, re
out=Path('/tmp/accelnet-common-default'); os.sched_setaffinity(0,{6})
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY')
rows=[];audit=[]
for compiler in ['gnu','nvhpc']:
 build=Path('/tmp/accelnet-other-descriptors')/f'build-{compiler}-serial'
 exe=build/'bin/accelnet-target-benchmark'
 flags=[(build/f'AccelNetPredictor/CMakeFiles/{t}.dir/flags.make').read_text() for t in ['AccelNetPredictor','AccelNetTarget','accelnet-target-benchmark']]
 symbols=subprocess.check_output(['nm','-u',str(exe)],text=True)
 assert not re.search(r'\b(?:GOMP_|__kmpc_|__nvomp_|omp_)[A-Za-z0-9_]*',symbols)
 assert not any(re.search(r'(?:^|\s)(?:-fopenmp|-mp)(?:=\S+)?(?:\s|$)',x) for x in flags)
 audit.append(dict(compiler=compiler,executable=str(exe),sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),flags=flags,no_openmp_symbols=True))
 cases=[(family,n,mode,1.7) for family in ['lj','g4-distinct','g4-series','g5','g5-series','behler','lj-behler','chebyshev'] for n in [8,512,4096] for mode in ([1,2] if family=='chebyshev' else [0] if family in ['g5','g5-series','behler','lj-behler'] else [1])]
 cases += [('g5-series',n,0,1.2) for n in [8,512]]
 for rnd in range(2):
  for family,n,mode,spacing in (cases if rnd==0 else cases[::-1]):
   command=[str(exe),str(n),'8','0.10',str(mode),str(spacing),'cpu-shared',family,'no-neighbors']
   p=subprocess.run(command,env=env,text=True,capture_output=True,timeout=600)
   name=f'{compiler}-{family}-{n}-mode{mode}-spacing{spacing}-round{rnd}'
   (out/(name+'.log')).write_text(p.stdout+p.stderr)
   if p.returncode:raise RuntimeError(name+' failed')
   samples=[list(map(float,l.split()[2:])) for l in p.stdout.splitlines() if l.startswith('TIMING ')]
   case=next(l.split() for l in p.stdout.splitlines() if l.startswith('CASE '))
   assert len(samples)==5
   row=dict(compiler=compiler,family=family,natoms=n,mode=mode,spacing=spacing,round=rnd,command=command,max_error=float(case[4]),neighbors=int(case[3]),samples=samples,reference_ms=1000*statistics.median(s[0] for s in samples),common_ms=1000*statistics.median(s[1] for s in samples),ratio=statistics.median(s[1]/s[0] for s in samples))
   rows.append(row); (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
   print(name, 'ratio',round(row['ratio'],3),'ms',round(row['common_ms'],4),flush=True)
(out/'serial-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
