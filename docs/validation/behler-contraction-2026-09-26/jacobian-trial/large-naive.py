from pathlib import Path
import json,os,subprocess,statistics,hashlib
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-jacobian-trial');os.sched_setaffinity(0,{6})
rows=[]
for backend,label,uuid in [('gpu','h100','GPU-2644154d-7268-af42-6631-59e1f3c6e7f3'),('gpu','blackwell','GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9')]:
 env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES=uuid)
 if backend=='host':
  variants=[('common',Path(f'/tmp/accelnet-other-descriptors/build-{label}-serial/bin/accelnet-target-benchmark')),('jacobian',out/f'build-{label}/bin/accelnet-target-benchmark')]
 else:
  variants=[('jacobian',out/'build-gpu/bin/accelnet-target-benchmark')]
 for family in ['g4','g5','behler']:
  for mode,exe in variants:
   name=f'large-{label}-{family}-{mode}'
   cmd=[str(exe),'4096','5','0.08','1','1.7',backend,family,'no-neighbors']
   p=subprocess.run(cmd,env=env,text=True,capture_output=True)
   (out/(name+'.log')).write_text(p.stdout+p.stderr)
   if p.returncode:raise SystemExit(f'FAIL {name}\n'+p.stdout[-2000:]+p.stderr[-2000:])
   samples=[list(map(float,l.split()[2:])) for l in p.stdout.splitlines() if l.startswith('TIMING ')]
   row=dict(backend=label,family=family,variant=mode,command=cmd,sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),samples=samples,
            ms=1000*statistics.median(s[1] for s in samples),ratio=statistics.median(s[1]/s[0] for s in samples))
   rows.append(row);(out/'large-naive.json').write_text(json.dumps(rows,indent=2)+'\n')
   print(name,'ms',round(row['ms'],4),'ratio',round(row['ratio'],3),flush=True)
