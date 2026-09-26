from pathlib import Path
import os, subprocess, statistics, json
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-generic-common')
os.sched_setaffinity(0,{6})
results=[]
for backend,compiler,uuid in [('host','gnu',''),('host','nvhpc',''),('gpu','h100','GPU-2644154d-7268-af42-6631-59e1f3c6e7f3'),('gpu','blackwell','GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9')]:
 before=out/('before-gpu' if backend=='gpu' else 'before-'+compiler)/'accelnet-target-benchmark'
 after=Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark') if backend=='gpu' else Path(f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/accelnet-target-benchmark')
 for n in [512,4096]:
  for family in ['lj','g4','g5','behler']:
   for round in range(2):
    for name,exe in ([('before',before),('after',after)] if round==0 else [('after',after),('before',before)]):
     label=f'{compiler}-{n}-{family}-{round}-{name}'
     cmd=[str(exe),str(n),'5','0.08','1','1.7',backend,family,'no-neighbors']
     env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES=uuid)
     p=subprocess.run(cmd,cwd=root,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
     (out/(label+'.log')).write_text(p.stdout)
     if p.returncode: raise SystemExit('FAIL '+label+'\n'+p.stdout[-3000:])
     samples=[list(map(float,l.split()[2:])) for l in p.stdout.splitlines() if l.startswith('TIMING ')]
     phases=[list(map(float,l.split()[3:])) for l in p.stdout.splitlines() if l.startswith('PROFILE ')]
     assert len(samples)==5
     r=dict(compiler=compiler,n=n,family=family,round=round,variant=name,cmd=cmd,samples=samples,phases=phases,
            target=statistics.median(x[1] for x in samples),ratio=statistics.median(x[1]/x[0] for x in samples))
     results.append(r);(out/'timings.json').write_text(json.dumps(results,indent=2)+'\n')
     print(label,'ms',f"{r['target']*1000:.4f}",'ratio',f"{r['ratio']:.3f}",flush=True)
