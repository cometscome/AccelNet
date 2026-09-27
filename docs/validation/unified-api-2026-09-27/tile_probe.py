import os,subprocess,statistics,json,sys
from pathlib import Path
os.sched_setaffinity(0,{6})
u=Path('/tmp/accelnet-unified-api'); out={}
for kind,exes in [('atomic-energy',{'legacy':u/'before-public-api','fixed4':u/'public-fixed4','fixed8':Path('/tmp/accelnet-n2p2-cpu-opt/build-gnu/bin/accelnet-public-api-benchmark')}),('gpu-moment',{'fixed4':u/'gpu-fixed4','fixed8':Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark')})]:
 samples={k:[] for k in exes}
 for r in range(4):
  for key in list(exes)[::1 if r%2==0 else -1]:
   args=['aenet','/home/nagai/AccelNetGPU/AccelNetTrainer.jl/test/data/fortran_predict','192','.2','atomic-energy'] if kind=='atomic-energy' else ['2048','8','.05','2','1.7','gpu','chebyshev','no-neighbors','64','candidate-only']
   result=subprocess.run([str(exes[key])]+args,capture_output=True,text=True,check=True)
   (u/f'fixed8-{kind}-{key}-{r}.log').write_text(result.stdout+result.stderr)
   samples[key]+=[float(l.split()[1 if kind=='atomic-energy' else 3]) for l in result.stdout.splitlines() if l.startswith('TIMING ')]
 out[kind]={'samples':samples,'medians':{k:statistics.median(v) for k,v in samples.items()}}
 print(kind,out[kind]['medians'],flush=True)
(u/'fixed8-probe.json').write_text(json.dumps(out,indent=2)+'\n')
