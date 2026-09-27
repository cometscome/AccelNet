import os,subprocess,statistics,json,hashlib
from pathlib import Path
os.sched_setaffinity(0,{6});u=Path('/tmp/accelnet-energy-common-opt');out=u/'gpu-final2';out.mkdir(exist_ok=True)
current=Path('/tmp/accelnet-gpu-research/build-nvhpc/bin')
model='/home/nagai/AccelNetGPU/AccelNetTrainer.jl/test/data/fortran_predict'
report={'gpu':os.environ['CUDA_VISIBLE_DEVICES'],'rows':[]}
for api,kind,n,mode in [('energy','chebyshev',2048,m) for m in (1,2)]+[('energy','aenet',192,m) for m in (0,1,2)]+[('force','chebyshev',2048,m) for m in (1,2)]+[('force','g5-series',1024,m) for m in (1,2)]:
 bins={'before':u/('before-gpu-energy' if api=='energy' else 'before-gpu'),'after':current/('accelnet-target-energy-benchmark' if api=='energy' else 'accelnet-target-benchmark')}
 row={'api':api,'model':kind,'atoms':n,'mode':mode,'samples':{k:[] for k in bins},'runs':[],'binaries':{k:{'path':str(v),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in bins.items()}}
 for r in range(4):
  for label in list(bins)[::1 if r%2==0 else -1]:
   if api=='energy':args=[str(n),'8','.15',str(mode),model if kind=='aenet' else '-']
   else:args=[str(n),'8','.1',str(mode),'1.7','gpu',kind,'no-neighbors','64','candidate-only']
   cmd=[str(bins[label]),*args];p=subprocess.run(cmd,capture_output=True,text=True,timeout=300,check=True)
   (out/f'{api}-{kind}-{n}-{mode}-{r}-{label}.log').write_text(p.stdout+p.stderr)
   row['samples'][label]+=[float(l.split()[1 if api=='energy' else 3]) for l in p.stdout.splitlines() if l.startswith('TIMING ')]
   row['runs'].append({'command':cmd,'profile_or_errors':[l for l in p.stdout.splitlines() if l.startswith(('PROFILE ','MAX_ENERGY_ERROR ','CASE '))]})
 row['seconds']={k:statistics.median(v) for k,v in row['samples'].items()};row['ratio']=row['seconds']['after']/row['seconds']['before'];report['rows'].append(row)
 print(api,kind,n,mode,row['seconds'],row['ratio'],flush=True)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
