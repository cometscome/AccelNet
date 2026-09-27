import json,os,subprocess,statistics
from pathlib import Path
out=Path('/tmp/accelnet-unified-api/gpu-performance-final');out.mkdir(exist_ok=True)
exes={'before':'/tmp/accelnet-unified-api/baseline-gpu/bin/accelnet-target-benchmark','after':'/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark'}
os.sched_setaffinity(0,{6})
report={'rows':[]}
for fam,n,o in [('chebyshev',2048,8),('g5-series',2048,4),('g5-high',256,16)]:
 for mode in [1,2 if fam=='chebyshev' else 3]:
  samples={k:[] for k in exes};runs=[]
  for r in range(2):
   for key in (['before','after'] if r%2==0 else ['after','before']):
    cmd=[exes[key],str(n),str(o),'.05',str(mode),'1.7','gpu',fam,'no-neighbors','64','candidate-only']
    x=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=300)
    (out/f'{fam}-{mode}-{r}-{key}.log').write_text(x.stdout+x.stderr)
    values=[float(l.split()[3]) for l in x.stdout.splitlines() if l.startswith('TIMING ')]
    assert len(values)==5
    samples[key]+=values;runs.append({'command':cmd,'errors':[l for l in x.stdout.splitlines() if l.startswith('CASE ')]})
  secs={k:statistics.median(v) for k,v in samples.items()};row=dict(family=fam,atoms=n,order=o,mode=mode,samples=samples,seconds=secs,ratio=secs['after']/secs['before'],runs=runs)
  report['rows'].append(row);print(fam,mode,secs,row['ratio'],flush=True)
  (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
report['passed']=True
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
