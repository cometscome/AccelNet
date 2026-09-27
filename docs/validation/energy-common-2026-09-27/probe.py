import json,subprocess,statistics,os,sys,hashlib
from pathlib import Path
os.sched_setaffinity(0,{6});u=Path('/tmp/accelnet-energy-common-opt');tag=sys.argv[1]
bins={'legacy':'/tmp/accelnet-unified-api/before-public-api','before':str(u/'before-public'),'after':'/tmp/accelnet-n2p2-cpu-opt/build-gnu/bin/accelnet-public-api-benchmark'}
report={'binaries':{k:{'path':v,'sha256':hashlib.sha256(Path(v).read_bytes()).hexdigest()} for k,v in bins.items()},'rows':[]}
for api in ['atomic-energy','energy','atomic','batch']:
 samples={k:[] for k in bins};reference=None
 for r in range(4):
  for label in list(bins)[::1 if r%2==0 else -1]:
   cmd=[bins[label],'aenet','/home/nagai/AccelNetGPU/AccelNetTrainer.jl/test/data/fortran_predict','192','.15',api]
   result=subprocess.run(cmd,capture_output=True,text=True,check=True)
   (u/f'{tag}-{api}-{label}-{r}.log').write_text(result.stdout+result.stderr)
   rows=result.stdout.splitlines();samples[label].append(float(next(l for l in rows if l.startswith('TIMING ')).split()[1]))
   vals=[float(x) for l in rows if l.startswith(('ENERGY ','FORCE ','VIRIAL ')) for x in l.split()[1:]]
   if reference is None:reference=vals
   assert len(vals)==len(reference) and all(abs(x-y)<=2e-10*(1+abs(y)) for x,y in zip(vals,reference)),(api,label)
 med={k:statistics.median(v) for k,v in samples.items()};entry={'api':api,'samples':samples,'seconds':med,'after_over_legacy':med['after']/med['legacy']};report['rows'].append(entry);print(api,med,entry['after_over_legacy'],flush=True)
 (u/f'{tag}.json').write_text(json.dumps(report,indent=2)+'\n')
