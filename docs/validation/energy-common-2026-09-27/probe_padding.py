import os,subprocess,statistics,json
from pathlib import Path
os.sched_setaffinity(0,{6});u=Path('/tmp/accelnet-energy-common-opt');r=[]
bins={'legacy':'/tmp/accelnet-unified-api/before-public-api','unpad':str(u/'before-padding-public'),'padded':'/tmp/accelnet-n2p2-cpu-opt/build-gnu/bin/accelnet-public-api-benchmark'}
for n in (192,512):
 for api in ('energy','atomic-energy','batch'):
  row={'n':n,'api':api,'samples':{k:[] for k in bins}};ref=None
  for i in range(4):
   for label in list(bins)[::1 if i%2==0 else -1]:
    p=subprocess.run([bins[label],'aenet','/home/nagai/AccelNetGPU/AccelNetTrainer.jl/test/data/fortran_predict',str(n),'.2',api],capture_output=True,text=True,check=True)
    (u/f'padding-{n}-{api}-{label}-{i}.log').write_text(p.stdout+p.stderr)
    row['samples'][label].append(float(next(l for l in p.stdout.splitlines() if l.startswith('TIMING ')).split()[1]))
    e=float(next(l for l in p.stdout.splitlines() if l.startswith('ENERGY ')).split()[1]);ref=e if ref is None else ref
    assert abs(e-ref)<=2e-10*(1+abs(ref))
  row['seconds']={k:statistics.median(v) for k,v in row['samples'].items()};r.append(row);print(n,api,row['seconds'],flush=True)
  (u/'padding.json').write_text(json.dumps(r,indent=2)+'\n')
