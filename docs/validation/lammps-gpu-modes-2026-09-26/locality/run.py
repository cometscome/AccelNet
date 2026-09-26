from pathlib import Path
import os,subprocess,re,json,argparse,statistics
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--replicates',type=int,nargs='+',default=[6,10]);p.add_argument('--variants',nargs='+',default=['baseline','tile1','tile32','tile128','tile512']);p.add_argument('--samples',type=int,default=2);p.add_argument('--output',required=True)
a=p.parse_args();base=Path('/tmp/accelnet-moment-analysis');out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3',OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',OPENBLAS_NUM_THREADS='1')
text=Path('/tmp/accelnet-lammps-modes-h100/24000-moment-0/in.benchmark').read_text()
report=[]
def snapshot(p):
 l=p.read_text().splitlines();i=next(i for i,x in enumerate(l) if x.startswith('ITEM: ATOMS'))+1
 return np.array([list(map(float,x.split())) for x in l[i:]])
def thermo(t):
 result=[]
 for line in t.splitlines():
  f=line.split()
  if len(f)==8 and re.fullmatch(r'\d+',f[0]):
   try:result.append(list(map(float,f)))
   except ValueError:pass
 return np.array(result)
for rep in a.replicates:
 reference=None;refthermo=None
 for sample in range(a.samples):
  order=a.variants if sample%2==0 else a.variants[::-1]
  for variant in order:
   d=out/f'{24*rep**3}-{variant}-{sample}';d.mkdir(exist_ok=True)
   (d/'in.benchmark').write_text(re.sub(r'replicate \d+ \d+ \d+',f'replicate {rep} {rep} {rep}',text))
   binary=Path('/home/nagai/AccelNetGPU/lammps-accelnet-gpu/lmp') if variant=='baseline' else base/variant/'lmp'
   with (d/'stdout.txt').open('w') as f:subprocess.run(['taskset','-c','6',str(binary),'-in','in.benchmark'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=240)
   log=(d/'stdout.txt').read_text();got=snapshot(d/'final.dump');th=thermo(log)
   assert np.isfinite(got).all() and np.isfinite(th).all()
   if reference is None:
    assert variant=='baseline';reference=got;refthermo=th
   np.testing.assert_allclose(got[:,:2],reference[:,:2],atol=0,rtol=0)
   np.testing.assert_allclose(got[:,2:],reference[:,2:],atol=2e-8,rtol=2e-9)
   np.testing.assert_allclose(th,refthermo,atol=2e-5,rtol=2e-9)
   loops=re.findall(r'Loop time of ([\d.eE+-]+)',log)
   item=dict(atoms=len(got),variant=variant,sample=sample,loop_ms=float(loops[1])*10,max_force_error=float(np.max(np.abs(got[:,5:8]-reference[:,5:8]))),max_atom_energy_error=float(np.max(np.abs(got[:,8]-reference[:,8]))),max_thermo_difference=float(np.max(np.abs(th-refthermo))))
   report.append(item);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(item),flush=True)
