from pathlib import Path
import subprocess,sqlite3,re,json,statistics
import numpy as np
base=Path('/tmp/accelnet-moment-analysis')
labels=['geometry_radial_powers','moment_build','moment_to_descriptor','network','edge_force','force_virial_scatter']
def snapshot(p):
 l=p.read_text().splitlines();start=next(i for i,x in enumerate(l) if x.startswith('ITEM: ATOMS'))+1
 return np.array([list(map(float,x.split())) for x in l[start:]])
report=[]
for d in sorted((base/'profile').iterdir()):
 if not (d/'trace.nsys-rep').exists():continue
 atoms,variant=d.name.split('-');atoms=int(atoms)
 with (d/'kernel-summary.csv').open('w') as f:subprocess.run(['nsys','stats','--report','cuda_gpu_kern_sum','--format','csv',str(d/'trace.nsys-rep')],stdout=f,check=True)
 c=sqlite3.connect(d/'trace.sqlite')
 rows=c.execute('select s.value,k.start,k.end,k.gridX,k.blockX,k.registersPerThread,k.localMemoryPerThread from CUPTI_ACTIVITY_KIND_KERNEL k join StringIds s on k.demangledName=s.id order by k.start').fetchall()
 groups={}
 for row in rows:
  if 'run_target_batch' in row[0]:groups.setdefault(row[0],[]).append(row)
 assert len(groups)==6
 phases={}
 for label,(name,events) in zip(labels,groups.items()):
  assert len(events)==123
  selected=events[22:122];ms=[(e[2]-e[1])/1e6 for e in selected]
  phases[label]=dict(mean_ms=statistics.mean(ms),median_ms=statistics.median(ms),min_ms=min(ms),max_ms=max(ms),kernel=name,gridX=selected[0][3],blockX=selected[0][4],registersPerThread=selected[0][5],localMemoryPerThread=selected[0][6])
 old=snapshot(Path(f'/tmp/accelnet-lammps-modes-h100/{atoms}-moment-0/final.dump'));got=snapshot(d/'final.dump');np.testing.assert_allclose(got,old,atol=2e-8,rtol=2e-9)
 loop=float(re.findall(r'Loop time of ([\d.eE+-]+)',(d/'stdout.txt').read_text())[1])*10
 item=dict(atoms=atoms,variant=variant,phases=phases,loop_ms=loop,max_force_error=float(np.max(np.abs(old[:,5:8]-got[:,5:8]))));report.append(item)
 print(atoms,variant,'loop',round(loop,3),' '.join(f'{k}={v["mean_ms"]:.3f}' for k,v in phases.items()))
(base/'profile/report.json').write_text(json.dumps(report,indent=2)+'\n')
