from pathlib import Path
import subprocess,sqlite3,re,json,statistics
import numpy as np
base=Path('/tmp/accelnet-modes-profile')
report=[]
labels={41:'geometry_radial_direct_or_powers',126:'moment_build',145:'moment_to_descriptor',179:'network_and_moment_adjoint',236:'edge_force',312:'force_virial_scatter'}
def snapshot(p):
 l=p.read_text().splitlines();start=next(i for i,x in enumerate(l) if x.startswith('ITEM: ATOMS'))+1
 return np.array([list(map(float,x.split())) for x in l[start:]])
for atoms in [192,5184,24000]:
 for mode in ['direct','moment']:
  d=base/f'{atoms}-{mode}'
  with (d/'kernel-summary.csv').open('w') as f:
   subprocess.run(['nsys','stats','--report','cuda_gpu_kern_sum','--format','csv',str(d/'trace.nsys-rep')],stdout=f,check=True)
  c=sqlite3.connect(d/'trace.sqlite')
  rows=c.execute('select s.value,k.start,k.end,k.gridX,k.blockX,k.registersPerThread,k.localMemoryPerThread from CUPTI_ACTIVITY_KIND_KERNEL k join StringIds s on k.demangledName=s.id order by k.start').fetchall()
  groups={}
  for row in rows:groups.setdefault(row[0],[]).append(row)
  phases={}
  for name,events in groups.items():
   if 'run_target_batch' not in name and 'import_lammps_inputs' not in name:continue
   assert len(events)==123,(name,len(events))
   # run 20: initial evaluation + 20 steps (indices 0..20).
   # run 100: initial evaluation (21) + measured steps (22..121).
   # final run 0: index 122. Exclude all initial evaluations/warmup.
   selected=events[22:122]
   line=int(re.search('F1L(\d+)_',name)[1])
   label=labels[line] if 'run_target_batch' in name else f'import_{line}'
   ms=[(e[2]-e[1])/1e6 for e in selected]
   phases[label]=dict(mean_ms=statistics.mean(ms),median_ms=statistics.median(ms),min_ms=min(ms),max_ms=max(ms),kernel=name,steps=len(ms),gridX=selected[0][3],blockX=selected[0][4],registersPerThread=selected[0][5],localMemoryPerThread=selected[0][6])
  old=snapshot(Path(f'/tmp/accelnet-lammps-modes-h100/{atoms}-{mode}-0/final.dump'))
  got=snapshot(d/'final.dump');np.testing.assert_allclose(got,old,atol=2e-8,rtol=2e-9)
  stdout=(d/'stdout.txt').read_text()
  loops=re.findall(r'Loop time of ([\d.eE+-]+)',stdout)
  item=dict(atoms=atoms,mode=mode,profiled_loop_ms_per_step=float(loops[1])*10,phases=phases,max_force_difference_from_unprofiled=float(np.max(np.abs(old[:,5:8]-got[:,5:8]))))
  report.append(item)
  print(atoms,mode,'loop',round(item['profiled_loop_ms_per_step'],3),' '.join(f'{k}={v["mean_ms"]:.3f}' for k,v in phases.items()))
(base/'report.json').write_text(json.dumps(report,indent=2)+'\n')
