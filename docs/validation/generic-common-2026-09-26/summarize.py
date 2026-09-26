from pathlib import Path
import json,statistics
out=Path('/tmp/accelnet-generic-common');rows=json.loads((out/'timings.json').read_text())
summary=[]
for compiler in ['gnu','nvhpc','h100','blackwell']:
 for n in [512,4096]:
  for family in ['lj','g4','g5','behler']:
   sub=[x for x in rows if x['compiler']==compiler and x['n']==n and x['family']==family]
   before=[x for x in sub if x['variant']=='before'];after=[x for x in sub if x['variant']=='after']
   if len(before)!=2 or len(after)!=2:continue
   med=lambda entries,key:statistics.median(x[key] for x in entries)
   r=dict(backend=compiler,n=n,family=family,before_ms=1000*med(before,'target'),after_ms=1000*med(after,'target'),
      after_before_reference_ratio=med(after,'ratio')/med(before,'ratio'),after_reference_ratio=med(after,'ratio'))
   r['speedup']=r['before_ms']/r['after_ms']
   for label,items in [('before',before),('after',after)]:
    r[label+'_phase_ms']=[1000*statistics.median(statistics.median(s[k] for s in x['phases']) for x in items) for k in range(8)]
   summary.append(r)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
for r in summary:
 if r['backend'] in ['gnu','nvhpc']:
  print(f"| {r['backend']} | {r['n']} | {r['family']} | {r['after_before_reference_ratio']:.3f} | {r['after_reference_ratio']:.3f} |")
 else:
  print(f"| {r['backend']} | {r['n']} | {r['family']} | {r['before_ms']:.4f} | {r['after_ms']:.4f} | {r['speedup']:.2f}× |")
