from pathlib import Path
import json,statistics,hashlib,shutil
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-common-default'); dest=root/'docs/validation/common-cpu-default-2026-09-27';dest.mkdir(exist_ok=True)
rows=json.loads((out/'results.json').read_text()); groups={}
for r in rows:groups.setdefault((r['compiler'],r['family'],r['natoms'],r['mode'],r['spacing']),[]).append(r)
summary=[]
for key,rs in groups.items():
 assert len(rs)==2
 compiler,family,n,mode,spacing=key
 summary.append(dict(compiler=compiler,family=family,natoms=n,mode=mode,spacing=spacing,reference_ms=statistics.median(r['reference_ms'] for r in rs),common_ms=statistics.median(r['common_ms'] for r in rs),ratio=statistics.median(r['ratio'] for r in rs)))
(dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
for pattern in ['*.log','*.json','*.py']:
 for p in out.glob(pattern):shutil.copy2(p,dest/p.name)
shutil.copytree(out/'chebyshev-control',dest/'chebyshev-control',dirs_exist_ok=True)
paths=['AccelNetPredictor/src/accelnet_batch.f90','AccelNetPredictor/src/accelnet_target_kernels.f90','AccelNetPredictor/src/accelnet_batch_target.f90','AccelNetPredictor/test/test_batch.f90','AccelNetPredictor/test/test_batch_target.f90','AccelNetPredictor/benchmark/target_benchmark.f90','AccelNetPredictor/CMakeLists.txt']
(dest/'source-sha256.json').write_text(json.dumps({p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths},indent=2)+'\n')
print('new family ratio ranges')
for c in ['gnu','nvhpc']:
 for family in sorted({r['family'] for r in summary}):
  rs=[r for r in summary if r['compiler']==c and r['family']==family]
  print(c,family,round(min(r['ratio'] for r in rs),3),round(max(r['ratio'] for r in rs),3))
print('max error',max(r['max_error'] for r in rows))
print('rows',len(rows))
