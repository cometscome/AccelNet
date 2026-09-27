from pathlib import Path
import subprocess,re,json,hashlib
out=Path('/tmp/accelnet-g4-improve');rows=[]
for label in ['gnu','nvhpc']:
 build=Path(f'/tmp/accelnet-other-descriptors/build-{label}-serial')
 for tag,exe in [('before',out/'before'/label),('after',build/'bin/accelnet-target-benchmark')]:
  symbols=subprocess.check_output(['nm','-u',str(exe)],text=True)
  omp=[s for s in symbols.splitlines() if re.search(r'\b(?:GOMP_|__kmpc_|__nvomp_|omp_)[A-Za-z0-9_]*',s)]
  used=[line for name in ['AccelNetPredictor','AccelNetTarget','accelnet-target-benchmark'] for line in (build/f'AccelNetPredictor/CMakeFiles/{name}.dir/flags.make').read_text().splitlines() if line.startswith('Fortran_FLAGS')]
  assert not omp and not any(re.search(r'(?:^|\s)(?:-fopenmp|-mp)(?:=\S+)?(?:\s|$)',line) for line in used)
  assert 'ACCELNET_TARGET_SERIAL:BOOL=ON' in (build/'CMakeCache.txt').read_text()
  rows.append(dict(compiler=label,variant=tag,path=str(exe),sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),target_serial=True,flags=used,undefined_openmp_symbols=omp,cpu_affinity=[6]))
(out/'serial-audit.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS: OpenMP-disabled builds and no direct OpenMP runtime references')
