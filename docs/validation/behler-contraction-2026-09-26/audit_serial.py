from pathlib import Path
import subprocess,re,json,hashlib
out=Path('/tmp/accelnet-angular-contraction');rows=[]
for label in ['gnu','nvhpc']:
 build=Path(f'/tmp/accelnet-other-descriptors/build-{label}-serial')
 for variant,exe in [('before',out/f'baseline-{label}/accelnet-target-benchmark'),('after',build/'bin/accelnet-target-benchmark')]:
  symbols=subprocess.check_output(['nm','-u',str(exe)],text=True)
  omp=[s for s in symbols.splitlines() if re.search(r'\b(?:GOMP_|__kmpc_|__nvomp_|omp_)[A-Za-z0-9_]*',s)]
  flags={name:(build/f'AccelNetPredictor/CMakeFiles/{name}.dir/flags.make').read_text() for name in ['AccelNetPredictor','AccelNetTarget','accelnet-target-benchmark']}
  used=[line for text in flags.values() for line in text.splitlines() if line.startswith('Fortran_FLAGS')]
  bad=[line for line in used if re.search(r'(?:^|\s)(?:-fopenmp|-mp)(?:=\S+)?(?:\s|$)',line)]
  assert not omp and not bad,(label,omp,bad)
  cache=(build/'CMakeCache.txt').read_text();assert 'ACCELNET_TARGET_SERIAL:BOOL=ON' in cache
  row=dict(compiler=label,variant=variant,path=str(exe),sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),target_serial=True,flags=used,undefined_openmp_symbols=omp,ldd=subprocess.check_output(['ldd',str(exe)],text=True),cpu_affinity_in_driver=[6])
  rows.append(row)
(out/'serial-audit.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS: both compilers, both candidates, OpenMP compilation disabled and no direct OpenMP runtime calls')
