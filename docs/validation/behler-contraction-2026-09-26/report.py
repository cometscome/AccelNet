from pathlib import Path
import json,statistics,shutil,hashlib,re
repo=Path('/home/nagai/AccelNetGPU/AccelNet-clone');tmp=Path('/tmp/accelnet-angular-contraction');out=repo/'docs/validation/behler-contraction-2026-09-26'
a=json.loads((tmp/'final-timings.json').read_text());b=json.loads((tmp/'adopted-timings.json').read_text())
assert len(a)==96 and len(b)==48
rows=[]
for r in b:
 old=next(x for x in a if x['variant']=='before' and (x['backend'],x['n'],x['family'])==(r['backend'],r['n'],r['family']))
 speedup=old['ratio']/r['ratio'] if r['backend'] in ['gnu','nvhpc'] else old['ms']/r['ms']
 rows.append(dict(backend=r['backend'],atoms=r['n'],family=r['family'],before_ms=old['ms'],after_ms=r['ms'],speedup=speedup,after_over_cpu_reference=r['ratio'],max_error=r['max_error']))
(out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
for name in ['adopted-timings.json','adopted-compare-progress.log','adopted-checks-progress.log','build-adopted-progress.log','serial-audit.json']:
 shutil.copy2(tmp/name,out/name)
for p in tmp.glob('adopted-*.log'):
 shutil.copy2(p,out/p.name)
checks=tmp/'adopted-checks'
(out/'checks').mkdir(exist_ok=True)
for p in checks.glob('*.log'):shutil.copy2(p,out/'checks'/p.name)
for p in checks.glob('chebyshev-*'):
 if p.is_dir():shutil.copytree(p,out/'checks'/p.name,dirs_exist_ok=True)
for p in checks.glob('lammps/*/report.json'):
 (out/'checks'/'lammps'/p.parent.name).mkdir(parents=True,exist_ok=True);shutil.copy2(p,out/'checks'/'lammps'/p.parent.name/p.name)
shutil.copy2(checks/'final-binaries.json',out/'final-binaries.json')
text='''\n## Adopted measurements\n\nThe tables use the final implementation with coefficients held in existing NN-gradient\nslots. CPU speedup is the before/after ratio of reference-normalized times; GPU\nspeedup is before/after absolute time. `Final / CPU reference` below one means the\ncommon serial implementation is faster than the established CPU implementation.\n\n### OpenMP-disabled single-core CPU\n\n| Compiler | Atoms | Model | Before (ms) | Final (ms) | Normalized common speedup | Final / CPU reference |\n|---|---:|---|---:|---:|---:|---:|\n'''
for r in rows:
 if r['backend'] in ['gnu','nvhpc']:text+=f"| {r['backend']} | {r['atoms']} | {r['family']} | {r['before_ms']:.3f} | {r['after_ms']:.3f} | {r['speedup']:.2f}x | {r['after_over_cpu_reference']:.3f} |\n"
text+='''\nThe independent reference can also vary across runs. For example, the GNU
4096-atom four-input G4 reference was slower in the final run than in the initial
campaign, so its normalized common speedup exceeds the ratio of raw common
times. Both are shown; this is not evidence that all of that difference comes
from this optimization. The per-run final/reference comparison remains explicit.

### GPU synchronous batches\n\n| GPU | Atoms | Model | Before (ms) | Final (ms) | Speedup |\n|---|---:|---|---:|---:|---:|\n'''
for r in rows:
 if r['backend'] not in ['gnu','nvhpc']:text+=f"| {r['backend']} | {r['atoms']} | {r['family']} | {r['before_ms']:.3f} | {r['after_ms']:.3f} | {r['speedup']:.2f}x |\n"
text+='''\n### Remaining CPU gap and interpretation\n\nThe four-output G4 cases still favor the established CPU implementation. Its\ndescriptor stage computes an unordered pair once, keeps both derivatives, and\nreuses radial/geometry factors across angular parameters. The shared no-Jacobian\nCPU path visits pairs in the value stage and again after NN backpropagation for\nforces. With few outputs, coefficient contraction saves too little work to fully\noffset that repeated geometry/radial evaluation. With a wider shared angular\nbasis, it avoids many output-specific gradient vectors and the common path closes\nthe G4 gap while substantially accelerating G5. These operation-count differences\nare visible in the source; the timings do not identify every compiler/cache cost.\n\nThe default non-Chebyshev CPU dispatch is therefore retained. This is a measured\nconditional performance result, not a general assertion that the CPU algorithm is\nunsuitable for GPUs or that shared Fortran requires CPU threading.\n\n## Validation\n\n- GNU and NVHPC OpenMP-disabled serial tests: 10 each.\n- H100: 31 GPU/host tests; Blackwell: 29 GPU tests.\n- GNU bounds/runtime-checked build: 9 tests.\n- Descriptor equivalence/finite-difference suite: **99 cases**, including sparse\n  powers, degree 16, noninteger and high powers, near-integer inputs, collinear\n  endpoint angles, species-map changes, radial-group split/merge, and coefficient\n  workspace reuse. This is contained in the suites above, not an additional test count.\n- Compute Sanitizer memcheck on H100: **0 errors** for those 99 cases (`--leak-check no`,\n  as in previous campaigns; this is not a device-pool leak claim).\n- LAMMPS: **63** LJ/Behler/n2p2 CPU/GPU comparisons, including 1/2 MPI ranks,\n  orthogonal/triclinic cells, GPU neighbor options, empty-rank coverage, and NVE.\n- Ordinary CPU: 42 numerical tests and 2 performance checks passed.\n- Chebyshev direct/moment before/after spot checks: both CPU compilers and H100.\n  An initial GNU 4096-atom direct result about 5% slower was not reproduced by\n  the follow-up with reversed executable order; see `chebyshev-repeat-gnu`.\n  Small single-campaign timing differences are not treated as universal changes.\n\nThe checks use the independent established CPU implementation, plus coordinate\nand strain finite differences. Agreement is numerical, not bitwise. Raw logs and\nphase times, including full-output checks and workspace-residency assertions, are\nretained. `source-sha256.json`, `implementation.diff`, `final-binaries.json`, and\n`baseline-library-sha256.json` identify the working-tree implementation and builds.\nThe Git base hash alone does not identify these uncommitted changes.\n\n## Reproduction\n\nThe preserved `before-src` files restore the revision-1.1 numeric backend in an\nisolated copy of this working tree; retain the current benchmark/model builder\nin both copies. Build CPU candidates with `ACCELNET_TARGET_SERIAL=ON` and without\nOpenMP compiler flags. The archived scripts record exact local paths and commands;\nuse separate build directories and substitute your own paths. For example:\n\n```sh\npython3 AccelNetPredictor/benchmark/compare_chebyshev_variants.py \\\n  --variant before host /path/to/before/accelnet-target-benchmark \\\n  --variant after host /path/to/after/accelnet-target-benchmark \\\n  --family g4-series --orders 8 --modes 1 --sizes 512 4096 \\\n  --rounds 2 --seconds 0.12 --cpu 6 --output comparison-g4-series\n```\n\nFor GPU measurements use `gpu` instead of `host` and select one GPU with\n`CUDA_VISIBLE_DEVICES`. Both candidates still verify all energy/force/virial\ncomponents against the established CPU reference. The old diagnostic source-edit\nexperiments are archived; use this binary-comparison driver for current kernels.\n'''
p=out/'README.md';s=p.read_text().split('\n## Adopted measurements')[0];p.write_text(s+text)
# Compact final measurements in the mathematical record.
p=repo/'speedupmethods.md';s=p.read_text().split('\n### 13.4 Final measurements')[0]
s+='''\n### 13.4 Final measurements\n\nThe table uses 4096 atoms and the 48-input integer-power models defined in the\n[report](docs/validation/behler-contraction-2026-09-26/README.md). CPU entries are\n**final common / established CPU** time with OpenMP compilation disabled (less\nthan one is faster). GPU entries are **previous common / final common** speedups;\nthey are not GPU-versus-CPU speedups or full MD measurements.\n\n| Model | GNU CPU ratio | NVHPC CPU ratio | H100 speedup | Blackwell speedup |\n|---|---:|---:|---:|---:|\n'''
for fam in ['g4-series','g5-series']:
 rs={r['backend']:r for r in rows if r['atoms']==4096 and r['family']==fam}
 s+=f"| {fam} | {rs['gnu']['after_over_cpu_reference']:.3f} | {rs['nvhpc']['after_over_cpu_reference']:.3f} | {rs['h100']['speedup']:.2f}x | {rs['blackwell']['speedup']:.2f}x |\n"
s+='''\nThe smaller four-input G4 models still favor the established CPU path; the report\nalso includes those cases, G5, mixed Behler, and LJ at 512 and 4096 atoms. Hence\nthese numbers support algebraic reuse for compatible bases, not a universal\nperformance guarantee or an unconditional CPU-dispatch replacement. The final\n99-case descriptor checks, device memcheck, and LAMMPS comparisons passed.\n''';p.write_text(s)
print(json.dumps([r for r in rows if r['atoms']==4096 and r['family'] in ['g4-series','g5-series']],indent=2))
