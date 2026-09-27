from pathlib import Path
import json,shutil,statistics,hashlib
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');temp=Path('/tmp/accelnet-g4-improve');out=root/'docs/validation/g4-common-2026-09-27';out.mkdir(exist_ok=True)
rows=[]
for label in ['gnu','nvhpc','h100','blackwell']:
 rows+=json.loads((temp/f'final-{label}.json').read_text())
assert len(rows)==40
(out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
for p in temp.glob('*.json'):shutil.copy2(p,out/p.name)
for p in temp.glob('*.log'):shutil.copy2(p,out/p.name)
for p in temp.glob('*.py'):shutil.copy2(p,out/p.name)
shutil.copytree(temp/'checks',out/'checks',dirs_exist_ok=True,ignore=shutil.ignore_patterns('lammps'))
for f in ['lj','behler','n2p2']:
 p=out/'checks/lammps'/f;p.mkdir(parents=True,exist_ok=True);shutil.copy2(temp/f'checks/lammps/{f}/report.json',p/'report.json')
lookup={(r['backend'],r['n'],r['family'],r['variant']):r for r in rows}
def table(labels):
 text='| Backend | Atoms | Model | Before common ms | After common ms | Speedup | After / legacy CPU time |\n|---|---:|---|---:|---:|---:|---:|\n'
 for label in labels:
  for n in [512,4096]:
   for family in ['g4','g4-series','g4-distinct']:
    a=lookup[(label,n,family,'after')];b=lookup.get((label,n,family,'before'))
    text+=f"| {label} | {n} | {family} | "+(f"{b['ms']:.3f}" if b else '—')+f" | {a['ms']:.3f} | "+(f"{b['ms']/a['ms']:.3f}x" if b else '—')+' | '+(f"{a['ratio']:.3f}" if label in ['gnu','nvhpc'] else '—')+' |\n'
 return text
s='''# G4 common-kernel follow-up — 2026-09-27

This is document revision **1.3**, AccelNet **1.0.1**, base Git commit
`c6631460a1bbb990c82e3e0ff5e73c36c52f6f9b` plus the uncommitted GPU/common-kernel
work. The immediate baseline is the adopted revision 1.2 implementation, including
NN-gradient-slot reuse, not its intermediate coefficient-array experiment.

## Implementation and cause

The two production files changed are `accelnet_target_kernels.f90` and
`accelnet_target_math.f90`. Both serial CPU and OpenMP GPU variants compile these
same numerical sources. No extra workspace, GPU launch, transfer, or CPU parallel
region was added. See [the math](../../../speedupmethods.md#14-reducing-g4-pair-work-in-the-common-kernel-revision-13),
[production diff](implementation.diff), [before sources](before-src/), and
[source hashes](source-sha256.json).

1. The single-power value kernel previously recomputed identical G4/G5 descriptor
   values despite combining their NN-gradient coefficients in the force stage.
   Representatives now compute once and broadcast unnormalized values. Separate
   normalization and NN weights remain intact; this is not NN-input removal.
2. G4 now rejects squared j-k distances outside the valid interval before sqrt,
   radial/exponential evaluation, or angular powers. The old value helper evaluated
   an angular power even for an already zero radial cutoff.
3. Pair helper integer arguments now have the packer's explicit extent (14).
   Assumed-shape array descriptor setup was unnecessary inside these hot loops.
4. Pair forces contract scalar coefficients before forming Cartesian components,
   reusing the j-k term and placing divisions outside vector expressions.

The last change contributes only approximately one percent in the three-round
alternating GNU probe (`ablation.json`); it does not explain the entire gain.
Initial incremental experiments are archived for diagnosis, while the final
before/after tables below use matched binaries and the final implementation.
The exploratory `before-now` run overlapped a build and is not used in these tables.
The first `duplicates-tests` invocation omitted `--host` and only reported GPU
unavailability; validation claims below refer to the final completed test suite.

The residual architectural difference remains: the established G4 CPU path
computes values and a descriptor Jacobian in one pair traversal, then contracts
that Jacobian. The common path computes values and contracted forces in separate
traversals, avoiding a global Jacobian but repeating pair geometry and nonlinear
functions. This tradeoff is visible especially for a small, distinct basis. The
previous Jacobian experiment was not adopted because of its measured GPU costs;
this revision does not introduce another CPU-only implementation.

## Measurement conditions

- GNU Fortran 11.4.0 (`-O3`) and NVHPC 25.3-0 (`-fast -O3`). Both CPU builds set
  `ACCELNET_TARGET_SERIAL=ON`, without `-fopenmp` or `-mp`. CPU affinity is core 6
  on Xeon Gold 6526Y. `serial-audit.json` records flags and confirms no direct
  OpenMP runtime references in either CPU benchmark. Merely setting one OpenMP
  thread is not the basis for the single-core claim.
- GPU: NVHPC 25.3, CUDA 12.8, driver 590.48.01, FP64, mandatory offload,
  `-mp=gpu -gpu=cc90,cc120`. Devices: H100 NVL (GPU 0, UUID ending
  `59e1f3c6e7f3`) and RTX PRO 6000 Blackwell Max-Q (GPU 2, UUID ending
  `139cedd981c9`). The other H100 was not used.
- Five timed samples per invocation, alternating legacy/reference and common
  execution order, with warmup and a 0.10 s minimum measurement interval.
  No benchmark in the final tables overlapped our builds or other tests.
- Atoms 512/4096, spacing 1.7, direct mode, two species, neighbor list fixed and
  built outside timing. Synchronous common timing includes transfers, descriptor
  evaluation, NN, forces, and virial. These are batch timings, not full MD rates.
- `g4`: four inputs, including two identical mixed-species angular descriptors.
  `g4-series`: 48 **distinct** inputs, integer powers 1–8, lambda ±1, three
  unordered species pairs, sharing radial parameters. `g4-distinct`: a new
  four-input control with the second mixed-species term's eta changed from 0.12
  to 0.18, so all four inputs differ. The saved before executable does not contain
  this new fixture; its before time is deliberately absent.
- CPU-reference code and compiler flags did not change. Absolute wall times can
  vary with CPU frequency/system conditions. Each raw log also contains its
  contemporaneous legacy CPU time and the eight stage measurements. CPU comparison
  ratios below are median per-sample ratios, not ratios of separately rounded medians.

## Final CPU measurements

Speedup means previous common / new common time. The last column compares the
new common kernel to the **established CPU algorithm**; below 1 is faster. These
are different comparisons. OpenMP is disabled for both algorithms.

'''+table(['gnu','nvhpc'])+'''
## Final GPU measurements

Speedup means previous common / new common GPU time, not GPU / single-core CPU.

'''+table(['h100','blackwell'])+'''
## Validation and adoption

- GNU serial: 10 batch/host tests passed; NVHPC serial: 10 passed.
- H100: 31 GPU/host tests passed; Blackwell: 29 GPU tests passed.
- Checked GNU build (`-fcheck=all -fbacktrace`): 9 batch/host tests passed.
- Descriptor suite: **130 cases**, including 10 cutoff types, G4 cutoff boundaries
  from both sides, integer/fractional/high/near-integer powers, reversed species,
  grouped and distinct radial parameters, collinear neighbors, finite differences,
  and workspace reuse. H100 memcheck: **0 errors**, maximum absolute E/F/W error
  about **7.77e-16** in that run. Agreement is within tolerances, not bitwise identity.
- LAMMPS: **63 comparisons** (21 each for LJ, Behler, n2p2), including neighbor
  modes, 1/2 MPI ranks, orthogonal/triclinic cells, empty partitions, and NVE.
- Ordinary CPU: 42 correctness tests and 2 performance checks passed.
- Raw checks and the LAMMPS reports are in `checks/`. Final G4 benchmark initial
  comparisons have maximum absolute E/F/W difference '''+f"{max(r['max_error'] for r in rows):.4g}"+'''; each timed sample also
  checks the E/F/W tolerance and verifies persistent allocation/upload counts.

The common kernels adopt these changes. **Default CPU LJ/Behler dispatch remains
unchanged.** A speedup relative to the previous common implementation does not
establish universal superiority over the legacy CPU path, particularly for the
four-distinct-input control. Chebyshev's existing default common dispatch is unchanged.

The validated LAMMPS binary was copied to
`/home/nagai/AccelNetGPU/lammps-accelnet-gpu/lmp` with the prior version retained as
`lmp-before-g4-pair-work`. See `deployment.json` for both SHA256 identities.

## Reproduction

`bench.py`, `build.py`, `checks.py`, and `audit.py` record this machine's exact
paths and commands. `binaries.json` records saved-before and final benchmark
hashes. Rebuild the before variant by restoring the two `before-src/` files in
an isolated copy of the full working tree and using the same compiler/options;
the base Git commit alone omits earlier GPU work and is not this baseline.

For example, with OpenMP-disabled before/after CPU binaries:

```sh
python3 AccelNetPredictor/benchmark/compare_chebyshev_variants.py \\
  --variant before host /path/to/before \\
  --variant after host /path/to/after \\
  --family g4 --sizes 512 4096 --orders 8 --modes 1 \\
  --rounds 3 --seconds 0.10 --cpu 6 --output /tmp/g4-recheck
```

For a GPU comparison select its UUID with `CUDA_VISIBLE_DEVICES`, use mandatory
OpenMP offload, and replace the two `host` backend arguments by `gpu`. Run timing
separately from builds/tests. `g4-distinct` requires binaries rebuilt with the new
benchmark fixture; archived revision-1.2 binaries cannot evaluate that fixture.
'''
(out/'README.md').write_text(s)
print('Wrote report',len(rows),'timing runs')

# Archive the independent distinct-basis control and the near-flat LJ recheck.
for name in ['control','lj-repeat']:
 shutil.copytree(temp/name,out/name,dirs_exist_ok=True)
shutil.copy2(temp/'control-before/build.json',out/'control-before-build.json')
controls=json.loads((temp/'control/report.json').read_text())
extra='''\n## Distinct-basis follow-up and other descriptors\n\nFor an actual before/after comparison of `g4-distinct`, the new benchmark fixture\nwas also linked against the frozen revision-1.2 GNU libraries in\n`/tmp/accelnet-angular-contraction/build-reuse-gnu`. Their numerical source files\nmatch the saved before sources byte for byte. `control-before-build.json` records\nthe command, library hashes, and binary identity. Only the fixture/driver was\nrecompiled; neither old numerical library nor the legacy CPU algorithm was edited.\n\nTwo rounds reverse the executable order (before/after, then after/before), with\nfive alternating legacy/common samples per invocation. These are OpenMP-disabled\nGNU runs pinned to CPU 6. The rows aggregate the two per-round medians.\n\n| Atoms | Old common ms | New common ms | Speedup | Old common / legacy | New common / legacy |\n|---|---:|---:|---:|---:|---:|\n'''
for n in [512,4096]:
 groups={label:[x for x in controls if x['natoms']==n and x['label']==label] for label in ['before','after']}
 times={label:1000*statistics.median(x['seconds'] for x in a) for label,a in groups.items()}
 ratios={label:statistics.median(x['ratio_to_reference'] for x in a) for label,a in groups.items()}
 extra+=f"| {n} | {times['before']:.3f} | {times['after']:.3f} | {times['before']/times['after']:.3f}x | {ratios['before']:.3f} | {ratios['after']:.3f} |\n"
extra+='''
Thus the common G4 improvement is not confined to duplicate inputs. The small
unique basis still favors legacy CPU substantially. Its GNU ratio to legacy in
this follow-up is higher than in the first table (about 1.58 vs 1.44 at 4096),
while common absolute time is slightly lower. Reference wall time changed between
campaigns; both results are retained rather than interpreting the difference as
an algorithm change. No hardware-counter claim about the cause of this timing
variation is made. `control/report.json` contains all samples and phase timings.

The shared-helper impact on G5, mixed Behler, and LJ was checked at 512 atoms:

| Backend | Model | Old common ms | New common ms | Speedup |\n|---|---|---:|---:|---:|\n'''
for label in ['gnu','nvhpc','h100','blackwell']:
 regress=json.loads((temp/f'regression-{label}.json').read_text())
 for family in ['g5','behler','lj']:
  a=next(x for x in regress if x['family']==family and x['variant']=='after')
  b=next(x for x in regress if x['family']==family and x['variant']=='before')
  extra+=f"| {label} | {family} | {b['ms']:.3f} | {a['ms']:.3f} | {b['ms']/a['ms']:.3f}x |\n"
lj=json.loads((temp/'lj-repeat/report.json').read_text())
t={label:statistics.median(x['seconds'] for x in lj if x['label']==label)*1000 for label in ['before','after']}
extra+=f'''
The initial GNU LJ sample suggested a small slowdown, so it was repeated in
three alternating-order rounds. Median common time was {t['before']:.4f} ms before
and {t['after']:.4f} ms after ({(t['after']/t['before']-1)*100:.2f}% longer), with overlapping
per-round time ranges. This is approximately flat, not a claimed LJ speedup;
`lj-repeat/report.json` preserves the samples. NVHPC CPU and both GPUs were also
approximately flat for LJ. The ordinary default CPU performance gate passed;
this revision does not route default LJ/Behler CPU evaluation to the common kernel.
'''
p=out/'README.md';p.write_text(p.read_text()+extra)
