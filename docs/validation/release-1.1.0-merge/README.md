# AccelNet 1.1.0 merge validation

Date: 2026-09-28.

Parents: main `d20bb0d` (1.0.2 grouped G4 improvements) and gpu
`af87a96` (1.1.0 version update). The unpublished local annotated tag
`1.1.0` is moved to the resulting merge commit.

The main G4 changes and fused cutoff helper are retained in the standalone
reference evaluator, including the two new descriptor tests. The production
CPU/GPU shared numerical kernels and n2p2 extensions from gpu are retained.
All package versions remain 1.1.0. Earlier validation archives are unchanged.

## Correctness

- GNU Fortran 11.4 release build: 54/54 non-GPU, non-performance tests passed.
- GNU bounds-check build: 24/24 selected descriptor, host-target, public API,
  batch and virial tests passed.
- NVHPC 25.3, H100, mandatory device offload, one host thread: 35/35 GPU tests
  passed, including multi-element extended n2p2 and high-order G5 comparisons.

See the corresponding `*-tests.log` files. LAMMPS was not rerun for this merge.

## Performance and remaining limitation

CPU benchmarks compile OpenMP out and use one CPU core. The performance suite
passed 3/4 checks: common CPU kernels, G5 moment advantage, and the order-scaling
smoke check. The aenet energy regression gate failed against the archived
original CPU public API executable (`/tmp/accelnet-unified-api/before-public-api`).
No acceptance threshold was changed.

The first run measured a 1.1114 after/before ratio for 512-atom atomic-energy
calls, above the 1.10 limit. A longer run (8 alternating rounds, 0.5 seconds per
sample, core 6) reproduced this at 1.1101. The longer run also reported 1.1512
for 64-atom structure energy, but its samples span a substantial timing change
mid-run; that aggregate is not a stable estimate. Raw samples are retained in
`aenet-first.json` and `aenet-long.json`.

To distinguish the merge from the original CPU comparison, the same benchmark
driver and support module were compiled with GNU Fortran -O3 against the
pre-merge 1.1.0 installed modules and static libraries. A direct comparison
(6 alternating rounds, 0.3 seconds per sample, core 6) gave:

| Atoms | Structure energy, merged/pre-merge | Atomic energy, merged/pre-merge |
|---:|---:|---:|
| 64 | 1.0012 | 0.9698 |
| 192 | 1.0008 | 0.9726 |
| 512 | 1.0028 | 0.9782 |

See `aenet-vs-gpu.json` for executable hashes, paths, commands, samples and
numerical comparison results. These measurements show no meaningful slowdown
from this merge in the tested cases. The roughly 11% difference from the
original CPU implementation remains unresolved; this release must not be
described as passing every performance gate.
