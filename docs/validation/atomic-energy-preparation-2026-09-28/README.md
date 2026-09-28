# Atomic energy preparation: AccelNet 1.1.0 / methods 1.14

Measured 2026-09-28 (JST). Starting revision: main `c38c563`, tagged `1.1.0`.
Independent original-CPU baseline: `543b180`, using the identical public API
driver and immutable Ti/O networks. The final source changes are in the working
tree; this report does not move the existing release tag.

## Cause and change

The common per-atom energy API previously created an intermediate environment
in `accelnet.f90`, zeroed its force accumulator, then called the general CSR
entry. That entry validated the generated CSR, copied its species, indices and
displacements again into `target_workspace`, and invoked the common kernels.
For 512 atoms this took about 11% longer than the original CPU evaluator even
though OpenMP was compiled out on both sides.

`evaluate_atomic_energy_target` constructs the same single-row CSR directly in
the workspace. It preserves model/shape/species/local-environment/version-10
checks, model reload handling and independent periodic image slots. The public
C/Fortran API retains its status codes. No per-neighbor caller force accumulator
is prepared for an energy-only call.

The new entry calls the existing `execute_workspace` pipeline, including
composite descriptors and the existing network. Numerical descriptor, moment,
network and force kernel sources are unchanged from `c38c563`. The atomic entry
itself is compiled for both serial CPU and OpenMP target. It does not restore a
legacy CPU evaluation path or add a CPU-specific numerical formula.

## Validation

- GNU CPU correctness: 54/54 passed.
- GNU array bounds/runtime checks: 22/22 selected tests passed.
- NVHPC 25.3, H100 NVL, mandatory offload: 35/35 GPU tests passed.
- Performance suite: 4/4 passed, including the unchanged aenet 1.10 limit.
- LAMMPS and other GPU vendors were not rerun for this API-preparation change.

Tests now compare directly prepared atomic energy to the independent reference
for both elements in each existing batch case, including direct/moment,
Chebyshev versions, generic descriptors, element mappings, composites, zero
neighbors, periodic images, reused workspaces and subsequent force calls.
Malformed model/species/coordinate inputs are also checked. Binary-exact
coordinate shifts preserve the hard-cutoff boundary fixture; the initial test
used decimal shifts that rounded an exactly-on-cutoff distance across its
discontinuity. The test was corrected, not its numerical tolerance.

## Timing method

GNU Fortran 11.4, release -O3; CPU core 6; no OpenMP references in either CPU
executable. Alternating before/after execution, with the same fixtures, driver,
per-result FP64 checks and unchanged 1.10 performance threshold. No benchmark
overlapped our builds or correctness suites. CPU frequency was not locked;
small residual differences are not claimed to be meaningful speedups.

The final original-CPU comparison used 8 alternating rounds and a minimum of
0.5 seconds per sample. Times below cover all per-atom energy calls for a
structure, with prepared neighbor environments; they exclude neighbor-list
construction. The initial reproduction used 6 rounds of 0.3 seconds. All
reported numerical checks passed.

| Atoms | Original CPU (ms) | Fixed common CPU (ms) | Fixed/original | Before fix/original |
|---:|---:|---:|---:|---:|
| 64 | 1.475 | 1.449 | 0.9827 | 1.0970 |
| 192 | 3.791 | 3.702 | 0.9765 | 1.0767 |
| 512 | 11.946 | 12.114 | 1.0141 | 1.1110 |

The original approximately 11% slowdown at 512 atoms is reduced to about 1.4%.
At 64 and 192 atoms the fixed entry is about 2% faster than the original CPU
in this run. This is near parity, not a claim of universal speedup.

Structure-energy calls keep their existing preparation path. Their final
192/512-atom ratios are 1.0307/1.0241. The 64-atom structure-energy samples cross
a timing/frequency regime change; do not interpret its aggregate ratio below
one as an improvement from this atomic-API change. The JSON preserves every
sample and binary hash. No OpenMP parallel speedup is included.

The initial reproduction is in `baseline/report.json`. Exploration also tested
radial-loop separation, scratch arrays, small tiles and skipped dispatch scans;
none was retained. The adopted change isolates environment-preparation overhead
without changing the arithmetic. This matters: a large numerical stage in an
absolute profile need not be the source of the *difference* from the reference.
