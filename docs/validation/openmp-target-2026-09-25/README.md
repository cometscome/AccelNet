# OpenMP target stages 2–3 validation — 2026-09-25

This report records the initial implementation. See the subsequent
[residency/moment validation](../gpu-residency-moments-2026-09-25/README.md)
for the optimized implementation and updated measurements.

Implemented a separate FP64 GPU backend for one Chebyshev component per element.
Descriptors, NN output/input gradients, direct derivative contraction, atomic
force scatter and virial accumulation all run on the device. No descriptor
Jacobian is stored. CPU inference source and compiler options are unchanged.
See [build/API instructions](../../openmp-target.md) for the supported scope.

## Environment and device identity

- Base commit: `c6631460a1bbb990c82e3e0ff5e73c36c52f6f9b`; implementation is the
  accompanying working-tree change. Source checksums are in [sources.json](sources.json).
- NVIDIA HPC SDK 25.3, Release (`-fast -O3`); GPU library/link flags
  `-mp=gpu -gpu=cc90,cc120`. Driver 590.48.01.
- NVIDIA H100 NVL: `GPU-2644154d-7268-af42-6631-59e1f3c6e7f3`, PCI `2A:00.0`.
- NVIDIA RTX PRO 6000 Blackwell Max-Q Workstation Edition:
  `GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9`, PCI `AB:00.0`.
- Final validation and benchmarks select these explicit UUIDs, not ordinal indices.
- Intel Xeon Gold 6526Y; timing runs pinned to logical CPU 6, `OMP_NUM_THREADS=1`.

**Identity correction:** CUDA's default ordinal 0 is Blackwell on this machine,
while `nvidia-smi` index 0 is H100. Early commentary and the stage-1 report used
`CUDA_VISIBLE_DEVICES=0` and incorrectly labeled it H100. Driver API PCI identity
inspection established the mapping. The stage-1 report was corrected and both
GPU suites/benchmarks were rerun with UUID selection. Numerical results are
unaffected by this labeling error.

## Correctness

| Check | Result | Evidence |
| --- | --- | --- |
| H100 GPU suite | 13/13 passed | [CTest log](h100-uuid-tests.log), [JUnit](h100-tests.xml) |
| Blackwell GPU suite | 13/13 passed | [CTest log](blackwell-uuid-tests.log), [JUnit](blackwell-tests.xml) |
| H100 Compute Sanitizer memcheck | 0 errors | [Log](h100-native-memcheck.log) |
| Blackwell Compute Sanitizer memcheck | 0 errors | [Log](blackwell-native-memcheck.log) |
| GNU 11.4 CPU Release suite | 44/44 passed | [Log](cpu-tests.log) |
| CPU performance gate against pre-change library | 12/12 passed | [Raw measurements](cpu-performance.json) |
| GNU 11.4 optional backend compile / negative tests | compiled; 4/4 passed | [Log](gnu-target-rejection-tests.log) |

An installation to a temporary prefix was also checked: a separate NVHPC
consumer found `AccelNetPredictor`, linked `AccelNet::Target`, imported the GPU
Fortran module and ran its model/workspace lifecycle calls successfully.

Each GPU suite includes 56 synthetic and 18 real Ti/O model/geometry/mode
comparisons, plus invalid-input rejection tests. Every atomic energy, force
component and all nine virial components are checked at
`2e-10 + 2e-10*abs(reference)`. Coordinate and nine-component strain finite
differences independently compare GPU energies with GPU forces/virial.

Maximum absolute CPU/GPU E/F/W differences:

- H100: synthetic `1.055e-15`, real Ti/O `3.340e-12`.
- Blackwell: synthetic `1.048e-15`, real Ti/O `3.354e-12`.

Tests cover Chebyshev versions 0/1/10, CPU direct/auto/moment references,
cutoffs 0–9, activations 0–11, distinct element network widths/depths, species
mapping and normalization, isolated atoms, periodic self images, triclinic cells,
subsets, additive force/virial, empty batches, model copies/reinitialization and
workspace growth/reuse. Memcheck uses the smaller `--quick` test with periodic
angular interactions, partitions and GPU coordinate/strain finite differences.
It does not cover every model/size in the full suite.

The CPU gate compares Chebyshev, LJ, G4 and G5 at 8/64/512 atoms against the
unchanged baseline library built from the base commit, with the same GNU driver
and flags. Every printed CPU E/F/W component is identical (maximum error 0).
Paired median time ratios range from 0.9921 to 1.0073, below the 1.10 noise gate.
This is measured preservation on those fixtures, not a guarantee for all models.

An initial memcheck on Blackwell with an H100-only `cc90` binary reported a
recoverable CUDA module-loading API error from the NVHPC runtime while numerical
checks passed. The [original log](blackwell-memcheck-runtime-probe.log) is retained.
The final binary includes native `cc90,cc120` code; both unfiltered memchecks
then passed with zero errors. No API errors were suppressed.

## Initial timings (not a tuned GPU implementation)

Times are milliseconds per call. Speedup is median of paired CPU/GPU times;
**greater than 1 means GPU faster**. Batch times include target allocations and
all transfers with an existing neighbor list. Full speedups include rebuilding
that list on the CPU each call. CPU references use NVHPC, the same model and
single-threaded execution. Both outputs are checked before timing and after
measuring each method. Five alternating-order samples, each at least 0.2 s,
follow warm-up. Initial model packing/device startup is excluded.

| GPU | Atoms | Order | CPU batch ms | GPU batch ms | Batch speedup | Full speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| h100 | 8 | 3 | 0.090 | 6.206 | 0.014 | 0.065 |
| h100 | 64 | 3 | 0.683 | 4.839 | 0.141 | 0.452 |
| h100 | 512 | 3 | 3.949 | 6.795 | 0.583 | 0.973 |
| h100 | 1024 | 3 | 7.139 | 6.125 | 1.166 | 0.998 |
| h100 | 512 | 8 | 16.091 | 7.627 | 2.107 | 1.069 |
| blackwell | 8 | 3 | 0.064 | 9.599 | 0.007 | 0.031 |
| blackwell | 64 | 3 | 0.492 | 7.418 | 0.066 | 0.251 |
| blackwell | 512 | 3 | 4.106 | 10.029 | 0.410 | 0.957 |
| blackwell | 1024 | 3 | 7.157 | 9.761 | 0.733 | 0.997 |
| blackwell | 512 | 8 | 16.110 | 12.455 | 1.294 | 1.026 |

Raw sample times, errors, device UUIDs and conditions:
[H100](h100-timings.json), [Blackwell](blackwell-timings.json).
The [earlier Blackwell cc90/PTX measurements](historical-blackwell-cc90-timings.json)
are retained separately with corrected device identity. Measurements on this
shared host can vary with clocks and load; small differences are not a reliable
speedup claim.

The GPU is slower for small/low-order examples. At 512 atoms/order 8, H100 batch
speedup is approximately 2.11, but full speedup is only 1.07 because CPU neighbor
construction dominates. Existing CPU APIs always keep using the CPU. There is
no automatic GPU routing, so small existing workloads incur no GPU overhead.

## Remaining scope

The GPU backend currently enumerates angular pairs. GPU moment kernels,
persistent device allocation/model data, within-atom parallelism, automatic
selection and tuning are stage-4 work. LJ, Behler/n2p2 and composite GPU models,
CLI/C/LAMMPS GPU integration, and AMD/Intel GPU validation are not implemented
by these stages. The tested Fortran GPU API is available through the optional
`AccelNet::Target` CMake target and is disabled in the default CPU build.
