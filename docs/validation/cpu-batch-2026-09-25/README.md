# CPU batch API validation — 2026-09-25

Stage 1 adds a CPU batch API and reusable outer workspace. AccelNet inference
is still CPU-only. The OpenMP target test below verifies device readiness only.

## Reference and environment

- Upstream `main` pulled to `c6631460a1bbb990c82e3e0ff5e73c36c52f6f9b`.
- Baseline library built from an archive of that commit. Both executables use
  identical benchmark driver/fixture sources from this change.
- Intel Xeon Gold 6526Y; GNU Fortran 11.4, CMake Release (`-O3 -DNDEBUG -O3`).
  Driver preprocessing enabled; single-threaded runs pinned to logical CPU 6.
- Seven alternating baseline/candidate pairs per case, at least 0.2 seconds per
  sample after warm-up. Reported ratio is median(candidate time / baseline time)
  over adjacent pairs; less than 1 means faster. Neighbor construction, energy,
  force and virial evaluation are included. Parsing/initialization are excluded.
- Correctness tolerance: `2e-10 + 2e-10*abs(reference)` for every energy, force
  and virial component. Timing gate: ratio <= 1.10, allowing measurement noise.

## Results

- GNU Release: **44/44 CTest tests passed**, including real Ti/O golden model
  comparisons, existing virial tests and the CPU performance gate.
- GNU Debug with `-fcheck=all -fbacktrace`: **7/7 batch tests passed**.
- NVIDIA HPC SDK 25.3 CPU build: **7/7 batch tests passed**.
- NVIDIA RTX PRO 6000 Blackwell Max-Q, NVIDIA HPC SDK 25.3, `-mp=gpu -gpu=cc90`,
  `CUDA_VISIBLE_DEVICES=0`: **1/1 device smoke test passed**. It rejects host
  fallback and checks array mapping plus FP64 atomic scatter on the GPU.
- Both performance comparisons: **12/12 cases passed**. The established CPU
  API produces identical printed FP64 results in every case (maximum error 0).
  Batch accumulation differs by at most approximately `6.2e-12`.
- Fault injection confirmed the comparison script exits nonzero for a doubled
  reported execution time and, independently, an energy changed by 1.

| Model | Atoms | Existing API / baseline time | Batch / baseline time | Batch max absolute E/F/W error |
| --- | ---: | ---: | ---: | ---: |
| chebyshev | 8 | 1.0012 | 1.0005 | 1.78e-15 |
| chebyshev | 64 | 0.9976 | 1.0040 | 1.49e-13 |
| chebyshev | 512 | 1.0009 | 1.0018 | 6.2e-12 |
| lj | 8 | 1.0169 | 0.9592 | 1.78e-15 |
| lj | 64 | 0.9990 | 1.0091 | 1.78e-13 |
| lj | 512 | 1.0020 | 0.9994 | 5.34e-12 |
| n2p2-g5 | 8 | 1.0066 | 1.0130 | 0 |
| n2p2-g5 | 64 | 0.9975 | 1.0193 | 0 |
| n2p2-g5 | 512 | 0.9998 | 0.9998 | 0 |
| n2p2-g4 | 8 | 1.0059 | 1.0080 | 2.22e-16 |
| n2p2-g4 | 64 | 0.9963 | 1.0008 | 2.84e-14 |
| n2p2-g4 | 512 | 0.9968 | 1.0002 | 1.66e-12 |

Raw final measurements: [existing API](cpu-existing-api.json),
[batch API](cpu-batch-api.json). Each includes individual times, paired ratios,
thresholds and errors. Paths in these files identify the temporary local builds
used for this run; rebuild instructions are in [the API guide](../../batch-api.md).

## Findings retained from development

An initial implementation selected direct derivative contraction separately per
species. The G4/G5 fixture became about 12% slower for 8/64 atoms. The final CPU
implementation uses the established structure evaluator's model-wide selection
policy. This removed that regression. The rejected implementation's
[raw measurements](historical-g4-before-fix.json) are retained.

A later run using ratios of separate medians reported a spurious 1.386 ratio
for 512-atom Chebyshev: individual baseline and candidate times switched between
roughly 0.19 and 0.136 seconds during that run, with adjacent ratios near 1.
[Those measurements](historical-unpaired-noisy.json) and a
[longer confirmation](historical-long-confirmation.json) (nine 1-second samples,
ratio approximately 1.001) are retained. The final gate takes the median of
paired ratios and both complete comparisons were rerun with that method.
This reduces sensitivity to host load/frequency changes; it cannot eliminate it.

## Limits

These are four small synthetic model families at 8/64/512 atoms, not a guarantee
for all production models or sizes. The existing CPU implementation is unchanged;
future GPU work must continue running this gate and broaden representative models.
No GPU inference timing is reported by this stage. AMD/Intel GPUs have not been
validated. **GPU identity correction during stages 2–3:** the initial report
labeled CUDA ordinal 0 as H100 based on `nvidia-smi` index 0. Driver API / PCI
identity inspection showed CUDA ordinal 0 is Blackwell (AB:00.0), while ordinal
2 is H100 (2A:00.0). Subsequent testing uses explicit GPU UUIDs. CPU measurements
and numerical results above are unaffected.
