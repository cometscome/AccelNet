# Common G5 moments — 2026-09-27

Methods document **revision 1.7**, AccelNet **1.0.1**. Baseline: `gpu`
checkpoint **1d6985d** (revision 1.6). `source-sha256.json` identifies the
implementation sources; `audit.json` records measured binaries, flags and GPU
UUIDs. The mathematical derivation is in [Section 18](../../../speedupmethods.md#18-common-g5-moments-without-losing-the-original-advantage-revision-17).

## What changed

G5 integer moments now use the same numerical source on CPU and GPU. Raw moments
are retained through the NN. Bilinear contractions are computed once per radial
group/species pair/degree. NN force coefficients are contracted by degree before
forming monomial adjoints, then differentiated multivariate Horner evaluates the
force polynomial. Pure-moment rows skip direct feature scans. These reuse the
Chebyshev optimization principles without approximating the G5 formula.

G5 auto again follows the original per-component 16-angular-neighbor threshold;
forcing G5 moments no longer falls back to the old CPU implementation. Exact
integer orders 1–10 are eligible; fractional/near-integer/higher orders remain
common direct. Chebyshev and G5 mode controls are independent. LAMMPS accepts
`pair_style accelnet/gpu auto H.nn.ascii O.nn.ascii g5 moment`.

## Controlled comparison

The earlier degree-8, cutoff-3.4 fixture cannot establish whether the original
moment advantage was preserved. This campaign uses the descriptor parameters and
local environments from `AccelNetDescriptors/benchmark/benchmark_g5_scaling.f90`:

- 54 descriptors: three eta values (0.000357, 0.028569, 0.089277), zeta 1/2/4,
  lambda ±1, three unordered species pairs, cutoff 6.5 and radial shift 0.35.
- Identical deterministic spherical environments with 16/32/64/128 neighbors per
  center, replicated over 64/512 centers. These are arbitrary local CSR
  environments, **not a physical periodic MD configuration**.
- All four methods use the same synthetic 54–8–4–1 network. Timings include the
  NN, forces and virial, whereas the original descriptor-only benchmark did not.
  CPU common timings use the public batch API and include model packing each
  call. GPU timings use a prepared model and include per-call transfers.
  Neighbor construction and initial allocations/model upload are excluded.
- Xeon Gold 6526Y, pinned CPU 6. GNU Fortran 11.4.0 `-O3`; NVHPC 25.3 `-fast -O3`.
  Both CPU executables compile without OpenMP, with no detected OpenMP runtime
  imports; this is not merely setting the thread count to one.
- NVHPC GPU compilation uses `-mp=gpu -gpu=cc90,cc120`, CUDA 12.8; H100 NVL and
  RTX PRO 6000 Blackwell Max-Q. Device offload is mandatory.
- Two rounds in reverse case order, five alternating reference/common samples
  per invocation, each at least 0.05 seconds. No builds/tests/other campaign
  measurements overlapped. Tables show medians of the two round medians in ms.
  These are descriptive measurements, without statistical-significance claims.
- All 128 invocations check energy, every force component and virial. Maximum
  initial absolute difference: **3.18e-14**. Mode 1 compares against retained
  CPU direct; mode 3 compares against retained CPU moment. Separate correctness
  suites compare moments directly with the original direct implementation.

`results.json` contains all samples, commands and phase timings; `summary.json`
contains the table values. GPU phase totals are not independent benchmarks.

### Single-core CPU, OpenMP compilation disabled

| Compiler | Centers | Neighbors | Old direct ms | Old moment ms | Common direct ms | Common moment ms |
|---|---:|---:|---:|---:|---:|---:|
| gnu | 64 | 16 | 1.859 | 3.768 | 1.986 | 1.031 |
| gnu | 64 | 32 | 7.145 | 8.457 | 6.854 | 1.703 |
| gnu | 64 | 64 | 23.660 | 18.179 | 21.410 | 2.561 |
| gnu | 64 | 128 | 94.351 | 53.909 | 82.577 | 4.901 |
| gnu | 512 | 16 | 12.438 | 25.171 | 13.374 | 8.472 |
| gnu | 512 | 32 | 47.755 | 56.451 | 45.965 | 13.016 |
| gnu | 512 | 64 | 189.614 | 145.568 | 171.331 | 22.418 |
| gnu | 512 | 128 | 765.812 | 431.084 | 660.719 | 41.189 |
| nvhpc | 64 | 16 | 1.417 | 3.573 | 1.457 | 0.793 |
| nvhpc | 64 | 32 | 5.289 | 7.669 | 5.090 | 1.260 |
| nvhpc | 64 | 64 | 20.687 | 18.648 | 19.059 | 2.214 |
| nvhpc | 64 | 128 | 82.050 | 51.460 | 72.700 | 4.163 |
| nvhpc | 512 | 16 | 11.420 | 28.501 | 11.620 | 6.979 |
| nvhpc | 512 | 32 | 42.462 | 61.183 | 40.556 | 10.898 |
| nvhpc | 512 | 64 | 165.276 | 149.391 | 151.821 | 18.684 |
| nvhpc | 512 | 128 | 657.410 | 410.708 | 582.071 | 34.489 |

Original moment beats original direct at 64/128 neighbors in this campaign.
Common moment retains and strengthens that advantage. At 512 centers/64 neighbors,
GNU common moment is **6.49x faster than old moment** and **7.64x faster than
common direct**. At 16/32 neighbors the original moment is slower than original
direct, while the common moment is faster than both in these measured cases.

### GPU, same model and local environments

| GPU | Centers | Neighbors | Common direct ms | Common moment ms | Direct/moment speedup |
|---|---:|---:|---:|---:|---:|
| h100 | 64 | 16 | 1.055 | 0.800 | 1.32x |
| h100 | 64 | 32 | 1.247 | 0.939 | 1.33x |
| h100 | 64 | 64 | 2.724 | 1.591 | 1.71x |
| h100 | 64 | 128 | 4.591 | 1.889 | 2.43x |
| h100 | 512 | 16 | 1.509 | 1.292 | 1.17x |
| h100 | 512 | 32 | 2.569 | 1.513 | 1.70x |
| h100 | 512 | 64 | 4.412 | 1.937 | 2.28x |
| h100 | 512 | 128 | 13.280 | 2.779 | 4.78x |
| blackwell | 64 | 16 | 1.542 | 0.827 | 1.86x |
| blackwell | 64 | 32 | 2.010 | 1.094 | 1.84x |
| blackwell | 64 | 64 | 3.971 | 1.694 | 2.34x |
| blackwell | 64 | 128 | 8.103 | 2.703 | 3.00x |
| blackwell | 512 | 16 | 2.209 | 1.391 | 1.59x |
| blackwell | 512 | 32 | 3.008 | 1.703 | 1.77x |
| blackwell | 512 | 64 | 6.667 | 2.384 | 2.80x |
| blackwell | 512 | 128 | 20.126 | 3.694 | 5.45x |

At 512 centers/64 neighbors H100 moment is **2.28x faster than H100 direct**,
**11.57x faster than GNU common single-core moment**, and **75.15x faster than
GNU old single-core moment**. The latter two use the separate OpenMP-disabled
GNU measurements, not the CPU section of a GPU-enabled executable. They include
the CPU/GPU model-packing scope difference described above. These are batch
kernel/API results, not LAMMPS timestep speedups.

## Controls and limits

A separate 512-atom degree-8 G5-series model (48 inputs, cutoff 3.4, spacing 1.7)
checks the earlier high-order case. This is deliberately separate from the
original scaling fixture. Two reversed rounds give:

| Backend | Old direct ms | Old moment ms | Common direct ms | Common moment ms |
|---|---:|---:|---:|---:|
| gnu | 59.681 | 135.102 | 21.248 | 15.468 |
| nvhpc | 58.930 | 165.916 | 18.118 | 13.042 |
| h100 | 58.892 | 165.287 | 1.345 | 1.361 |

Here the original moment is slower than original direct. The optimized common
CPU moment is faster than common direct, but H100 moment is approximately 1%
slower than direct. Thus moment is not universally faster on GPU; degree, radial
groups, neighbor count and batch size determine the crossover. Auto preserves
the original neighbor threshold, rather than silently replacing moments with
direct to hide a regression. More accurate automatic cost selection remains a
possible future change.

An independent before/after control uses immutable revision-1.6 binaries and
mode-1 direct evaluation at 512 atoms. Median common-path after/before ratios:

| Compiler | LJ | G4 distinct | G5 series direct |
|---|---:|---:|---:|
| gnu | 1.033 | 1.002 | 1.013 |
| nvhpc | 1.038 | 1.005 | 1.009 |

LJ shows a small approximately 3–4% increase in this fixed-neighbor public-batch
control; increased packed metadata/dispatch overhead is a plausible cause, not
an isolated profiler finding. No large direct/G4 regression is observed. The
end-to-end default-CPU regression gate below passes for both compilers.

## Validation

- GNU serial: 11 selected batch/host tests passed; NVHPC serial: 11 passed;
  checked GNU: 10 passed; ordinary non-target CPU build: 42 tests passed.
- H100: all 33 selected GPU/target-host tests passed. Blackwell: four selected
  suites passed. The new G5 suite has 185 cases, alongside 137 generic descriptor
  and 78 Chebyshev cases. GNU G5 maximum absolute E/F/W error is 2.22e-16.
- G5 coverage includes ten cutoffs, all modes, finite differences, exact order
  10, ineligible orders 11/16, fractional and near-integer powers, lambda ±1/0/0.7,
  collinear/zero-coordinate cases, periodic self images, isolated atoms, mixed
  components/central elements and per-component angular cutoffs, buffer reuse
  and mode changes. Default CPU batch results are checked as well as explicit
  target/host evaluation. C API tests cover independent mode selection and
  invalid modes.
- H100 Compute Sanitizer memcheck: zero errors in the forced-moment scaling
  fixture. This is a selected memory check, not exhaustive race detection.
- 63 LAMMPS comparisons passed: LJ/Behler/converted n2p2, CPU/GPU/hybrid neighbor
  modes, one/two MPI ranks, orthogonal/triclinic/empty-rank cases and short NVE.
  GPU requests G5 moment; the independent CPU reference explicitly uses direct.
- Existing `predictor_common_cpu_performance`: all 12 conditions passed on both
  compilers, with seven paired end-to-end samples per condition. Common/reference
  ratios are 0.850–1.017 GNU and 0.879–1.014 NVHPC (limit 1.10).
- New `predictor_g5_moment_performance`: passed on both compilers. At its fixed
  64-center/64-neighbor fixture, common/old moment ratios are 0.140 GNU and 0.119
  NVHPC; common moment/direct ratios are 0.119 and 0.116. Three alternating rounds
  protect both the original moment speed and its advantage over common direct.

## Reproduction and artifacts

`bench.py`, `checks.py`, and `final_checks.py` record the campaign, commands and
local build paths. Adapt those paths to equivalent builds. For the new opt-in
performance test:

```sh
cmake -S . -B /path/to/serial-build \
  -DACCELNET_BUILD_OPENMP_TARGET=ON -DACCELNET_TARGET_SERIAL=ON \
  -DACCELNET_PREDICTOR_BUILD_TESTING=ON \
  -DACCELNET_TEST_G5_MOMENT_PERFORMANCE=ON
cmake --build /path/to/serial-build --target accelnet-target-benchmark
ctest --test-dir /path/to/serial-build -R '^predictor_g5_moment_performance$' -V
```

Do not add OpenMP compiler flags to the serial build. The script rejects known
OpenMP runtime imports; inspect compiler flags as well, as the campaign does.
Run performance tests without competing jobs.

A single fixed-neighbor evaluation comparison is:

```sh
/path/to/accelnet-target-benchmark 512 4 0.05 3 1.7 cpu-shared g5-scaling no-neighbors 64
```

Change mode 3 to 1 for direct; use `gpu` with a GPU build for offload. Every
invocation also measures and checks the matching independent CPU reference.

Validated LAMMPS executable:
`/home/nagai/AccelNetGPU/lammps-accelnet-gpu/lmp-g5-moments`.
`lammps-binary.json` records its SHA256; older executables remain available.
`implementation.diff` records source changes from revision 1.6 with zero context.
Archived logs have trailing whitespace trimmed. Performance-test
registration was added after numerical builds; numerical sources and benchmark
executables were unchanged during the final measurements. `checks/`, `timings/`
and `lammps/` retain the validation logs and comparison reports.
