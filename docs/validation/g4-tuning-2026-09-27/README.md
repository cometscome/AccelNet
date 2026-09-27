# G4 GPU tuning with shared CPU/GPU arithmetic — 2026-09-27

Methods document **revision 1.5**, AccelNet **1.0.1**. The baseline for the source
diff is `gpu` checkpoint **2bcc603** (revision 1.4). The older two-traversal
checkpoint **cd00106** (revision 1.3) is also measured. Final source and binary
identities are recorded in `source-sha256.json`, `implementation.diff`, and
`binaries.json`. The candidate is `lmp-g4-optimized`; existing `lmp` remains v1.3
and `lmp-g4-fused` remains v1.4 for comparison.

## Final change

- G4 still computes values and both neighbor derivatives together, saves a full
  Jacobian, runs the NN and contracts that Jacobian without a second pair pass.
- Scalar cutoff/exponential/radial/angular caches replace per-pair global scratch
  and scans of unrelated descriptor groups. Only matching species-pair lists run.
- Packed fields 19–21 flatten these lists. Descriptor ownership is fixed by the
  global packed position modulo the owner count, in both initialization and every
  pair update. This avoids G4 value/Jacobian atomics and pair-loop barriers.
- A single flat center/owner GPU loop replaces nested parallel regions. The CPU
  compiled without OpenMP has one owner and executes the same numerical loop.
  GPU owners replicate pair geometry; this is a scheduling tradeoff, not a second
  CPU/GPU mathematical implementation. Per-descriptor pair order is preserved.
- GPU owner count is a power of two capped at 32, targeting the largest matching
  list and at least 16384 center/owner tasks. It is 32 for both 512-atom fixtures,
  4 for the 4096-atom four-input fixture, and 16 for the 4096-atom 48-input fixture.
  The 32-thread limit and this heuristic were tested on these two NVIDIA GPUs;
  they do not establish the optimal layout on other hardware or models.
- The CSR Jacobian layout remains `(Cartesian, descriptor, edge)`. Removing
  `g4_scratch` saves 18 MiB on device and the same host reservation for 4096 atoms
  and 48 inputs. The 130.46 MiB Jacobian remains, as do small species-pair heads.
  There are no per-evaluation Jacobian transfers or device-kernel allocations.
- Default CPU LJ/Behler dispatch is unchanged. CPU timings explicitly select the
  common `host` backend compiled without OpenMP; the established CPU evaluator is
  separately timed in the same executable.

Equations, ownership and tradeoffs are in
[Section 16](../../../speedupmethods.md#16-g4-scalar-caches-and-flat-descriptor-ownership-revision-15).

## Numerical and integration checks

All final suites run after the implementation was selected. Every exercised test
executable was rebuilt/relinked against the final libraries.

- GNU and NVHPC serial: 10 tests each.
- H100 GPU/host: 31 tests; Blackwell GPU: 29 tests.
- GNU runtime/bounds-checked build: 9 tests.
- Descriptor suite: **137 cases**, including the existing cutoff, fractional/
  near-integer power, mixed-family, reload/residency, and finite-difference cases.
  New cases cover interleaved cache groups and a 72-input G4 model wider than the
  owner count. The latter assigns multiple columns to an owner across its lists.
- H100 Compute Sanitizer memcheck: three selected 8-atom benchmark cases
  (`g4-distinct`, `g4-series`, `behler`), zero errors. Full numerical suites ran
  separately; the complete 137-case suite was not run under memcheck.
- LAMMPS: 63 comparisons, 21 each for LJ/Behler/n2p2, including 1/2 ranks,
  neighbor modes, empty ranks, NVE and orthogonal/triclinic cells.
- Ordinary CPU: 42 correctness tests and two performance gates.

An initial broad `all` build also attempted unrelated standalone benchmarks:
GNU 11.4 hit an internal compiler error in `benchmark-g4-derivative` and a line
length error in `benchmark-g5-scaling`. Building extra NVHPC C API tests exposed
an existing PIE link mismatch. These diagnostics are retained. They are not
claimed to pass or attributed to this G4 change. The final build commands select
and relink all executables used by the stated suites; the NVHPC target C API uses
its configured Fortran linker. The ordinary CPU C API tests use the GNU build.
The checked GNU build also caught an overlong new OpenMP clause; it was wrapped
within the standard line limit before the final rebuild and validation.

## Formal benchmark conditions

- Intel Xeon Gold 6526Y, affinity CPU 6. GNU Fortran 11.4.0 `-O3`; NVHPC 25.3
  `-fast -O3`. CPU builds have `ACCELNET_TARGET_SERIAL=ON`, no OpenMP compiler
  flags or direct OpenMP runtime references (`serial-audit.json`). Revision 1.3
  CPU identities/audit are also retained in the preceding validation archive.
- H100 NVL (UUID ending `59e1f3c6e7f3`) and RTX PRO 6000 Blackwell Max-Q (UUID
  ending `139cedd981c9`), NVHPC 25.3, CUDA 12.8, FP64,
  `-mp=gpu -gpu=cc90,cc120`, mandatory offload. The other occupied H100 was unused.
- `g4-distinct`: four distinct inputs. `g4-series`: 48 distinct inputs, powers
  1–8, lambda ±1, three unordered species pairs. 512/4096 atoms, spacing 1.7.
- Fixed neighbor lists outside timing. Common timings include synchronous data
  transfers, descriptors, NN, force and virial; initial allocation/model upload
  is excluded. These are batch timings, not full LAMMPS MD throughput.
- Two rounds reverse the order of v1.3/v1.4/v1.5. Each invocation alternates five
  legacy/common samples of at least 0.10 s, after warmup: **96 invocations**.
  No formal benchmark overlaps our builds, tests or profiling.
- Tables are medians of the two round medians. CPU/common-to-legacy ratios use
  per-sample ratios. Every sample checks E/F/virial and stable residency counts.
  Profiled runs and short exploratory runs are separate from the formal results.

## CPU results

A speedup column above 1 means v1.5 is faster than that common-code revision.
The last column is different: below 1 means v1.5 beats the established CPU
evaluator. Every CPU algorithm here is single-core with OpenMP compilation off.

| Backend | Atoms | Inputs | v1.3 ms | v1.4 ms | v1.5 ms | v1.4 / v1.5 | v1.3 / v1.5 | v1.5 / legacy CPU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| gnu | 512 | 4 | 24.425 | 12.912 | 9.562 | 1.35x | 2.55x | 0.613 |
| gnu | 512 | 48 | 30.899 | 32.531 | 24.528 | 1.33x | 1.26x | 0.779 |
| gnu | 4096 | 4 | 197.279 | 104.175 | 77.746 | 1.34x | 2.54x | 0.621 |
| gnu | 4096 | 48 | 247.330 | 296.399 | 213.038 | 1.39x | 1.16x | 0.807 |
| nvhpc | 512 | 4 | 17.453 | 10.003 | 7.340 | 1.36x | 2.38x | 0.592 |
| nvhpc | 512 | 48 | 24.680 | 29.073 | 23.390 | 1.24x | 1.06x | 0.792 |
| nvhpc | 4096 | 4 | 139.729 | 80.511 | 59.906 | 1.34x | 2.33x | 0.601 |
| nvhpc | 4096 | 48 | 197.924 | 259.779 | 204.053 | 1.27x | 0.97x | 0.853 |

## GPU results

These speedups compare GPU implementations, not GPU versus CPU.

| Backend | Atoms | Inputs | v1.3 ms | v1.4 ms | v1.5 ms | v1.4 / v1.5 | v1.3 / v1.5 | v1.5 / legacy CPU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| h100 | 512 | 4 | 1.464 | 5.965 | 1.336 | 4.47x | 1.10x | — |
| h100 | 512 | 48 | 1.501 | 26.364 | 1.496 | 17.62x | 1.00x | — |
| h100 | 4096 | 4 | 2.933 | 6.840 | 3.420 | 2.00x | 0.86x | — |
| h100 | 4096 | 48 | 5.574 | 30.262 | 4.620 | 6.55x | 1.21x | — |
| blackwell | 512 | 4 | 3.820 | 13.664 | 2.998 | 4.56x | 1.27x | — |
| blackwell | 512 | 48 | 2.513 | 34.166 | 2.823 | 12.10x | 0.89x | — |
| blackwell | 4096 | 4 | 9.962 | 14.224 | 7.680 | 1.85x | 1.30x | — |
| blackwell | 4096 | 48 | 11.346 | 34.911 | 10.837 | 3.22x | 1.05x | — |

## Measured outcome

At 4096 atoms, the following table compares the starting one-traversal v1.4
against the final v1.5. Times are milliseconds per complete batch evaluation.

| Backend | 4 distinct inputs: v1.4 → v1.5 ms | 48 inputs: v1.4 → v1.5 ms |
|---|---:|---:|
| gnu | 104.175 → 77.746 | 296.399 → 213.038 |
| nvhpc | 80.511 → 59.906 | 259.779 → 204.053 |
| h100 | 6.840 → 3.420 | 30.262 → 4.620 |
| blackwell | 14.224 → 7.680 | 34.911 → 10.837 |

Every measured CPU/GPU case improved relative to v1.4. Final common CPU time
relative to the established CPU evaluator ranges from **0.592 to
0.853** across the tested compilers, sizes and models.
GPU speedups over v1.4 range from **1.85x to 17.62x**.

The older two-traversal v1.3 comparison is less uniform. The following measured
cases are still slower than v1.3; the formal tables above also retain the cases
that improved. These short campaigns do not establish statistical significance
for differences of only a few percent.

| Backend | Atoms | Inputs | v1.5 / v1.3 time |
|---|---:|---:|---:|
| nvhpc | 4096 | 48 | 1.031 |
| h100 | 4096 | 4 | 1.166 |
| blackwell | 512 | 48 | 1.123 |

The 4096-atom GPU phase breakdown below compares v1.4 → v1.5, in milliseconds.
The descriptor stage includes geometry/radial preparation, species-pair heads
and value/Jacobian accumulation. Other stages and transfers are part of total
time but not all are shown in this table.

| GPU | Inputs | Descriptor stage | NN | Force stage |
|---|---:|---:|---:|---:|
| h100 | 4 | 5.831 → 2.411 | 0.041 → 0.042 | 0.108 → 0.106 |
| h100 | 48 | 28.989 → 3.348 | 0.135 → 0.136 | 0.268 → 0.267 |
| blackwell | 4 | 13.227 → 6.689 | 0.046 → 0.046 | 0.098 → 0.097 |
| blackwell | 48 | 33.566 → 9.540 | 0.185 → 0.185 | 0.296 → 0.263 |

The final shared implementation is retained. Large private Jacobians, atomically
accumulated pair parallelism, interleaved Jacobian storage and nested descriptor
worksharing are not selected at runtime. The only CPU/GPU G4 scheduling difference
is the owner count; there is no device-specific value/derivative formula.
Further tuning could avoid geometry for owners without a matching descriptor
and reduce the full Jacobian's memory traffic. This report does not claim those
additional changes were implemented or measured.

## Exploratory trials

Short pilot runs use five samples with at least 0.03 s per timing interval.
They guided selection; they are not substitutes for the alternating final runs.
Times below are milliseconds at 4096 atoms. A dash means that trial was not run.
The eight-center layout reserved extra padded memory and is not in the final code.
The fixed private Jacobian prototype had a capacity fallback, also not retained.

| Trial | H100 4 | H100 48 | Blackwell 4 | Blackwell 48 |
|---|---:|---:|---:|---:|
| before | 6.837 | 30.483 | 14.550 | 35.077 |
| threads32 | 5.678 | 23.466 | 12.337 | 28.428 |
| threads64 | 6.158 | 25.805 | 12.842 | 31.929 |
| threads128 | 6.777 | 30.818 | 14.232 | 35.072 |
| scalars | 4.156 | 15.797 | 9.651 | 21.351 |
| local-jac | 4.247 | 28.000 | 10.023 | 34.111 |
| pair-parallel | 8.437 | 36.708 | 13.622 | 44.050 |
| tile8 | 4.265 | 27.232 | 9.715 | 27.602 |
| feature-parallel | 29.912 | 46.357 | 50.545 | 58.406 |
| feature-teams | 9.155 | 14.025 | 17.189 | 20.536 |
| flat | 3.601 | 4.600 | 8.156 | 11.004 |
| flat32 | 4.682 | 5.327 | 12.567 | 11.154 |
| flat4 | 3.389 | — | 7.669 | — |

`flat` used 2/16 owners for the four/48-input models; `flat32` used 32 for both;
`flat4` tested four owners on the four-input 4096-atom case. The final launch
heuristic combines the small-system and large-system findings in one rule.

## Profiling evidence and limits

`cuobjdump --dump-resource-usage` reports roughly 99 KB of stack per thread for
the private-Jacobian prototype versus 128 bytes for scalar reuse on sm_90.
Nested pair/descriptor parallelism reports roughly 1.5 KB of stack and 255
registers. These are compiler resource allocations, not measured occupancy.
Nsight Systems on H100, 4096 atoms/48 inputs, records:

| Trial | CUDA grid x | CUDA block x |
|---|---:|---:|
| scalars | 128 | 32 |
| nested feature worksharing | 132 | 64 |
| nested with explicit team count | 4096 | 64 |
| final flat value/Jacobian kernel | 2048 | 32 |
| final species-pair head preparation | 128 | 32 |

Thus a source-level 32-thread OpenMP limit did not imply a 32-thread physical
CUDA block for nested execution. Flat ownership avoids that nested structure.
The final sm_90 value/Jacobian kernel reports 255 registers and 144 bytes of
stack per thread; the small head-preparation kernel reports 32 registers and
zero stack. Register allocation remains high, so this is not a claim that all
GPU resource limits have been removed.
Nsight Compute returned `ERR_NVGPUCTRPERM`; the failed attempt is retained, and
no hardware bandwidth/cache/occupancy counters are claimed. Atomics, scheduling,
geometry duplication and memory access all changed between some trials; the
timings do not isolate every individual hardware contribution.

## Reproduction

`run.py`, `final_bench.py`, `checks.py` and `audit.py` record actual commands and
paths on this machine. The saved before executables and the final ones are
identified by SHA256; formal JSON also records the exact command for every run.
Rebuild v1.3 from `cd00106`, v1.4 from `2bcc603`, and v1.5 from the archived
implementation with identical compiler options. The unchanged public benchmark
driver can also compare supplied executable paths:

```sh
python3 AccelNetPredictor/benchmark/compare_chebyshev_variants.py \
  --variant before host /path/to/v14 \
  --variant after host /path/to/v15 \
  --family g4-distinct --sizes 512 4096 --orders 8 --modes 1 \
  --rounds 2 --seconds 0.10 --cpu 6 --output /tmp/g4-tuning-comparison
```

Use `gpu` for both backends and select the intended GPU UUID for GPU comparison.
Prototype snapshots and launch traces are diagnostic artifacts; the final
source diff and hashes identify the production change. The maximum initial
E/F/virial absolute difference over the 96 formal runs is 2.22e-16.
