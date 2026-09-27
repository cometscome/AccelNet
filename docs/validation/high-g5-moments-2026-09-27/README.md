# Exact high-order G5 moments — revision 1.10

Baseline: `gpu` **2a96dcc**, AccelNet **1.0.1**. The common batch kernel now
supports exact integer G5/n2p2 type-9 orders **11--16** in explicit moment modes.
Chebyshev and lower-order G5 moments were already implemented. No approximation
of compact angular windows or fractional powers is introduced. Mathematical
scope and formulas are in [speedupmethods.md, Section 21](../../../speedupmethods.md#21-exact-high-order-g5-moments-revision-110).

## Code sharing and selection

CPU batch, LAMMPS CPU batch, and GPU target execution use the same Fortran
numerical source. The retained atomic CPU implementation is an independent
reference and compatibility API; its moment limit remains 10, so it evaluates
higher orders directly. The repository therefore still contains that older
implementation, but the new high-order algorithm has one CPU/GPU source.

Modes: `0` preserves automatic integer orders 1--10; `1` forces direct;
`2` permits integer orders 1--16 with the existing 16-neighbor threshold;
`3` forces eligible moments. Orders above 16, fractional and near-integer
orders remain direct. LAMMPS `g5 moment` selects mode 3. Explicit high-order
moments are not guaranteed to be faster, particularly on a single CPU core.

## Measurement protocol

* CPU: Intel Xeon Gold 6526Y, **one core (CPU 6)**, GNU 11.4 `-O3`, FP64,
  **OpenMP compiled out**, verified by runtime-symbol inspection. The timed
  public CPU batch path includes metadata packing.
* GPU: **H100 NVL**, NVIDIA HPC SDK 25.3, `-fast -O3`, OpenMP target
  `-mp=gpu -gpu=cc90,cc120`. Wall time includes uploads/downloads, descriptors,
  NN, and force/virial assembly; model loading and warmed allocation are excluded.
* Same deterministic two-element NN/model and CSR environment for direct and
  moment. All supplied edges are inside Rc=3.4. `g5-high` has six descriptors:
  one order, three unordered species pairs, both lambda signs. `g5-high-series`
  has 96 descriptors: orders 1--16 sharing their radial moment group.
* 256/2048 centers, 64 neighbors per center for isolated high orders; 256
  neighbors for the 256-center shared-order family. These are synthetic inference
  benchmarks, **not LAMMPS step times** or equally accurate trained potentials.
* Two reversed process rounds, five alternating warmed samples per process,
  at least 0.05 seconds per sample (or one full evaluation if longer). Results
  are medians of ten samples. The core warms for 20 seconds; its frequency is
  not locked. Builds and other tests do not run concurrently with timings.
* Every energy, force and virial component is compared against the retained
  direct reference before and during measurement (`2e-10 + 2e-10*abs(ref)`).
  Independent n2p2 accuracy tests are separate from the timed fixture.

Increasing the number of centers with the neighbor count held fixed does not
change the per-center pair count: both methods scale linearly with the number
of centers. Changes in speed ratio across those rows reflect hardware
utilization and storage effects, not a change from linear to quadratic scaling
in the number of centers. The relevant pair-to-moment reduction is in the
**number of neighbors per center**.

## High-order direct versus moment

Milliseconds per fixed-CSR inference, median of ten samples. GPU speedup is
GPU direct time divided by GPU moment time; a value below one is a slowdown.

| Centers | G5 orders | SFs | Neighbors | CPU direct | CPU moment | H100 direct | H100 moment | GPU moment speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 256 | 11 | 6 | 64 | 34.109 | 28.402 | 2.423 | 1.347 | 1.80x |
| 256 | 12 | 6 | 64 | 33.587 | 37.432 | 2.481 | 1.422 | 1.75x |
| 256 | 13 | 6 | 64 | 33.844 | 50.650 | 2.430 | 1.524 | 1.60x |
| 256 | 14 | 6 | 64 | 33.424 | 68.671 | 2.471 | 1.661 | 1.49x |
| 256 | 15 | 6 | 64 | 33.685 | 88.079 | 2.437 | 1.795 | 1.36x |
| 256 | 16 | 6 | 64 | 33.708 | 108.173 | 2.484 | 1.950 | 1.27x |
| 2048 | 11 | 6 | 64 | 272.277 | 333.227 | 4.378 | 2.889 | 1.52x |
| 2048 | 12 | 6 | 64 | 266.291 | 443.075 | 4.400 | 3.129 | 1.41x |
| 2048 | 13 | 6 | 64 | 272.986 | 564.107 | 4.392 | 3.349 | 1.31x |
| 2048 | 14 | 6 | 64 | 267.181 | 711.789 | 4.486 | 3.694 | 1.21x |
| 2048 | 15 | 6 | 64 | 269.201 | 880.058 | 4.391 | 4.091 | 1.07x |
| 2048 | 16 | 6 | 64 | 268.576 | 1062.185 | 4.352 | 4.551 | 0.96x |
| 256 | 1--16 | 96 | 256 | 770.896 | 400.068 | 8.242 | 6.353 | 1.30x |

At 256 centers, H100 moments are 1.27--1.80x faster for the six-SF isolated-order
models. At 2048 centers, order 16 is **4.6% slower** than GPU direct. On one CPU
core the isolated high-order transform can be up to **3.96x slower**, whereas
the 96-SF shared-order model is **1.93x faster** on CPU and **1.30x faster** on
H100. These results support keeping high-order moments opt-in. The shared-order
case changes both descriptor reuse and neighbor count; it does not isolate the
effect of either variable. The larger 2048-center shared-order timing was
omitted to bound redundant legacy-reference work; no result is claimed for it.

The 969-term basis at order 16 increases arithmetic and scratch storage. The
CPU slowdown and varying GPU ratio are consistent with these costs and reuse,
but these timings alone do not isolate cache bandwidth or occupancy.

Raw samples and commands: [high-g5-performance.json](high-g5-performance.json).

## CPU OpenMP OFF versus ON

Both builds execute the prepared common `host` path with GNU 11.4 `-O3` and
FP64. OFF has no OpenMP runtime calls; ON adds `-fopenmp -foffload=disable
-ffree-line-length-none`. Cores 6--13 are eight distinct physical cores on
socket/NUMA node 0 (no SMT). Runs explicitly bind 1/2/4/8 workers to the first
1/2/4/8 cores, disable dynamic/nested parallelism and GPU offload, and use
`OMP_WAIT_POLICY=PASSIVE`, `GOMP_SPINCOUNT=300000`. Runtime affinity logs verify
actual worker placement. The reference is evaluated once; five candidate
samples per process are all checked against it. Two reversed process rounds
produce ten samples per path. All cases below have 64 supplied CSR neighbors;
only the high-order family rescales every neighbor inside its cutoff.

Milliseconds per inference. The last column uses **OFF/1**, not ON/1, as baseline.

| Family | Centers | Method | OFF/1 | ON/1 | ON/2 | ON/4 | ON/8 | OFF/1 over ON/8 |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| Chebyshev degree 8 | 2048 | direct | 122.892 | 129.587 | 101.191 | 69.893 | 51.832 | 2.37x |
| Chebyshev degree 8 | 2048 | moment | 74.334 | 80.066 | 77.727 | 63.083 | 48.202 | 1.54x |
| G5 orders 1--4 | 2048 | direct | 208.350 | 439.972 | 264.938 | 153.474 | 162.978 | 1.28x |
| G5 orders 1--4 | 2048 | moment | 31.897 | 34.552 | 46.116 | 42.561 | 34.049 | 0.94x |
| G5 order 16 | 256 | direct | 34.186 | 46.073 | 50.476 | 33.995 | 22.844 | 1.50x |
| G5 order 16 | 256 | moment | 107.830 | 104.049 | 92.410 | 61.930 | 43.009 | 2.51x |
| G5 order 16 | 2048 | direct | 273.922 | 370.180 | 263.833 | 136.536 | 121.457 | 2.26x |
| G5 order 16 | 2048 | moment | 1067.091 | 1035.348 | 625.825 | 309.917 | 246.727 | 4.32x |

The low-order G5 direct case has substantial process-to-process variation:
OFF/1 round medians are 242.27 and 173.67 ms, ON/1 medians 440.09 and 341.50 ms,
and ON/8 medians 143.16 and 171.02 ms. Its aggregate ratio is not a precise
hardware-independent number. CPU frequency was not locked; these measurements
do not isolate the cause of this variation. The other complete samples and
stage timings are retained in [openmp-performance.json](openmp-performance.json).

### Why enabling OpenMP is not free

The initial thread audit caught a real control issue: unqualified
`if(device /= omp_get_initial_device())` on combined target/parallel directives
also disabled host parallelism. The correction is `if(target:...)`, selecting
only offload placement according to the [OpenMP if-clause rules](https://www.openmp.org/spec-html/5.2/openmpse17.html).
No scalar formula or direct/moment selection is changed. Setting
`OMP_NUM_THREADS` on the **ordinary serial CPU API still has no effect**;
this test explicitly uses the optional target API on the host. It is not a
LAMMPS CPU multithread benchmark.

G5 direct retains the host pair-once force traversal. With OpenMP OFF its
atomic directives disappear; with OpenMP ON they compile to locked updates,
even for one worker. The compiled force loop contains `lock cmpxchg` at its
pair and own-edge accumulation sites (see
[g5-force-openmp-disassembly.txt](g5-force-openmp-disassembly.txt)). For orders
1--4, the force-stage median rises from 144.49 ms (OFF/1) to 335.28 ms (ON/1);
the descriptor stage also rises from 61.53 to 101.44 ms. Thus force accumulation
is the largest observed source of the ON/1 regression; atomic instructions
are a concrete additional cost, though this experiment does not attribute
every extra millisecond exclusively to them.

Low-order G5 moment is already inexpensive in serial. At ON/8 its descriptor
stage falls from 22.38 to 14.37 ms, but its force stage rises from 5.96 to
16.11 ms, leaving total time 6.7% above OFF/1. Higher-order moment has enough
work to gain 4.32x at 2048 centers and eight threads. The force assembly and
parallel-region overhead therefore need improvement before making threaded
host execution the default. A next optimization is host ownership by center
so its pair-once force accumulation needs no inter-thread atomic updates,
combined with coarser parallel regions. The numerical routines can remain shared.

## Validation

* 49 existing GNU functional CTests passed; the new independent
  `n2p2_high_g5_moments` CTest also passed.
* Shared G5 tests: **209 checks** on GNU CPU, GNU with bounds checks, and H100;
  maximum E/F/virial difference from the retained reference **2.23e-16** (rounded
  upper bound). Every order 11--16 is exercised, with finite differences,
  collinear/zero-coordinate cases, multiple species, cutoffs, mixed modes and
  workspace resizing. Eligibility assertions guard against accidentally testing
  a direct fallback; auto/direct and unsupported-degree fallbacks are checked.
* Independent n2p2 comparisons passed on GNU CPU, H100 and Blackwell for a
  two-element high-G5 model and a four-element/all-descriptor mixture. Orders
  include 1/10/11--16, above-cap 17, fractional and near-integer values. Modes
  0--3 and converted native models agree. Maximum independent force finite-
  difference error is **5.85e-12**, and strain error **1.76e-11**, atomic units.
* LAMMPS high-order model: **21 CPU/GPU comparisons**, 1/2 MPI ranks,
  orthogonal/triclinic/empty-rank cases and GPU neighbor modes no/yes/hybrid,
  with `g5 moment` on GPU against CPU direct.
* OpenMP host: 209 G5 checks and 137 other-descriptor checks passed at both
  2 and 8 threads after the target-condition fix; the independent n2p2
  two-/four-element mixtures passed at 8 threads. H100 G5 and independent
  mixtures were rechecked after that fix.
* H100 Compute Sanitizer, four-element/all-descriptor mixture, forced moments:
  **0 errors**.

* Existing common-CPU and low-order G5 performance guards and the order-scaling
  smoke test passed. The focused direct regression against the previous revision
  and n2p2 passed for types 2/3/21/22 (n2p2 slowdown cap 1.10).

## Reproduction

Build with the same flags as the previous [direct optimization report](../n2p2-cpu-2026-09-27/README.md).
The independent oracle is the unmodified n2p2 v2.3.0 reference driver documented
there. `N2P2_REFERENCE_BENCHMARK` enables the new CPU CTest;
`ACCELNET_TEST_N2P2_GPU=ON` additionally enables its GPU counterpart.

```sh
CUDA_VISIBLE_DEVICES="$H100_UUID" \
python3 AccelNetPredictor/benchmark/compare_high_moments.py \
  --cpu "$GNU_TARGET_BENCHMARK" --gpu "$NVHPC_TARGET_BENCHMARK" \
  --affinity 6 --rounds 2 --seconds .05 --output /tmp/high-g5-timings

OMP_TARGET_OFFLOAD=MANDATORY CUDA_VISIBLE_DEVICES="$GPU_UUID" \
python3 AccelNetPredictor/test/check_n2p2_high_moments.py \
  --reference "$N2P2_REFERENCE" --candidate "$NVHPC_N2P2_BENCHMARK" \
  --converter "$NVHPC_CONVERTER" --backend gpu --output /tmp/high-g5-accuracy
```

The benchmark JSON retains commands, executable hashes, all samples and
reference errors. The GPU speed changes from the preceding direct optimization
are separately recorded in the [actual LAMMPS table](../n2p2-cpu-2026-09-27/README.md#actual-lammps-step-time);
they should not be conflated with the new synthetic high-order moment timings.

## Archived build identity

[Source hashes and measurement stages](audit.json) identify the measured scope.
High-G5 GPU timings, Blackwell, LAMMPS and sanitizer results precede the host
`if(target:...)` correction; the GPU condition remains true in both versions.
The final host-parallel and H100 accuracy checks cover the correction. Compiler
flags, raw performance logs, CTest logs and numerical summaries are retained
alongside this report. Reproduction of the thread sweep is documented in
[openmp-target.md](../../openmp-target.md#cpu-openmp-thread-comparison).
