# G4 shared value/Jacobian traversal — 2026-09-27

This report accompanies **speedupmethods.md revision 1.4**, AccelNet **1.0.1**.
The before baseline is `gpu` checkpoint **cd00106**, including revision 1.3's G4
pair-work optimizations. The after implementation is the working-tree change
identified by `source-sha256.json` and `implementation.diff`.

## What is aligned with the established CPU implementation

G4 now visits center → unordered neighbor pair (`j < k`) → descriptor, computing
values and **both** neighbor derivatives in the same traversal. Pair geometry is
computed once and cutoff, exponential, radial, and angular factors are reused.
A species-pair linked list avoids scanning all outputs for each contribution.
The NN runs afterwards and forces contract the saved Jacobian, without another
G4 geometry/radial/angular evaluation. Duplicate NN gradients remain separate
for G4 so that each saved derivative is weighted exactly once.

This is the same numerical loop for OpenMP-disabled CPU and OpenMP target GPU.
There is no device-specific G4 algorithm selection. GPU work is distributed over
centers; each center owns its CSR edges and needs no Jacobian atomics. G5 and
Chebyshev retain their previous force algorithms.

This is not a claim of identical complete pipelines: the established CPU API
handles one center's descriptors, NN and derivative contraction together. The
shared batch backend still has batch-wide descriptor/NN/force stages and a
batch-resident Jacobian. The old private value/Jacobian prototype also differed:
it iterated outputs before pairs. This revision shares pair work across outputs
in the same order as the established CPU implementation.

See [Section 15](../../../speedupmethods.md#15-g4-valuejacobian-traversal-shared-by-cpu-and-gpu-revision-14)
for equations and storage. Internal packed metadata grows from 14 to 19 integer
fields. This does not change the public C batch API. The existing default CPU
LJ/Behler dispatch is not switched by this experiment; `host` timing explicitly
selects the common OpenMP-disabled backend.

## Validation before timing

- GNU and NVHPC serial: 10 tests each passed.
- H100 GPU/host: 31 tests passed. Blackwell GPU: 29 tests passed.
- GNU runtime/bounds-checked build: 9 tests passed.
- Descriptor suite: **135 cases**, including finite differences, cutoff types
  and boundaries, integer/fractional/near-integer powers, collinear geometries,
  and mixed families. New cases verify same-capacity G5→G4 Jacobian growth/reuse
  and independent grouping across G4 components with different cutoffs/alpha.
- H100 memcheck: **0 errors** for the 4-input distinct G4, 48-input G4, and mixed
  Behler benchmark smoke cases. Full 135-case numerical suites ran separately;
  the entire descriptor suite was not run under memcheck in this revision.
- LAMMPS: **63 comparisons**, 21 each for LJ, Behler and n2p2, including neighbor
  modes, 1/2 ranks, empty ranks, NVE, and orthogonal/triclinic boxes.
- Ordinary CPU: 42 correctness and 2 performance tests passed.

The initial variable-size `private` array prototype produced invalid GPU reads
at kernel entry. Explicit bound captures did not resolve them. Moving pair caches
and species heads into explicitly mapped persistent workspace arrays removed
the failure; the final tests and timings use that implementation. The diagnostic
logs are retained. No conclusion about GPU speed is drawn from those failed runs.

## Benchmark conditions

- GNU Fortran 11.4.0 (`-O3`); NVHPC 25.3-0 (`-fast -O3`). Both CPU builds have
  `ACCELNET_TARGET_SERIAL=ON`, without `-fopenmp` or `-mp`; no direct OpenMP calls
  in the serial benchmark binaries. See `serial-audit.json`.
- CPU: Xeon Gold 6526Y, affinity core 6. This is single-core versus single-core,
  not an OpenMP thread-count setting on a parallel CPU build.
- GPU: H100 NVL (UUID ending `59e1f3c6e7f3`) and RTX PRO 6000 Blackwell Max-Q
  (UUID ending `139cedd981c9`), NVHPC 25.3, CUDA 12.8, FP64,
  `-mp=gpu -gpu=cc90,cc120`, `OMP_TARGET_OFFLOAD=MANDATORY`.
- 512 and 4096 atoms, spacing 1.7, direct mode. `g4` has four inputs including
  two duplicate mixed-species terms; `g4-distinct` has four distinct inputs;
  `g4-series` has 48 distinct inputs (powers 1–8, lambda ±1, three species pairs).
- Two rounds reverse executable order, with five alternating legacy/common
  samples per invocation (96 invocations total). Timing intervals are at least
  0.10 s, with warmup. No final benchmark overlaps our builds or tests.
- Fixed neighbor lists, built outside timing. Common timing is synchronous and
  includes input/output transfers, descriptors, NN, forces and virial. Initial
  allocations/model upload are outside steady-state timing. These are batch
  timings, not full LAMMPS MD rates.
- Tables show medians of the two per-round medians. CPU ratios use the per-sample
  common/reference ratio. All raw samples, phase timings and commands are saved.
  Changes in reference wall time across campaigns are not algorithm changes.

## CPU results

Before/fused above 1 means the new common implementation is faster. Fused/legacy
below 1 means it beats the established CPU implementation; these are separate
comparisons. Both CPU algorithms run without OpenMP compilation.

| Backend | Atoms | Model | Before common ms | Fused common ms | Before / fused | Fused / legacy CPU |
|---|---:|---|---:|---:|---:|---:|
| gnu | 512 | g4 | 16.941 | 11.829 | 1.432x | 0.825 |
| gnu | 512 | g4-distinct | 24.431 | 12.910 | 1.892x | 0.825 |
| gnu | 512 | g4-series | 30.881 | 32.752 | 0.943x | 1.038 |
| gnu | 4096 | g4 | 136.040 | 95.212 | 1.429x | 0.826 |
| gnu | 4096 | g4-distinct | 197.500 | 104.196 | 1.895x | 0.831 |
| gnu | 4096 | g4-series | 247.388 | 296.201 | 0.835x | 1.169 |
| nvhpc | 512 | g4 | 12.889 | 9.047 | 1.425x | 0.813 |
| nvhpc | 512 | g4-distinct | 17.453 | 10.000 | 1.745x | 0.804 |
| nvhpc | 512 | g4-series | 24.706 | 29.050 | 0.850x | 0.985 |
| nvhpc | 4096 | g4 | 102.951 | 72.622 | 1.418x | 0.815 |
| nvhpc | 4096 | g4-distinct | 139.868 | 80.665 | 1.734x | 0.812 |
| nvhpc | 4096 | g4-series | 197.842 | 260.596 | 0.759x | 1.095 |

## GPU results

Before/fused above 1 is a GPU improvement; below 1 is a GPU regression. This
column is not a GPU-versus-CPU speedup.

| Backend | Atoms | Model | Before common ms | Fused common ms | Before / fused | Fused / legacy CPU |
|---|---:|---|---:|---:|---:|---:|
| h100 | 512 | g4 | 1.317 | 5.522 | 0.238x | — |
| h100 | 512 | g4-distinct | 1.477 | 5.967 | 0.248x | — |
| h100 | 512 | g4-series | 1.500 | 26.896 | 0.056x | — |
| h100 | 4096 | g4 | 2.533 | 6.417 | 0.395x | — |
| h100 | 4096 | g4-distinct | 2.919 | 6.839 | 0.427x | — |
| h100 | 4096 | g4-series | 5.586 | 30.406 | 0.184x | — |
| blackwell | 512 | g4 | 3.369 | 12.445 | 0.271x | — |
| blackwell | 512 | g4-distinct | 3.799 | 13.639 | 0.279x | — |
| blackwell | 512 | g4-series | 2.514 | 33.957 | 0.074x | — |
| blackwell | 4096 | g4 | 8.014 | 13.140 | 0.610x | — |
| blackwell | 4096 | g4-distinct | 9.943 | 14.309 | 0.695x | — |
| blackwell | 4096 | g4-series | 11.330 | 34.829 | 0.325x | — |


## Interpretation and measured phase costs

For the four-input models, the fused common CPU code is faster than both the
previous common code and the established CPU evaluator. This confirms that the
earlier extra traversal was an avoidable cost in those CPU cases. The result is
not universal: the 48-input case loses performance against the previous common
code, and the table above separately records whether it beats the established
CPU evaluator with each compiler and size.

Both GPUs regress for all measured G4 cases. The force stage becomes cheaper,
but descriptor/Jacobian construction grows much more expensive. The following
timings are synchronized phase wall times at 4096 atoms, in milliseconds; totals
also include geometry, transfers and output handling. Each entry is before → fused.

| GPU | Model | Descriptor stage ms | NN ms | Force stage ms | Total ms |
|---|---|---:|---:|---:|---:|
| h100 | g4-distinct | 0.999 → 5.827 | 0.042 → 0.041 | 1.032 → 0.109 | 2.919 → 6.839 |
| h100 | g4-series | 3.395 → 29.124 | 0.146 → 0.135 | 1.182 → 0.268 | 5.586 → 30.406 |
| blackwell | g4-distinct | 3.973 → 13.304 | 0.046 → 0.046 | 5.087 → 0.098 | 9.943 → 14.309 |
| blackwell | g4-series | 6.109 → 33.483 | 0.153 → 0.185 | 4.226 → 0.296 | 11.330 → 34.829 |

The code now distributes G4 work only over centers and writes each descriptor's
two derivative contributions into a batch-resident Jacobian. The old common code
distributed work over centers/outputs and directed edges and contracted NN
gradients without storing that full Jacobian. These structural differences are
established by the source. The phase data identify the descriptor/Jacobian stage
as the source of the regression; they do not distinguish limited parallelism,
register pressure, memory transactions or cache behavior. No hardware-counter
measurement was made in this campaign. Thus these results apply to this shared
implementation, not to every possible one-traversal GPU implementation.

The established CPU evaluator also consumes a center's Jacobian immediately,
whereas this batch pipeline retains derivatives for all centers across the NN
stage. Exploring that storage lifetime and the GPU work distribution remains
useful; this comparison does not establish that separate CPU/GPU math is needed.
The measured fused implementation is retained in source, without automatic
device-based selection of the old G4 algorithm.

## LAMMPS candidate

The validated candidate is `/home/nagai/AccelNetGPU/lammps-accelnet-gpu/lmp-g4-fused`.
The existing `lmp` executable remains revision 1.3 for comparison; it does not
contain this experiment. `binaries.json` identifies the new build; the candidate
copy has the same SHA256. No full-MD performance improvement is claimed here.

## Additional workspace

For these fresh workspaces, Jacobian bytes are `8 * 3 * D * E`, pair-cache bytes
are `8 * 12 * D * N`, and head-table bytes are `4 * 2 * 2 * N`, where `D` is the
input dimension and `E` the CSR edge count. These arrays also have host allocations,
stay resident on the GPU, and are not transferred each evaluation. Capacities are
retained on reuse. Other work/model arrays are not included in this table.

| Atoms | Inputs | Edges | Jacobian MiB | Pair caches + heads MiB | Additional total MiB |
|---|---:|---:|---:|---:|---:|
| 512 | 4 | 14824 | 1.357 | 0.195 | 1.552 |
| 512 | 48 | 14824 | 16.286 | 2.258 | 18.544 |
| 4096 | 4 | 118750 | 10.872 | 1.562 | 12.434 |
| 4096 | 48 | 118750 | 130.463 | 18.062 | 148.525 |

## Reproduction and identities

`bench.py`, `build.py`, `checks.py`, and `audit.py` contain the actual commands and
paths used on this machine. `before-binaries.json` and `binaries.json` identify the
saved before and final binaries; `before-src/` stores the relevant baseline sources.
Rebuild before from `cd00106` and after from this implementation with the same
compiler and options. Public benchmark fixtures are identical between both.

For example, substitute the appropriate executable paths:

```sh
python3 AccelNetPredictor/benchmark/compare_chebyshev_variants.py \
  --variant before host /path/to/before \
  --variant after host /path/to/after \
  --family g4-distinct --sizes 512 4096 --orders 8 --modes 1 \
  --rounds 2 --seconds 0.10 --cpu 6 --output /tmp/g4-fused-comparison
```

Use `gpu` for both backend arguments and select the intended GPU UUID for the GPU
comparison. The benchmark checks energy, force and virial tolerances after every
sample and asserts that the allocation/upload counts remain stable. Initial
benchmark E/F/W comparisons have maximum absolute difference 2.22e-16.
