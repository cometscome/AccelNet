# Behler contraction and CPU-algorithm experiments — 2026-09-26

This report accompanies document revision 1.2 of
[speedupmethods.md](../../../speedupmethods.md), on AccelNet 1.0.1 plus the
uncommitted changes based on `c6631460a1bbb990c82e3e0ff5e73c36c52f6f9b`.
It distinguishes the adopted common implementation from the archived Jacobian
and force-accumulation experiments.

## Adopted algorithm

- Preserve the radial cache and fused LJ work from the
  [preceding revision](../generic-common-2026-09-26/README.md).
- Partition G4/G5 angular outputs by descriptor kind, unordered species pair,
  radial/cutoff parameters, and lambda. Group compatible positive integer powers
  up to 16. Fractional and higher powers retain the exact original power formula;
  only equal exponents can share a scalar coefficient in that fallback.
- Aggregate negative NN input gradients into angular coefficients at the end of
  the existing NN kernel. There is no extra launch or host transfer. Evaluate
  the polynomial and derivative together by Horner in `t=(1+lambda*cos(theta))/2`.
- Groups with a single distinct power retain the original direct value kernel.
  Models with multi-power groups share pair geometry and radial factors while
  accumulating local power sums. GPU teams reduce neighbor contributions;
  OpenMP-disabled CPU builds execute the same numerical loop serially.
- Use the same pair formula for both execution backends. CPU evaluation visits
  each unordered force pair once and updates both edges. GPU threads own directed
  edges, avoiding the force atomics that slowed the measured two-sided variant.
- Reuse existing NN-gradient slots for angular coefficients. Equal-power terms
  sum into their first slot; distinct powers stay in place. Gather a bounded local
  coefficient vector before Horner. There is no extra global coefficient array or
  edge-by-descriptor Jacobian. Pure radial/LJ CPU batches also avoid the unordered
  force path's unnecessary zeroing/additions.
- Exactly integer zeta is required for polynomial grouping. Near-integer inputs
  retain the original power helper and derivative prefactor. Degree changes within
  existing model dimensions reuse the workspace without growing a coefficient array.

The ordinary CPU LJ/Behler dispatch is unchanged. A faster common implementation
is not sufficient evidence to replace the established CPU path for every model,
especially when a default CPU G5 call could choose its existing moment algorithm.
The benchmark explicitly forces direct G5 in both candidates.

## Measurement protocol

- Intel Xeon Gold 6526Y, CPU affinity 6. GNU Fortran 11.4.0 `-O3`; NVHPC 25.3
  `-fast -O3`. Both CPU candidates are compiled without `-fopenmp` or `-mp`.
  `serial-audit.json` records flags and absence of direct OpenMP runtime calls;
  `OMP_NUM_THREADS=1` is an additional setting, not the basis for this claim.
- NVHPC 25.3 OpenMP target, CUDA FP64, H100 NVL and RTX PRO 6000 Blackwell Max-Q.
  See the root mathematical document for exact toolchain/hardware versions.
- Synchronous batches with neighbor construction and initial packing/allocation
  excluded. Input checks, copies/transfers, synchronization, energies, forces,
  and virial are included. These timings are not full LAMMPS MD timings.
- Before/after libraries use the same benchmark/model builder. The baseline
  driver is relinked against preserved revision-1.1 libraries so both binaries
  support the additional dense-basis models; numerical library code is unchanged.
- Each run alternates five established-CPU-reference/common-candidate samples
  and checks all outputs. CPU speedups use ratios to the same-run reference to
  reduce sensitivity to clock/load variation. GPU speedups use absolute times.
  The first before/after campaign is archived as candidate A. The adopted implementation
  is remeasured after replacing the extra coefficient array with existing gradient
  slots. CPU comparisons retain same-run references; GPU baseline times come from
  the initial campaign on the same devices. Small differences are not conclusive.
- `g4` and `g5` each have four inputs, `behler` has fourteen G1–G5 inputs, and
  `lj` has four inputs. `g4-series`/`g5-series` have **48 distinct angular inputs**:
  zeta 1–8, lambda ±1, and species pairs (1,1), (1,2), (2,2), sharing one radial
  function. They exercise reuse across a basis without duplicate features.
- Two species, input–8–4–1 NN, spacing 1.7, 512/4096 atoms. Independent timing
  jobs and builds from this task do not run concurrently.

## Earlier experiments

The `jacobian-trial` directory preserves commands, raw timing/phase logs, numerical
checks, serial-build evidence, and diffs against the revision-1.1 source. Its
`pilot.json` uses 512 atoms and `large.json` uses 4096 atoms. The scalar-per-feature
GPU implementation and neighbor-team implementation are separate variants;
`large-naive.json` adds the former at 4096 atoms.

The Jacobian strategy is valid and considerably improves serial CPU execution.
For example, GNU 4096-atom G5 changes from 185.75 to 94.79 ms, and its time relative
to the established CPU reference changes from 1.296 to 0.665. This was measured
with OpenMP compilation disabled. There is no demonstrated general prohibition
against using the old CPU algorithm on a GPU.

The neighbor-team Jacobian trial improved H100 G4 at 512 atoms (1.506 to 0.993 ms),
but regressed at 4096 atoms (3.293 to 5.048 ms). Blackwell improved in both of those
cases (4.242 to 1.791 ms and 11.357 to 8.283 ms). Phase timings locate the H100
large-case cost in descriptor/Jacobian construction. They do not prove a specific
memory-bandwidth or occupancy explanation without hardware-counter measurements.
This implementation was therefore not adopted universally.

The contraction experiments also compare directed force ownership with evaluating
both force directions once and atomically updating the two edge slots on GPUs.
The latter helps the OpenMP-disabled CPU, but was slower in the measured GPU
cases. Local value sums remove repeated writes from the inner pair loop; applying
that team layout indiscriminately to tiny bases still regressed large H100 cases.
The final implementation keeps the original value kernel for single-power groups.

See `experiments` for the successive contraction source diffs and pilot logs.

## Adopted measurements

The tables use the final implementation with coefficients held in existing NN-gradient
slots. CPU speedup is the before/after ratio of reference-normalized times; GPU
speedup is before/after absolute time. `Final / CPU reference` below one means the
common serial implementation is faster than the established CPU implementation.

### OpenMP-disabled single-core CPU

| Compiler | Atoms | Model | Before (ms) | Final (ms) | Normalized common speedup | Final / CPU reference |
|---|---:|---|---:|---:|---:|---:|
| gnu | 512 | lj | 1.366 | 1.417 | 0.97x | 1.085 |
| gnu | 512 | g4 | 36.098 | 21.659 | 1.65x | 1.523 |
| gnu | 512 | g5 | 22.371 | 13.711 | 1.63x | 0.801 |
| gnu | 512 | behler | 61.569 | 38.029 | 1.61x | 0.755 |
| gnu | 512 | g4-series | 315.910 | 32.813 | 9.60x | 1.045 |
| gnu | 512 | g5-series | 164.990 | 21.862 | 7.44x | 0.367 |
| gnu | 4096 | lj | 11.022 | 11.299 | 0.98x | 1.078 |
| gnu | 4096 | g4 | 292.336 | 174.569 | 2.12x | 1.191 |
| gnu | 4096 | g5 | 178.512 | 108.780 | 1.64x | 0.793 |
| gnu | 4096 | behler | 496.441 | 305.444 | 1.66x | 0.737 |
| gnu | 4096 | g4-series | 2553.721 | 263.400 | 9.60x | 1.048 |
| gnu | 4096 | g5-series | 1324.610 | 174.881 | 7.48x | 0.364 |
| nvhpc | 512 | lj | 0.972 | 0.997 | 0.97x | 0.756 |
| nvhpc | 512 | g4 | 26.607 | 16.532 | 1.61x | 1.493 |
| nvhpc | 512 | g5 | 22.660 | 14.183 | 1.60x | 0.778 |
| nvhpc | 512 | behler | 50.936 | 32.430 | 1.57x | 0.727 |
| nvhpc | 512 | g4-series | 262.638 | 28.273 | 9.42x | 0.958 |
| nvhpc | 512 | g5-series | 155.569 | 19.857 | 7.90x | 0.332 |
| nvhpc | 4096 | lj | 7.829 | 7.981 | 0.98x | 0.756 |
| nvhpc | 4096 | g4 | 213.825 | 132.527 | 1.61x | 1.490 |
| nvhpc | 4096 | g5 | 180.453 | 112.604 | 1.60x | 0.762 |
| nvhpc | 4096 | behler | 408.257 | 259.069 | 1.61x | 0.710 |
| nvhpc | 4096 | g4-series | 2116.955 | 226.072 | 9.38x | 0.957 |
| nvhpc | 4096 | g5-series | 1246.897 | 159.004 | 7.96x | 0.334 |

The independent reference can also vary across runs. For example, the GNU
4096-atom four-input G4 reference was slower in the final run than in the initial
campaign, so its normalized common speedup exceeds the ratio of raw common
times. Both are shown; this is not evidence that all of that difference comes
from this optimization. The per-run final/reference comparison remains explicit.

### GPU synchronous batches

| GPU | Atoms | Model | Before (ms) | Final (ms) | Speedup |
|---|---:|---|---:|---:|---:|
| h100 | 512 | lj | 0.342 | 0.344 | 0.99x |
| h100 | 512 | g4 | 1.530 | 1.403 | 1.09x |
| h100 | 512 | g5 | 1.144 | 1.087 | 1.05x |
| h100 | 512 | behler | 1.832 | 1.638 | 1.12x |
| h100 | 512 | g4-series | 3.934 | 1.573 | 2.50x |
| h100 | 512 | g5-series | 2.415 | 1.367 | 1.77x |
| h100 | 4096 | lj | 1.121 | 1.114 | 1.01x |
| h100 | 4096 | g4 | 3.246 | 2.808 | 1.16x |
| h100 | 4096 | g5 | 3.047 | 2.675 | 1.14x |
| h100 | 4096 | behler | 4.733 | 4.196 | 1.13x |
| h100 | 4096 | g4-series | 13.732 | 6.584 | 2.09x |
| h100 | 4096 | g5-series | 11.935 | 5.610 | 2.13x |
| blackwell | 512 | lj | 0.407 | 0.407 | 1.00x |
| blackwell | 512 | g4 | 4.212 | 3.879 | 1.09x |
| blackwell | 512 | g5 | 2.692 | 2.467 | 1.09x |
| blackwell | 512 | behler | 4.986 | 4.388 | 1.14x |
| blackwell | 512 | g4-series | 11.576 | 2.689 | 4.31x |
| blackwell | 512 | g5-series | 5.087 | 1.893 | 2.69x |
| blackwell | 4096 | lj | 1.196 | 1.224 | 0.98x |
| blackwell | 4096 | g4 | 11.211 | 8.966 | 1.25x |
| blackwell | 4096 | g5 | 6.896 | 5.656 | 1.22x |
| blackwell | 4096 | behler | 16.274 | 13.002 | 1.25x |
| blackwell | 4096 | g4-series | 50.243 | 11.994 | 4.19x |
| blackwell | 4096 | g5-series | 19.600 | 7.901 | 2.48x |

### Remaining CPU gap and interpretation

The four-output G4 cases still favor the established CPU implementation. Its
descriptor stage computes an unordered pair once, keeps both derivatives, and
reuses radial/geometry factors across angular parameters. The shared no-Jacobian
CPU path visits pairs in the value stage and again after NN backpropagation for
forces. With few outputs, coefficient contraction saves too little work to fully
offset that repeated geometry/radial evaluation. With a wider shared angular
basis, it avoids many output-specific gradient vectors and the common path closes
the G4 gap while substantially accelerating G5. These operation-count differences
are visible in the source; the timings do not identify every compiler/cache cost.

The default non-Chebyshev CPU dispatch is therefore retained. This is a measured
conditional performance result, not a general assertion that the CPU algorithm is
unsuitable for GPUs or that shared Fortran requires CPU threading.

## Validation

- GNU and NVHPC OpenMP-disabled serial tests: 10 each.
- H100: 31 GPU/host tests; Blackwell: 29 GPU tests.
- GNU bounds/runtime-checked build: 9 tests.
- Descriptor equivalence/finite-difference suite: **99 cases**, including sparse
  powers, degree 16, noninteger and high powers, near-integer inputs, collinear
  endpoint angles, species-map changes, radial-group split/merge, and coefficient
  workspace reuse. This is contained in the suites above, not an additional test count.
- Compute Sanitizer memcheck on H100: **0 errors** for those 99 cases (`--leak-check no`,
  as in previous campaigns; this is not a device-pool leak claim).
- LAMMPS: **63** LJ/Behler/n2p2 CPU/GPU comparisons, including 1/2 MPI ranks,
  orthogonal/triclinic cells, GPU neighbor options, empty-rank coverage, and NVE.
- Ordinary CPU: 42 numerical tests and 2 performance checks passed.
- Chebyshev direct/moment before/after spot checks: both CPU compilers and H100.
  An initial GNU 4096-atom direct result about 5% slower was not reproduced by
  the follow-up with reversed executable order; see `chebyshev-repeat-gnu`.
  Small single-campaign timing differences are not treated as universal changes.

The largest initial E/F/W absolute difference across the 48 final benchmark runs
was `4.275e-15`; all benchmark sample checks passed.

The checks use the independent established CPU implementation, plus coordinate
and strain finite differences. Agreement is numerical, not bitwise. Raw logs and
phase times, including full-output checks and workspace-residency assertions, are
retained. `source-sha256.json`, `implementation.diff`, `final-binaries.json`, and
`baseline-library-sha256.json` identify the working-tree implementation and builds.
The Git base hash alone does not identify these uncommitted changes.

## Reproduction

The preserved `before-src` files restore the revision-1.1 numeric backend in an
isolated copy of this working tree; retain the current benchmark/model builder
in both copies. Build CPU candidates with `ACCELNET_TARGET_SERIAL=ON` and without
OpenMP compiler flags. The archived scripts record exact local paths and commands;
use separate build directories and substitute your own paths. For example:

```sh
python3 AccelNetPredictor/benchmark/compare_chebyshev_variants.py \
  --variant before host /path/to/before/accelnet-target-benchmark \
  --variant after host /path/to/after/accelnet-target-benchmark \
  --family g4-series --orders 8 --modes 1 --sizes 512 4096 \
  --rounds 2 --seconds 0.12 --cpu 6 --output comparison-g4-series
```

For GPU measurements use `gpu` instead of `host` and select one GPU with
`CUDA_VISIBLE_DEVICES`. Both candidates still verify all energy/force/virial
components against the established CPU reference. The old diagnostic source-edit
experiments are archived; use this binary-comparison driver for current kernels.

## Local LAMMPS executable

The validated executable is installed at `/home/nagai/AccelNetGPU/lammps-accelnet-gpu/lmp`.
The previous Horner build is preserved as `lmp-before-generic-common`. The deployed
SHA256 is `9a7200fba689020f9f53381de9ed6dd22e973d3e9c24bb440ed9af4700f8b5b5`; `deployment.json` records both identities.
