# Batch evaluation and CPU performance regression checks

This batch API executes on the CPU. Supported Chebyshev, LJ and Behler G1–G5
models use the same numerical source as the
[OpenMP target backend](openmp-target.md), compiled without OpenMP directives.
Chebyshev and G5 support direct/moment; G4 uses direct pairs. Multiple/mixed
Chebyshev components use the common component pipeline. Production inference
has no legacy fallback. G5 auto retains its per-component 16-neighbor threshold
and exact integer orders 1--10; explicit moment modes also permit orders 11--16.
Fractional/near-integer and orders above 16 remain direct.
Structure/file, Fortran/C atomic, CLI and ordinary LAMMPS CPU calls use the
common serial kernels. No OpenMP runtime or GPU compiler is required for a
normal CPU build. Old code is in [legacy/cpu-reference](../legacy/cpu-reference/README.md).

## API

```fortran
use accelnet_batch, only: batch_workspace, evaluate_batch
type(batch_workspace) :: work

! model is a loaded predictor_model; neighbors is a current neighbor_data.
! Allocate centers(natoms), energies(natoms), forces(3,natoms) in the caller.
centers = [(i, i=1,structure%natoms)]
forces = 0.0_real64
virial = 0.0_real64
call evaluate_batch(model, structure%species, centers, neighbors%offsets, &
    neighbors%atom_indices, neighbors%displacements, energies, forces, work, virial)
total_energy = sum(energies)
```

- All indices are **one-based**. Row `r` evaluates central atom `centers(r)`;
  edges are `offsets(r):offsets(r+1)-1`. Offsets must be nondecreasing.
- `species` covers all atoms, including any ghost atoms. `indices` addresses
  this array and the columns of `forces`. Species use the model's global order.
- `displacements(:,edge)` is **neighbor image position minus central position**.
  The caller supplies the current geometry, all required neighbors, and periodic
  image translations. This routine does not build/reuse a Verlet list, apply a
  minimum-image convention, or validate geometric completeness.
- Repeated targets, including periodic self images, are valid and must not be
  deduplicated. Virial is accumulated from image displacements before force
  contributions are folded onto atom indices.
- `energies(r)` is overwritten with the physical atomic energy, including the
  energy normalization, shift and atomic reference. Each center is evaluated
  once per occurrence; duplicate centers intentionally count more than once.
- `forces` and optional `virial` are **added to**, like the atomic API. Initialize
  them once before evaluating one or several batches. Virial has the existing
  `W(a,b) = sum dr(a)*F_neighbor(b)` convention, in energy units, with no factor
  1/2 or volume division.
- Rows can be partitioned without copying/rebasing the complete edge arrays:
  pass `centers(first:last)` and `offsets(first:last+1)`. Forces on atoms outside
  the central subset are still included. This also supports multiple independent
  structures packed into a single index space (without cross-structure edges).
- Zero rows and zero neighbors are supported. Invalid CSR sizes, indices and
  species stop with an error, consistent with the object API.

Keep one `batch_workspace` per independently executing caller. It retains
descriptor, gradient, species and force buffers and grows as necessary. Model
metadata is refreshed at each call, so a workspace can be reused after reloads.
`work%reserve(model, maximum_neighbors)` can preallocate the reference-path buffers;
shared buffers grow on the first evaluation, when the CSR size is known.
`work%allocations()` counts **growth events**, not individual allocator calls;
`work%release()` frees the buffers and resets that counter. It does not own the
model or neighbor list and is not safe for simultaneous use by multiple threads.

On the common path, G4 computes values and derivatives in one pair traversal,
then contracts the saved Jacobian after the NN evaluation. CPU and GPU use the
same formulas with different descriptor-owner counts. Other common descriptors
use cached geometry and direct derivative contraction. G5 moments share cached radial groups, reuse raw
moments after the NN, and evaluate contracted force polynomials with differentiated
Horner. Auto follows the established neighbor threshold; this policy is not a
guarantee of the fastest choice for every model or neighbor density.
The separate GPU API uses packed arrays, persistent device buffers and model
parameters, and GPU direct/moment kernels. See the [GPU API](openmp-target.md)
for lifetime rules and phase profiling.

`evaluate_batch_reference` has the same arguments as `evaluate_batch` and
retains the former CPU implementation for independent correctness and speed
comparisons. The default batch path repacks model metadata each call
to honor edits/reloads, while retaining scratch buffers. A prepared
`target_model` in an `ACCELNET_TARGET_SERIAL=ON` build can avoid repeated packing
when the caller explicitly manages the model snapshot.

CMake generates CPU module instances from the target source, changing module
names and removing OpenMP directives; there is no separately maintained CPU
copy of these kernels. This also prevents OpenMP region startup when the user
supplies global OpenMP flags. GPU builds can link both instances in one program.
See the [Chebyshev validation report](validation/chebyshev-common-2026-09-26/README.md)
for the earlier Chebyshev switch and direct arithmetic optimizations. The
[LJ/Behler default-dispatch report](validation/common-cpu-default-2026-09-27/README.md)
records the later switch, including model packing and G5 auto-mode comparisons.

## Correctness tests

Build `test_batch` with `BUILD_TESTING=ON`, then run:

```sh
ctest --test-dir build -R predictor_batch --output-on-failure
```

The tests compare energies, every force component and all nine virial components
against the existing structure API. They exercise Chebyshev versions 0/1/10,
direct/moment modes, LJ, composite models, n2p2 G4/G5, different per-element NN
widths/depths, normalization, isolated atoms, orthogonal/triclinic periodic cells,
repeated images, empty batches, central subsets, reordered rows, additive output,
workspace growth/reuse/release, and model changes. Independent coordinate finite
differences check forces. The standard suite retains the existing virial strain
finite-difference tests. Run a bounds-checked GNU build as well.

## CPU regression gate

To guard the default common CPU dispatch against the retained structure
evaluator, configure `-DACCELNET_TEST_COMMON_CPU_PERFORMANCE=ON`, build
`accelnet-cpu-regression-benchmark`, and run:

```sh
ctest --test-dir build -R '^predictor_common_cpu_performance$' --output-on-failure
```

This opt-in test uses the same executable/compiler for both paths, includes
neighbor construction, checks all energy/force/virial components, and compares
seven alternating single-core samples. The default time-ratio limit is 1.10,
configurable with `ACCELNET_CPU_MAX_SLOWDOWN`. Use a build without OpenMP flags
and an otherwise idle machine. It covers Chebyshev, LJ and real n2p2 G4/G5
fixtures at 8, 64 and 512 atoms. For fixed-neighbor measurements including model
packing, `accelnet-target-benchmark ... cpu-shared FAMILY no-neighbors` compares
the public default batch API against `evaluate_batch_reference`; generic mode 0
also includes the reference's automatic G5 moment selection.

`accelnet-cpu-regression-benchmark` is built when predictor testing is enabled.
Its synthetic fixtures require no external training corpus. Benchmark cases use
Chebyshev, LJ, n2p2 G5 and G4 models with 8, 64 and 512 atoms. Both evaluation
modes include neighbor construction, energies, forces and virial. File parsing
and initialization are outside the timed region. This is an end-to-end CPU
check, not a measurement of GPU kernels or the batch evaluator alone.

Create an immutable reference executable **before** changing CPU evaluation.
The same current fixture/driver sources can be linked against the original
`c663146` library because the baseline driver uses only its existing API:

```sh
# Run from the repository root. Choose a new, empty directory.
baseline_dir=$(mktemp -d /tmp/accelnet-baseline.XXXXXX)
mkdir "$baseline_dir/source"
git archive c663146 | tar -x -C "$baseline_dir/source"
cmake -S "$baseline_dir/source" -B "$baseline_dir/build" \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build "$baseline_dir/build" --parallel --target AccelNetPredictor
gfortran -O3 -cpp -J "$baseline_dir" \
  -I "$baseline_dir/build/AccelNetPredictor/modules" \
  -I "$baseline_dir/build/AccelNetDescriptors/modules" \
  AccelNetPredictor/test/batch_test_support.f90 \
  AccelNetPredictor/benchmark/cpu_regression_benchmark.F90 \
  "$baseline_dir/build/lib/libaccelnet.a" \
  "$baseline_dir/build/lib/libAccelNetDescriptors.a" \
  -o "$baseline_dir/cpu-benchmark"

cmake -S . -B build -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DBUILD_SHARED_LIBS=OFF \
  -DACCELNET_CPU_BASELINE_EXECUTABLE="$baseline_dir/cpu-benchmark"
cmake --build build --parallel --target accelnet-cpu-regression-benchmark
ctest --test-dir build -R predictor_cpu_performance_regression --output-on-failure
```

Use the same compiler/version/options for both libraries and drivers. Rebuild
the reference driver if fixtures change. For other compilers, use that compiler
and its equivalent preprocessing/module flags throughout.

The opt-in CTest test uses one pinned CPU, seven samples of at least 0.2 seconds
per case, warm-up, alternating baseline/candidate execution, and the median of
adjacent candidate/baseline wall-time ratios. Pairing reduces bias when host
load/frequency changes during a run. It fails if this ratio exceeds `ACCELNET_CPU_MAX_SLOWDOWN` (default
1.10) in any case, or any energy/force/virial component changes beyond
`2e-10 + 2e-10*abs(reference)`. Raw timings and errors are saved to
`build/AccelNetPredictor/cpu-performance.json`. Performance tests are serial and
are opt-in because shared CI runners are noisy; correctness runs in normal CI.
The 10% threshold is a noise allowance, not a promised acceptable slowdown.

To also compare the batch path against the original structure path:

```sh
python3 AccelNetPredictor/benchmark/compare_cpu_performance.py \
  --baseline "$baseline_dir/cpu-benchmark" \
  --candidate build/bin/accelnet-cpu-regression-benchmark \
  --candidate-mode batch --output batch-performance.json
```

`--cpu`, `--samples`, `--seconds`, `--sizes`, and `--families` control the run.
Check host load/frequency and repeat a noisy result before attributing it to code.
Keep raw failed measurements too. These four small synthetic models establish
a reproducible initial gate; representative production models and larger scales
must be added as the GPU implementation develops.

## Optional OpenMP target environment check

This checks actual device execution, array mapping and FP64 atomic scatter. It
fails on CPU fallback. It does **not** turn on GPU evaluation of AccelNet and its
compiler flags apply only to the smoke-test executable, not the CPU library.
For NVIDIA HPC SDK and an H100, for example:

```sh
cmake -S . -B build-nvhpc -DCMAKE_Fortran_COMPILER=nvfortran \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
  -DACCELNET_TEST_OPENMP_TARGET=ON \
  -DACCELNET_OPENMP_TARGET_TEST_FLAGS="-mp=gpu -gpu=cc90"
cmake --build build-nvhpc --target test_openmp_target
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=0 \
  ctest --test-dir build-nvhpc -R openmp_target_device_smoke --output-on-failure
```

Use suitable flags for other compilers/devices. The environment must allow
access to the GPU driver and device nodes.
