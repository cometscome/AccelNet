# OpenMP target GPU batch evaluation

The optional `AccelNet::Target` library evaluates **Chebyshev, LJ, and Behler
G1–G5 descriptors, the neural network and its input gradient, force contraction,
and force/virial scatter in FP64**. Existing CPU entry points retain their
compiler options and do not acquire an OpenMP/GPU runtime dependency.

Supported configurations include one Chebyshev component per element (versions
0, 1 and 10), LJ/Behler components and their combinations, and different families
for different elements. All ten cutoff types, element-specific network shapes,
and all twelve activation codes are supported. G4/G5 use direct pair evaluation;
GPU G5 moments are not implemented and forced G5 moment modes are rejected.
Mixed/multiple Chebyshev components within a single element remain unsupported.
Angular parameters require finite values, |lambda| <= 1 and zeta >= 1.

The same target kernels can explicitly execute on the CPU:
`call packed%initialize(model, use_host=.true.)`. This option bypasses neither
validation nor numerical checks: initialization verifies the actual execution
location, and every target/data/update directive selects the requested backend.
The default still requires a GPU and rejects accidental fallback. In NVHPC 25.3,
`device(omp_get_initial_device())` alone did not select the CPU in our probe;
conditional `target if(...)` is therefore used as well.

Cutoffs, angular powers, and activation formulas have one source in
`AccelNetDescriptors/src/shared/*.inc`, included by the CPU and device modules.
The default CPU CSR batch API now uses a serial compilation of the same kernels
for single-component Chebyshev models, in both direct and moment modes. Other
CPU descriptor families keep their measured existing implementation.
`evaluate_batch_reference` retains the old batch path for independent comparisons;
the object/atomic APIs remain available as additional independent references.

Direct force contraction uses saved unit directions and a differentiated
Clenshaw recurrence for the Chebyshev series. This reduces normalization and
polynomial work without adding scratch arrays or GPU launches. It is shared
between the CPU and GPU instances. Direct still evaluates angular pairs from
both directed edges; computing each pair once is a separate parallel-layout
optimization, not part of this change.

The [initial validation report](validation/openmp-target-2026-09-25/README.md)
records the first GPU implementation. The
[residency/moment report](validation/gpu-residency-moments-2026-09-25/README.md)
records the optimized implementation, numerical/memory tests, and phase timings.
The [moment construction/force report](validation/gpu-moment-force-2026-09-26/README.md)
records the subsequent per-atom work ordering and force coefficient/gradient
contraction optimization. Raw moment scratch is reused for force coefficients
after descriptor evaluation and rebuilt on the next call; no new resident
buffers are needed. Direct/moment species mixtures are covered by finite
differences and workspace reuse tests.

The [descriptor/CPU comparison](validation/gpu-descriptors-cpu-comparison-2026-09-26/README.md)
records the earlier OpenMP-disabled measurements and G4/G5/LJ diagnostic experiments.
The [Chebyshev shared-kernel report](validation/chebyshev-common-2026-09-26/README.md)
records the default CPU batch switch, direct arithmetic optimizations, and
before/after CPU, H100, and LAMMPS timings.

The [generic common-kernel optimization report](validation/generic-common-2026-09-26/README.md)
records persistent grouped G4/G5 radial caching, value-only helpers, and fused LJ
components. These improvements apply to both serial CPU and GPU instances; the
default non-Chebyshev CPU batch still uses the independent established CPU path. The
[Behler contraction report](validation/behler-contraction-2026-09-26/README.md)
extends this with angular coefficient aggregation fused into the NN stage,
differentiated Horner for shared integer powers, and local grouped value sums.
Singleton powers retain their direct evaluation. CPU builds remain OpenMP-disabled;
GPU directed-edge force ownership avoids atomic updates used by the tested
unordered-pair variant. Both paths use the same scalar formulas.
The [G4 follow-up report](validation/g4-common-2026-09-27/README.md) adds
equal-power value reuse, squared-distance rejection before expensive pair work,
explicit packed-argument extents, and scalar force contraction. These changes
remain in the common source; the report includes a four-input model with no
duplicate descriptors as a control. The default CPU LJ/Behler dispatch is unchanged.
The subsequent [G4 value/Jacobian comparison](validation/g4-fused-2026-09-27/README.md)
changes common G4 evaluation to the established CPU traversal: each unordered
pair contributes values and both derivatives once, and forces read the saved
Jacobian after the NN. CPU and GPU compile that same numerical loop. G5 keeps
coefficient contraction; earlier reports describe the former G4 implementation.

## Build and test on H100

```sh
cmake -S . -B build-target \
  -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/bin/nvfortran \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
  -DACCELNET_BUILD_OPENMP_TARGET=ON \
  -DACCELNET_OPENMP_TARGET_FLAGS="-mp=gpu -gpu=cc90,cc120"
cmake --build build-target --parallel --target test_batch_target accelnet-target-benchmark
# List GPU UUIDs, then choose the intended physical device explicitly.
nvidia-smi --query-gpu=name,uuid --format=csv
export CUDA_VISIBLE_DEVICES=GPU-2644154d-7268-af42-6631-59e1f3c6e7f3
ctest --test-dir build-target -R predictor_target --output-on-failure
```

The target flags apply to the optional GPU library and its consumers' link
steps only. Use an appropriate compiler and offload flags on other platforms.
The example embeds H100 (`cc90`) and Blackwell (`cc120`) code; `-gpu=cc90`
suffices for H100 alone. Replace the example H100 UUID with your device's UUID.
CUDA ordinal order can differ from `nvidia-smi` index order on mixed-GPU hosts.
GNU Fortran 11 compiles this source with `-fopenmp`, but that alone does not
enable GPU execution: a working offload compiler/runtime and device are required.
The default backend rejects missing GPUs, invalid device IDs and CPU fallback.
Explicit `use_host=.true.` also works with a GNU `-fopenmp` build without a GPU.

The default self-contained test checks synthetic models. Set
`ACCELNET_PREDICTOR_GOLDEN_DIR` to the directory containing `Ti.nn.ascii`,
`O.nn.ascii`, `structure0001.xsf` and `structure2935.xsf` to also enable the real
Ti/O model test. This directory has the same meaning as in the CPU suite.

## API and lifetime

```fortran
use accelnet_batch_target, only: target_model, target_workspace, evaluate_batch_target
type(target_model) :: gpu_model
type(target_workspace) :: gpu_work

! model is a loaded CPU predictor_model. Pack a snapshot once after loading.
call gpu_model%initialize(model)          ! optional device= selects OpenMP device ID
! Optional mode=: Chebyshev 0=auto (default config mode), 1=direct, 2=moment.
! For an explicit comparison on CPU, initialize with use_host=.true..
forces = 0.0_real64
virial = 0.0_real64
call evaluate_batch_target(gpu_model, species, centers, offsets, indices, &
    displacements, energies, forces, gpu_work, virial)
! The call is synchronous: outputs are ready on return.
call gpu_work%release()
call gpu_model%release()
```

The CSR indexing, periodic-image displacements, ghost atom indices, row subsets,
and additive force/virial semantics are those of the [CPU batch API](batch-api.md).
Atomic energies are overwritten per row. Virial is optional and uses the same
`W(a,b) = sum dr(a)*F_neighbor(b)` convention. No image deduplication occurs.
Empty batches are no-ops; isolated atoms are supported. For multi-element version
10, the historical lookup of the center species through the neighbor array is
preserved, including its minimum-neighbor-count requirement.

`target_model` owns an independent packed **snapshot**, including normalization,
species mappings and the CPU-generated moment coefficients. Call `initialize`
again after changing/reloading the CPU model. Copies own their arrays independently.
Initialization probes actual GPU execution once; ordinary evaluations avoid that
extra probe kernel.

`target_workspace` owns **persistent GPU mappings**, a private model cache and
host buffers. The first evaluation uploads the model and allocates buffers.
Subsequent calls compare packed model values with that cache, so copied models
and same-shape reloads cannot silently reuse stale weights. Only a changed model
is uploaded again. `uploads()` counts model uploads; `allocations()` counts
buffer growth events. Model comparison/validation cost is included in profiling.
Changing the OpenMP device releases the old device's mappings first.

Each evaluation transfers current species/CSR/displacements and downloads atomic
energies plus force/virial contributions. The API adds contributions to the
caller's accumulators on the CPU. Row subsets pack and transfer only their used
edge range, with rebased offsets; ghost force indices still address the full
species array. Intermediate descriptors, moments, NN activations and gradients
stay on the GPU. Current geometry is always uploaded; there is no assumption of
an unchanged neighbor list or coordinates.

Buffers grow as needed and do not shrink automatically. `release()` deletes GPU
mappings before deallocating host memory and resets counters. Scope finalization
also releases mappings. **Assigning a workspace creates an empty scratch cache**,
not a copy of GPU ownership; its next evaluation allocates independently. Use
one workspace per independently executing caller. Releasing a `target_model`
does not release a workspace's cached copy: release the workspace to reclaim its
GPU memory. No process-global model cache is introduced.

The evaluation stages are:

1. Compute neighbor geometry/radial descriptors and choose the angular method.
   For G4/G5 rows, fill shared radial-group values/derivatives in the same stage.
2. In moment mode, construct moments in parallel over both atoms and monomials,
   then form angular descriptors using the CPU model's polynomial coefficients.
3. For G4, accumulate descriptor values and edge Jacobians together in one
   center/pair/descriptor traversal, reusing cutoff, exponential and angular factors.
4. Evaluate NN outputs and input gradients; transform G5/Chebyshev angular coefficients.
5. Contract derivatives in parallel over neighbor edges, reusing cached geometry,
   and moments. In moment mode, nested x/y/z Horner recurrences evaluate the
   contracted polynomial and all three Cartesian derivatives together. Store
   three force components per edge. G4 reads saved Jacobians without another pair traversal.
6. Aggregate edge contributions by center and scatter forces/virial with FP64 atomics.

G4 adds a persistent `(3, descriptor, edge)` Jacobian and per-center pair caches.
They remain on the device and are allocated/reused with the workspace. Other
families do not materialize this Jacobian. Main scratch storage also scales
with rows × network size, rows × number of moments, and edges × angular order
(for cached coordinate powers). Generic angular rows also use a persistent
edge × radial-group × 2 cache for values and derivatives. G5 angular force coefficients
reuse existing NN-gradient slots, with no new global coefficient buffer. Polynomial
degree is bounded at 16; higher/noninteger powers keep their original evaluation.
CSR batching bounds this memory. Chebyshev row kernels use
32-thread teams, while G4 distributes centers and monomial/edge kernels expose more parallel work to the
runtime. These are portable OpenMP constructs, with performance validated here
on NVIDIA GPUs only.

Moment metadata is packed in lexicographic x/y/z order so the force recurrence
reads coefficients consecutively without a per-term lookup. This retains the
original x/y summation order within each total degree when forming descriptors.
No extra coefficient buffer or target launch is introduced. Cached direction
powers are still used to construct the moments, but the force recurrence no
longer reads them. CPU serial and GPU builds use the same source. See the
[Horner measurements](validation/moment-horner-2026-09-26/README.md) for the
initial indirect prototype, CPU comparisons, and GPU/LAMMPS results.

`initialize(..., mode=1)` forces pair enumeration and `mode=2` forces moments.
Auto mode chooses moments when `angular_neighbors * (angular_order+1)` is at
least the number of monomials, `(p+1)*(p+2)*(p+3)/6`, with `p=angular_order`.
This initial work estimate was checked against direct/moment timing sweeps; it
is not a guarantee of the fastest method on every GPU/model. With no override,
the CPU config's `evaluation_mode` is retained (normally auto). Forced modes let
callers compare methods on their hardware. Original CPU kernels and mode policy
are unchanged.

## Correctness and performance checks

The device tests compare every atomic energy, force component and all nine
virial components with CPU batch evaluation, using
`2e-10 + 2e-10*abs(reference)`. Separate coordinate and strain finite differences
use GPU energy evaluations to check GPU forces and virial. Tests cover:

- versions 0/1/10, CPU direct/auto/moment references, order zero and higher orders;
- all cutoffs and activations, different element network topologies, normalization;
- nonperiodic, orthogonal and triclinic cells, repeated periodic self images;
- single/multiple species, empty/isolated systems, subsets, reordered centers;
- additive accumulators, omitted virial, workspace growth/reuse, model reinitialization;
- invalid inputs, unsupported models and rejection of CPU fallback.

Run the [CPU performance regression gate](batch-api.md#cpu-regression-gate) too.
Enabling the separate GPU library does not redirect existing CPU calls.

```sh
# NATOMS, polynomial ORDER, minimum seconds, optional MODE and lattice SPACING
OMP_NUM_THREADS=1 \
  taskset -c 6 build-target/bin/accelnet-target-benchmark 512 8 0.2 0 1.7
```

This benchmark verifies all results before timing and after measuring each
method, warms up each method, then alternates execution order across five
samples. It prints `CASE` (atoms, order, edge count, maximum absolute error) and
`TIMING` rows with four seconds per evaluation: CPU batch, GPU batch, CPU with
neighbor construction, GPU with neighbor construction. It also fails if model
upload/allocation counters change during steady-state measurement.

`PROFILE sample method` reports eight seconds-per-call fields: CPU neighbor
construction, host preparation/model comparison, input upload, GPU descriptors,
GPU NN, GPU force contraction/aggregation, output download/addition, and total
GPU API time. Descriptor/force phases include their multiple kernel launches.
These are synchronous host wall times including launch/synchronization overhead,
not isolated CUDA event timings. Other GPU processes can affect them. Initial
packing, first device allocation and model upload are excluded from warmed
measurements; data transfers remain included.

Profiling is also available directly through the optional `profile=` argument
of type `target_profile`, with fields `prepare`, `upload`, `descriptors`,
`network`, `forces`, `download`, `total` (all seconds). Full timings rebuild the
CPU neighbor list every call. Comparisons use a single CPU thread. Inspect
GPU identity and current utilization before drawing performance conclusions;
CUDA ordinal order and `nvidia-smi` index order can differ.

Run `test_batch_target --switch-device` with two visible devices to exercise
0 → 1 → 0 migration of one workspace. `--quick` is a smaller periodic case for
Compute Sanitizer; the full direct/moment suites also check copies, finalization
and repeated buffer growth. GPU data residency does not switch existing CPU,
CLI or LAMMPS callers to the GPU. Small systems may still favor the CPU.
AMD and Intel GPU execution must be validated on those devices before claiming
support.

## Reproduce the CPU/common-kernel comparison

Build `accelnet-target-benchmark` and use the same executable for both paths:

```sh
OMP_NUM_THREADS=1 python3 AccelNetPredictor/benchmark/compare_target_backends.py \
  --benchmark build-target/bin/accelnet-target-benchmark \
  --output comparison --cpu 6 --sizes 64 512 4096
```

The script alternates five warmed samples, pins CPU execution to one core,
compares every energy/force/virial component, and reports target-time/CPU-time.
It times the established CPU batch path and the identical target source on the
CPU and GPU, with resident CSR neighbors; GPU timings include transfers.
Neighbor construction and initialization are excluded. `--backends host` also
works with GNU Fortran plus `-fopenmp`. Measurements from different compilers
must be reported separately. `--descriptors` tests all new families, all ten
cutoffs, finite differences, descriptor/species-map reloads, and mixed families.
Tests near discontinuous hard cutoffs use distances away from roundoff ambiguity;
an exactly representable boundary and both sides are checked separately.

### Compile the identical kernels with OpenMP completely disabled

To separate OpenMP overhead from the numerical algorithm, configure a second
build with **no `-fopenmp` / `-mp` compiler or link flag**:

```sh
cmake -S . -B build-target-serial -DCMAKE_BUILD_TYPE=Release \
  -DACCELNET_BUILD_OPENMP_TARGET=ON -DACCELNET_TARGET_SERIAL=ON
cmake --build build-target-serial --target test_batch_target accelnet-target-benchmark
ctest --test-dir build-target-serial -L target-host --output-on-failure
python3 AccelNetPredictor/benchmark/compare_target_backends.py \
  --benchmark build-target-serial/bin/accelnet-target-benchmark \
  --backends host --output serial-comparison --cpu 6
```

The numerical `.f90` files are unchanged; all `!$omp` directives are ignored.
A small runtime facade replaces timing/device queries with `system_clock` and
host constants. GNU's resulting executable does not link an OpenMP runtime.
NVHPC also links `libnvomp` through its standard Fortran runtime, but the kernel
objects contain no OpenMP parallel/offload calls. GPU-required tests are disabled
in this serial build. Use the same compiler and optimization flags for each pair.
Timings exclude one-time model initialization, and **include per-step region
entry, synchronization, workspace checks, and input/output copies**.

### Diagnose rather than assume why a path is slow

The recipe below is historical: it applies to the source before the persistent
generic radial cache and fused LJ implementation. The diagnostic script rejects
the current kernel ABI rather than mixing incompatible sources. Archived sources
and measurements are in the [original diagnostic report](validation/gpu-descriptors-cpu-comparison-2026-09-26/README.md).
For the current common code, compare saved binaries with the same compiler:

```sh
python3 AccelNetPredictor/benchmark/compare_chebyshev_variants.py \
  --variant before host /path/to/saved-before \
  --variant after host build-target-serial/bin/accelnet-target-benchmark \
  --family g4 --modes 1 --sizes 512 4096 --rounds 2 --output comparison-g4
```

Use `--family lj`, `g5`, or `behler` for other families. For GPU measurements use
`gpu` instead of `host`, a GPU build, and the intended `CUDA_VISIBLE_DEVICES` UUID.
The [new report](validation/generic-common-2026-09-26/README.md) records matched
CPU and GPU results. Historical experiment recipe:

The phase profile separates descriptors, NN, and forces. Controlled CPU-only
experiments preserve all energy/force/virial checks and keep production sources
unchanged:

```sh
python3 AccelNetPredictor/benchmark/diagnose_target_cpu.py \
  --build build-target-serial --output experiment-cached --variant cached
OMP_NUM_THREADS=1 experiment-cached/test_batch_target --host --descriptors
python3 AccelNetPredictor/benchmark/compare_target_backends.py \
  --benchmark experiment-cached/accelnet-target-benchmark --backends host \
  --sizes 512 --families g4 g5 behler --output experiment-cached/timing
```

`cached` reuses angular radial values/derivatives before the pair loops.
`lj-powers` additionally replaces LJ integer powers with multiplications;
`lj-fused` instead calculates both LJ features together and shares their cutoff
and powers. These are **serial diagnostic implementations**, not device kernels
or changes to the default CPU path. The output contains the exact modified
sources, library, numerical test and timing executable.
