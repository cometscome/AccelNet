# AccelNetPredictor

Fortran inference library for supported ænet and n2p2 potentials. Energies,
analytic forces and virials use the same maintained CPU/GPU numerical kernels.
Structure/file, atomic Fortran/C, ordinary CSR batch and ænet-compatible SFB
calls use their serial compilation, with OpenMP removed. The optional target
library runs the common code on an explicitly selected GPU or threaded host.

See [current descriptor/moment coverage](../docs/implementation-status.md),
[model compatibility](../docs/model-compatibility.md), the
[CSR batch API](../docs/batch-api.md), and the
[OpenMP target API](../docs/openmp-target.md). Multiple Chebyshev blocks and
mixed descriptor families are supported. The independent former evaluator is
available as `evaluate_batch_reference` in
[legacy/cpu-reference](../legacy/cpu-reference/README.md); production inference
has no legacy fallback. Standalone low-level descriptor/NN compatibility APIs
still use retained reference implementations; the common-code claim applies
to production potential inference, not every utility or API wrapper.

## CPU build and command-line prediction

Commands in this README assume the working directory is `AccelNetPredictor/`.
The standalone build writes executables directly under `build/`. A build from
the repository root writes them under `build/bin/`; see the
[root quick start](../README.md#cpu-quick-start). Build the descriptors with the
same Fortran compiler as the predictor.

The default is serial CPU evaluation. Requirements are CMake 3.20+, a Fortran
2008 compiler and a C compiler; GNU Fortran 11/13 are covered by CI. No GPU,
CUDA, MPI or BLAS installation is needed for this library build. Some tests
use Python; optional upstream n2p2 C++ reference tests need a C++ compiler.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
  -DACCELNET_BUILD_OPENMP_TARGET=OFF -DN2P2_SCALING_EXECUTABLE=
cmake --build build --parallel
ctest --test-dir build -LE "gpu|performance" --output-on-failure
```

Set `N2P2_SCALING_EXECUTABLE` to a compatible local `nnp-scaling` executable
to enable the optional upstream descriptor comparison.

The ordinary CLI and object/atomic/SFB/batch APIs remain serial even when the
optional target library is built. `OMP_NUM_THREADS` does not parallelize those
entry points. See the separate [CPU OpenMP](#cpu-openmp-threading) and
[GPU](#gpu-build-and-execution) sections for explicit backend selection.

Energy-only CLI:

```sh
build/accelnet-predict 2 Ti.fingerprint.stp O.fingerprint.stp \
  Ti.nn.ascii O.nn.ascii structure.xsf
```

Standard n2p2 2G-HDNNP model directories can also be loaded directly. The
directory must contain `input.nn` and the conventional `weights.%03d.data`
files (where the number is the element atomic number), plus `scaling.data`
when required by the selected scaling mode:

```sh
build/accelnet-predict --n2p2 /path/to/model structure.xsf
```

From Fortran, use `load_predictor_from_n2p2(model_directory, model)` or reload
an existing object with `call model%reload_from_n2p2(model_directory)`.
Imported models currently support 2G symmetry-function types 2/3/9/12/13/20--25,
all n2p2 `cutoff_type` values 0 through 8 and AccelNet's fractional extension
as type 9, with `0 <= cutoff_alpha < 1` (`cutoff_alpha > 0` for type 9),
all n2p2 activation functions,
the standard scaling modes, data-set energy normalization, and atomic energy
offsets, per-element network topologies, and `normalize_nodes`.
4G/charge models are rejected with an error rather than
evaluated with different semantics. XSF coordinates and
the parameters in `input.nn` must use the same physical length unit; returned
energies and forces use the model's physical energy and length units.
The loader accepts `nnp_type` values `2G`, `2G-HDNNP`, and numeric `2`.
Normalized models must provide `mean_energy`, `conv_energy`, and `conv_length`
together. Type-3/type-9 angular radial shifts are retained in row 7 of
AccelNet's extended descriptor metadata.

For example, global topology defaults can be overridden for one element:

```text
global_hidden_layers_short 2
global_nodes_short 20 20
global_activation_short t t l

element_hidden_layers_short O 2
element_nodes_short O 30 20
element_activation_short O p t l

normalize_nodes
```

`normalize_nodes` divides the weighted sum and bias feeding a layer by the
number of nodes in the preceding layer. The loader applies the equivalent
transformation to that layer's weights and biases when the model is loaded.

Optional integration tests also run n2p2 v2.3.0 `nnp-scaling` and compare its raw
symmetry-function output with AccelNet atom by atom. They cover types 2, 3,
and 9 for both elements and every n2p2 cutoff type from 0 through 8. A separate
reference fixture checks per-element widths and activation functions together
with `normalize_nodes`, comparing the energy and every force component with
upstream n2p2 v2.3.0 values. Another fixture covers an element-specific hidden
layer count. Weighted/compact and multi-element comparisons are recorded in
the [n2p2 extension report](../docs/validation/n2p2-extensions-2026-09-27/README.md).
Types 13/21/24 have exact direct implementations; lack of a finite moment
representation is not missing descriptor support.

The aenet-style atomic API can load the same directory without conversion:

```fortran
use accelnet
call accelnet_init_n2p2("/path/to/model", stat)
```

The C equivalent is `accelnet_init_n2p2(directory, &stat)`. Applications that
already called `accelnet_init` with the model's canonical species order can
instead call `accelnet_load_n2p2`. Both entry points publish the usual
`accelnet_Rc_max`, atom conversion, atomic-energy, and force APIs.

n2p2 structure files with one or more `begin`/`end` blocks are accepted by:

```sh
build/accelnet-predict --n2p2-data /path/to/model input.data
```

The public Fortran reader is `read_n2p2_data`. It reads atom coordinates,
elements, and either zero or three lattice vectors. Reference energies,
charges, stored forces, and comments are not part of `atomic_structure` and
are ignored.

The original AccelNet `predict.in` sections are accepted directly. Descriptor
setups are reconstructed from the fingerprint metadata embedded in each NN.
Both `*.nn.ascii` and AccelNet native sequential-unformatted binary NN files
are accepted:

```sh
build/accelnet-predict predict.in
```

For Chebyshev networks, `predict.in` can select the compatibility convention
with a new optional `VERSION` tag. The section form is recommended:

```text
NETWORKS
Ti Ti.nn.ascii
O  O.nn.ascii

VERSION
0

FILES
1
structure.xsf
```

`VERSION 1`, `VERSION=10`, and the more explicit `CHEBYSHEV_VERSION` spelling
are also accepted. Valid values are 0, 1, and 10. If the tag is omitted, the
default is `version=0`, matching the newer ænet v2.04 Chebyshev convention.
The tag has no effect on LJ or Behler2011 networks.

Append `--forces` only to the explicit `NSPECIES SETUP... NETWORK... XSF`
form to print analytic Cartesian forces. `predict.in`, `--n2p2` and
`--n2p2-data` already print forces and do not accept that extra flag.
A warm-up-separated whole-inference benchmark is also installed:

```sh
build/accelnet-predict-benchmark 2 Ti.fingerprint.stp O.fingerprint.stp \
  Ti.nn.ascii O.nn.ascii structure.xsf 1000
```

The benchmark can also load an n2p2 model directory directly. This compares
the Fortran evaluation using exactly the same `input.nn`, `scaling.data`, and
weights as n2p2:

```sh
build/accelnet-predict-benchmark --n2p2 /path/to/model structure.xsf 1000
```

For a controlled type-9 comparison in the same binary, append `direct`,
`moment`, or `auto`:

```sh
build/accelnet-predict-benchmark --n2p2 /path/to/model structure.xsf 1000 direct
build/accelnet-predict-benchmark --n2p2 /path/to/model structure.xsf 1000 moment
```

Model loading and XSF parsing are outside the reported in-memory API timings.

## CPU OpenMP threading

Build the optional `AccelNet::Target` library with GNU Fortran host OpenMP.
These commands still start in `AccelNetPredictor/`, so test executables are
directly inside the build directory:

```sh
cmake -S . -B build-openmp -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
  -DACCELNET_BUILD_OPENMP_TARGET=ON -DACCELNET_TARGET_SERIAL=OFF \
  -DACCELNET_OPENMP_TARGET_FLAGS="-fopenmp -foffload=disable -ffree-line-length-none" \
  -DN2P2_SCALING_EXECUTABLE=
cmake --build build-openmp --parallel
OMP_NUM_THREADS=1 build-openmp/test_batch_target --host --quick
OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close \
  build-openmp/test_batch_target --host --quick
OMP_NUM_THREADS=1 ctest --test-dir build-openmp -L target-host --output-on-failure
```

Link `AccelNet::Target`, use `accelnet_batch_target`, and initialize the packed
model with `call packed%initialize(model, use_host=.true.)` before calling
`evaluate_batch_target`. See [declarations and lifetime](../docs/openmp-target.md#api-and-lifetime).
This selects threaded CPU execution without requiring a GPU. The default
initialization requires a GPU instead. `ACCELNET_TARGET_SERIAL=ON` disables
OpenMP and is for serial benchmarking, not threaded execution.

There is no threaded CLI switch or LAMMPS `accelnet/omp` pair style.
LAMMPS `accelnet` uses one serial common-kernel instance per MPI rank; follow
the [LAMMPS CPU guide](../interfaces/lammps/README.md#cpu-build-and-run).

## GPU build and execution

The validated NVIDIA configuration uses **NVIDIA HPC SDK 25.3 `nvfortran`**,
the SDK's **CUDA 12.8** runtime/toolkit, and a compatible NVIDIA driver.
H100 NVL and RTX PRO 6000 Blackwell have been tested. `nvcc` alone cannot
compile these Fortran OpenMP target kernels; ordinary GNU `-fopenmp` alone
is not this NVIDIA offload configuration. AMD/Intel GPU execution and other
Fortran offload toolchains are not validated here.

For an H100, from `AccelNetPredictor/`:

```sh
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3
cmake -S . -B build-gpu -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER="$NVHPC_ROOT/compilers/bin/nvfortran" \
  -DCMAKE_C_COMPILER=gcc \
  -DACCELNET_BUILD_OPENMP_TARGET=ON -DACCELNET_TARGET_SERIAL=OFF \
  -DACCELNET_OPENMP_TARGET_FLAGS="-mp=gpu -gpu=cc90" \
  -DACCELNET_DESCRIPTORS_BUILD_N2P2_REFERENCE=OFF \
  -DN2P2_SCALING_EXECUTABLE=
cmake --build build-gpu --parallel
nvidia-smi --query-gpu=name,uuid --format=csv
# Replace with the UUID of the intended device:
export CUDA_VISIBLE_DEVICES=GPU-REPLACE-WITH-YOUR-DEVICE-UUID
export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_THREADS=1
build-gpu/test_batch_target --quick
ctest --test-dir build-gpu -L gpu --output-on-failure
```

Use `cc120` for the tested Blackwell GPU or `cc90,cc120` for both architectures.
The recipe disables optional external n2p2 C++ reference tests, avoiding a
mixed host-C++/NVHPC link; bundled correctness tests remain enabled.
Adjust SDK paths and architecture flags to your hardware. Keep compiler
families in separate build directories and use matching Fortran module files.

Applications link `AccelNet::Target` and use the Fortran target API or
`accelnet_target.h`. `call packed%initialize(model)` selects GPU execution
and rejects unintended CPU fallback. The [target guide](../docs/openmp-target.md)
describes device selection, CSR input, workspaces and energy-only evaluation.
The usual `accelnet-predict` CLI has no GPU option.

The target C API also loads n2p2 directories directly with
`accelnet_target_create_n2p2`, and selects Chebyshev convention 0/1/10 with
`accelnet_target_create_versioned`. Existing constructors retain version 0.
`accelnet_target_get_species` exposes the model's species order. See
[C model loading](../docs/openmp-target.md#c-model-loading) for signatures,
error behavior and explicit CPU execution with `ACCELNET_TARGET_HOST`.
Callers using a different Fortran compiler can use `ISO_C_BINDING` at this
C ABI boundary without importing the GPU library's compiler-specific modules.
For shared libraries, link the C-only `AccelNet::TargetC` target from
`find_package(AccelNetC CONFIG REQUIRED)`. This target does not propagate
NVHPC/OpenMP flags to the caller. The [C-only build and installation guide](../docs/openmp-target.md#c-only-cmake-package)
also covers serial CPU builds without a GPU SDK. Embedded network paths may
be supplied in any order, and malformed n2p2 models return diagnostics to the
C caller without terminating its process.

LAMMPS GPU use additionally requires **NVHPC `nvc++`**, the CUDA toolkit and
`PKG_GPU=ON`, `GPU_API=cuda`, `GPU_PREC=double`. The
[LAMMPS GPU guide](../interfaces/lammps/README.md#gpu-build-and-run) gives the
separate build and `pair_style accelnet/gpu` input; enabling GPU support does
not redirect `pair_style accelnet` to GPU.

## Object API

The Fortran library also accepts an already-loaded structure, avoiding XSF
I/O. The force array is caller-owned and can be reused between calls:

```fortran
use iso_fortran_env, only: real64
use accelnet_descriptors, only: atomic_structure, read_xsf
use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2

type(predictor_model) :: model
type(atomic_structure) :: structure
real(real64) :: energy
real(real64), allocatable :: forces(:, :)

call load_predictor_from_n2p2("/path/to/model", model)
call read_xsf("structure.xsf", model%species_names, structure)
allocate(forces(3, structure%natoms))
call model%predict_energy(structure, energy)
call model%predict_energy_forces(structure, energy, forces)
```

Pass an optional fourth output to obtain the full configurational virial:

```fortran
real(real64) :: virial(3,3)
call model%predict_energy_forces(structure, energy, forces, virial)
! The XSF filename overload accepts the same optional output.
```

The structure/filename APIs overwrite `forces` and `virial` on each call.
The virial convention, units, and periodic-image handling are described below.

## aenet-style atomic API

The `accelnet` Fortran module provides the complete aenet-style global API for
codes that already evaluate one local atomic environment at a time. It
includes initialization/finalization, ASCII and binary potential loading,
load-state and cutoff metadata, atom-type conversion, free-atom energies,
atomic energies and additive forces, a built-in neighbor list, and the
standalone structural-fingerprint basis routines. The C ABI is declared by
`include/accelnet.h`.

Native ænet 2.0.4 activation codes 0 through 4 are evaluated as linear, tanh,
sigmoid, modified tanh, and twist respectively. Imported n2p2 softplus uses
the non-conflicting AccelNet extension code 11; legacy AccelNet-generated n2p2
ASCII networks that used code 3 are recognized and upgraded while loading.
Reference tests compare the values and first derivatives of codes 0--4 with
ænet's `ff_activate`, and all supported n2p2 activations with n2p2's
`NeuralNetwork` implementation over negative, zero, and positive inputs.

```fortran
use accelnet

character(len=2) :: species(2) = ["Ti", "O"]
integer :: stat

call accelnet_init(species, stat)
call accelnet_set_chebyshev_evaluation(ACCELNET_CHEBYSHEV_MOMENT, stat)
call accelnet_load_potential(1, "Ti.nn.ascii", stat, is_ascii=.true.)
call accelnet_load_potential(2, "O.nn.ascii", stat, is_ascii=.true.)
call accelnet_atomic_energy(coo_i, type_i, n_j, coo_j, type_j, energy_i, stat)
call accelnet_final(stat)
```

The atomic coordinates are Cartesian, species and atom indices are one-based,
and `accelnet_atomic_energy_and_forces` adds its contribution to the caller's
force array. Thus the force array must be zeroed once before looping over
central atoms. The neighbor-list API returns Cartesian positions of periodic
images and the original one-based atom index for each image. Select Chebyshev
compatibility version 0, 1, or 10 with
`accelnet_set_chebyshev_version` before loading the last potential; version 0
is the default.

Choose the Chebyshev angular evaluation algorithm at runtime with
`accelnet_set_chebyshev_evaluation`. The accepted constants are
`ACCELNET_CHEBYSHEV_AUTO`, `ACCELNET_CHEBYSHEV_DIRECT`, and
`ACCELNET_CHEBYSHEV_MOMENT`. The selection may be made before or after model
loading and can be queried with `accelnet_get_chebyshev_evaluation()`. `AUTO`
is the default. The common kernel selects moments when
`N_angular * (angular_order + 1) >= M`, where `M` is the number of Cartesian
moments. This estimates operation count rather than timing the hardware.
`DIRECT` and `MOMENT` force the corresponding path regardless of neighbor count.

G5 selection is independent: `accelnet_set_g5_evaluation(mode, stat)` accepts
`ACCELNET_G5_AUTO`, `ACCELNET_G5_DIRECT`, `ACCELNET_G5_MOMENT`, and
`ACCELNET_G5_MOMENT_FORCE`. Auto uses integer powers 1--10 and at least 16
angular neighbors. The explicit moment modes support powers 1--16;
`MOMENT_FORCE` also bypasses the neighbor-count threshold. Ineligible powers
remain direct. See [the mode table](../docs/implementation-status.md#descriptor-coverage).

The C ABI provides the same Chebyshev interface:

```c
accelnet_set_chebyshev_evaluation(ACCELNET_CHEBYSHEV_DIRECT, &stat);
int mode = accelnet_get_chebyshev_evaluation();
```

The exported limits `accelnet_nsf_max`, `accelnet_nnb_max`,
`accelnet_Rc_min`, and `accelnet_Rc_max` have the same purpose as their aenet
counterparts. API failures are reported through `ACCELNET_OK` and the
`ACCELNET_ERR_*` status values. The neighbor-list procedures intentionally
retain aenet's status-free signatures.

### Virial

`accelnet_atomic_energy_and_forces_virial` returns an atomic energy and adds
the corresponding forces and full 3-by-3 virial to caller-owned accumulators:

```fortran
! Zero these once before the central-atom loop, not inside it.
forces = 0.0_real64
virial = 0.0_real64
! For each central atom i, supply its full neighbor list:
call accelnet_atomic_energy_and_forces_virial(coo_i, type_i, index_i, n_j, &
    coo_j, type_j, index_j, natoms, energy_i, forces, virial, stat)
```

For each central atomic energy, the contribution is
`W(a,b) = sum_j (coo_j(a,j) - coo_i(a)) * F_j(b)`, where `F_j` is that
atomic energy's force contribution to neighbor image `j`. Equivalently,
`W(a,b) = -dE/d epsilon(b,a)` for a homogeneous deformation
`r' = (I + epsilon) r` applied to both positions and lattice vectors.
This matches the displacement-times-force convention of the PIMD ænet
interface; all nine components are returned, rather than only one triangle.

The tensor has energy units: eV for a model using eV and Angstrom.
It is not divided by the cell volume, has no kinetic contribution, and needs
no factor of one half when summed over central atoms. Convert to a pressure
or stress convention explicitly in the caller; tensile-positive
configurational stress is `-W / volume` for these rotationally invariant
potentials. Input/output units otherwise follow the loaded model, as for
the energy/force API.

Supply Cartesian coordinates of the actual periodic images in `coo_j`.
Different images can share an original atom index in `index_j`, including
the central atom's index. The implementation forms each image's virial
before its force contributions are combined under those atom indices. This
also handles cells small enough to contain multiple images within the cutoff.
Reconstructing the periodic virial from folded atom forces and wrapped
coordinates alone is not equivalent. Nonperiodic structures are supported too.

The C symbol has the same name. Arrays use Fortran column-major layout:
`virial[a + 3*b]` holds `W(a+1,b+1)` for zero-based C indices.
Both accumulators must be initialized by the caller, just as with the
existing additive force API. Detected argument/initialization errors leave
them unchanged and report a nonzero `stat`.
The original `accelnet_atomic_energy_and_forces` signature and C ABI remain
unchanged, and do not calculate the additional tensor.

The virial tests sweep 13 dimensionless strain steps from `1e-2` to `1e-8`
using energy-only central differences of all nine tensor components. A pass
requires two adjacent steps to satisfy `abs(error) <= 2e-6 + 2e-7*abs(W)`
for every component; the tolerance is in model energy units and does not
scale with the large constant atomic reference energy. Synthetic Chebyshev
(AUTO/DIRECT/MOMENT), LJ, and bundled n2p2 models run without external data.
CSV files named `virial-convergence*.csv` are written under the predictor
build directory. If `ACCELNET_PREDICTOR_GOLDEN_DIR` points to the Ti/O corpus,
the tests also check `structure0001.xsf` and `structure2935.xsf`, their sheared
cells, and all three Chebyshev evaluation modes. Very small steps can amplify
floating-point cancellation and are recorded rather than required to improve
monotonically. See the [validation report](../docs/validation/virial-2026-09-25/README.md).

Additional self-contained tests cover nonlinear, per-element n2p2 networks
and a synthetic H/O model with G2/G4/G5 descriptors and nontrivial scaling.
The latter combines direct contraction and saved descriptor Jacobians
in the common kernel. Rotation covariance (`W' = R W R^T`), atom/species
permutations, independent lattice translations of individual atoms, and image
metadata are checked in molecular, orthogonal, and triclinic geometries.
G5 DIRECT/MOMENT/MOMENT_FORCE modes also undergo energy-only strain checks.
The C API tests invalid indices/types/counts, unloaded/finalized states,
accumulator preservation on errors, reloading, and the smooth cutoff boundary.
See [additional regression results](../docs/validation/virial-2026-09-25/additional-tests.md).

An existing model can be replaced explicitly without changing the variable or
reallocating the caller-owned structure/force objects. Old network and setup
allocations are released as part of the reload:

```fortran
! Reconstruct descriptors from the replacement NN metadata.
call model%reload(replacement_networks, chebyshev_version=0)

! Or reload explicit setup files together with the networks.
call model%reload(replacement_setups, replacement_networks)
```

The equivalent module procedures are
`reload_predictor_from_networks(model, networks, version)` and
`reload_predictor(model, setups, networks)`. Reloading is intentionally not
part of the timed inference path. If the replacement changes the species order,
structures must be read or remapped again using the new `model%species_names`.

## Validation and benchmarks

`accelnet-predict-benchmark` reports the in-memory energy API, in-memory
energy+force API, and the XSF energy API separately. For a direct comparison
against the untouched original implementation, first build AccelNet outside
its source tree, then configure with `ACCELNET_ORIGINAL_BUILD_DIR`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DACCELNET_ORIGINAL_BUILD_DIR=/path/to/original-accelnet-build
cmake --build build -j

build/accelnet-original-predict-benchmark \
  Ti.nn O.nn structure.xsf 1000
```

The original comparison calls AccelNet's
`AccelNet_atomic_energy_and_forces_novirial` once per atom and includes its
neighbor-list construction, matching the structure-level inference path while
excluding model loading and file I/O.

An existing ænet build can be included in the same direct-API comparison:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DAENET_BUILD_DIR=/path/to/aenet/build
cmake --build build -j

build/aenet-predict-benchmark \
  Ti.nn.ascii O.nn.ascii structure.xsf 1000
```

This harness calls ænet's public `aenet_atomic_energy` and
`aenet_atomic_energy_and_forces` routines. Both models and the XSF structure
are loaded once before timing; each measured call includes neighbor-list
construction but excludes model and file loading.

To measure only the atomic inference APIs when the caller has already supplied
the neighbor lists, use the direct comparison benchmark:

```sh
build/accelnet-atomic-api-benchmark \
  Ti.nn.ascii O.nn.ascii structure.xsf 1000
```

Append three replication counts to benchmark a periodic supercell without
creating another XSF file. For example, a 24-atom input becomes 1,152 atoms:

```sh
build/accelnet-atomic-api-benchmark \
  Ti.nn.ascii O.nn.ascii structure.xsf 20 4 4 3
```

This executable constructs every local environment once and passes exactly the
same Cartesian neighbor-image coordinates, species IDs, and original atom
indices to AccelNet and ænet. Model loading, XSF parsing, and neighbor-list
construction are outside the timed regions. It reports energy-only and
energy-plus-force time per complete structure, their speedups, and the maximum
energy/force disagreement before timing. The target is available when CMake is
configured with `AENET_BUILD_DIR`.

The standalone descriptor-order benchmark measures the current moment kernel
at angular orders 0 through 16 while holding the radial order and local atomic
environment fixed:

```sh
build/accelnet-order-scaling-benchmark structure.xsf 2000
```

For a corpus of sequentially named XSF files (`structure0001.xsf`, ...), the
corpus harnesses load every structure before timing:

```sh
build/accelnet-corpus-benchmark 2 Ti.fingerprint.stp O.fingerprint.stp \
  Ti.nn.ascii O.nn.ascii /path/to/TiO2-xsf 1000 1
build/aenet-corpus-benchmark Ti.nn.ascii O.nn.ascii \
  /path/to/TiO2-xsf 1000 1
```

The final two arguments are the number of structures and full-corpus repeats.
Both programs exclude model and XSF loading while retaining neighbor-list,
descriptor, network, energy, and force work.

To verify every energy and Cartesian force component against aenet rather than
only comparing aggregate checksums, use:

```sh
build/chebyshev-corpus-validate 2 Ti.fingerprint.stp O.fingerprint.stp \
  Ti.nn.ascii O.nn.ascii /path/to/TiO2-xsf 7815
```

The validator reports maximum absolute errors and their structure/atom/component
locations, and exits unsuccessfully if an energy or force error exceeds `1e-8`.

The test suite compares energies and analytic forces with the archived
AccelNet Fortran outputs, verifies `predict.in`, and converts the same Ti/O
ASCII networks to AccelNet's binary record layout before checking that binary
and ASCII inference are identical. The original AccelNet source and fixtures
are read-only inputs and are not modified.

Chebyshev, LJ, and Behler2011 fingerprints embedded in NN files are
reconstructed by the predictor. Standard ænet NN metadata does not store the
Chebyshev implementation version, so `predict.in` carries this additional
tag while remaining backward compatible with existing files that omit it.
Extended networks store `cutoff_type` and `cutoff_alpha` in descriptor
parameter rows 5 and 6 for all three descriptor families. Older networks with
fewer rows retain their historical defaults: cosine for Chebyshev and
Behler2011, and hard truncation for LJ.
For Behler2011 networks, embedded functions are placed in ænet's canonical
species/function ordering rather than their textual metadata order, matching
the input order used by the trained neural network and scaling arrays.

### CPU energy regression and GPU checks

The [revision 1.14 report](../docs/validation/atomic-energy-preparation-2026-09-28/README.md)
records the latest atomic-energy preparation optimization and regression checks.
The [revision 1.13 report](../docs/validation/energy-common-2026-09-27/README.md)
records the preceding shared-kernel energy/force optimization.
CPU comparisons compile OpenMP out on both sides; setting `OMP_NUM_THREADS=1`
alone does not remove runtime overhead. To enable the Ti/O energy gate:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DACCELNET_PREDICTOR_GOLDEN_DIR=/path/to/fortran_predict \
  -DACCELNET_PUBLIC_API_BASELINE_EXECUTABLE=/path/to/archived/accelnet-public-api-benchmark
cmake --build build --parallel
ctest --test-dir build -R predictor_aenet_energy_performance --output-on-failure
```

The archived and current executables must use the same driver, compiler and
optimization flags. The test verifies energies at 64/192/512 atoms for both
structure-energy and per-atom-energy APIs, checks the absence of OpenMP runtime
symbols, and applies `ACCELNET_CPU_MAX_SLOWDOWN` (default 1.10). Use the Python
comparison script directly to select another CPU core or measurement duration.

GPU builds and `use_host=.true.` threading checks are described in the
[target guide](../docs/openmp-target.md). Build all test dependencies before
running a full CTest label. Performance tests should run separately, without
concurrent builds or correctness tests. Energy-only target calls preserve the
caller's force/virial accumulators; the suite tests reuse after force calls.
