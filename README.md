# AccelNet

**New in 1.1.2:** Fix shared C API visibility with GNU Fortran 16.2 and add macOS export, linking, and relocation regression tests — [details](changes.md).

AccelNet is a Fortran library and command-line toolkit for evaluating
machine-learning interatomic potentials from **ænet and n2p2**. It computes
energies, analytic forces and configurational virials, provides Fortran/C APIs,
and integrates with LAMMPS on CPU and GPU. It is an inference package; training
remains in the upstream tools.

The method is described in the [AccelNet paper](https://arxiv.org/abs/2608.03280).

## Start here

The default build and the usual command-line/API workflow remain **serial CPU**.
GPU support and CPU OpenMP threading are optional, explicit backends.

| Task | Guide |
|---|---|
| Build and run with the usual CPU workflow | [CPU quick start](#cpu-quick-start), [install and link](#install-and-link) |
| Use multiple CPU threads | [CPU OpenMP threading](#cpu-openmp-threading) |
| Run LAMMPS on CPU | [LAMMPS on CPU](#lammps-on-cpu), [complete CPU instructions](interfaces/lammps/README.md#cpu-build-and-run) |
| Use a GPU, including LAMMPS | [GPU build and execution](#gpu-build-and-execution) |
| Check models and common-code coverage | [Supported models and methods](#supported-models-and-methods), [implementation status](docs/implementation-status.md) |
| Embed AccelNet | [Predictor APIs](AccelNetPredictor/README.md), [CSR batch API](docs/batch-api.md) |
| Convert model formats | [Fortran and Julia converters](AccelNetModelConverter/README.md) |
| Inspect equations and performance | [Optimization methods](speedupmethods.md), [validation](#validation-and-performance) |

## CPU quick start

Requirements: CMake 3.20+, a Fortran 2008 compiler, and a C compiler. GNU Fortran
11/13 are covered by CI; the recent CPU measurements use GNU Fortran 11.4.
The core libraries require no BLAS, LAPACK, MPI, ænet or n2p2 installation.
Python is used by some tests; optional upstream n2p2 reference tests also need
a C++ compiler. External reference tools and models are optional.

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
  -DACCELNET_BUILD_OPENMP_TARGET=OFF -DN2P2_SCALING_EXECUTABLE=
cmake --build build --parallel
ctest --test-dir build -LE "gpu|performance" --output-on-failure
```

Executables are in `build/bin`; libraries are in `build/lib`. The default build
is static and CPU-only. GPU support is opt-in. Use a separate build directory
when changing compiler families. The command disables auto-discovery of the
optional external `nnp-scaling` executable; set `N2P2_SCALING_EXECUTABLE` to a
working local executable to enable that upstream comparison.

### Run a prediction

With ænet/AccelNet setup files, networks and an XSF structure:

```sh
build/bin/accelnet-predict 2 \
  Ti.fingerprint.stp O.fingerprint.stp \
  Ti.nn.ascii O.nn.ascii structure.xsf --forces

# Or use an existing prediction input file:
build/bin/accelnet-predict predict.in
```

A bundled n2p2 fixture provides a runnable smoke example from the repository:

```sh
build/bin/accelnet-predict --n2p2-data \
  AccelNetPredictor/test/data/n2p2 \
  AccelNetPredictor/test/data/n2p2/input.data
```

For your own model, replace the directory and structure paths. An n2p2 directory
contains `input.nn`, `weights.%03d.data`, and `scaling.data` when required by
its scaling mode. Multiple n2p2 `input.data` structures can also be evaluated:

```sh
build/bin/accelnet-predict --n2p2-data /path/to/model input.data
```

Input coordinates and returned energies/forces use the model's physical units.
`predict.in`, `--n2p2` and `--n2p2-data` already output forces; do not append
`--forces` to those forms. These commands use the ordinary serial CPU path, including
in a GPU-enabled build. Setting `OMP_NUM_THREADS` does not parallelize this CLI.

## Install and link

```sh
cmake --install build --prefix /path/to/install
```

The install contains libraries, Fortran module files, C headers, CLI programs
and CMake package files. Consumer projects can use:

```cmake
find_package(AccelNetPredictor CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE AccelNet::AccelNet)
# The optional threaded/offload API instead links AccelNet::Target;
# see the corresponding backend section below.
```

Use the same Fortran compiler family for library and application module files.
Static C consumers also need the matching Fortran runtime. Headers and examples
are in [AccelNetPredictor](AccelNetPredictor/README.md).

| CMake option | Default | Purpose |
|---|---|---|
| `BUILD_SHARED_LIBS` | `OFF` | Build shared libraries |
| `BUILD_TESTING` | `ON` | Build correctness tests and benchmark drivers |
| `ACCELNET_BUILD_OPENMP_TARGET` | `OFF` | Build the optional target library |
| `ACCELNET_OPENMP_TARGET_FLAGS` | Empty | Compiler/link flags for that library |
| `ACCELNET_TARGET_SERIAL` | `OFF` | Compile the target API without OpenMP for serial comparisons |
| `ACCELNET_BUILD_REFERENCE_TESTS` | `OFF` | Compare with an external historical AccelNet tree |
| `ACCELNET_PREDICTOR_GOLDEN_DIR` | Optional sibling data directory | Enable real Ti/O-model tests if the files exist |
| `ACCELNET_DESCRIPTORS_BUILD_N2P2_REFERENCE` | `ON` | Build upstream n2p2 C++ reference tests if external sources exist |
| `N2P2_SCALING_EXECUTABLE` | Optional sibling executable | External n2p2 descriptor comparison; set empty to disable |
| `ACCELNET_PUBLIC_API_BASELINE_EXECUTABLE` | Empty | Enable the archived-baseline CPU energy performance test |

## CPU OpenMP threading

This is an optional **CPU multithreaded** build of the common target API; no GPU
or CUDA installation is needed. GNU Fortran 11.4 has been tested with the flags
below. Keep this build separate from the ordinary serial build:

```sh
cmake -S . -B build-openmp -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
  -DACCELNET_BUILD_OPENMP_TARGET=ON -DACCELNET_TARGET_SERIAL=OFF \
  -DACCELNET_OPENMP_TARGET_FLAGS="-fopenmp -foffload=disable -ffree-line-length-none" \
  -DN2P2_SCALING_EXECUTABLE=
cmake --build build-openmp --parallel

OMP_NUM_THREADS=1 build-openmp/bin/test_batch_target --host --quick
OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close \
  build-openmp/bin/test_batch_target --host --quick
OMP_NUM_THREADS=1 ctest --test-dir build-openmp -L target-host --output-on-failure
```

Use physical cores available to your process and avoid oversubscribing them.
The `--host` option above belongs to the test driver, not `accelnet-predict`.
Do not run the GPU-labelled tests for this host-only build.

In a Fortran application, link `AccelNet::Target`, use the `accelnet_batch_target`
module and explicitly select the host when packing an already loaded model:

```fortran
call packed%initialize(model, use_host=.true.)
call evaluate_batch_target(packed, species, centers, offsets, indices, &
    displacements, energies, forces, work, virial)
```

The full declarations and CSR/lifetime contract are in the
[target API guide](docs/openmp-target.md#api-and-lifetime).
`OMP_NUM_THREADS=8` controls this target-host execution. It does **not** change
ordinary CLI, object, atomic, SFB or `evaluate_batch` calls: these still use the
serial module compiled from the same source. There is no threaded CLI switch.
`ACCELNET_TARGET_SERIAL=ON` is a benchmarking option that disables OpenMP even
for the target API; leave it **OFF** when requesting host threads.

LAMMPS `pair_style accelnet` is also serial inside each MPI rank; there is no
`accelnet/omp` style. Use MPI ranks for its CPU parallelism. LAMMPS's `BUILD_OMP`
or `-pk omp` options do not turn this adapter into a threaded pair style.

## LAMMPS on CPU

The usual CPU interface requires the serial AccelNet libraries, a C++ compiler
and, for MPI runs, MPI. It needs no CUDA or GPU compiler.

Use **LAMMPS 22 Jul 2025 Update 6**, the latest stable release verified on
2026-09-28. The full instructions pin `stable_22Jul2025_update6` and use
`python3 interfaces/lammps/install.py /path/to/lammps` for automatic installation.

| Release | CPU build | Pair style |
|---|---|---|
| **22 Jul 2025 Update 6 (recommended)** | CMake, `PKG_ACCELNET=ON`, `PKG_GPU=OFF` | `accelnet` |
| 29 Aug 2024 Update 4 (compatibility) | Same installer and CMake options | `accelnet` |
| 4 Feb 2020 (legacy) | Traditional make, GNU Fortran libraries | `accelnet` |

The [compatibility report](docs/validation/lammps-current-2026-09-28/README.md)
records CPU/GPU checks for the recommended release. The newer 2 Sep 2026
release candidate requires a GPU API update and is not an installer target.

Follow the [CPU installation and run instructions](interfaces/lammps/README.md#cpu-build-and-run).
For network files, after creating the simulation box:

```lammps
pair_style accelnet auto Ti.nn.ascii O.nn.ascii
pair_coeff * *
```

A CPU run can also load an n2p2 directory directly. List elements in LAMMPS
atom-type order:

```lammps
pair_style accelnet n2p2 /path/to/model Ti O
pair_coeff * *
```

Run with `lmp -in in.cpu` for one process or `mpirun -np 4 lmp -in in.cpu`
for four MPI ranks. AccelNet remains serial within each rank. GPU installation
and input commands are collected in the [GPU section below](#lammps-on-gpu).

## GPU build and execution

### Required compiler and runtime

For the validated NVIDIA path, install **NVIDIA HPC SDK** and use its
**`nvfortran`** compiler. The tested combination is **NVHPC 25.3**, the SDK's
**CUDA 12.8** toolkit, a compatible NVIDIA driver, and H100 NVL or RTX PRO 6000
Blackwell GPUs. CMake 3.20+ and a C compiler are also required. `nvidia-smi`
should detect the intended GPU before running an offload test.

`nvcc` alone cannot compile the Fortran OpenMP target kernels. An ordinary
`gfortran -fopenmp` build provides host threading, not this validated NVIDIA
GPU configuration. Other Fortran offload toolchains and AMD/Intel GPU execution
have not been validated here; portable OpenMP source does not establish that
every compiler/device combination works.

For LAMMPS GPU integration there is an additional requirement: **NVHPC `nvc++`**
and the **CUDA toolkit (`nvcc`, headers and libraries)**. The adapter's CMake
checks enforce NVHPC C++ and Fortran, CUDA, and double precision. An MPI
installation is required when building LAMMPS with `BUILD_MPI=ON`.

### Build AccelNet and verify offload

Set `NVHPC_ROOT` to the actual SDK version directory on your machine. From the
AccelNet repository root, using H100 as the example:

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
# Replace this placeholder with a UUID printed above:
export CUDA_VISIBLE_DEVICES=GPU-REPLACE-WITH-YOUR-DEVICE-UUID
export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_THREADS=1
build-gpu/bin/test_batch_target --quick
ctest --test-dir build-gpu -L gpu --output-on-failure
```

Use `-gpu=cc120` for the tested Blackwell GPU, or `-gpu=cc90,cc120` to embed both
architectures, with a toolchain that supports them. Select architecture flags
for your actual hardware; the example does not cover all NVIDIA GPU models.
The GPU recipe disables optional external n2p2 C++ reference tests to avoid
mixing the host C++ linker with NVHPC Fortran objects; bundled GPU correctness
tests remain enabled. Run those optional upstream tests in the GNU CPU build.
GPU and CPU compiler families must use separate build directories and matching
Fortran module files. The target library's compile/link flags propagate to its
CMake consumers.

### Choose a GPU entry point

Building GPU support does not redirect ordinary prediction commands to a GPU.
Use the Fortran target API (`AccelNet::Target`) or the C API in
`accelnet_target.h`, or use the LAMMPS GPU adapter below. The ordinary
`accelnet-predict` CLI has no GPU option.

```fortran
! model is already loaded; default initialization requires actual GPU execution.
call packed%initialize(model)  ! optional device= selects an OpenMP device ID
call evaluate_batch_target(packed, species, centers, offsets, indices, &
    displacements, energies, forces, work, virial)
```

The [target API guide](docs/openmp-target.md) gives declarations, C handles,
model snapshots, resident workspaces, CSR data and energy-only calls. Default
initialization rejects unintended CPU fallback. `use_host=.true.` explicitly
selects the CPU OpenMP backend instead.

### LAMMPS on GPU

Use **22 Jul 2025 Update 6** and `pair_style accelnet/gpu`; the same adapter
also supports **29 Aug 2024 Update 4**. It combines the LAMMPS **GPU package**
(`PKG_GPU=ON`, `GPU_API=cuda`,
`GPU_PREC=double`) with AccelNet's Fortran OpenMP target library. It is not a
Kokkos or LAMMPS OPENMP pair style. The complete compiler, installation, build
and run commands are in the [LAMMPS GPU section](interfaces/lammps/README.md#gpu-build-and-run).

```lammps
# Before read_data/create_box:
package gpu 1 neigh yes newton on split 1
# After creating the box, with types matching the embedded global species order:
pair_style accelnet/gpu auto Ti.nn.ascii O.nn.ascii
pair_coeff * *
```

Select GPU visibility as above, and run `lmp -in in.gpu` or an MPI launch of that
GPU-enabled executable. This explicit `pair_style`/`package` example does not
require `-sf gpu` or another `-pk gpu` command-line override.
`neigh no` still evaluates the potential on the GPU; it only selects CPU
neighbor-list construction. `neigh yes/hybrid`, global virial and per-atom
energy are supported. `newton on` and `split 1` are required.

GPU input uses embedded network files. Direct n2p2-directory loading belongs
to the CPU pair style; convert n2p2 models with the
[Fortran converter](AccelNetModelConverter/README.md) before using them on GPU.
Per-atom stress and other adapter restrictions, device assignment and data
transfers are described in [the GPU integration guide](docs/lammps-gpu.md).

## Supported models and methods

All supported production potential-inference entry points use the **same
maintained Fortran numerical kernels**. The ordinary CPU library compiles them
with OpenMP directives removed. The optional target library compiles them for
OpenMP host threads or GPU offload. Descriptor evaluation, the neural network,
force contraction and virial evaluation run on the selected backend in FP64.

This means common **production inference mathematics**, not identical API
wrappers, device transfers or compiled binaries. Independent old evaluators
remain in [legacy/cpu-reference](legacy/cpu-reference/README.md) for explicit
reference tests and standalone low-level descriptor/NN compatibility APIs.
Those lower-level utilities are not all migrated to the common inference
pipeline; production potential evaluation has no automatic legacy fallback.

| Descriptor family | Common CPU/GPU evaluation | Exact moment support |
|---|---|---|
| ænet Chebyshev, versions 0/1/10 | Direct and moment | Yes |
| ænet Behler G1/G2/G3; n2p2 type 2 | Single-neighbor sums | No separate moment method needed |
| ænet Behler G4; n2p2 type 3 | Direct angular pairs | No general finite factorization |
| ænet Behler G5; n2p2 type 9 | Direct and moment | Integer angular powers 1--16; eligibility depends on mode |
| n2p2 types 12/20/23 | Weighted/compact radial sums | No separate moment method needed |
| n2p2 types 13/21/24 | Exact direct angular pairs | No general finite single-neighbor factorization |
| n2p2 types 22/25 | Exact direct angular pairs | No approximate moment expansion is used |
| AccelNet LJ extension | Radial sums | No separate moment method needed |

**Types 13/21/24 are implemented**, including forces and virial. Lack of a
finite moment representation does not mean lack of model support. Third-distance
and angular window terms prevent the general finite polynomial factorization
used for G5. See [the exactness boundary](docs/implementation-status.md#descriptor-coverage).

Multi-element models, multiple Chebyshev basis blocks per element, and mixed
Chebyshev/LJ/Behler blocks are supported. Each basis block may have its own
orders and cutoffs; there is no one-Chebyshev-block-per-element restriction.

Method selection is explicit and reproducible:

- **Chebyshev:** `auto`, `direct`, or `moment`. Auto uses angular neighbor count,
  angular order and moment count to estimate work.
- **G5:** auto uses exact integer powers 1--10 with at least 16 angular neighbors.
  Explicit moment modes also allow powers 11--16. Fractional, near-integer and
  powers above 16 use direct evaluation. LAMMPS `g5 moment` bypasses the
  neighbor-count threshold for eligible powers.
- Auto is an operation-count/eligibility rule, not a device-specific timing
  tuner. Moment is not always faster than direct.

### Model compatibility

The compatibility targets are ænet 2.0.4 and n2p2 2.3.0:

- ænet/AccelNet ASCII networks and compatible native binary networks;
- n2p2 short-range **2G-HDNNP** directories with SF types **2/3/9/12/13/20--25**;
- per-element network topology, activation functions, descriptor scaling,
  energy normalization/reference energies and n2p2 `normalize_nodes`;
- n2p2 cutoff types 0--8 and the AccelNet fractional-cutoff extension type 9;
- XSF structures, AccelNet `predict.in`, and n2p2 `input.data` structures.

Native binary files depend on the Fortran record representation; ASCII is the
more portable interchange format. Standard ænet network metadata does not
encode the Chebyshev version: AccelNet defaults to version 0; select 1 or 10
explicitly when required. AccelNet extended ASCII metadata is not guaranteed
to be readable by upstream ænet.

Training, n2p2 4G/Q charge/electrostatic models, and general ASE input formats
are not implemented. Unsupported settings are rejected. Accepted syntax,
unit conventions, conversion restrictions and endpoint behavior are specified
in [model compatibility](docs/model-compatibility.md).

## Validation and performance

[CI](.github/workflows/tests.yml) checks GNU 11/13 Release builds and GNU 13
Debug/shared builds with runtime checks. It uses bundled/synthetic fixtures;
it does not run on a GPU. Optional reference comparisons and real-model tests
are described in the component READMEs.

Release-specific changes, validation results and performance measurements are
collected in [changes.md](changes.md).

Run performance gates separately from correctness tests and competing work:

```sh
ctest --test-dir build -L performance --output-on-failure
```

To enable the energy regression gate, configure the same serial build with
`ACCELNET_PUBLIC_API_BASELINE_EXECUTABLE` pointing to an archived
`accelnet-public-api-benchmark`, and `ACCELNET_PREDICTOR_GOLDEN_DIR` pointing to
the Ti/O models. The test rejects linked OpenMP runtimes and uses the default
1.10 time-ratio limit. The [reproduction guide](docs/validation/energy-common-2026-09-27/README.md#reproduction)
explains compiler matching, fixtures and paired measurements. Historical reports
retain their original scope; the [implementation status](docs/implementation-status.md)
is the current coverage reference.

## Citation

If you use AccelNet, please cite:

- **AccelNet:** Y. Nagai, “AccelNet: Exact backward-compatible acceleration of
  polynomial angular descriptors through Cartesian moment factorization,”
  arXiv:2608.03280 [cond-mat.mtrl-sci] (2026),
  [doi:10.48550/arXiv.2608.03280](https://doi.org/10.48550/arXiv.2608.03280).

Please also cite the publications relevant to the upstream model and
descriptor used:

- **ænet:** N. Artrith and A. Urban, “An implementation of artificial
  neural-network potentials for atomistic materials simulations: Performance
  for TiO2,” *Computational Materials Science* **114**, 135--150 (2016),
  [doi:10.1016/j.commatsci.2015.11.047](https://doi.org/10.1016/j.commatsci.2015.11.047).
- **Behler--Parrinello HDNNP method:** J. Behler and M. Parrinello,
  “Generalized neural-network representation of high-dimensional
  potential-energy surfaces,” *Physical Review Letters* **98**, 146401 (2007),
  [doi:10.1103/PhysRevLett.98.146401](https://doi.org/10.1103/PhysRevLett.98.146401).
- **ænet Chebyshev descriptors, when used:** N. Artrith, A. Urban, and
  G. Ceder, “Efficient and accurate machine-learning interpolation of atomic
  energies in compositions with many species,” *Physical Review B* **96**,
  014112 (2017),
  [doi:10.1103/PhysRevB.96.014112](https://doi.org/10.1103/PhysRevB.96.014112).
- **n2p2:** A. Singraber, J. Behler, and C. Dellago, “Library-Based LAMMPS
  Implementation of High-Dimensional Neural Network Potentials,” *Journal of
  Chemical Theory and Computation* **15**, 1827--1840 (2019),
  [doi:10.1021/acs.jctc.8b00770](https://doi.org/10.1021/acs.jctc.8b00770),
  together with the [n2p2 software archive](https://doi.org/10.5281/zenodo.1344446).
- **Fractional cutoff (cutoff type 9), when used:** H. Mori, T. Tsuru,
  M. Okumura, D. Matsunaka, Y. Shiihara, and M. Itakura, “Dynamic interaction
  between dislocations and obstacles in bcc iron based on atomic potentials
  derived using neural networks,” *Physical Review Materials* **7**, 063605
  (2023), Appendix A, Eqs. (A4)--(A6),
  [doi:10.1103/PhysRevMaterials.7.063605](https://doi.org/10.1103/PhysRevMaterials.7.063605).
