# LAMMPS interfaces

Use the **CPU workflow below by default**. GPU installation and execution are
collected in the separate [GPU section](#gpu-build-and-run).
Commands below start from the AccelNet repository root unless noted otherwise.

**Recommended release: LAMMPS 22 Jul 2025 Update 6**, the latest stable release
verified on 2026-09-28. Pin the tag `stable_22Jul2025_update6` to reproduce the
validated build; do not substitute the moving `develop` or `stable` branch.

| LAMMPS release | CPU | GPU |
|---|---|---|
| **22 Jul 2025 Update 6 (recommended)** | `pair_style accelnet`, CMake | `pair_style accelnet/gpu`, CUDA GPU package + Fortran OpenMP target |
| 29 Aug 2024 Update 4 (compatibility) | Same CMake adapter and installer | Same GPU adapter |
| 4 Feb 2020 (legacy) | `pair_style accelnet`, traditional make | No adapter provided |

The [compatibility validation](../../docs/validation/lammps-current-2026-09-28/README.md)
records CPU and H100 GPU checks for the recommended release. The 2 Sep 2026
stable release candidate is not a supported installer target: its GPU API
changed. Use the pinned stable release below.

Both adapters call the common Fortran **numerical inference kernels**. CPU/GPU
input packing, device transfers and MPI integration differ. Independent legacy
low-level/reference evaluators remain in
[legacy/cpu-reference](../../legacy/cpu-reference/README.md), but production
LAMMPS potential evaluation has no legacy fallback.

## CPU build and run

### Requirements and AccelNet libraries

Use CMake 3.20+, GNU Fortran (11/13 in CI), a C compiler, a C++ compiler, and MPI
for parallel LAMMPS runs. GNU Fortran 11.4 is used in the recorded CPU timings.
No GPU, CUDA toolkit or NVIDIA compiler is needed. Use the same Fortran compiler
family for AccelNet and the LAMMPS link step. The legacy 2020 install script
specifically configures GNU Fortran runtime libraries.

```sh
# From the AccelNet repository root:
export ACCELNET_SOURCE="$PWD"
cmake -S "$ACCELNET_SOURCE" -B "$ACCELNET_SOURCE/build" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_C_COMPILER=gcc -DBUILD_SHARED_LIBS=OFF \
  -DACCELNET_BUILD_OPENMP_TARGET=OFF -DN2P2_SCALING_EXECUTABLE=
cmake --build "$ACCELNET_SOURCE/build" --parallel
```

The adapter needs `build/lib/libaccelnet.a` and
`build/lib/libAccelNetDescriptors.a`.

### Install into the recommended LAMMPS release

Download the pinned release and run the shared installer. If it has already
been downloaded, set `LAMMPS_SOURCE` to that tree and skip `git clone`.
Python 3 is required for installation. The same command installs the files
needed by either CPU or GPU builds; the CMake options select the backend.

```sh
export LAMMPS_SOURCE="$ACCELNET_SOURCE/../lammps-22Jul2025-update6"
git clone --depth 1 --branch stable_22Jul2025_update6 \
  https://github.com/lammps/lammps.git "$LAMMPS_SOURCE"
python3 "$ACCELNET_SOURCE/interfaces/lammps/install.py" "$LAMMPS_SOURCE"

cmake -S "$LAMMPS_SOURCE/cmake" -B "$LAMMPS_SOURCE/build-accelnet-cpu" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DBUILD_MPI=ON -DBUILD_OMP=OFF -DPKG_ACCELNET=ON -DPKG_GPU=OFF \
  -DACCELNET_DIR="$ACCELNET_SOURCE/build"
cmake --build "$LAMMPS_SOURCE/build-accelnet-cpu" --parallel
```

The installer copies the current canonical C headers, registers ACCELNET with
CMake and applies the triclinic GPU sorting fix. Re-running it updates the
adapter without duplicating registration or patches. There is no manual
`cp`/`patch` step. It checks the exact release/update and rejects unsupported
versions before modifying them.

For an existing **29 Aug 2024 Update 4** tree, run the same installer and build
commands with its path in `LAMMPS_SOURCE`. Both versions use one adapter source;
its historical `29Aug2024/` directory name does not limit supported versions.
The old `29Aug2024/install.py` command remains a compatibility entry point.

`ACCELNET_DIR` points to the build directory containing `lib/`; the CMake module
also accepts a source tree with the default `build/lib/`. It links the matching
Fortran runtime. Set `BUILD_MPI=OFF` for a serial LAMMPS build without MPI.
If needed, give CMake `-DMPI_CXX_COMPILER=/path/to/mpicxx` for MPI discovery.
Keep `PKG_GPU=OFF` for this CPU recipe; enabling both GPU and ACCELNET invokes
the GPU adapter's additional compiler/library requirements.

### CPU input and launch

Example `in.cpu` (run from a directory containing the data/model files):

```lammps
units metal
atom_style atomic
newton on
read_data TiO2.data
pair_style accelnet auto Ti.nn.ascii O.nn.ascii
pair_coeff * *
neighbor 0.6 bin
neigh_modify every 10 delay 0 check yes
fix integrate all nve
run 100
```

`units metal` is an example for models/data expressed in eV and Angstrom;
choose LAMMPS units consistent with your actual model. For embedded networks,
LAMMPS atom types must follow the model's global species order.

```sh
# One MPI rank/process:
"$LAMMPS_SOURCE/build-accelnet-cpu/lmp" -in in.cpu
# Four MPI ranks:
OMP_NUM_THREADS=1 mpirun -np 4 "$LAMMPS_SOURCE/build-accelnet-cpu/lmp" -in in.cpu
```

There is no `package gpu` line or GPU suffix in this CPU input.

### CPU OpenMP versus MPI

`pair_style accelnet` is **serial within each MPI rank**. It uses the common
serial library with OpenMP directives removed. `OMP_NUM_THREADS=8`, LAMMPS
`BUILD_OMP=ON`, `-pk omp` or `-sf omp` does not make this adapter multithreaded;
no `accelnet/omp` style is provided. Use multiple MPI ranks for LAMMPS CPU
parallelism. AccelNet's explicit target-host API supports CPU OpenMP outside
this adapter; see the [host-threading instructions](../../README.md#cpu-openmp-threading).

### LAMMPS 4 Feb 2020: traditional CPU installation

After the GNU Fortran AccelNet build above, set `LAMMPS_SOURCE` to a separate
4 Feb 2020 source tree. From the AccelNet repository root:

```sh
export LAMMPS_SOURCE=/path/to/lammps-4Feb2020
cp -R interfaces/lammps/4Feb2020/USER-ACCELNET "$LAMMPS_SOURCE/src/"
mkdir -p "$LAMMPS_SOURCE/lib/accelnet/include" "$LAMMPS_SOURCE/lib/accelnet/lib"
cp AccelNetPredictor/include/accelnet.h "$LAMMPS_SOURCE/lib/accelnet/include/"
cp build/lib/libaccelnet.a build/lib/libAccelNetDescriptors.a "$LAMMPS_SOURCE/lib/accelnet/lib/"
cd "$LAMMPS_SOURCE/src"
make yes-user-accelnet
make mpi -j8
mpirun -np 4 ./lmp_mpi -in /path/to/in.cpu
```

The input syntax below applies to both CPU versions except where noted. The
2020 adapter does not supply GPU or OpenMP pair styles.

## CPU input and model mapping

For both versions, the LAMMPS input syntax is:

```lammps
pair_style accelnet H.ann O.ann
pair_coeff * *
```

An n2p2 2G model directory can be loaded directly. Element names after the
directory map LAMMPS atom types to model elements and must be given in LAMMPS
type order; the interface converts them to n2p2's internal element order.

```lammps
# LAMMPS type 1 = Ti, type 2 = O
pair_style accelnet n2p2 /path/to/model Ti O
pair_coeff * *
```

All MPI ranks read the same `input.nn`, `weights.%03d.data`, and optional
`scaling.data` files. The optional leading Chebyshev mode is also accepted
with the n2p2 form. The modern CMake interface additionally accepts the trailing
`g5 MODE` option after the element list.

No additional pair-style option is needed for n2p2 per-element network
topologies or `normalize_nodes`; both are handled while the shared AccelNet
model loader reads `input.nn`.

The 4Feb2020 interface can force the Chebyshev angular algorithm by placing a
mode before the potential files. Omitting it selects `auto`.

```lammps
pair_style accelnet auto H.ann O.ann
pair_style accelnet direct H.ann O.ann
pair_style accelnet moment H.ann O.ann
```

The modern CMake interface also accepts an optional G5 evaluation mode after all
potential files:

```lammps
pair_style accelnet O.nn Ti.nn g5 direct
pair_style accelnet O.nn Ti.nn g5 moment
```

`direct` disables the integer-zeta G5 moment path. `moment` forces that path
without applying the automatic neighbor-count threshold. Omitting the option,
or selecting `g5 auto`, preserves the default automatic behavior.


## GPU build and run

### Required NVIDIA toolchain

Use **LAMMPS 22 Jul 2025 Update 6** with its **GPU package**, not Kokkos or
the LAMMPS OPENMP package. The same adapter retains 29 Aug 2024 Update 4
compatibility. The validated stack is:

| Component | Requirement / tested configuration |
|---|---|
| Fortran compiler | NVIDIA HPC SDK `nvfortran`, tested NVHPC 25.3 |
| C++ compiler for LAMMPS | Same SDK's `nvc++`; enforced by the adapter CMake module |
| CUDA | Toolkit with `nvcc`, headers, libraries and `bin2c`; tested SDK CUDA 12.8 |
| GPU and driver | Compatible NVIDIA driver; tested H100 NVL and RTX PRO 6000 Blackwell |
| Other build tools | CMake 3.20+, C compiler, Python 3 for installation, MPI for `BUILD_MPI=ON` |

`nvcc` alone cannot compile the Fortran kernels. GNU `-fopenmp` alone is not
this GPU offload build. AMD/HIP/OpenCL interoperability is not implemented in
the LAMMPS adapter. Use separate CPU/GPU build directories; do not mix GNU
Fortran libraries or module files with an NVHPC build.

### Build the GPU library and install the adapter

Start in the AccelNet repository root. The commands below target H100. Adjust
SDK and source paths to your installation. Download the pinned release below
only if it is not already available; otherwise skip `git clone`:

```sh
export ACCELNET_SOURCE="$PWD"
export LAMMPS_SOURCE="$ACCELNET_SOURCE/../lammps-22Jul2025-update6"
git clone --depth 1 --branch stable_22Jul2025_update6 \
  https://github.com/lammps/lammps.git "$LAMMPS_SOURCE"
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3
export ACCELNET_CUDA_ROOT="$NVHPC_ROOT/cuda/12.8"
export ACCELNET_GPU_FLAGS="-mp=gpu -gpu=cc90"

cmake -S "$ACCELNET_SOURCE" -B "$ACCELNET_SOURCE/build-gpu" \
  -DCMAKE_Fortran_COMPILER="$NVHPC_ROOT/compilers/bin/nvfortran" \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF \
  -DACCELNET_BUILD_OPENMP_TARGET=ON -DACCELNET_TARGET_SERIAL=OFF \
  -DACCELNET_OPENMP_TARGET_FLAGS="$ACCELNET_GPU_FLAGS"
cmake --build "$ACCELNET_SOURCE/build-gpu" --parallel

python3 "$ACCELNET_SOURCE/interfaces/lammps/install.py" "$LAMMPS_SOURCE"
cmake -S "$LAMMPS_SOURCE/cmake" -B "$LAMMPS_SOURCE/build-accelnet-gpu" \
  -DCMAKE_CXX_COMPILER="$NVHPC_ROOT/compilers/bin/nvc++" \
  -DCMAKE_Fortran_COMPILER="$NVHPC_ROOT/compilers/bin/nvfortran" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG" \
  -DBUILD_MPI=ON -DBUILD_OMP=OFF \
  -DPKG_ACCELNET=ON -DPKG_GPU=ON -DGPU_API=cuda -DGPU_PREC=double \
  -DGPU_ARCH=sm_90 -DCUDA_BUILD_MULTIARCH=OFF \
  -DCUDA_TOOLKIT_ROOT_DIR="$ACCELNET_CUDA_ROOT" \
  -DBIN2C="$ACCELNET_CUDA_ROOT/bin/bin2c" \
  -DACCELNET_DIR="$ACCELNET_SOURCE/build-gpu" \
  -DACCELNET_TARGET_FLAGS="$ACCELNET_GPU_FLAGS"
cmake --build "$LAMMPS_SOURCE/build-accelnet-gpu" --parallel
```

The recommended 2025 release was checked on H100; Blackwell validation was
performed with the retained 2024 release. For that Blackwell device, use
`cc120` and `GPU_ARCH=sm_120`. For a combined H100/Blackwell binary, use `cc90,cc120`, retain `GPU_ARCH=sm_90`, and add
`-DCUDA_NVCC_FLAGS="-gencode=arch=compute_120,code=sm_120"` to the LAMMPS
configuration, as in the archived validation build. Both the Fortran and CUDA
architecture choices must cover the device you run on.

The installer copies the current CPU/GPU adapters and canonical C headers,
registers the package and applies the required triclinic sorting fix. It can
also update a tree with the CPU adapter already installed. Copying only the
CPU package is insufficient for GPU use. If necessary, specify
`-DMPI_CXX_COMPILER=/path/to/mpicxx` for the LAMMPS configure command.

The supplied CMake module applies `-fortranlibs` and the target link flags so
NVHPC orders its own Fortran/offload runtimes. Do not replace this with a guessed
manual runtime-library list. `BUILD_OMP=OFF` refers to LAMMPS's host OpenMP
option; Fortran GPU offload is independently enabled by `-mp=gpu`.

### GPU input and launch

Example `in.gpu`, using the same units/data as the CPU example:

```lammps
# Must precede read_data/create_box:
package gpu 1 neigh yes newton on split 1
units metal
atom_style atomic
read_data TiO2.data
pair_style accelnet/gpu auto Ti.nn.ascii O.nn.ascii
pair_coeff * *
neighbor 0.6 bin
neigh_modify every 10 delay 0 check yes
fix integrate all nve
run 100
```

```sh
nvidia-smi --query-gpu=name,uuid --format=csv
# Choose a UUID from the output:
export CUDA_VISIBLE_DEVICES=GPU-REPLACE-WITH-YOUR-DEVICE-UUID
export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_THREADS=1
"$LAMMPS_SOURCE/build-accelnet-gpu/lmp" -in in.gpu
# Two MPI ranks sharing the one visible GPU requested by package gpu 1:
mpirun -np 2 "$LAMMPS_SOURCE/build-accelnet-gpu/lmp" -in in.gpu
```

This uses an explicit GPU pair style and package command, so no `-sf gpu` or
additional `-pk gpu` override is needed. `OMP_NUM_THREADS=1` sets host thread
count; it does not limit GPU threads. To use multiple GPUs, expose the intended
devices and adjust the `package gpu` GPU count/rank layout accordingly.

`neigh no`, `neigh yes` and `neigh hybrid` are supported. **`neigh no` still
runs the potential on GPU**; it constructs neighbors on CPU. `neigh yes/hybrid`
use the GPU package's neighbor arrays. `newton on` and `split 1` are required.
The GPU-enabled executable can also run the CPU input using `pair_style
accelnet`; GPU visibility does not redirect that CPU pair style.

### GPU model files and supported observables

Use embedded `.nn`/`.nn.ascii` models, with LAMMPS types following the embedded
global species order. Direct `n2p2 /path/to/model ...` input is supported by the
**CPU** style only. For GPU, first convert with the Fortran converter (built
by the full AccelNet build above):

```sh
"$ACCELNET_SOURCE/build-gpu/bin/accelnet-model-converter-fortran" \
  n2p2-to-accelnet /path/to/n2p2-model converted
```

Then use the actual output filenames, for example:

```lammps
pair_style accelnet/gpu auto converted/H.nn.ascii converted/O.nn.ascii g5 moment
pair_coeff * *
```

The converter preserves supported n2p2 types 2/3/9/12/13/20--25, scaling and
network weights. Chebyshev and G5 selectors are independent: the leading mode
selects Chebyshev, the trailing `g5 MODE` selects G5. G5 auto uses eligible exact
integer powers 1--10 and at least 16 angular neighbors; `g5 moment` permits
powers 1--16 without that neighbor threshold. Other powers remain direct.
Types 13/21/24 have exact direct evaluation. Auto is not a hardware benchmark,
and moment is not universally faster.

Global energy/virial and per-atom energy are supported. Per-atom stress,
pair hybrid, r-RESPA, molecular topology, neighbor exclusions/include groups
and model serialization into restart files are not supported. See the
[GPU data-flow and validation guide](../../docs/lammps-gpu.md) for device
assignment, transfers, lifecycle checks and benchmark commands.

The [archived shared-kernel validation](../../docs/validation/energy-common-2026-09-27/README.md)
includes 21 LAMMPS CPU/GPU comparisons across 1/2 MPI ranks, neighbor modes,
orthogonal/triclinic cells and empty ranks. These results refer to the revisions
and compiler settings recorded there.
