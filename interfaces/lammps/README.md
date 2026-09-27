# LAMMPS interfaces

This directory contains pair styles, build integration and tests for using
AccelNet in LAMMPS. Supply your trained models separately. Numerical results
and timings are recorded in the [validation reports](../../docs/implementation-status.md#validation-and-performance-scope).

Two LAMMPS generations are supported:

| Directory | LAMMPS version | Build system | Pair style |
|---|---|---|---|
| `4Feb2020` | 4 Feb 2020 | traditional make | `accelnet` |
| `29Aug2024` | 29 Aug 2024 Update 4 | CMake | `accelnet`, optional `accelnet/gpu` |

Both pair styles use the common Fortran inference kernels. `accelnet` is
serial within each MPI rank; it does not switch to GPU when OpenMP threads or
GPU visibility variables are set. `accelnet/gpu` uses the LAMMPS CUDA GPU
package together with Fortran OpenMP target. See the
[GPU build and execution guide](../../docs/lammps-gpu.md).

## CPU build

First build the static AccelNet libraries from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF
cmake --build build --parallel
```

The required outputs are `build/lib/libaccelnet.a` and
`build/lib/libAccelNetDescriptors.a`.

## LAMMPS 29Aug2024 Update 4

Starting in the AccelNet repository root, copy the package and CMake module into
a clean LAMMPS source tree and apply the package-registration patch:

```sh
cp -R interfaces/lammps/29Aug2024/ACCELNET /path/to/lammps/src/
cp interfaces/lammps/29Aug2024/cmake/ACCELNET.cmake \
  /path/to/lammps/cmake/Modules/Packages/
patch -d /path/to/lammps -p1 \
  < interfaces/lammps/29Aug2024/lammps-cmake.patch
```

Configure and build:

```sh
cmake -S /path/to/lammps/cmake -B /path/to/lammps/build-accelnet \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=ON \
  -D PKG_ACCELNET=ON \
  -D ACCELNET_DIR="$PWD"
cmake --build /path/to/lammps/build-accelnet --parallel
```

`ACCELNET.cmake` enables Fortran and propagates its runtime libraries and
search paths. Use the same Fortran compiler family as the AccelNet build.
`ACCELNET_DIR` may instead point to a build directory containing `lib/`.

## LAMMPS 4Feb2020

Copy the legacy package into the LAMMPS source tree:

```sh
cp -R interfaces/lammps/4Feb2020/USER-ACCELNET /path/to/lammps/src/
mkdir -p /path/to/lammps/lib/accelnet/include
mkdir -p /path/to/lammps/lib/accelnet/lib
cp AccelNetPredictor/include/accelnet.h \
  /path/to/lammps/lib/accelnet/include/
cp build/lib/libaccelnet.a build/lib/libAccelNetDescriptors.a \
  /path/to/lammps/lib/accelnet/lib/
```

Then use the traditional LAMMPS build:

```sh
cd /path/to/lammps/src
make yes-user-accelnet
make mpi -j8
```

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
with the n2p2 form. The 29Aug2024 interface additionally accepts the trailing
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

The 29Aug2024 interface also accepts an optional G5 evaluation mode after all
potential files:

```lammps
pair_style accelnet O.nn Ti.nn g5 direct
pair_style accelnet O.nn Ti.nn g5 moment
```

`direct` disables the integer-zeta G5 moment path. `moment` forces that path
without applying the automatic neighbor-count threshold. Omitting the option,
or selecting `g5 auto`, preserves the default automatic behavior.


## GPU package interface (29Aug2024 Update 4)

Use the [GPU installation guide](../../docs/lammps-gpu.md) to build
`libaccelnet_target` and install the adapter with `29Aug2024/install.py`.
The validated stack uses NVHPC 25.3, CUDA, `GPU_PREC=double`, and H100 or
Blackwell GPUs. The installation helper also applies the required triclinic
sorting patch; copying only the CPU package above is insufficient for GPU use.

```lammps
# Before read_data/create_box:
package gpu 1 neigh yes newton on split 1
# After creating the box, with atom types in the network's species order:
pair_style accelnet/gpu auto Ti.nn.ascii O.nn.ascii
pair_coeff * *
```

`neigh no`, `neigh yes` and `neigh hybrid` are supported. Direct n2p2-directory
loading is a CPU pair-style feature; GPU input uses embedded network files.
The [Fortran converter](../../AccelNetModelConverter/README.md) preserves the
supported weighted/compact types as well as types 2/3/9. Follow the converted
model's embedded global species order when assigning LAMMPS atom types.

Chebyshev and G5 selectors are independent in the 2024 CPU/GPU interfaces:

```lammps
pair_style accelnet/gpu moment Ti.nn.ascii O.nn.ascii g5 auto
```

Chebyshev auto estimates work from angular order, moment count and neighbors.
G5 auto uses exact integer powers 1--10 and at least 16 angular neighbors;
`g5 moment` allows powers 1--16 without that neighbor-count threshold. Other
powers remain direct. Types 13/21/24 are implemented with exact direct pairs.
Neither selector benchmarks the hardware, and moment is not always faster.

The adapter supports global energy/virial and per-atom energies. Per-atom
stress, pair hybrid and the other unsupported configurations are listed in the
[GPU scope](../../docs/lammps-gpu.md). The
[latest shared-kernel validation](../../docs/validation/energy-common-2026-09-27/README.md)
includes 21 LAMMPS CPU/GPU comparisons across 1/2 MPI ranks, neighbor modes,
orthogonal/triclinic cells and empty ranks. Commands for numerical and MD-loop
performance checks are in the GPU guide.
