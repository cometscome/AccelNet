# LAMMPS interfaces

This directory contains only the source files and build integration needed to
compile AccelNet into LAMMPS.  It intentionally excludes trained potentials,
atomic structures, trajectories, benchmark output and validation logs.

Two LAMMPS generations are supported:

| Directory | LAMMPS version | Build system | Pair style |
|---|---|---|---|
| `4Feb2020` | 4 Feb 2020 | traditional make | `accelnet` |
| `29Aug2024` | 29 Aug 2024 Update 4 | CMake | `accelnet` |

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

`ACCELNET.cmake` enables the Fortran language so that the correct GNU Fortran
runtime libraries and search paths are propagated to the final LAMMPS link.

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

For both versions, the LAMMPS input syntax is:

```lammps
pair_style accelnet H.ann O.ann
pair_coeff * *
```

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
