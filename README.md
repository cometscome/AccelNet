# AccelNet

AccelNet is a single-release Fortran monorepo containing the descriptor,
prediction, and model-conversion components.

## Build everything

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Executables are written to `build/bin` and libraries to `build/lib`. The
default static build produces `libaccelnet.a`, which contains the Fortran
`accelnet` module and the C ABI declared in `AccelNetPredictor/include/accelnet.h`.
Build shared libraries instead with `-DBUILD_SHARED_LIBS=ON`.

The API follows aenet's atomic-environment calling style, with `accelnet_*`
procedure names. See `AccelNetPredictor/README.md` for Fortran and C examples.
Chebyshev angular evaluation is runtime-selectable through
`accelnet_set_chebyshev_evaluation`: `AUTO` chooses by neighbor count, while
`DIRECT` and `MOMENT` force the corresponding exact algorithm.

## Install

```sh
cmake --install build --prefix /path/to/install
```

Installed CMake consumers can use either target name after calling
`find_package(AccelNetPredictor CONFIG REQUIRED)`:

```cmake
target_link_libraries(my_program PRIVATE AccelNet::AccelNet)
```

`AccelNetPredictor::AccelNetPredictor` remains available for compatibility.

Reference comparisons against the historical source tree are disabled by
default. Enable them with `-DACCELNET_BUILD_REFERENCE_TESTS=ON` and, when
needed, set `ACCELNET_ORIGINAL_DIR` explicitly.

## LAMMPS

Compile-only interfaces for LAMMPS 4Feb2020 and 29Aug2024 Update 4 are provided
in [`interfaces/lammps`](interfaces/lammps/README.md).  Both expose AccelNet as
`pair_style accelnet`; trained potentials and calculation data are not included
in the interface directories.
