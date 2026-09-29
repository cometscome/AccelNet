# Target C model loading validation (2026-09-28)

Source: changes on top of AccelNet 1.1.0,
`015128c66f45f528972b9ebf8a8ccd74eaab2366`.

The new C constructors select Chebyshev version 0/1/10 or load n2p2 model
directories directly. Embedded model files can be supplied in any order.
Existing constructor signatures retain version 0.
All constructors support explicit host execution; nonnegative device IDs
continue to require actual GPU execution. The species query lets a caller map
its atom types to the model's one-based ordering.

## Results

| Configuration | Checks | Result |
|---|---|---|
| GNU Fortran 11.4, Debug/shared, `-fcheck=all -fbacktrace`, OpenMP host | Non-GPU/non-performance tests | All 41 passed, including the follow-up below |
| GNU Fortran 11.4, Release/shared, target serial (OpenMP compiled out) | Complete `target-host` label | 9/9 passed |
| NVHPC 25.3, Release/shared, `-mp=gpu -gpu=cc90,cc120`, H100 NVL | Complete self-contained `gpu` and `target-host` labels | 39/39 passed |
| Same NVHPC library, RTX PRO 6000 Blackwell | New C loading test using GNU-generated reference fixtures | Passed |
| GNU Fortran caller + NVHPC shared library, H100 NVL | Direct n2p2 load and atomic energies/forces/all-nine-component virial | Passed; maximum absolute error `8.3266726846886741e-17` |

Counts overlap and should not be added. Two GNU tests initially could not run
because their executables had not been built; after building `test_neighbor_images`
and `accelnet-n2p2-benchmark`, both passed. The final C-loading/error/package tests
were also rerun successfully. Optional external golden/reference fixtures were
disabled. The standalone mixed-compiler check is included as
[gnu_client.f90](gnu_client.f90); it imports only `ISO_C_BINDING` and uses
`bind(C)` declarations, with no NVHPC module files. The PIMD adapter itself is
not changed or tested by this report.

The new regression test compares against `evaluate_batch_reference`, the
independent legacy CPU evaluator. It covers all three Chebyshev conventions
and evaluation modes, four n2p2 fixtures (including different element network
depths and angular/scaling/normalization cases), all four G5 modes, resident
workspace reuse, row subsets, independently live handles, species queries,
and invalid constructor/query arguments. The checked version-0 and version-1
fixtures produce different energies, so ignoring the new version parameter
does not pass. The GNU CI matrix now runs the host case automatically.

The shared-library tests also cover 36 malformed n2p2 cases, repeating each
failure four times and loading a valid model after each case while retaining
another live handle. They check unsupported models, malformed topology,
descriptor parameters, nonfinite values, missing/invalid weights, and
missing/invalid scaling (including out-of-range indices). Errors return a
diagnostic and clear outputs without terminating the caller.

`predictor_target_cmake_host` and `predictor_target_cmake_gpu` install only the
`TargetC` component, move the installation, and configure a C-only client with
`find_package(AccelNetC)`. They then enable GNU Fortran and run the Fortran C ABI
client. Both link `AccelNet::TargetC`, with no native Fortran module paths or
NVHPC flags. The in-tree C client also uses that target and links with `cc`.
The serial GNU library's `ldd` output contains neither OpenMP nor GPU runtimes.

## Reproduction

From the repository root, configure a host build:

```sh
cmake -S . -B /tmp/accelnet-c-loading-host \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_Fortran_FLAGS="-fcheck=all -fbacktrace" -DBUILD_SHARED_LIBS=ON \
  -DACCELNET_BUILD_OPENMP_TARGET=ON \
  -DACCELNET_OPENMP_TARGET_FLAGS="-fopenmp -foffload=disable -ffree-line-length-none" \
  -DACCELNET_DESCRIPTORS_BUILD_N2P2_REFERENCE=OFF \
  -DACCELNET_DESCRIPTORS_BUILD_AENET_REFERENCE=OFF \
  -DACCELNET_PREDICTOR_GOLDEN_DIR= -DN2P2_SCALING_EXECUTABLE=
cmake --build /tmp/accelnet-c-loading-host --parallel 8 --target \
  test_target_c_loading write_target_c_loading_fixtures test_batch_target
ctest --test-dir /tmp/accelnet-c-loading-host -L target-host --output-on-failure
```

For the GPU build, use a separate directory `/tmp/accelnet-c-loading-gpu`,
Release, NVHPC's `nvfortran`, shared libraries, and target flags
`-mp=gpu -gpu=cc90,cc120`. Keep the external-reference options above.
Build the three targets above plus `test_target_c` and `write_target_fixtures`.
Select the intended GPU with `CUDA_VISIBLE_DEVICES` (prefer its UUID), then run:

```sh
OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 \
  ctest --test-dir /tmp/accelnet-c-loading-gpu -L gpu --output-on-failure

# Compile and link the caller with GNU Fortran, not nvfortran:
gfortran -Wall -fcheck=all \
  docs/validation/target-c-loading-2026-09-28/gnu_client.f90 \
  -L/tmp/accelnet-c-loading-gpu/lib \
  -Wl,-rpath,/tmp/accelnet-c-loading-gpu/lib -laccelnet_target \
  -o /tmp/accelnet-gnu-target-client
OMP_TARGET_OFFLOAD=MANDATORY OMP_NUM_THREADS=1 /tmp/accelnet-gnu-target-client \
  AccelNetPredictor/test/data/n2p2-virial-angular \
  /tmp/accelnet-c-loading-host/AccelNetPredictor/target-c-loading-fixtures/n2p2-virial-angular.ref
```

The NVHPC runtimes must be available to the dynamic loader. GPU model loading
does not require conversion of the n2p2 directory. Invalid model contents return
status/message errors to C callers. Existing Fortran callers still get fatal
errors unless they request the optional status output; see the
[API contract](../../openmp-target.md#c-model-loading) and
[C-only installation guide](../../openmp-target.md#c-only-cmake-package).
