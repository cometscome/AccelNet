# AccelNet

AccelNet is a Fortran library and command-line toolkit for evaluating
machine-learning interatomic potentials. It provides a common descriptor and
inference implementation for supported ænet and n2p2 models, an aenet-style
Fortran/C atomic API, model-conversion tools, and a LAMMPS `pair_style`.

AccelNet is an inference package. It does not train potentials or replace the
training and data-preparation tools provided by ænet, ænet-PyTorch, or n2p2.

## Features

- Read ænet/AccelNet ASCII and compatible native-binary neural networks.
- Read supported n2p2 2G-HDNNP model directories without conversion.
- Evaluate energies and analytic Cartesian forces from XSF structures.
- Read molecular and periodic structures from n2p2 `input.data` files.
- Provide object-based Fortran interfaces and an aenet-style Fortran/C atomic
  API.
- Convert supported models between n2p2 and AccelNet ASCII representations.
- Run supported ænet and n2p2 models from LAMMPS with
  `pair_style accelnet`.

The descriptor definitions implemented in AccelNet are also used by
`AccelNet.jl`, a training package that has not yet been publicly released.
`AccelNet.jl` and its training functionality are not included in this
repository.

## Compatibility overview

The following table summarizes compatibility with ænet 2.0.4 and n2p2 2.3.0.
"Conditional" means that only the model families and settings listed in the
Notes column are accepted; unsupported settings are rejected with an error.

| Feature | ænet / AccelNet | n2p2 | Notes |
|---|---|---|---|
| Model input | Supported | Conditional | ænet/AccelNet ASCII and compatible native binary; n2p2 short-range 2G model directories |
| Descriptors | Supported | Conditional | ænet Chebyshev and Behler G1--G5, plus AccelNet LJ; n2p2 SF types 2, 3, and 9 |
| Activation functions | Supported | Supported | ænet native codes 0--4; all n2p2 2.3.0 activation characters |
| Cutoff functions | Supported | Supported | Standard ænet metadata; n2p2 cutoff types 0--8 and the [Mori *et al.* fractional cutoff](https://doi.org/10.1103/PhysRevMaterials.7.063605) as extension type 9 |
| Descriptor scaling | Supported | Supported | Affine ænet scaling; n2p2 scale, center, scale+center, and sigma modes |
| Energy normalization and atomic reference energies | Supported | Supported | Applied consistently to energies and analytic forces |
| Energy and force inference | Supported | Supported | XSF files and in-memory structures |
| Structure input | XSF, `predict.in` | XSF, `input.data` | n2p2 `input.data` supports multiple molecular or periodic structures |
| Fortran and C atomic API | Supported | Supported | n2p2 directories can be loaded directly without conversion |
| LAMMPS `pair_style accelnet` | Supported | Supported | n2p2 directories can be loaded directly with an explicit atom-type mapping |
| Model conversion | Conditional | Conditional | Only features representable by both formats are converted |
| Potential training | Not supported | Not supported | Use ænet, ænet-PyTorch, or n2p2 for training |
| 4G/Q, charge, weighted, or compact models | Not applicable | Not supported | n2p2 4G/Q and SF types 12, 13, and 20--25 are rejected |
| Per-element topology and `normalize_nodes` | Not applicable | Supported | n2p2 global defaults and per-element overrides are accepted; node normalization is folded into weights and biases |

The exact accepted syntax, formulas, ordering conversions, and tested cases are
documented in
[`docs/model-compatibility.md`](docs/model-compatibility.md).

## Requirements

The standard build requires:

- CMake 3.20 or newer;
- a Fortran 2008-compatible compiler;
- a C compiler.

GNU Fortran and Intel Fortran/IntelLLVM compiler flags are configured by the
build system. GNU Fortran 15 or newer is recommended for performance. Use the
same Fortran compiler family for AccelNet and Fortran applications that consume
its module files. The C API still requires the corresponding Fortran runtime
when statically linked.

The core library does not require BLAS, LAPACK, MPI, ænet, or n2p2. Those
packages are needed only for selected reference comparisons or external
integration builds. The optional Julia model converter requires Julia 1.10.

## Build and test

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Executables are written to `build/bin` and libraries to `build/lib`. The
default static build produces:

- `libAccelNetDescriptors.a`;
- `libaccelnet.a`;
- `accelnet-descriptor`;
- `accelnet-setup-descriptor`;
- `accelnet-predict`;
- `accelnet-model-converter-fortran`.

Useful CMake options are:

| Option | Default | Purpose |
|---|---:|---|
| `BUILD_SHARED_LIBS` | `OFF` | Build shared instead of static libraries |
| `BUILD_TESTING` | `ON` | Build the standard test suite |
| `ACCELNET_BUILD_REFERENCE_TESTS` | `OFF` | Compare against an external historical AccelNet tree |

For a shared-library build:

```sh
cmake -S . -B build-shared -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON
cmake --build build-shared --parallel
```

Reference tests against separate ænet and n2p2 source trees are documented in
[`AccelNetDescriptors/README.md`](AccelNetDescriptors/README.md) and
[`AccelNetPredictor/README.md`](AccelNetPredictor/README.md). They are not
required for normal use.

## Run predictions

### ænet/AccelNet networks

An energy calculation from descriptor setup files, neural networks, and an
XSF structure can be run with:

```sh
build/bin/accelnet-predict 2 \
  Ti.fingerprint.stp O.fingerprint.stp \
  Ti.nn.ascii O.nn.ascii structure.xsf
```

Original AccelNet `predict.in` files are also accepted:

```sh
build/bin/accelnet-predict predict.in
```

Append `--forces` to an applicable prediction command to print analytic
Cartesian forces. Native binary networks use Fortran sequential-unformatted
records and are therefore less portable than ASCII networks.

### n2p2 model directories

A supported n2p2 directory contains `input.nn`, `weights.%03d.data`, and,
when required by the model, `scaling.data`. It can be evaluated directly:

```sh
build/bin/accelnet-predict --n2p2 /path/to/model structure.xsf
```

One or more structures in n2p2 `input.data` format can be read with:

```sh
build/bin/accelnet-predict --n2p2-data /path/to/model input.data
```

AccelNet reads coordinates, element names, and zero or three lattice vectors.
Reference energies, stored forces, charges, and comments in `input.data` are
not used for inference.

The n2p2 loader accepts global network defaults together with
`element_hidden_layers_short`, `element_nodes_short`, and
`element_activation_short` overrides. It also supports `normalize_nodes`.
Normalization is folded into each layer's weights and biases during model
loading. The same support is available through the Fortran/C atomic APIs and
LAMMPS direct directory loading.

## Install and link

Install the libraries, module files, C header, executables, and CMake package
files with:

```sh
cmake --install build --prefix /path/to/install
```

An installed CMake project can link the predictor with:

```cmake
find_package(AccelNetPredictor CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE AccelNet::AccelNet)
```

`AccelNetPredictor::AccelNetPredictor` is also available. The C interface is
declared in
[`AccelNetPredictor/include/accelnet.h`](AccelNetPredictor/include/accelnet.h),
and Fortran/C examples are provided in
[`AccelNetPredictor/README.md`](AccelNetPredictor/README.md).

## LAMMPS

Interfaces are provided for these LAMMPS releases:

| LAMMPS release | Integration |
|---|---|
| 4Feb2020 | traditional make package |
| 29Aug2024 Update 4 | CMake package |

First build static AccelNet libraries:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF
cmake --build build --parallel
```

For LAMMPS 29Aug2024 Update 4, copy the package and CMake module into a clean
LAMMPS source tree and apply the supplied registration patch:

```sh
cp -R interfaces/lammps/29Aug2024/ACCELNET /path/to/lammps/src/
cp interfaces/lammps/29Aug2024/cmake/ACCELNET.cmake \
  /path/to/lammps/cmake/Modules/Packages/
patch -d /path/to/lammps -p1 \
  < interfaces/lammps/29Aug2024/lammps-cmake.patch

cmake -S /path/to/lammps/cmake -B /path/to/lammps/build-accelnet \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=ON \
  -D PKG_ACCELNET=ON \
  -D ACCELNET_DIR="$PWD"
cmake --build /path/to/lammps/build-accelnet --parallel
```

An ænet/AccelNet model is selected with:

```lammps
pair_style accelnet Ti.nn.ascii O.nn.ascii
pair_coeff * *
```

An n2p2 directory can be loaded without conversion. Element names are given in
LAMMPS atom-type order:

```lammps
# LAMMPS type 1 = Ti, type 2 = O
pair_style accelnet n2p2 /path/to/model Ti O
pair_coeff * *
```

Every MPI rank must be able to read the same model directory. Full instructions
for both LAMMPS releases, compiler-runtime linking, and Chebyshev/G5 evaluation
modes are in [`interfaces/lammps/README.md`](interfaces/lammps/README.md).

## Supported models and limitations

AccelNet supports ænet Chebyshev and Behler-style descriptors used by the
documented network formats. Its n2p2 loader supports short-range 2G models
with symmetry-function types 2, 3, and 9, supported cutoff functions, network
activations, scaling, energy normalization, and atomic reference energies.

The following are not currently supported:

- training of neural-network potentials;
- n2p2 4G/charge models;
- weighted and compact n2p2 symmetry functions;
- direct input of general ASE formats or ænet training-set files.

Unsupported model settings are rejected rather than silently approximated.
See [`docs/model-compatibility.md`](docs/model-compatibility.md) for the exact
accepted syntax, descriptor formulas, activation and cutoff mappings,
conversion restrictions, and reference-test coverage.

## Model conversion

The Fortran converter is built with the main project. For example:

```sh
# n2p2 to AccelNet ASCII
build/bin/accelnet-model-converter-fortran n2p2-to-accelnet \
  /path/to/n2p2-model /path/to/accelnet-output

# AccelNet ASCII or compatible native binary to n2p2
build/bin/accelnet-model-converter-fortran accelnet-to-n2p2 \
  /path/to/n2p2-output Ti.nn O.nn
```

Conversion is limited to model features representable by both formats. See
[`AccelNetModelConverter/README.md`](AccelNetModelConverter/README.md) for the
Fortran and Julia interfaces and their current restrictions. The Fortran
converter preserves per-element topology in both directions and absorbs
`normalize_nodes` exactly when importing n2p2 models.

## Citation

A software citation and archival DOI will be added with the first public
release. Until then, please cite the publications relevant to the model and
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

## License

Original AccelNet code is released under the [MIT License](LICENSE). The following
third-party-derived files retain their original licenses and copyright
notices:

- `AccelNetDescriptors/src/accelnet_legacy_lcl.f90`: Mozilla Public License
  2.0;
- `interfaces/lammps/**/pair_accelnet.cpp` and `pair_accelnet.h`: GNU General
  Public License version 2.

The linked-cell file is the only ænet-derived source file retained under the
MPL-2.0 in the core library. A LAMMPS executable built with the supplied pair
style remains subject to the LAMMPS GPL terms. Full license texts and provenance
are provided in [`LICENSES/`](LICENSES/) and
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
