# AccelNet

AccelNet is a Fortran library and command-line toolkit for evaluating
machine-learning interatomic potentials from **ænet and n2p2**. It computes
energies, analytic forces and configurational virials, provides Fortran/C APIs,
and integrates with LAMMPS on CPU and GPU. It is an inference package; training
remains in the upstream tools.

**Library version: 1.1.0.** The CPU/GPU methods and validation are documented in
[speedupmethods.md](speedupmethods.md), revision 1.13. The method is described
in the [AccelNet paper](https://arxiv.org/abs/2608.03280).

## Start here

| Task | Guide |
|---|---|
| Build and run on CPU | [Quick start](#cpu-quick-start) |
| Use an NVIDIA GPU | [GPU build](#gpu-build-and-execution), [target API](docs/openmp-target.md) |
| Run molecular dynamics | [LAMMPS interfaces](interfaces/lammps/README.md), [GPU package guide](docs/lammps-gpu.md) |
| Check supported models and direct/moment methods | [Coverage below](#supported-models-and-methods), [implementation status](docs/implementation-status.md) |
| Embed AccelNet | [Predictor APIs](AccelNetPredictor/README.md), [CSR batch API](docs/batch-api.md) |
| Convert model formats | [Fortran and Julia converters](AccelNetModelConverter/README.md) |
| Inspect equations and measured performance | [Optimization methods](speedupmethods.md), [validation](#validation-and-performance) |

## Supported models and methods

All supported production potential-inference entry points use the **same
maintained Fortran numerical kernels**. The ordinary CPU library compiles them
with OpenMP directives removed. The optional target library compiles them for
OpenMP host threads or GPU offload. Descriptor evaluation, the neural network,
force contraction and virial evaluation run on the selected backend in FP64.

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

## CPU quick start

Requirements: CMake 3.20+, a Fortran 2008 compiler, and a C compiler. GNU Fortran
11/13 are covered by CI; the recent CPU measurements use GNU Fortran 11.4.
The core libraries require no BLAS, LAPACK, MPI, ænet or n2p2 installation.
Python is used by some tests; external reference tools and models are optional.

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DN2P2_SCALING_EXECUTABLE=
cmake --build build --parallel
ctest --test-dir build -LE performance --output-on-failure
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
`--forces` to those forms. These commands use the ordinary CPU path, including
in a GPU-enabled build.

## GPU build and execution

The tested offload stack is **NVHPC 25.3 with NVIDIA H100 and RTX PRO 6000
Blackwell GPUs**. AMD/Intel GPU execution has not been validated. Building with
GNU `-fopenmp` alone enables host threading, not NVIDIA offload.

With `nvfortran` on `PATH`, build the optional target library and its smoke test:

```sh
cmake -S . -B build-gpu -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DACCELNET_BUILD_OPENMP_TARGET=ON \
  -DACCELNET_OPENMP_TARGET_FLAGS="-mp=gpu -gpu=cc90,cc120"
cmake --build build-gpu --parallel --target test_batch_target

# Inspect UUIDs and select an available device from this output:
nvidia-smi --query-gpu=name,uuid --format=csv
export CUDA_VISIBLE_DEVICES=GPU-REPLACE-WITH-YOUR-DEVICE-UUID
export OMP_TARGET_OFFLOAD=MANDATORY
export OMP_NUM_THREADS=1
build-gpu/bin/test_batch_target --quick
```

Use `-gpu=cc90` for H100 alone. To run the configured GPU suite, first build all
its executables with `cmake --build build-gpu --parallel`, then run
`ctest --test-dir build-gpu -L gpu --output-on-failure` with the same environment.
The target API rejects unintended CPU fallback.

| Entry point | Execution |
|---|---|
| CLI, structure/file and atomic Fortran/C APIs, ænet-compatible SFB, ordinary batch API | Common serial CPU kernels; OpenMP compiled out |
| Target API initialized with `use_host=.true.` | Common kernels with OpenMP host threads |
| Target API initialized for a GPU | Common kernels with OpenMP target offload |
| LAMMPS `accelnet` / `accelnet/gpu` | Common serial CPU / GPU kernels, respectively |

The [target API guide](docs/openmp-target.md) covers model snapshots, persistent
workspaces, CSR data, C handles and energy-only calls. A packed model must be
reinitialized after its source model changes. The independent former evaluator
is retained in [legacy/cpu-reference](legacy/cpu-reference/README.md) for explicit
reference/low-level compatibility calls; production inference has no legacy
fallback. Host force accumulation avoids fine-grained atomics; GPU force
scatter still uses atomic additions.

## LAMMPS

| Release | CPU | GPU |
|---|---|---|
| 4 Feb 2020 | `pair_style accelnet`, traditional make | Not provided |
| 29 Aug 2024 Update 4 | `pair_style accelnet`, CMake | `pair_style accelnet/gpu`, CUDA GPU package + Fortran OpenMP target |

CPU inputs can load an n2p2 directory directly, with elements in LAMMPS type order:

```lammps
pair_style accelnet n2p2 /path/to/model Ti O
pair_coeff * *
```

GPU inputs use embedded network files, including converted n2p2 models:

```lammps
# Before read_data/create_box:
package gpu 1 neigh yes newton on split 1
# After creating the box, with types matching the network species order:
pair_style accelnet/gpu auto Ti.nn.ascii O.nn.ascii
pair_coeff * *
```

The GPU adapter supports `neigh no/yes/hybrid`, global virial and per-atom energy.
It currently requires CUDA, NVHPC, double precision, `newton on` and `split 1`.
Direct n2p2-directory loading is available in the CPU pair style; for GPU use
[the Fortran converter](AccelNetModelConverter/README.md) first. Per-atom stress,
charge models and the other adapter restrictions are listed in the
[GPU integration guide](docs/lammps-gpu.md). See the
[LAMMPS README](interfaces/lammps/README.md) for installation and mode selection.

## Install and link

```sh
cmake --install build --prefix /path/to/install
```

The install contains libraries, Fortran module files, C headers, CLI programs
and CMake package files. Consumer projects can use:

```cmake
find_package(AccelNetPredictor CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE AccelNet::AccelNet)
# For an installation built with the optional target library:
# target_link_libraries(my_program PRIVATE AccelNet::Target)
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
| `N2P2_SCALING_EXECUTABLE` | Optional sibling executable | External n2p2 descriptor comparison; set empty to disable |
| `ACCELNET_PUBLIC_API_BASELINE_EXECUTABLE` | Empty | Enable the archived-baseline CPU energy performance test |

## Validation and performance

[CI](.github/workflows/tests.yml) checks GNU 11/13 Release builds and GNU 13
Debug/shared builds with runtime checks. It uses bundled/synthetic fixtures;
it does not run on a GPU. Optional reference comparisons and real-model tests
are described in the component READMEs.

The [revision 1.13 report](docs/validation/energy-common-2026-09-27/README.md)
records the `f655fb0` implementation (2026-09-27): 52 CPU tests, 22 bounds/runtime-check tests, 35 GPU tests, 21 LAMMPS
comparisons, H100 memcheck with zero errors, and host 2/8-thread checks. These
counts describe that configuration and available external fixtures, not every
fresh checkout.

CPU comparisons use **OpenMP compiled out on both sides**, identical drivers
and a fixed physical core. The latest real Ti/O energy-only sweep at 64/192/512
atoms is within about 0.7--3.1% of the former CPU evaluator. On Blackwell, the
same model's GPU energy-only auto/moment time is 13--14% shorter than the prior
common implementation; the synthetic Chebyshev moment force case is 4.5% slower.
These are separate baselines and workloads, not a universal CPU/GPU speedup.

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
