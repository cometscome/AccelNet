# AccelNetPredictor

Fortran energy/force predictor built on the single canonical
`AccelNetDescriptors` implementation. No descriptor source is copied into this
package.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Energy-only CLI:

```sh
build/accelnet-predict 2 Ti.fingerprint.stp O.fingerprint.stp \
  Ti.nn.ascii O.nn.ascii structure.xsf
```

Standard n2p2 2G-HDNNP model directories can also be loaded directly. The
directory must contain `input.nn`, `scaling.data`, and the conventional
`weights.%03d.data` files (where the number is the element atomic number):

```sh
build/accelnet-predict --n2p2 /path/to/model structure.xsf
```

From Fortran, use `load_predictor_from_n2p2(model_directory, model)` or reload
an existing object with `call model%reload_from_n2p2(model_directory)`.
Imported models currently support 2G symmetry-function types 2, 3, and 9,
cosine `cutoff_type 1` with `cutoff_alpha 0`, all n2p2 activation functions,
the standard scaling modes, data-set energy normalization, and atomic energy
offsets. Per-element network topologies, `normalize_nodes`, weighted/compact
symmetry functions, non-cosine cutoffs, and 4G/charge models are rejected with
an error rather than evaluated with different semantics. XSF coordinates and
the parameters in `input.nn` must use the same physical length unit; returned
energies and forces use the model's physical energy and length units.

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

Append `--forces` to print analytic Cartesian forces. A warm-up-separated
whole-inference benchmark is also installed:

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

The Fortran library also accepts an already-loaded structure, avoiding XSF
I/O. The force array is caller-owned and can be reused between calls:

```fortran
use iso_fortran_env, only: real64
use accelnet_descriptors, only: atomic_structure, read_xsf
use accelnet_predictor, only: predictor_model

type(atomic_structure) :: structure
real(real64) :: energy
real(real64), allocatable :: forces(:, :)

call read_xsf("structure.xsf", model%species_names, structure)
allocate(forces(3, structure%natoms))
call model%predict_energy(structure, energy)
call model%predict_energy_forces(structure, energy, forces)
```

## aenet-style atomic API

The `accelnet` Fortran module provides the complete aenet-style global API for
codes that already evaluate one local atomic environment at a time. It
includes initialization/finalization, ASCII and binary potential loading,
load-state and cutoff metadata, atom-type conversion, free-atom energies,
atomic energies and additive forces, a built-in neighbor list, and the
standalone structural-fingerprint basis routines. The C ABI is declared by
`include/accelnet.h`.

```fortran
use accelnet

character(len=2) :: species(2) = ["Ti", "O"]
integer :: stat

call accelnet_init(species, stat)
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

The exported limits `accelnet_nsf_max`, `accelnet_nnb_max`,
`accelnet_Rc_min`, and `accelnet_Rc_max` have the same purpose as their aenet
counterparts. API failures are reported through `ACCELNET_OK` and the
`ACCELNET_ERR_*` status values. The neighbor-list procedures intentionally
retain aenet's status-free signatures.

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
For Behler2011 networks, embedded functions are placed in ænet's canonical
species/function ordering rather than their textual metadata order, matching
the input order used by the trained neural network and scaling arrays.
# AccelNetPredictor
