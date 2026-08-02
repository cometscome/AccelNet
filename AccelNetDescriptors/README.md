# AccelNetDescriptors

`AccelNetDescriptors` is a standalone Fortran implementation of the
Chebyshev, Behler 2011, and Lennard-Jones atomic-environment descriptors
used by AccelNet. It can read AccelNet `.fingerprint.stp` setup files
directly.

The package supports all Chebyshev compatibility modes accepted by AccelNet:

- Chebyshev `version=0`, `version=1`, and `version=10`
- Lennard-Jones (`sum(r^-6)`, `sum(r^-12)`) values and analytical derivatives
- Behler 2011 G1, G2, G3, G4, and G5 values and analytical derivatives
- a statically dispatched descriptor model that can concatenate Chebyshev,
  Behler, and LJ blocks
- one or more chemical species
- descriptor values and analytical coordinate derivatives
- periodic and non-periodic structures
- periodic images in cells smaller than the cutoff
- XSF input for the command-line regression driver

The original AccelNet source tree is never modified. When
`ACCELNET_DESCRIPTORS_BUILD_REFERENCE=ON`, the original Chebyshev and linked
cell-list sources are compiled in the build directory as read-only reference
libraries. The tests cover two distinct reproducibility levels:

- with identical neighbor coordinates and ordering, the value-only kernels
  produce identical Float64 results;
- through each implementation's own neighbor-list path, descriptor files
  are compared coefficient by coefficient with a scaled tolerance of
  `5e-13`. The tolerance accounts only for summation-order roundoff.

The reference build never copies generated files or module files into the
original AccelNet directory.

## Build and test

```sh
cmake -S . -B build \
  -DACCELNET_ORIGINAL_DIR=/path/to/AccelNet-a76ec3a83c8b0246328e9cca4817cc0afc84c925
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The standalone Fortran G4 derivative-layout benchmark is available as
`build/benchmark-g4-derivative`. It checks numerical equivalence before timing
the array-expression, contiguous fused-scalar, and member-first SoA kernels.

If the sibling `TiO2-xsf` corpus is present, CMake automatically adds a
six-structure cross-implementation regression and a 95-atom performance
gate. Without that corpus, the self-contained 24-atom fixture is still
tested.

## Generate a descriptor file

The setup-aware CLI takes one AccelNet setup file per central species:

```sh
build/accelnet-setup-descriptor descriptors.hex 2 \
  Ti.fingerprint.stp O.fingerprint.stp \
  structure0001.xsf structure0002.xsf
```

The setup-file order defines the global species indices used while reading
the XSF files. In this example, `Ti=1` and `O=2`. Each setup's `ENV` order is
mapped independently to the local species indices expected by its
descriptor model. Setup parsing happens once before evaluating structures,
so it is outside the performance-critical atom and neighbor loops.

The older fixed-configuration CLI is retained for regression testing. It
uses the Ti/O reference setup from AccelNet:
`radial_Rc=6.5`, `radial_N=20`, `angular_Rc=5.0`, `angular_N=6`,
`version=0`.

```sh
build/accelnet-descriptor descriptors.hex structure0001.xsf
```

Floating-point values are written as hexadecimal IEEE-754 bit patterns.
This format is for exact regression testing, not as a training-data format.

## AccelNet setup-file support

The `accelnet_setup` module reads `DESCR`, `ATOM`, `ENV`, `RMIN`, `BASIS`,
`SYMMFUNC`, and `FUNCTIONS` sections. It currently supports:

- `BASIS type=Chebyshev` with `version=0`, `version=1`, or `version=10`,
  including AccelNet's omitted-key defaults;
- `BASIS type=LJ`;
- `SYMMFUNC type=Behler2011` and `FUNCTIONS type=Behler2011`, with G1
  through G5;
- `BASIS type=multi` containing Chebyshev and LJ blocks.

Behler functions are reordered exactly as AccelNet orders its output
coefficients, rather than simply retaining their order in the setup file.
Unsupported `multi` members such as Spline and Spherical are rejected
explicitly; they will be added as their kernels are implemented and
validated.

All setup descriptor families accept the same cutoff selection as n2p2. Both
keys are optional; omitting them preserves the historical default (cosine for
Behler and Chebyshev, hard truncation for LJ):

```text
FUNCTIONS type=Behler2011 cutoff_type=7 cutoff_alpha=0.2
BASIS type=Chebyshev cutoff_type=7 cutoff_alpha=0.2
BASIS type=LJ cutoff_type=7 cutoff_alpha=0.2
```

`cutoff_type` supports all n2p2 cutoff functions plus the fractional cutoff
from Mori et al. for Behler, Chebyshev, and LJ:
`0` hard, `1` cosine,
`2` unnormalized tanh, `3` normalized tanh, `4` exponential, and polynomial
orders `1` through `4` as types `5` through `8`. Type `9` is the fractional
cutoff $X^2/(1+X^2)$ with $X=(r-R_c)/(\alpha R_c)$, so
`cutoff_alpha=h/Rc` and must be positive. Otherwise `cutoff_alpha` must satisfy
`0 <= alpha < 1`. As in n2p2, alpha controls the inner flat region for cosine,
exponential, and polynomial cutoffs and is ignored by hard and tanh cutoffs.
The selection applies to Behler G1--G5 values and analytical derivatives,
including both the direct and Cartesian-moment G5 implementations. Chebyshev
and LJ retain their descriptor-specific native cutoff conventions.

```fortran
use accelnet_setup

type(descriptor_setup), allocatable :: setups(:)
character(len=16), allocatable :: species(:)
character(len=256) :: files(2)

files = [character(len=256) :: "Ti.fingerprint.stp", "O.fingerprint.stp"]
call read_accelnet_setup_set(files, setups, species)
```

Each returned `descriptor_setup` contains the central species, description,
minimum distance, environment species, global-to-local species mapping, and
a ready-to-evaluate `descriptor_model`. Use
`setup%map_species(global_neighbors, local_neighbors)` before passing
neighbors to the model.

## Descriptor model API

`accelnet_descriptor_models` provides the extensible public model. Dispatch
occurs once per descriptor block, outside all neighbor and neighbor-pair
loops, so the original specialized Chebyshev kernel remains the fast path.

```fortran
use accelnet_descriptors
use accelnet_lj
use accelnet_behler
use accelnet_descriptor_models

type(descriptor_config) :: chebyshev
type(lj_config) :: lj
type(behler_config) :: behler
type(descriptor_model) :: model

call initialize_config(chebyshev, 2, 6.5d0, 20, 5.0d0, 6, version=0)
call initialize_lj_config(lj, 2, 6.5d0)
call initialize_behler_config(behler, 2)
call add_g1(behler, 1, 6.0d0)
call add_g2(behler, 2, 6.0d0, 0.0d0, 0.5d0)
call add_g4(behler, 1, 2, 5.0d0, 1.0d0, 2.0d0, 0.1d0)
call add_chebyshev(model, chebyshev)
call add_behler(model, behler)
call add_lj(model, lj)
call evaluate_model_values(model, displacements, neighbor_species, values)
```

For a manually constructed multi-species `version=10` configuration, also
pass `central_type_index`. The setup reader supplies it automatically:

```fortran
call initialize_config(chebyshev, 2, 6.5d0, 20, 5.0d0, 6, &
                       version=10, central_type_index=1)
```

AccelNet's version 10 source obtains the weighted center channel by indexing
the mapped neighbor-type array with this central type ID. This package
preserves that observable behavior for compatibility. Consequently, exact
version 10 reproduction requires the same neighbor ordering.
`build_neighbor_list` therefore uses the AccelNet/AENET linked-cell ordering
for periodic structures. The vendored implementation retains its MPL-2.0
license header and is namespaced so it does not collide with an AccelNet
reference build. Version 1 and version 10 setup-driven descriptor files are
compared against the untouched AccelNet implementation coefficient by
coefficient. Version 0 and version 1 are otherwise independent of neighbor
ordering except for floating-point summation roundoff.

The order in which `add_g1` through `add_g5` are called defines the Behler
coefficient order for manually constructed models. The setup reader creates
one model per central species and applies AccelNet's coefficient ordering
automatically.

Chebyshev angular evaluation has three runtime modes. `AUTO` uses direct pair
enumeration below 16 angular neighbors and exact Cartesian moments otherwise;
`DIRECT` and `MOMENT` force one path for every environment. Use
`set_chebyshev_evaluation(config, mode)` with
`CHEBYSHEV_EVALUATION_AUTO`, `CHEBYSHEV_EVALUATION_DIRECT`, or
`CHEBYSHEV_EVALUATION_MOMENT`.

Integer-exponent G5/type-9 functions support exact Cartesian-moment
evaluation, including nonzero radial shifts. Selection deliberately uses only
two fixed bounds: moment evaluation is used with at least 16 neighbors and for
integer angular orders `zeta <= 10`. Fewer neighbors and orders `zeta >= 11` use
direct pairs. A model may therefore evaluate its low-order G5 functions with
moments and its high-order functions directly in the same call. Encountering
a high-order function emits one warning per configuration. The default order
limit of 10 is intended to cover broader descriptor compositions such as TiO2
and H2O; the precise direct/moment crossover depends on the number of species
and G5 functions. `G5_EVALUATION_MOMENT_FORCE` bypasses only the neighbor-count
bound and exists for tests and benchmarking; it does not bypass the maximum
order.
Measure the crossover for a descriptor configuration with:

```sh
build/benchmark-g5-scaling 500
```

The LJ implementation deliberately initializes every atom's output and
provides derivatives. The AccelNet LJ reference currently does neither
reliably; reference tests therefore zero its value buffer before each atom
and compare the intended mathematical result.

## Performance comparison

```sh
build/compare_performance /path/to/TiO2-xsf/structure7513.xsf 20
build/compare_performance /path/to/TiO2-xsf/structure7513.xsf 20 1
build/compare_performance /path/to/TiO2-xsf/structure7513.xsf 20 10
```

The command measures both the value-only kernel with identical neighbor data
and the complete neighbor-list-plus-descriptor path. It exits unsuccessfully
if either new path is more than 5% slower than the untouched AccelNet
reference under the same compiler and optimization flags.

`compare_extension_performance` separately gates the common-model overhead,
LJ kernel, and Behler kernel. This prevents later descriptor additions from
silently slowing the existing Chebyshev fast path.
# AccelNetDescriptors
