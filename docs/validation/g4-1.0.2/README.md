# G4 validation for AccelNet 1.0.2

Version 1.0.2 automatically groups zero-shift G4 descriptors by neighbor
species pair and cutoff radius. Multiple cutoffs and shifted descriptors
can coexist in one model. Nonzero shifts retain the general kernel.
Both value-only and derivative entry points use the same G4 value
arithmetic; the exact-equality failure in the diagnostic prototype is fixed.
No input flag, model conversion or retraining is required.

The release also includes active parameter-group evaluation, fused cutoff
values/derivatives, and reuse of derivative factors. No fast-math flags or
approximations to the descriptor definition are introduced.

## Correctness

- All 23 self-contained CTests passed in Release/static and Debug/shared
  builds with GNU Fortran 11.4.0. Debug uses `-fcheck=all -fbacktrace`.
- `descriptor_cutoff_pair` and `descriptor_g4_groups` are labeled `g4` and
  included in GitHub Actions. The hosted matrix has not yet run for this release.
- The G4 test checks a scalar definition and numerical coordinate derivatives
  for all ten cutoff types, mixed/all-zero/all-nonzero shifts, multiple
  fast-path cutoffs in one species pair, negative shifts, interleaved outputs,
  fractional/integer angular powers, cutoff boundaries and near-collinear
  geometry. Value-only and derivative G4 values must match exactly. Energy-only
  calls are also checked for unused singular angular derivative evaluation.
- All eight PIMD `accelnet` CTest groups passed, including MPI/XMPI, numerical
  forces/cell derivatives, and real aenet/n2p2 model comparisons.
- On the displaced 36-atom ice model, all 108 force components and all nine
  cell derivatives passed central differences at `fdiff=1e-4` bohr. Maximum
  discrepancies in printed derivatives were 3.3e-9 and 2.0e-8, respectively.
  PIMD prints cell derivatives to eight decimal places.
- On initial and displaced ice structures, maximum differences against n2p2
  were 4.45e-15 Ha in energy, 3.79e-15 Ha/bohr in force and 1.20e-13 Ha in virial.

Optional standalone external descriptor-wrapper tests were disabled; the
PIMD comparisons against both actual aenet and n2p2 backends were enabled.

## PIMD ice timing

36 atoms, four beads, one MPI rank and one thread. Twenty warmup steps precede
200 measured steps. Each variant was run three times sequentially in rotating
order, without concurrent builds or validation runs. Startup and model loading
are excluded. These are complete integration wall times, measured from PIMD's
timestamps at 0.01 s resolution.

| Variant | Three times (s) | Median (s) |
|---|---|---:|
| AccelNet 1.0.1 | 12.02, 11.99, 11.99 | 11.99 |
| Previous local G4 optimizations | 11.95, 9.97, 9.95 | 9.97 |
| **AccelNet 1.0.2 grouped fast path** | **9.06, 9.03, 9.01** | **9.03** |
| n2p2 | 7.87, 7.84, 7.82 | 7.84 |

The grouped path takes 9.4% less time than the previous local optimization,
and 24.7% less than 1.0.1 (1.33x throughput). n2p2 remains faster on this model.
The first run of the previous kernel was slower; all measurements are retained
and medians are reported. This table measures the final implementation with
independently grouped cutoffs and both value/derivative entry points.

All printed thermodynamic rows agreed. The maximum final restart difference
between the grouped path and previous local kernel was 8.76e-11, within the
absolute-plus-relative tolerance of 2e-8.

Hardware: Intel Xeon Gold 6526Y. Open MPI 5.0.7; AccelNet uses GNU Fortran 11.4.0
and `-O3`. The reference is PIMD-patched n2p2 2.2.0 with grouping/cache enabled
and `-O3 -march=native`. The PIMD integration working tree is `AccelNet_v1`,
based on `ea2e7541`, with its trained-model examples and backend fixes.

## Reproduction and data

[results.json](results.json) contains individual timings, executable/source
hashes, numerical checks and aenet/n2p2 comparison errors. The PIMD driver is
`examples/10_validation/AccelNet/ice_profile/benchmark_optimization.py` in the
integration working tree. Supply saved executables as `--variant NAME BACKEND
EXECUTABLE`, and use the same example/model for each variant.

For the self-contained G4 tests in an existing configured build:

```sh
cmake --build build --target test_cutoff_pair test_g4_groups test_behler
ctest --test-dir build -L g4 --output-on-failure
ctest --test-dir build -R '^descriptor_behler$' --output-on-failure
```

For PIMD, build with the local AccelNet checkout and both reference backends:
set `ACCELNET=ON`, `AENET=ON`, `N2P2=ON`,
`FETCHCONTENT_SOURCE_DIR_ACCELNET=/path/to/AccelNet`, and
`ACCELNET_DESCRIPTORS_SOURCE_DIR=/path/to/AccelNet/AccelNetDescriptors` in a fresh
build directory, along with the reference dependency paths for the installation.
Then run `ctest --test-dir build -L accelnet --output-on-failure`.
