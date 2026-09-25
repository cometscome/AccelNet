# Additional regression tests

The expanded suite passes **34/34 tests** with the local Ti/O corpus.
The GitHub Actions configure/build/test commands pass **21/21 tests** in
Release/static and **21/21** in Debug/shared with `-fcheck=all -fbacktrace`,
both reproduced locally with GNU Fortran 11.4.0. The configured GNU 13
GitHub runners remain unexecuted until these local changes are pushed.

## Added coverage

| Test | Failure it detects |
|---|---|
| Multispecies n2p2 strain sweep | Wrong element mapping, network topology, or descriptor normalization |
| H/O G2/G4/G5 strain sweep | Missing angular contributions or inconsistent derivative paths |
| Atomic vs structure API with mixed derivative paths | Inconsistent per-species contraction/Jacobian results |
| G5 DIRECT/MOMENT/MOMENT_FORCE strain derivatives | Mode-specific force or virial errors |
| Atom/species permutation | Index/order-dependent energies, forces, or virials |
| Rotation covariance | Wrong tensor components; verifies `W' = R W R^T` and `F' = R F` |
| Independent periodic translations of atoms | Lost images when individual atoms lie multiple cells outside the primary cell |
| Neighbor image metadata | Inconsistency between image shifts, Cartesian image coordinates, and displacements |
| Invalid C API types, indices, counts, and states | Wrong status or modification of the caller's accumulators on failure |
| Finalize/reinitialize/reload | Stale global model state |
| Smooth cutoff: below, exactly at, and above | Spurious forces/virials or inconsistent energy derivatives at the cutoff |

The new multispecies strain sweeps use 13 steps and all nine components, with
the same per-component tolerances as the original validation. The driver now
uses three molecular atoms so angular terms also contribute without periodicity.
The nonlinear fixture is the existing `n2p2-per-element` model; the mixed
G4/G5 fixture is the new hand-authored `n2p2-virial-angular` model. Neither
requires a private corpus or a reference n2p2 executable.

## Bug detected and fixed

The previous triclinic brute-force neighbor search enumerated a bounded set
of lattice translations about each input coordinate. Moving individual atoms
by several lattice vectors placed valid neighbors outside that search range.
Such a reimaging must leave the physical structure unchanged, but the tests
failed before the fix. The largest virial discrepancy for the angular fixture
was about **0.459 model energy units**; its energy discrepancy was **0.193**.

The search now recenters the relative fractional position of each atom pair
before enumerating images, and includes that offset in the returned image
metadata. Both counting and population passes use the same translation.
Pairs already inside the primary cell retain the previous enumeration order.
The former failure and the complete existing test suite pass after the fix;
the invariance comparison tolerance is `2e-9` in model units.

## Reproduction

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --target test_virial_invariants test_virial_convergence test_virial_c --parallel 2
ctest --test-dir build -R 'virial_(invariants|convergence_n2p2-)|predictor_c_virial' --output-on-failure
```

Both added test executables/targets are included in the GitHub Actions
workflow. The C API additions run under the existing `predictor_c_virial`
test. These checks establish the tested numerical invariants; they do not
replace the pending PIMD-level integration and `TESTVIRIAL` checks.
