# Common LJ/Behler kernel optimization — 2026-09-26

This change improves the same numerical source compiled for an OpenMP-disabled
CPU and for Fortran OpenMP target on GPUs. It does not switch the default LJ/Behler
CPU batch to this implementation: the established CPU path remains faster in some
conditions and remains the independent reference.

## Implementation

1. Group G4/G5 features by exactly equal `(Rc, eta, Rs, cutoff type, cutoff alpha)`.
   Species pairs, lambda, zeta, and whether the descriptor is G4 or G5 do not affect
   the center–neighbor radial function. Pack the group number and representative
   feature once during model initialization.
2. Store `q(r)` and `q'(r)` per edge and radial group in a persistent workspace
   array. Fill it during the existing geometry stage, without another GPU launch.
   Reuse it in both descriptor values and force contraction. The buffer is device
   resident; steady-state evaluations do not allocate it or transfer it to the host.
3. Use value-only Behler radial/pair helpers in the value stage. G4's
   neighbor–neighbor radial factor still depends on the pair and is computed there.
4. Evaluate adjacent LJ6/LJ12 outputs together, sharing cutoff and inverse powers;
   contract radial force coefficients before multiplying the unit direction.

The cache uses two FP64 scalars per edge per radial group, rather than per feature
or a descriptor Jacobian. The current allocator reserves at least one group slot,
including for models without angular Behler terms. The mapping lifecycle, copying,
growth, and same-shape reload checks cover this additional array.

The new regression cases make G4 and G5 share one radial group, split one parameter
without changing model dimensions, and merge it again. They check allocation growth,
reuse, numerical equivalence, and coordinate/strain finite differences.

The formulas are documented in [speedupmethods.md](../../../speedupmethods.md).

## Measurement protocol

- Xeon Gold 6526Y, pinned to CPU 6; GNU Fortran 11.4.0 `-O3` and NVHPC 25.3
  `-fast -O3`, **no OpenMP compiler flags** in CPU comparisons.
- NVHPC 25.3 OpenMP target, H100 NVL and RTX PRO 6000 Blackwell Max-Q, FP64.
- Two species, input–8–4–1 NN, spacing 1.7, 512/4096 atoms. LJ has four outputs;
  G4/G5 each have four; mixed Behler has fourteen G1–G5 outputs.
- Both candidates use direct G5 evaluation. The independent CPU reference is
  unchanged. Each executable alternates five reference/candidate samples and
  verifies all energy/force/virial components at absolute/relative `2e-10`.
- Two rounds, reversing before/after binary order. No other benchmark from this
  task runs concurrently. CSR construction and first allocation/model packing are
  excluded; normal validation, copies, transfers, synchronization, and all outputs
  are included. These are synchronous batch timings, not full MD timings.
- CPU changes are compared using same-run reference ratios. GPU times are the
  median of the two round medians. Small changes are not treated as significant.

Raw samples and phase timings are recorded in `timings.json` and individual logs.
`summary.json` records aggregated results; `compare.py` records the exact local
commands. The tables below summarize the completed before/after measurements.

## OpenMP-disabled CPU results

After/before uses the ratio of same-run reference-normalized timings; values below
one indicate improvement. Final/reference compares with the independent established
CPU implementation, not with the previous common code.

| Compiler | Atoms | Family | After/before | Final/reference |
|---|---:|---|---:|---:|
| gnu | 512 | lj | 0.609 | 1.041 |
| gnu | 512 | g4 | 0.471 | 2.524 |
| gnu | 512 | g5 | 0.362 | 1.299 |
| gnu | 512 | behler | 0.425 | 1.217 |
| gnu | 4096 | lj | 0.600 | 1.041 |
| gnu | 4096 | g4 | 0.466 | 2.559 |
| gnu | 4096 | g5 | 0.359 | 1.292 |
| gnu | 4096 | behler | 0.428 | 1.215 |
| nvhpc | 512 | lj | 0.538 | 0.731 |
| nvhpc | 512 | g4 | 0.593 | 2.404 |
| nvhpc | 512 | g5 | 0.550 | 1.242 |
| nvhpc | 512 | behler | 0.579 | 1.139 |
| nvhpc | 4096 | lj | 0.554 | 0.783 |
| nvhpc | 4096 | g4 | 0.592 | 2.425 |
| nvhpc | 4096 | g5 | 0.537 | 1.224 |
| nvhpc | 4096 | behler | 0.564 | 1.115 |

The common CPU kernel improves in every measured case. G4/G5 still trail the old
CPU implementation. G4's pair-dependent radial factor and the three visits per
unordered pair remain; caching center–neighbor radial factors does not remove
those differences. No hardware cache-miss or bandwidth counters were collected.
GNU LJ is close to the old CPU but not consistently faster, so the default CPU
selection remains unchanged for these families.

## Synchronous GPU batch results

| GPU | Atoms | Family | Before [ms] | After [ms] | Speedup |
|---|---:|---|---:|---:|---:|
| h100 | 512 | lj | 0.3719 | 0.3531 | 1.05× |
| h100 | 512 | g4 | 2.7287 | 1.5398 | 1.77× |
| h100 | 512 | g5 | 2.2120 | 1.1571 | 1.91× |
| h100 | 512 | behler | 3.2230 | 1.7959 | 1.79× |
| h100 | 4096 | lj | 1.1974 | 1.1965 | 1.00× |
| h100 | 4096 | g4 | 5.0921 | 3.3483 | 1.52× |
| h100 | 4096 | g5 | 4.0135 | 3.1701 | 1.27× |
| h100 | 4096 | behler | 9.2494 | 4.8719 | 1.90× |
| blackwell | 512 | lj | 0.4779 | 0.4045 | 1.18× |
| blackwell | 512 | g4 | 8.9111 | 4.2246 | 2.11× |
| blackwell | 512 | g5 | 7.3671 | 2.6692 | 2.76× |
| blackwell | 512 | behler | 10.5241 | 4.9849 | 2.11× |
| blackwell | 4096 | lj | 1.3663 | 1.2345 | 1.11× |
| blackwell | 4096 | g4 | 21.6252 | 11.4307 | 1.89× |
| blackwell | 4096 | g5 | 17.4718 | 7.1130 | 2.46× |
| blackwell | 4096 | behler | 34.2400 | 16.5011 | 2.08× |

H100 LJ at 4096 atoms is effectively unchanged. Its descriptor phase decreased
from approximately 0.1245 to 0.1102 ms and force/scatter from 0.1262 to 0.1202 ms,
but upload alone is approximately 0.66–0.68 ms; small transfer-time variation can
hide this amount of arithmetic improvement in the complete call. These phase
measurements include the synchronous stage overhead, not just CUDA event times.
G4/G5 improvements are much larger and occur on both GPUs and both CPU compilers.
These synthetic descriptor benchmarks do not imply the same full-MD speedups for
arbitrary trained models.

## Numerical checks and reproducibility

All 128 before/after benchmark executions completed their full E/F/virial checks.
The largest initial-case absolute difference in the optimized benchmark runs was
`1.124e-15`; each timed sample was also checked component by component.
The acceptance tolerances were not relaxed.

The optimized shared kernels retain the old independent CPU implementation for
validation. The current descriptor suite has 75 cases, including the three new
radial-group split/merge cases. Full-suite validation results are recorded below.

`implementation.diff` and `before-src/` record the numerical changes;
`source-sha256.json` identifies the final source files. `measured-binaries.json`
identifies timing executables, and `final-binaries.json` identifies the final
validation/build executables. Final rebuilds include source formatting/OpenMP line
wrapping and the added regression cases; the arithmetic is the measured algorithm.

For another machine, use the commands in [OpenMP target instructions](../../openmp-target.md)
and the reusable comparison driver (the filename is retained for compatibility):

```sh
python3 AccelNetPredictor/benchmark/compare_chebyshev_variants.py \
  --variant before host /path/to/saved-before \
  --variant after host build-target-serial/bin/accelnet-target-benchmark \
  --family g4 --modes 1 --sizes 512 4096 --rounds 2 --seconds 0.08 \
  --output comparison-g4
```

Repeat for `lj`, `g5`, and `behler`; use GPU executables and `gpu` as the backend
with the intended GPU UUID for device comparisons. Original CPU diagnostic variants
are historical; `diagnose_target_cpu.py` now rejects the changed kernel ABI and
points to their archived results rather than generating incompatible code.

### Completed validation of this revision

- GNU/NVHPC OpenMP-disabled CPU suites: 10 tests each, passed.
- H100: 29 GPU + 2 host tests; Blackwell: 29 GPU tests, passed.
- GNU OpenMP with `-fcheck=all -fbacktrace`: 9 tests, passed.
- H100 Compute Sanitizer 2025.4.1.0 (CUDA 13.1 distribution): all 75 descriptor
  cases, `ERROR SUMMARY: 0 errors`; `--leak-check no`.
- LAMMPS: 63 LJ/Behler/converted-n2p2 comparisons, covering CPU/GPU neighbor
  modes, 1/2 ranks, orthogonal/triclinic/empty-rank cases and short trajectories.
- Existing CPU numerical suite: 42 tests; performance gates: 2 tests, passed.
- Chebyshev direct/moment comparison: one round at 512/4096 atoms on both CPU
  compilers and H100; raw results in `chebyshev-*`. Small timing differences in
  this spot check are not claims of a Chebyshev algorithm improvement.

The new arithmetic tests and final rebuilds preserve all numerical tolerances.
The compiled LAMMPS hash is recorded in `final-binaries.json`. This revision's
binary was validated before the subsequent CPU-style Jacobian experiment.
