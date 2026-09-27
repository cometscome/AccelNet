# Default common CPU batch dispatch — 2026-09-27

Methods document **revision 1.6**, AccelNet **1.0.1**. Source baseline: `gpu`
checkpoint **b5e2fcd** (revision 1.5). `source-sha256.json` identifies the measured
sources; `serial-audit.json` records executable hashes and compilation flags.

## Change and scope

`evaluate_batch` now tries the common serial implementation for every model
accepted by the target-model packer. Previously only single-component Chebyshev
models selected it. This adopts LJ and Behler G1–G5, including supported combined
LJ/Behler models. The numerical kernels are unchanged from revision 1.5.

The serial modules are generated from the GPU source with OpenMP directives
removed. CPU execution does not start OpenMP regions or require a GPU compiler.
Model metadata is repacked each call to honor public model edits; workspace
buffers persist. These timings include that packing cost, unlike the preceding
prepared-`target_model` comparisons. The full G4 Jacobian remains a memory cost
proportional to descriptor count times batch edge count.

Unsupported packing retains `evaluate_batch_reference`: mixed/multiple
Chebyshev components within an element, forced G5 moment modes, and other
unsupported configurations keep the established CPU behavior. G5 auto uses
common direct pairs; the reference's auto policy can use moments. Generic
benchmark mode 0 now measures this policy difference; mode 1 still compares
direct algorithms. No algorithm change is made to GPU execution.

This is the **CSR batch API** default. Object/atomic APIs, CLI and ordinary
LAMMPS `pair_style accelnet` retain their established evaluators. The independent
reference remains available for numerical and performance comparisons.

## Fixed-neighbor measurements

- Xeon Gold 6526Y; pinned CPU 6; GNU Fortran 11.4.0 `-O3` and NVHPC 25.3
  `-fast -O3`. Both builds have OpenMP compilation disabled and no direct
  OpenMP runtime references. This is not merely `OMP_NUM_THREADS=1`.
- Eight model families; 8, 512 and 4096 atoms, spacing 1.7. G4-series and
  G5-series use 48 distinct inputs; G4-distinct uses four. G5/Behler use auto
  mode, so the independent CPU may select moments. Chebyshev order 8 is
  measured in both direct and moment modes.
- Additional G5-series auto cases at spacing 1.2, 8/512 atoms, probe denser
  neighborhoods. These are examples, not a density-independent guarantee.
- Two rounds reverse the case order. Each invocation alternates five
  reference/default samples of at least 0.10 seconds. All 116 invocations check
  energies, all force components, and the virial. Initial allocations are
  warmed up; neighbor lists stay outside timing. No builds/tests/profile runs
  overlap this measurement campaign.
- `results.json` retains every sample and command; `summary.json` takes the
  median of the two round medians. Ratios are medians of paired default/reference
  times, not quotients of independently rounded table entries.
- Maximum initial absolute energy/force/virial difference: **4.61e-15**.

The table gives the range of **default/reference CPU time ratios** across the
measured sizes (and the extra density for G5-series). Below 1 means faster.

| Family | GNU ratio range | NVHPC ratio range |
|---|---:|---:|
| LJ | 0.952–1.053 | 0.959–1.193 |
| G4, four distinct inputs | 0.498–0.579 | 0.598–0.627 |
| G4, 48 distinct inputs | 0.628–0.662 | 0.958–1.029 |
| G5 | 0.373–0.437 | 0.314–0.380 |
| G5, 48 distinct inputs | 0.174–0.320 | 0.127–0.237 |
| Behler G1–G5 | 0.355–0.428 | 0.355–0.425 |
| LJ + Behler | 0.366–0.436 | 0.361–0.428 |

The largest slowdown among newly switched fixtures is **19.3%**, NVHPC LJ with
8 atoms: approximately **22.64 → 26.99 microseconds**, an increase of 4.35 us.
At 512/4096 atoms, NVHPC LJ is slightly faster and GNU LJ is about 5% slower.
NVHPC G4-series is approximately unchanged; it is not assigned the speedup
previously measured with a prepared target model. Small percentage differences
are descriptive measurements, without a statistical-significance claim.

## Chebyshev control

Chebyshev was already routed through the common code. Order-8 direct execution
at 512/4096 atoms is slower than the independent old CPU evaluator in this
campaign. It is retained in the measurements rather than hidden among the new
families. `chebyshev-control/` separately alternates the immutable pre-change
binary and the current binary, both using `cpu-shared`, to distinguish this
pre-existing difference from the present dispatch change. The after/before
common-path time ratios are 1.014/1.002 for GNU at 512/4096 atoms and
0.996/1.000 for NVHPC, consistent with no material change in this control.
The common/reference direct ratios remain approximately 1.47–1.68 in these
high-order cases; this adoption does not resolve that earlier limitation.

## Numerical checks and end-to-end regression gate

- GNU serial: 10 selected batch/host tests passed.
- NVHPC serial: the same 10 tests passed.
- GNU runtime/bounds-checked build: 9 selected tests passed.
- H100: target equivalence and the 137-case generic descriptor suite passed,
  including the newly added default-CPU-versus-independent-reference checks.
- Ordinary non-target build: all 42 correctness tests passed after rebuilding.
- The expanded batch suite covers generic models in G5 modes 0–3, shared/fallback
  transitions, model reloads, finite differences, partitioned and reordered rows,
  additive force/virial outputs, empty batches and workspace reuse.

The new opt-in `predictor_common_cpu_performance` CTest passed with **both
compilers**. Each run measures 12 conditions: Chebyshev, LJ, real n2p2 G4 and
G5 models, at 8/64/512 atoms. It uses seven alternating samples of at least
0.2 seconds, includes neighbor construction in both the independent structure
API and default batch API, and checks every energy/force/virial component.
Default/reference ratios range **0.895–1.012 for GNU** and **0.863–1.023 for
NVHPC**, below the configured 1.10 limit. Maximum absolute differences are
6.20e-12 and 5.85e-12, respectively, within the existing absolute/relative
2e-10 tolerances. Raw values and timings are in `gate-gnu.json` and
`gate-nvhpc.json`. This larger end-to-end timing scope is distinct from the
fixed-neighbor table above.

The initial GPU test attempt inside the restricted environment could not see a
GPU; that diagnostic remains in `tests-gpu.log`. The same two device suites were
rerun with GPU access and passed (`tests-gpu-device.log`). CPU host/batch suites
were also successful in that GPU-enabled build. No unresolved GPU test failure
is being treated as a pass.

## Reproduction

`bench.py` records the exact fixed-neighbor campaign; `chebyshev_control.py`
records the before/after control. `finish_checks.py` records the ordinary build
and end-to-end performance gate. Local paths identify the tested builds and
can be replaced with equivalent ones. The general driver can measure the
public API, for example:

```sh
python3 AccelNetPredictor/benchmark/compare_chebyshev_variants.py \
  --variant default cpu-shared /path/to/accelnet-target-benchmark \
  --family g5-series --sizes 8 512 4096 --orders 8 --modes 0 \
  --rounds 2 --seconds 0.10 --cpu 6 --output /tmp/common-cpu-check
```

The benchmark's other column is always the independent retained CPU reference.
`host` selects a prepared target model instead and is a different timing scope.
