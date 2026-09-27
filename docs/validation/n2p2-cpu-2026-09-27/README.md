# Direct CPU inference and LAMMPS comparison — revision 1.9

Baseline: `gpu` commit **18d1ac3**, AccelNet **1.0.1**. The implementation in
this revision groups extended angular descriptors, contracts wide angular
coefficients before Cartesian arithmetic, and connects the LAMMPS CPU adapter
to the same serial numerical source used by GPU inference. Mathematical details
are in [speedupmethods.md, Section 20](../../../speedupmethods.md#20-grouped-extended-angular-evaluation-and-lammps-cpu-batches-revision-19).

## Measurement conditions

* Intel Xeon Gold **6526Y**, one process pinned to CPU **6**. The governor is
  `schedutil`; frequency is not locked. Timings use reversed process order and
  repeated samples after a 20-second pinned-core warmup. Benchmarks and builds
  do not run concurrently. Final timing runs use the host execution context.
  Short exploratory runs showed large frequency transitions; their first-case
  aggregate medians are not used as the final comparison.
* GNU Fortran/C++ **11.4.0**, `-O3`, default ISA, FP64. Both CPU implementations
  are compiled **without OpenMP**; undefined symbols are checked for OpenMP
  runtimes. Thread-count environment variables are also set to 1. The serial Predictor and Descriptors
  libraries use `-fPIC` so NVHPC objects can link into C PIE
  executables; the original baseline did not require that flag.
* n2p2 **v2.3.0**, unmodified grouped libnnp, Eigen **3.4.0** with
  `EIGEN_DONT_PARALLELIZE`. Source identity is retained in the
  [extension archive](../n2p2-extensions-2026-09-27/audit.json). The serial
  LAMMPS interface uses `N2P2_NO_MPI`; the inference code is unchanged.
* LAMMPS **29 Aug 2024 Update 4**. CPU: GNU, `BUILD_MPI=OFF`,
  `BUILD_OMP=OFF`, ACCELNET and ML-HDNNP packages. GPU: NVIDIA HPC SDK **25.3**,
  OpenMPI, C++ `-O2 -DNDEBUG`, CUDA GPU package in double precision, Fortran OpenMP target
  `-mp=gpu -gpu=cc90,cc120`. GPU measurements use one MPI rank, one pinned CPU
  core, and an **H100 NVL**. CPU and GPU therefore have different compilers;
  the n2p2/AccelNet CPU comparison uses the same GNU toolchain.
* Published models and example geometries are used without retraining:
  Ethylbenzene (H/C, 288 atoms), Anisole (H/C/O, 256), DMABN (H/C/N, 21),
  and water (H/O, 1080). Model hashes accompany the results; model weights are
  not redistributed in this archive. The first three contain types 20/22;
  water contains types 2/3. All angular work in these models is **direct**.

## Actual LAMMPS step time

Both CPU pair styles use the same input geometry, units `electron`, neighbor
skin 0.6, and `neigh_modify every 5 delay 0 check no`. A seeded 100 K NVE
trajectory uses a 0.1 fs timestep. Two warmup steps precede three timed runs;
backend order reverses in the second process round. Each result is the median
of six run averages. Runs use 3 steps for dense compact models, 100 for the
isolated molecule, and 10 for water. Loading, setup and dumps are outside Loop
time; integration and in-loop neighbor rebuilding remain inside it.

Initial/final positions, every force component, energy and all six virial
pressure components are checked against n2p2. Snapshot tolerance is
`2e-8 + 2e-9*abs(reference)`; virial-pressure tolerance is
`0.05 Pa + 2e-8*abs(reference)`. Actual errors are retained in the JSON results.
GPU timings include LAMMPS GPU neighbor construction and host/device transfer.

The standalone n2p2 oracle uses `Mode::calculateForces`, which searches
neighboring centers' neighbor lists. Its LAMMPS interface instead scatters
stored derivatives directly through `InterfaceLammps::getForces`. This is why
standalone parity is insufficient evidence of LAMMPS parity.

Final LAMMPS results (**milliseconds/MD step**, median of six runs):

| Model | Atoms | n2p2 CPU | AccelNet CPU direct | n2p2 / AccelNet CPU | Previous H100 | Current H100 |
|---|---:|---:|---:|---:|---:|---:|
| Ethylbenzene_SCAN | 288 | 1267.387 | 1164.287 | 1.089× | 112.485 | 87.094 |
| Anisole_SCAN | 256 | 1170.568 | 1090.090 | 1.074× | 124.025 | 82.977 |
| DMABN_SCAN | 21 | 3.487 | 2.102 | 1.659× | 3.971 | 4.145 |
| H2O_RPBE-D3 | 1080 | 227.590 | 233.477 | 0.975× | 18.463 | 14.954 |

All four CPU models passed the **1.10× n2p2 time** gate. Three are faster;
water takes **2.6% longer**. The direct path therefore reaches near-parity on
these published models in the actual LAMMPS integration loop. CPU and GPU use
the common numerical source; no moment replacement is used in this table.

Relative to the previous H100 binary, the three periodic systems improve by
approximately **1.29×, 1.49×, and 1.23×**. The 21-atom DMABN GPU run is **4.4%
slower** than the previous GPU binary and about twice the new CPU time; its
small workload does not amortize GPU overhead. GPU speedup is not universal.

Across the new CPU/GPU runs, maximum differences from n2p2 were below
**2.92e-11 Ha** in total energy, **2.78e-14 Ha/Bohr** in force components,
**7.11e-15 Bohr** in final positions, and **0.00168 Pa** in virial pressure.
See `lammps-direct.json` for every sample and error.

## Standalone descriptor-family comparison

`compare_n2p2_cpu.py --direct` also forces G5/type 9 to use the direct method.
**Public CPU** is `evaluate_batch`, including per-call model packing;
**prepared CPU** retains a serial packed model. Both execute the same shared
Fortran numerical source with OpenMP compiled out. Fixed-neighbor timing
isolates inference. The optional full scope rebuilds the standalone neighbor
list on every evaluation; LAMMPS uses its own neighbor builder.

Each process reports five warmed samples lasting at least 0.15 s, or one
complete evaluation if longer. Two reversed process rounds give ten samples.
Every process's energy and all forces must match n2p2 within
`2e-10 + 2e-9*abs(reference)`. AccelNet additionally computes virial internally;
n2p2's final output-unit conversion is outside its timing.

Synthetic tests use 512 atoms, two elements, and six functions per central
species for each SF type (three powers for type 9). They isolate descriptor
families and do not represent equally accurate fitted potentials. The baseline
type-9 result uses automatic moments, so it is not a before/after direct timing.
The separate G5 moment performance gate protects the existing moment advantage.

Both benchmark scripts accept `--max-n2p2-slowdown 1.10` to reject a measured
CPU/n2p2 time ratio above 1.10. This is a configurable performance check, not a
universal accuracy or speed guarantee for every model and machine.

Final fixed-neighbor **direct** measurements (milliseconds/evaluation):

| Model / n2p2 type | n2p2 | Public CPU | Prepared CPU | n2p2 / public CPU |
|---|---:|---:|---:|---:|
| Ethylbenzene_SCAN | 1261.466 | 1255.546 | 1250.898 | 1.005× |
| Anisole_SCAN | 1159.658 | 1213.565 | 1208.464 | 0.956× |
| DMABN_SCAN | 3.969 | 2.335 | 1.930 | 1.700× |
| H2O_RPBE-D3 | 243.156 | 236.742 | 237.051 | 1.027× |
| type2 | 5.938 | 4.502 | 4.322 | 1.319× |
| type3 | 44.791 | 41.387 | 41.018 | 1.082× |
| type9 | 28.831 | 26.866 | 26.423 | 1.073× |
| type12 | 8.887 | 6.273 | 6.532 | 1.417× |
| type13 | 83.022 | 76.986 | 77.198 | 1.078× |
| type20 | 3.919 | 3.125 | 2.995 | 1.254× |
| type21 | 27.750 | 30.236 | 30.311 | 0.918× |
| type22 | 42.257 | 46.165 | 46.268 | 0.915× |
| type23 | 19.100 | 3.119 | 2.983 | 6.124× |
| type24 | 97.026 | 49.727 | 49.634 | 1.951× |
| type25 | 172.677 | 69.231 | 69.262 | 2.494× |

All 15 cases passed the numerical, historical CPU regression, and **1.10×
n2p2 time** gates. Twelve are faster than n2p2 in the public CPU path; Anisole,
type 21, and type 22 take approximately 4.6%, 9.0%, and 9.2% longer.
This establishes near-parity on these fixtures, not a universal speed guarantee.

## Interpretation and rejected experiments

The final implementation uses a single pair pass for descriptor values and
Jacobians, as in the established G4 implementation. It groups matching chemical
pairs and angular windows, shares radial/cutoff factors, and expands Cartesian
coefficients only after their scalar contraction. These changes carry over the
Chebyshev/G4 lessons without changing the descriptor definitions or NN inputs.

Small synthetic models expose overhead that dense fitted models can hide:
G2 caches must skip unused chemical species, and compact angular windows must
not introduce a small function call for every pair/member. The window
polynomial is a shared include used by the radial and angular evaluators;
CPU and GPU compile the same arithmetic. Profiling the intermediate type-22
kernel found approximately 1.6 million compact-window calls per evaluation.

Several plausible changes were measured and rejected: an array cache of radial
pair products added indexing/memory traffic; combining the three G4 Gaussian
exponents did not provide a consistent benefit; explicitly fixing more array
extents slowed type 22. The retained explicit Cartesian Jacobian extent is 3.
These observations are specific to the measured builds, not general compiler
rules. Workspace padding has no separately established speedup claim.

## Numerical and regression validation

* GNU CPU: **49 functional CTests passed**, including the additive C batch API,
  cache reloads, native conversion, legacy fixtures and periodic neighbor images.
* CPU performance CTests: **3 passed**. Across the common-code regression
  fixtures, the largest common/previous time ratio was **1.021**. On the G5
  moment fixture, common moment/direct time was **0.11145** (about **8.97×**
  faster); common moment/previous moment time was **0.1313**. These are the
  dedicated fixture results, not the n2p2 direct comparison.
* GNU bounds-checking build: **72 independent n2p2 cases passed** with
  `-fcheck=all -fbacktrace -g`. H100 and Blackwell: **72 cases each passed**.
  These include all compact subtypes, weighted cutoff variants, conversion and
  round trips, mixed descriptors, four elements, sparse/collinear cases, and
  direct/automatic/moment G5 selection.
* Independent finite differences: maximum force error **1.05e-11** and strain
  derivative error **6.44e-12** (rounded upper bounds, atomic units).
* H100 existing batch tests: **78 checks** (maximum E/F/virial error
  **7.98e-16**) and G5 moment tests: **185 checks** (**2.23e-16** upper bound).
* LAMMPS mixed-model regression: **21 comparisons** spanning 1/2 MPI ranks,
  orthogonal/triclinic/empty-rank configurations and CPU/GPU-neighbor modes
  `no`, `yes`, `hybrid`. Maximum total-energy error **2.23e-16**; maximum
  virial-pressure difference **1.39e-16** in that test's metal units.
* H100 Compute Sanitizer memcheck on the four-element/all-type model with
  G5 moments forced: **0 errors**.

Reports and logs are stored alongside this file. Timing checks are separate
from numerical checks and run without concurrent builds or other test jobs.

## Reproduction

Build AccelNet and n2p2 in serial release mode (GNU `-O3`, no OpenMP), and
build the LAMMPS ACCELNET/ML-HDNNP comparison binary with the same libraries.
Use `audit.json` for exact flags and source identities. The baseline n2p2
standalone reference driver and its build are documented in the
[revision-1.8 archive](../n2p2-extensions-2026-09-27/README.md).
The following placeholders refer to the resulting executables and the upstream
`examples/nnp-predict` directory:

```sh
python3 AccelNetPredictor/benchmark/compare_n2p2_cpu.py \
  --reference "$N2P2_REFERENCE" --candidate "$ACCELNET_CANDIDATE" \
  --examples "$N2P2_EXAMPLES" --suite both --scopes fixed --direct \
  --seconds .15 --rounds 2 --affinity 6 --max-n2p2-slowdown 1.10 \
  --output /tmp/accelnet-direct-cpu

CUDA_VISIBLE_DEVICES="$H100_UUID" \
python3 interfaces/lammps/29Aug2024/tests/benchmark_n2p2.py \
  --before "$LAMMPS_N2P2" --after "$LAMMPS_ACCELNET_CPU" \
  --gpu-before "$LAMMPS_PREVIOUS_GPU" --gpu "$LAMMPS_ACCELNET_GPU" \
  --converter "$ACCELNET_CONVERTER" --examples "$N2P2_EXAMPLES" \
  --variants n2p2 after gpu_before gpu --rounds 2 --samples 3 \
  --affinity 6 --max-n2p2-slowdown 1.10 --output /tmp/accelnet-direct-lammps
```

`--before` supplies the n2p2 pair style; it need not contain the old AccelNet
implementation unless the `before` variant is also requested. GPU arguments
can be omitted for `--variants n2p2 after`. `--types 2 21 22 --suite synthetic`
selects a focused standalone check. The CPU script also accepts a prior JSON
report through `--baseline` to guard against old CPU regressions; the historical
G5 automatic-moment timing is explicitly excluded from a direct-mode baseline
comparison. The final baseline gate used `cpu-prior-session.json` from this
machine, whose separate measurement session is a limitation of that comparison.

`performance-commands.json` retains the exact final invocations and exit codes.
The JSON timing reports contain binary/model hashes, raw timing samples,
numerical errors, and per-backend commands. No trained weights are copied here.
