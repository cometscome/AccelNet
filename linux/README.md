# Linux benchmark kit: H2O and TiO2

This directory builds one LAMMPS executable containing both the original
`pair_style aenet` and `pair_style accelnet`, verifies energy/force agreement,
and benchmarks fixed-coordinate and dynamic MD workloads for H2O and TiO2.

The default benchmark is deliberately short.  Run the smoke test first, then
select paper-length step counts and system sizes explicitly.

For the three-way comparison of original aenet, forced-direct AccelNet, and
forced-moment AccelNet—including detached background execution—see
[`VARIANT_BENCHMARK.md`](VARIANT_BENCHMARK.md).

## Included benchmark systems

| system | base cell | small | medium | large | model order |
|---|---:|---:|---:|---:|---|
| H2O | 192 atoms | 1,152 | 11,520 | 96,768 | H, O |
| TiO2 | 6 atoms | 1,152 | 8,064 | 98,304 | Ti, O |

- H2O uses the binary aenet networks in `examples/water`.
- TiO2 uses the trained Chebyshev networks from aenet's
  `example-chebyshev/03-predict/set001` and relaxed `structure6349.xsf` from the
  7,815-structure TiO2 reference set.  The converted LAMMPS data and network
  files are checked into `examples/tio2`.
- `small` uses 1,152 atoms for both systems, which makes the chemistry
  comparison especially simple.

The classical-water inputs under `benchmark/` remain available as a separate
reference.  This Linux driver first targets the exact, same-model comparison:
original aenet versus AccelNet.

The dynamic inputs in this kit are throughput workloads, not production
trajectories for comparing physical properties.  In particular, the TiO2
input uses a strongly coupled CSVR thermostat to keep the supplied reference
cell close to 300 K during timing.  For structure, diffusion, or thermodynamic
results, equilibrate each potential independently and document the production
ensemble separately.

## 1. Install Linux dependencies

Ubuntu/Debian example:

```sh
sudo apt update
sudo apt install -y \
  build-essential cmake gfortran python3 \
  openmpi-bin libopenmpi-dev \
  libopenblas-dev liblapack-dev
```

On a managed cluster, load equivalent compiler, Open MPI, CMake, OpenBLAS, and
Python modules.  Do not mix Fortran runtimes from different compiler stacks.
The automated build currently targets GNU Fortran and OpenBLAS.

## 2. Place the source trees together

The scripts use this layout by default:

```text
AccelNetGPU/
├── AccelNet/
├── AccelNet-lammps/
├── aenet-master/
├── aenet-lammps/
└── lammps-4Feb20/
```

For example, transfer clean source trees from the development machine with
build products excluded:

```sh
rsync -a --exclude 'build*' --exclude 'Obj_*' --exclude 'results/' \
  AccelNet AccelNet-lammps aenet-master aenet-lammps lammps-4Feb20 \
  user@linux-host:/path/to/AccelNetGPU/
```

The locations can instead be supplied through `ACCELNET_SOURCE`,
`AENET_SOURCE`, `AENET_LAMMPS_SOURCE`, and `LAMMPS_SOURCE`.

## 3. Build the combined LAMMPS executable

From the Linux machine:

```sh
cd /path/to/AccelNetGPU/AccelNet-lammps
./linux/build_linux.sh
```

The script performs the following operations:

1. Release/static build of AccelNet.
2. Release/static build of original aenet with OpenBLAS.
3. Installation of both libraries into the LAMMPS tree.
4. Installation of `USER-AENET` and `USER-ACCELNET`.
5. Enabling KSPACE, MOLECULE, and RIGID for the classical-water inputs.
6. MPI build of LAMMPS and verification that both pair styles are present.

Useful overrides:

```sh
BENCH_BUILD_JOBS=16 BENCH_RUN_TESTS=1 ./linux/build_linux.sh
```

The build compiles only the libraries needed by LAMMPS plus AccelNet's small,
self-contained unit tests.  Standalone benchmark executables and tests that
require an external original-aenet reference tree are deliberately excluded.
This also keeps older GNU Fortran versions from compiling unrelated benchmark
sources.  Set `BENCH_RUN_TESTS=0` only as a temporary build diagnostic; the
LAMMPS smoke test is still mandatory before benchmarking.
The expected executable is `lammps-4Feb20/src/lmp_mpi`.

## 4. Run correctness smoke tests

```sh
./linux/smoke_test.sh
```

This evaluates 1,152 atoms for each system with both pair styles and compares
the total energy and every Cartesian force component.  The command must end in
two `PASS` lines before timing measurements are accepted.

The smoke test uses bare model filenames from each model directory.  This is
intentional: the historical aenet-LAMMPS interface has a short filename buffer,
whereas AccelNet-LAMMPS accepts long paths.

## 5. Run a short end-to-end benchmark

The defaults are suitable for checking a new machine:

```sh
./linux/run_benchmarks.sh
```

Default conditions:

- systems: H2O and TiO2
- size: 1,152 atoms
- modes: fixed coordinates and dynamic NVT
- pair styles: aenet and AccelNet
- MPI ranks: 1 and 2
- repeats: 3
- warm-up: 100 steps
- measured: 500 steps
- OpenMP/BLAS threads: 1
- MPI binding: one rank per core

Runs are interleaved as aenet/AccelNet pairs within each trial to reduce drift
from thermal throttling or changing machine load.

## 6. Paper-length examples

### Main 1,152-atom comparison

```sh
BENCH_SIZES='small' \
BENCH_MPI_RANKS='1 2 4 8' \
BENCH_TRIALS=5 \
BENCH_WARMUP_STEPS=500 \
BENCH_MEASURE_STEPS=5000 \
./linux/run_benchmarks.sh
```

### Size scaling

Large original-aenet runs can take many hours.  Measure each size separately
and choose a step count that gives at least 10--30 seconds per trial before
launching the full matrix.

```sh
BENCH_SIZES='medium' BENCH_MPI_RANKS='1 2 4 8' \
BENCH_TRIALS=5 BENCH_WARMUP_STEPS=100 BENCH_MEASURE_STEPS=500 \
./linux/run_benchmarks.sh

BENCH_SIZES='large' BENCH_MPI_RANKS='4 8' \
BENCH_TRIALS=5 BENCH_WARMUP_STEPS=20 BENCH_MEASURE_STEPS=100 \
./linux/run_benchmarks.sh
```

Do a one-trial pilot before increasing these values.  Do not combine results
from different node types, compiler stacks, CPU governors, or binding policies.

### Fixed-coordinate only

```sh
BENCH_MODES='fixed' ./linux/run_benchmarks.sh
```

### Dynamic MD only

```sh
BENCH_MODES='dynamic' ./linux/run_benchmarks.sh
```

### SLURM allocation

Inside an allocation, use `srun` instead of `mpirun`:

```sh
BENCH_LAUNCHER=srun BENCH_MPI_RANKS='1 2 4 8' ./linux/run_benchmarks.sh
```

The script passes `--cpu-bind=cores` to `srun`.  If the cluster supplies its
own binding policy, set `BENCH_BIND_CORES=0` and record that policy separately.

## 7. Configuration variables

| variable | default | meaning |
|---|---|---|
| `LAMMPS_BINARY` | sibling `lammps-4Feb20/src/lmp_mpi` | executable |
| `BENCH_SYSTEMS` | `water tio2` | selected systems |
| `BENCH_SIZES` | `small` | `small`, `medium`, `large` |
| `BENCH_MODES` | `fixed dynamic` | selected workloads |
| `BENCH_PAIR_STYLES` | `aenet accelnet` | implementations |
| `BENCH_MPI_RANKS` | `1 2` | rank counts |
| `BENCH_TRIALS` | `3` | independent processes per condition |
| `BENCH_WARMUP_STEPS` | `100` | untimed warm-up steps |
| `BENCH_MEASURE_STEPS` | `500` | timed steps |
| `BENCH_LAUNCHER` | `mpirun` | `mpirun` or `srun` |
| `BENCH_BIND_CORES` | `1` | enable explicit core binding |
| `BENCH_OMP_THREADS` | `1` | OpenMP threads per rank |
| `BENCH_BLAS_THREADS` | `1` | OpenBLAS/MKL/BLIS threads per rank |
| `BENCH_RESULT_DIR` | timestamped directory | output location |

## 8. Outputs

Each run produces a timestamped directory under `linux/results`:

```text
run-YYYYMMDD-HHMMSS/
├── system-info.txt
├── summary.tsv
└── logs/
    └── SYSTEM-SIZE-ATOMS-MODE-PAIR-rRANKS-tTRIAL.log
```

`system-info.txt` records hardware, compilers, MPI, selected environment,
available Git revisions, and SHA-256 hashes of the executable and model files.
Every raw LAMMPS log contains the standard timing breakdown (`Pair`, `Neigh`,
`Comm`, `Modify`, `Kspace`, and `Other`).

`summary.tsv` reports median and IQR loop times, steps/s, atom-step/s, and
median timing components.  Recreate it at any time with:

```sh
python3 linux/summarize_benchmarks.py /path/to/run/logs > summary.tsv
```

Keep `system-info.txt`, `summary.tsv`, and all raw logs together when moving
results back from the Linux machine.

## 9. Acceptance criteria

Before using a timing set in a paper, verify:

- both smoke tests report energy/force differences below `1e-8`;
- every log reports `Dangerous builds = 0`;
- each condition has the requested number of completed trials;
- measured sections are long enough to suppress launch and clock noise;
- no unrelated jobs share the selected cores;
- `summary.tsv` was generated without parser errors;
- aenet and AccelNet used the same data, model files, rank count, and binding.
