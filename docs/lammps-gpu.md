# LAMMPS GPU package + Fortran OpenMP target

The standard CPU workflow remains the default; see the
[LAMMPS CPU instructions](../interfaces/lammps/README.md#cpu-build-and-run).
This guide describes the optional `pair_style accelnet/gpu` adapter. LAMMPS's
GPU package assigns devices and supplies neighbor lists; the common Fortran
OpenMP target kernels evaluate descriptors, the NN, forces and virials.

## Build requirements and instructions

Use **LAMMPS 22 Jul 2025 Update 6** (`stable_22Jul2025_update6`), **NVHPC `nvc++` and `nvfortran`**, and the
**CUDA toolkit** with `nvcc`, headers, libraries and `bin2c`. The tested stack is
NVHPC 25.3 / SDK CUDA 12.8, with a compatible NVIDIA driver on H100 NVL and RTX
PRO 6000 Blackwell. The 2025 release was checked on H100; Blackwell results
refer to the retained 29 Aug 2024 Update 4 release. CMake 3.20+, Python 3 for
the installer, a C compiler and
MPI for MPI-enabled builds are also needed. `nvcc` alone cannot compile the
Fortran kernels. The adapter enforces NVHPC C++/Fortran, `GPU_API=cuda` and
`GPU_PREC=double`; AMD/HIP/OpenCL interoperability is not implemented.

The [complete GPU build and run recipe](../interfaces/lammps/README.md#gpu-build-and-run)
includes separate CPU/GPU build directories, architecture flags (`cc90`/`sm_90`
for H100, `cc120`/`sm_120` for the tested Blackwell device), installation, MPI
configuration and n2p2 conversion. Use that recipe as the build reference.
Install with `python3 interfaces/lammps/install.py /path/to/lammps`. The same
installer supports 29 Aug 2024 Update 4 and copies the shared adapter sources
from their historical `29Aug2024/` directory. The 4 Feb 2020 adapter is CPU-only.
The 2 Sep 2026 release candidate is not supported by the installer because its
GPU API changed; see the [compatibility report](validation/lammps-current-2026-09-28/README.md).

The installer copies the CPU/GPU pair styles, lib/gpu adapter, canonical C
headers and CMake module; it also registers the package and patches triclinic
sorting. It updates an existing CPU installation as needed. The supplied CMake
module uses NVHPC `-fortranlibs` so the driver orders its own Fortran/offload
runtimes: an earlier manual runtime ordering caused `omp_get_num_devices()`
to return zero despite a working GPU. `BUILD_OMP=OFF` in the LAMMPS recipe
controls host OpenMP; Fortran GPU offload is enabled separately by `-mp=gpu`.

## Input and supported models

```lammps
# Before read_data/create_box:
package gpu 1 neigh yes newton on split 1
units metal
atom_style atomic
read_data TiO2.data
pair_style accelnet/gpu auto Ti.nn.ascii O.nn.ascii
pair_coeff * *
neighbor 0.6 bin
neigh_modify every 10 delay 0 check yes
fix integrate all nve
run 100
```

Use units consistent with the model; `metal` assumes eV and Angstrom here.
Select devices with `CUDA_VISIBLE_DEVICES` and use `OMP_TARGET_OFFLOAD=MANDATORY`
to require offload. LAMMPS passes each MPI rank's assigned device to Fortran.
`OMP_NUM_THREADS=1` controls host threads, not GPU parallelism. With the explicit
pair style and package command above, no `-sf gpu` or extra `-pk gpu` override
is needed. Changing to `neigh no` changes neighbor construction, **not** the
potential's GPU execution. The GPU-enabled binary can still run the ordinary
CPU input with `pair_style accelnet`.

Supported models use embedded `.nn`/`.nn.ascii` metadata: Chebyshev (the LAMMPS
interface uses version 0), LJ, Behler G1--G5 and converted n2p2 types
2/3/9/12/13/20--25. Multiple Chebyshev blocks and mixed families are supported.
LAMMPS atom types must match the embedded global species order. Direct n2p2
model-directory input is CPU-only; for GPU use the Fortran converter:

```sh
# From the AccelNet repository root, after a full build:
build/bin/accelnet-model-converter-fortran n2p2-to-accelnet n2p2-model converted
```

Use the actual output files, for example:

```lammps
pair_style accelnet/gpu auto converted/H.nn.ascii converted/O.nn.ascii g5 moment
pair_coeff * *
```

The leading mode selects Chebyshev `auto`/`direct`/`moment`; the independent
trailing `g5 MODE` selects G5. G5 auto admits exact integer powers 1--10 with
at least 16 angular neighbors. Explicit `g5 moment` admits powers 1--16 and
removes that neighbor threshold. Fractional, near-integer and powers above 16
remain direct. Types 13/21/24 are implemented with exact direct evaluation;
no approximate moment expansion is substituted. Auto estimates operation
counts and eligibility; it does not benchmark the hardware.

Local and ghost atoms, periodic boundaries, orthogonal and restricted triclinic
cells, global energy/forces/virial and per-atom energy are supported. `newton on`
and `split 1` are required. Per-atom stress, pair hybrid, r-RESPA, molecular
topology, neighbor exclusions/include groups and model serialization into
restart files are unsupported. Small systems can favor CPU execution.

## Data movement, communication and lifetime

With `neigh no`, the adapter packs a full CPU neighbor list into CSR and uploads
it in a batch. With `neigh yes/hybrid`, it borrows lib/gpu device coordinates and
neighbor arrays and converts them to CSR/displacements on GPU. It does not
download and re-upload the neighbor arrays; one integer CSR edge count is
returned to the host. Model, NN and descriptor/force buffers remain resident
and grow only when capacity is insufficient.

In this LAMMPS release's standard CUDA configuration, even `neigh yes` uses
host binning inside lib/gpu and GPU candidate-neighbor search. `neigh hybrid`
can share that implementation; this is not a claim that all binning runs on GPU.

Each `compute()` completes GPU evaluation and downloads local+ghost forces
before LAMMPS performs reverse communication. AccelNet does not also enqueue
the same forces in FixGPU's post-force answer queue. The GPU package and OpenMP
backend synchronize before exchanging external device pointers.

The supplied triclinic patch avoids reading uninitialized orthogonal sub-box
bounds during GPU atom sorting, using the normal bounding-box sort instead;
sorting remains enabled. Explicit scalar/rank-1 workspace finalizers avoid an
NVHPC 25.3 elemental-finalizer issue. Each C handle releases its model and
workspace. The GPU's primary CUDA context remains alive until process exit
because OpenMP runtime caches outlive LAMMPS `clear`. Reloading after `clear`
has dedicated validation.

## Correctness and performance checks

The scripts under `interfaces/lammps/29Aug2024/tests/` require Python, NumPy
and an MPI launcher. Run these from the AccelNet root, after selecting GPU
visibility. Replace model/executable paths and the MPI launcher as needed:

```sh
export OMP_NUM_THREADS=1
export OMP_TARGET_OFFLOAD=MANDATORY
python3 interfaces/lammps/29Aug2024/tests/check_gpu.py \
  --lammps /path/to/lmp --golden /path/to/fortran_predict \
  --mpiexec /path/to/mpiexec --output /tmp/gpu-check
python3 interfaces/lammps/29Aug2024/tests/check_gpu_errors.py \
  --lammps /path/to/lmp --input /tmp/gpu-check/orthogonal-auto-no-1rank/in.test \
  --output /tmp/gpu-errors
python3 interfaces/lammps/29Aug2024/tests/benchmark_gpu.py \
  --lammps /path/to/lmp --input /tmp/gpu-check/orthogonal-auto-cpu-1rank/in.test \
  --mpiexec /path/to/mpiexec --output /tmp/gpu-benchmark
python3 interfaces/lammps/29Aug2024/tests/benchmark_gpu_modes.py \
  --lammps /path/to/lmp --input /tmp/gpu-check/orthogonal-auto-cpu-1rank/in.test \
  --cpu 6 --output /tmp/gpu-modes
```

`check_gpu.py` compares against CPU/one MPI rank: total/per-atom energy, every
force, six pressure/virial components and a short NVE trajectory. Defaults
cover 1/2/4 ranks, all neighbor modes, orthogonal/triclinic/empty-rank cases,
and neighbor rebuilding/sorting each step. `--gpus 2 --ranks 2 4` exercises
multiple GPUs; `--cases migration --rebuild-every 5 --modes direct moment`
checks rank migration and neighbor reuse. `check_gpu_lifecycle.py` checks
capacity growth and reload after `clear`. `check_gpu_descriptors.py` adds
LJ, G1--G5 and converted n2p2 fixtures across CPU/GPU neighbors and MPI layouts.

The benchmark records warmed LAMMPS Loop/Pair times and repeated-sample medians;
logs retain Neigh/Comm times. Its defaults compare CPU 1/4 ranks and GPU 1 rank
at 24/192/5184/24000 atoms. `benchmark_gpu.py` currently pins single-rank runs
to CPU 6; edit that script if CPU 6 is unavailable. The mode benchmark exposes
`--cpu` for this choice and alternates direct/moment/auto across three samples,
with 20 warmup + 100 timed steps. Post-timing energy, force, virial and coordinate
checks are outside the timed region. Avoid competing workloads during timing.

Historical reports retain their measured revisions and workload scope:

- [Initial GPU integration](validation/lammps-gpu-2026-09-26/README.md).
- [Direct/moment comparison](validation/lammps-gpu-modes-2026-09-26/README.md).
- [Moment-force optimization](validation/gpu-moment-force-2026-09-26/README.md):
  24000-atom Ti/O changed from 25.07 to 10.81 ms/step on H100 and 54.25 to
  17.19 ms/step on Blackwell in that experiment.
- [Shared-kernel validation, revision 1.13](validation/energy-common-2026-09-27/README.md):
  21 LAMMPS CPU/GPU comparisons plus library checks.

For current supported descriptors and common-code scope, use the
[implementation status](implementation-status.md). These historical timings
are not a universal direct/moment or CPU/GPU speedup guarantee.
