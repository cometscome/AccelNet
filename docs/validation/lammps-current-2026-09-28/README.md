# Current LAMMPS compatibility probe — 2026-09-28

**Installer follow-up:** 22 Jul 2025 Update 6 is now the recommended release,
and the shared installer supports it directly. See the [installer validation](#shared-installer-follow-up)
and [current build instructions](../../../interfaces/lammps/README.md). The
initial manual probe below is retained as the compatibility evidence.

The existing 29Aug2024 AccelNet adapter can run on **LAMMPS 22 Jul 2025
Update 6 on CPU and GPU after manual installation adjustments**. The **2 Sep
2026 release candidate runs on CPU**, but its GPU adapter does not compile
unchanged. The initial probe is a compatibility investigation, not a claim
that all models/platforms have been validated.

## Releases and installation status

At the time of this check, the official
[latest stable release](https://github.com/lammps/lammps/releases/tag/stable_22Jul2025_update6)
is 22 Jul 2025 Update 6. The newer
[2 Sep 2026 release](https://github.com/lammps/lammps/releases/tag/patch_2Sep2026)
is a stable release candidate.

| Source | Exact commit | CPU | NVIDIA GPU |
|---|---|---|---|
| `stable_22Jul2025_update6` | `9c5ab448c78a14fd534619622162ba418d6a1fb1` | Build/run/comparisons passed | Build/run/21 comparisons passed |
| `patch_2Sep2026` | `d71abe6102c44577442ba7f03b7378a83166b9fd` | Build/run/comparisons passed | Adapter compile failure: removed `particle_split()` |

**At the time of the initial probe, neither version was accepted by the old `install.py`.**
Running it against each unmodified source tree exits with the explicit
29 Aug 2024 version error. Simply deleting that check is insufficient for
2026: its CMake package include loop starts with `GRAPHICS KSPACE`, whereas
the installer searches for a loop starting with `KSPACE`.

The probe used separate temporary source/build directories. The current
CPU/GPU pair styles, lib/gpu adapter, canonical public headers and ACCELNET
CMake module were copied unchanged into each tree. `ACCELNET` was added to
both CMake package lists using the actual loop in that version. The existing
triclinic atom-sorting guard was applied to both matching `atom.cpp` sites.
No descriptor, NN, force, virial or pair-style algorithm was changed.
The initial probe left that installer restriction unchanged. The follow-up
below replaces it with explicit support for the tested stable releases.

## Build configuration

AccelNet source base: `c38c563` (1.1.0), with the uncommitted methods-1.14
atomic-energy preparation optimization already present. The libraries were
built during the preceding README verification. [Source hashes](source-sha256.json)
identify the adapter and changed inference sources used here.

- CPU: GNU C++/Fortran 11.4, Release, `BUILD_MPI=ON`, `BUILD_OMP=OFF`,
  `PKG_ACCELNET=ON`, `PKG_GPU=OFF`; serial common AccelNet library.
- GPU: NVHPC 25.3 `nvc++`/`nvfortran`, CUDA 12.8, `-O2 -DNDEBUG` for C++,
  `BUILD_MPI=ON`, `BUILD_OMP=OFF`, `PKG_ACCELNET=ON`, `PKG_GPU=ON`,
  `GPU_API=cuda`, `GPU_PREC=double`, `GPU_ARCH=sm_90`,
  `CUDA_BUILD_MULTIARCH=OFF`, `ACCELNET_TARGET_FLAGS="-mp=gpu -gpu=cc90"`.
- GPU execution: H100 NVL, UUID
  `GPU-2644154d-7268-af42-6631-59e1f3c6e7f3`, `OMP_NUM_THREADS=1`,
  `OMP_TARGET_OFFLOAD=MANDATORY`.
- MPI wrapper/launcher: `/opt/ompi-cuda/bin/mpicxx` and `mpiexec`.

Build flags otherwise follow the [LAMMPS README](../../../interfaces/lammps/README.md).
Separate GNU and NVHPC AccelNet libraries were used. Configuration and build
logs, staging/check scripts, inputs and snapshots are in [raw logs](raw-logs.tar.gz).
The archived staging script is a one-off probe with local absolute paths,
not a replacement production installer.

## CPU checks

The real Ti/O Chebyshev fixture contains 24 atoms. Each release ran five NVE
steps with neighbor rebuilding/sorting each step, testing:

- `auto`, `direct` and `moment`;
- orthogonal, triclinic and empty-rank geometries;
- one and two MPI ranks.

Both releases passed all **nine two-rank versus one-rank comparisons**,
covering total/per-atom energy, forces, virial pressure and coordinates.
The same nine one-rank cases per release were also compared with the existing
29 Aug 2024 Update 4 CPU pair style (in the prior NVHPC executable): **18
cross-version comparisons passed**. This latter control uses a different
compiler family, not a CPU performance baseline.

Cross-version maximum total-energy difference was zero; maximum force
component difference was `4.671263376110346e-14` eV/Angstrom. The tolerances
were `atol=2e-5, rtol=2e-9` for thermo values and `atol=2e-8, rtol=2e-9`
for dump coordinates/forces/per-atom energies; atom identities/types matched
exactly.

Results: [stable MPI](stable-cpu-mpi.json), [candidate MPI](rc-cpu-mpi.json),
[cross-version comparisons](cpu-cross-version.json).

## Stable-release GPU checks

The unchanged `check_gpu.py` ran on the staged stable-release executable:

```sh
CUDA_VISIBLE_DEVICES=GPU-2644154d-7268-af42-6631-59e1f3c6e7f3 \
OMP_NUM_THREADS=1 OMP_TARGET_OFFLOAD=MANDATORY \
python3 interfaces/lammps/29Aug2024/tests/check_gpu.py \
  --lammps /tmp/accelnet-lammps-current-stable-gpu/lmp \
  --golden /home/nagai/AccelNetGPU/AccelNetTrainer.jl/test/data/fortran_predict \
  --output /tmp/accelnet-lammps-current-stable-gpu-check --ranks 1 2 --steps 5
```

All **21 comparisons** passed: three geometries, auto mode, CPU/two ranks
plus GPU `neigh no/yes/hybrid` on one/two ranks, relative to CPU/one rank.
Eighteen of these comparisons execute the GPU pair style.

Maximum observed differences were:

| Quantity | Maximum absolute difference |
|---|---:|
| Total energy | `7.275957614183426e-12` eV |
| Force component | `5.684341886080802e-14` eV/Angstrom |
| Per-atom energy | `5.684341886080802e-14` eV |
| Virial pressure component | `1.0331859812140465e-9` bar |

See [all results](stable-gpu.json). This run does not establish performance,
Blackwell compatibility, multi-GPU behavior, lifecycle/restart behavior or
coverage of every n2p2 descriptor on this newer LAMMPS release.

## Release-candidate GPU incompatibility

CMake configuration and `pair_accelnet_gpu.cpp` compilation passed, but
`lal_accelnet_ext.cpp` failed:

```text
class "LAMMPS_AL::Device<double, double>" has no member "particle_split"
```

The [compiler log](rc-gpu-adapter-build.log) captures the failure. The 2026
GPU package removed host/device particle splitting; its `fix_gpu.cpp` accepts
`split` only as a deprecated, ignored option. AccelNet's old check of
`global_device.particle_split()` therefore cannot compile. The stable 2025
version still has this method and [compiled the adapter](stable-gpu-adapter-build.log).

A supported 2026 port needs version-aware package installation and GPU API
compatibility handling, followed by the same numerical, neighbor/MPI and
lifecycle checks. Removing the obsolete split check may resolve this immediate
compile error; it does not by itself prove GPU runtime compatibility. No 2026
GPU executable was linked or run in this investigation.


## Shared installer follow-up

The recommended release is now **22 Jul 2025 Update 6**. Use:

```sh
python3 interfaces/lammps/install.py /path/to/lammps
```

The installer explicitly accepts that release/update and **29 Aug 2024
Update 4**. It registers CMake packages, copies the current canonical headers
and the same CPU/GPU adapter sources, and applies the triclinic sorting guard.
The historical `29Aug2024/install.py` path delegates to this shared installer.
No version-specific numerical kernel or duplicate pair style was added.
The README recipes now pin `stable_22Jul2025_update6` and use this entry point
for both CPU and GPU builds.

Validation after this change:

- All **five installer tests passed**, covering both releases, repeat
  installation, canonical headers, rejection without modification, malformed
  source registration/sorting and the old command path. CI runs these tests
  with `python3 -m unittest discover -s interfaces/lammps/tests -v`.
- Installed twice into a fresh archive of the stable tag. All ten installed
  or patched files matched the updated existing tree byte for byte; see
  [file hashes](installer-file-sha256.json). Fresh CPU CMake configuration passed.
- The real 2024 Update 4 tree accepted the compatibility entry point.
- Rebuilt the stable CPU and GPU executables after running the production
  installer. **Nine CPU comparisons** (three modes, three geometries, MPI 1/2)
  and **21 H100 comparisons** (auto, three geometries, CPU/GPU neighbor modes,
  MPI 1/2) passed again: [CPU results](installer-cpu.json),
  [GPU results](installer-gpu.json). These retain the initial probe's models,
  five-step trajectories, checks and tolerances.

[Follow-up logs](installer-logs.tar.gz) contain configure/build output and
numerical-test inputs, snapshots and logs. The 2026 candidate remains rejected;
its GPU incompatibility is unchanged. Blackwell was not rerun on the 2025
release, and no performance claim is made by this installer change.
