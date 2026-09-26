# LAMMPS GPU package validation — 2026-09-26

## Implementation and environment

- LAMMPS 29 Aug 2024 Update 4, optional ACCELNET and GPU packages.
- NVHPC 25.3 Fortran (`-fast -O3 -mp=gpu -gpu=cc90,cc120` for target code);
  NVHPC C++ (`-O2 -DNDEBUG`); CUDA 12.8 GPU library, FP64, sm_90 + sm_120.
- Open MPI 5.0.7 (`/opt/ompi-cuda`), Xeon Gold 6526Y.
- H100 NVL `GPU-2644154d-7268-af42-6631-59e1f3c6e7f3`.
- Two RTX PRO 6000 Blackwell Max-Q GPUs:
  `GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9` and
  `GPU-ebe7b5e0-03cd-5b5f-55e3-d6a16ed996f5`.
- Immutable Ti/O networks and `structure0001.xsf` from
  `AccelNetTrainer.jl/test/data/fortran_predict`. LAMMPS type 1=Ti, type 2=O.
- Only isolated LAMMPS source/build copies under `/tmp/accelnet-lammps-gpu`
  were used; the pre-existing LAMMPS source/build was not overwritten.

Usage and reproducible build/test commands: [lammps-gpu.md](../../lammps-gpu.md).

## Numerical results

All comparisons include total energy, six configurational pressure components,
per-atom energy, all local-atom forces after ghost reverse communication,
positions, and a 20-step NVE trajectory. Floating-point summation order may differ;
bitwise identity is not required.

| Test | Comparisons | Settings |
| --- | ---: | --- |
| [H100](h100-correctness.json) | 33 | CPU and GPU, 1/2/4 ranks; orthogonal/triclinic/empty ranks; neigh no/yes/hybrid |
| [Two Blackwell GPUs](blackwell-2gpu-correctness.json) | 24 | CPU and GPU, 2/4 ranks; same geometry/neighbor cases |
| [Migration and neighbor reuse](migration-reuse-correctness.json) | 22 | CPU and GPU, 1/2/4 ranks; direct/moment; neighbor rebuilding every 5 steps; imposed drift crosses rank and periodic boundaries |

Across these 79 comparisons: maximum force difference `6.38e-14 eV/Angstrom`,
maximum per-atom energy difference `2.28e-13 eV`, maximum total energy difference
`1.46e-11 eV`, maximum configurational pressure difference `1.76e-9 bar`.
The test thresholds are `2e-8 + 2e-9*abs(reference)` for atom data and
`2e-5 + 2e-9*abs(reference)` for thermodynamic data.

The existing CPU pair mistakenly tallied every atom's energy to atom 0.
The CPU interface now uses the actual center index for this tally. Global energy
and forces are unaffected; per-atom comparisons above use the corrected tally.

## Regressions and error handling

- [CPU CTest: 44/44 passed](cpu-tests.log), including the original energy/force/virial,
  binary/ASCII, n2p2, and CPU performance tests. Required test targets were rebuilt.
- [CPU performance gate](cpu-performance.json): 12 model/size combinations, all passed;
  candidate/baseline time ratio approximately `0.998–1.015`, numerical error zero.
  The gate permits at most 10% slowdown and uses repeated, pinned runs.
- [Fortran GPU CTest: 19/19 passed](gpu-tests.log). Includes new C API coverage for
  independent instances, mode agreement, malformed input/CSR, and repeated creation/destruction.
- [Six explicit rejection cases with two MPI ranks](mpi-error-tests.log): newton off,
  missing/truncated network, reversed model order, per-atom stress, neighbor exclusions.
  All terminate with the intended error rather than a crash or timeout.
- GNU 11.4 can still compile the optional target library with the host OpenMP test
  configuration; this is a compilation check, not a claim of GNU GPU execution.
- Installer applied twice to an isolated source subset: identical file contents,
  both triclinic GPU-sort guards installed.

## Runtime interoperability and memory checks

`neigh yes/hybrid` uses the GPU package's real device neighbor matrix. A Fortran
GPU kernel computes CSR offsets, indices and displacements from borrowed CUDA
pointers. Neighbor entries are not downloaded and re-uploaded. The edge count
(integer scalar) is copied to the host. Forces are returned before MPI reverse
communication; FixGPU's answer queue remains empty for AccelNet.

Four integration issues were reproduced and resolved:

1. Explicitly ordering NVHPC runtime libraries before its offload runtime caused
   `omp_get_num_devices()` to return zero. The final link uses `-fortranlibs` and
   avoids LAMMPS appending that runtime list again.
2. NVHPC 25.3 invoked the workspace assignment routine during elemental
   finalization through a C handle. Explicit scalar/rank-1 finalizers and the
   repeated handle-lifetime test cover the fix.
3. LAMMPS 2024 GPU auto-sorting read orthogonal subdomain coordinates in a triclinic
   box, producing uninitialized bin extents and crashes. The included patch uses
   the normal bbox-based sorter for triclinic boxes. Sorting remains enabled.
4. After `clear`, lib/gpu released the primary CUDA context while NVHPC still
   cached its kernel/pool state. The adapter now retains one reference per used
   device for the process lifetime. Model/workspace allocations are still released
   per instance; CUDA reclaims the retained context when the process exits.

[Compute Sanitizer memcheck](sanitizer-host-mpi.log) on the complete Blackwell
LAMMPS triclinic `neigh yes` trajectory reports **0 errors**. The run used
`OMPI_MCA_accelerator=null OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm,tcp`, which selects
host MPI buffers, as used by this pair style. No sanitizer suppressions or disabled
memory-access checks were used. `--leak-check no` was explicit; this is not a
process-exit leak-clean claim.

With default CUDA-aware MPI, [the initial sanitizer run](sanitizer-default-mpi.log)
reported 1,774 CUDA API errors; the printed diagnostics came from MPI/UCX host-pointer classification/probing
(`cuMemRetainAllocationHandle`, `cuPointerGetAttribute`, `cuCtxSetFlags`). Selecting
host MPI eliminated those diagnostics. Ordinary multi-rank correctness tests used
the default MPI configuration. Known NVHPC runtime pooling limitations from the
[standalone GPU validation](../gpu-residency-moments-2026-09-25/README.md) still apply.

## LAMMPS timings (H100)

Median of 3 independent runs, 10 warmup + 50 measured NVE steps, timestep
0.0001 ps. Neighbor check every 10 steps with skin 0.6 Angstrom; atom sorting
every 100 steps. No other task GPU benchmarks ran concurrently. Single-rank
runs were pinned to CPU 6; CPU 4-rank runs used MPI core binding. These are
Loop times, including communication and integration, not just GPU kernel time.

| Atoms | CPU 1 rank, ms/step | CPU 4 ranks | H100, neigh no | H100, neigh yes | CPU1/GPU yes | CPU4/GPU yes |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 24 | 1.462 | 1.019 | 1.551 | 1.573 | 0.93x | 0.65x |
| 192 | 11.673 | 3.164 | 1.881 | 1.603 | 7.28x | 1.97x |
| 5,184 | 226.642 | 79.845 | 12.626 | 5.114 | 44.32x | 15.61x |
| 24,000 | 1053.766 | 537.142 | 74.359 | 25.053 | 42.06x | 21.44x |

[Samples and Pair times](h100-performance.json); representative per-case logs
are in [benchmark-logs](benchmark-logs). The GPU package neighbor path is
2.47x faster than the host-neighbor GPU path at 5,184 atoms, and 2.97x at 24,000
atoms. This includes eliminating host CSR/displacement packing and edge-data
uploads; it is not a measurement of neighbor building alone.

The benchmark's total energy and six virial-pressure trajectories were also
[checked against CPU 1 rank](benchmark-thermo-check.json), for all 48 runs,
including 24,000 atoms.

These results concern this Ti/O network and geometry. They do not establish
a universal GPU speedup. Small systems can favor CPU execution; backend
selection remains explicit. CPU 4-rank timing at 24,000 atoms has visibly more
variation (0.523–0.578 s/step); retain the raw samples when comparing systems.

## Final lifetime and larger-system checks

- [Two-rank lifecycle test](lifecycle.log): 24 atoms → 192 atoms while retaining
  the model, then `clear` and reload at 24 atoms; E/F/positions agree with CPU.
- [H100 lifecycle memcheck](sanitizer-lifecycle.log): **0 errors**, including the
  complete growth/clear/reload sequence, with host MPI and leak checking disabled.
- [5,184-atom Blackwell comparison](large-blackwell-correctness.json): full force,
  per-atom energy, total energy and pressure match CPU for neigh no/yes/hybrid.
  Maximum force difference is below `4.0e-14 eV/Angstrom`; total-energy difference
  is below `3.9e-8 eV` for this larger extensive energy.
- [Final two-GPU triclinic rerun](final-multigpu-correctness.json): passed after the
  process-lifetime context fix.

The timing samples above precede the final context-retention change, which only
runs at pair initialization. The force computation and neighbor-update paths are
unchanged. The final executable and source hashes are recorded in
[provenance.json](provenance.json).
