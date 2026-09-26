# GPU residency and Chebyshev moments — 2026-09-25

This change retains model parameters and scratch buffers on the GPU, adds angular
moment descriptors and their force contraction, and exposes phase timings. Both
monomial construction and neighbor-edge contraction now expose parallel work
within each atom. The original CPU inference kernels and compiler flags remain
unchanged. See [API/build instructions](../../openmp-target.md).

## Environment and method

Base commit: `c6631460a1bbb990c82e3e0ff5e73c36c52f6f9b`; changes are in the working
tree. [Source SHA-256 values](sources.json) identify the tested GPU code.
NVHPC 25.3 Release (`-fast -O3`), `-mp=gpu -gpu=cc90,cc120`; GNU 11.4 for the CPU
regression gate and optional-backend syntax/negative checks.

- H100 NVL: `GPU-2644154d-7268-af42-6631-59e1f3c6e7f3`.
- RTX PRO 6000 Blackwell Max-Q: `GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9`.
- Xeon Gold 6526Y; GPU timing driver pinned to CPU 6, `OMP_NUM_THREADS=1`.

The H100 was sharing another process's workload during this work. Numerical tests
run on both GPUs; performance comparisons below use the available Blackwell.
No isolated H100 performance claim is made. CPU/GPU comparisons are single CPU
thread versus one GPU, not an all-core CPU comparison. Shared-host clocks/load
can change timings.

## Correctness and CPU preservation

| Check | Result | Evidence |
| --- | --- | --- |
| H100 GPU CTest suite | 17/17 pass | [Log](h100-final-tests.log) |
| Blackwell GPU CTest suite | 17/17 pass | [Log](blackwell-final-tests.log), [details](blackwell-detailed-tests.log) |
| One workspace, H100 → Blackwell → H100 | pass, max E/F/W error 1.11e-16 | [Log](device-switch.log) |
| GNU CPU suite including speed gate | 44/44 pass | [Log](cpu-final-tests.log) |
| GNU optional backend compile and rejection tests | 4/4 pass | [Log](gnu-final-tests.log) |
| Existing CPU API performance vs baseline | 12/12 pass, identical printed outputs | [Samples](cpu-performance.json) |

The GPU suite runs auto, forced-direct and forced-moment variants, each with 59
synthetic and 18 real Ti/O model/geometry/reference-mode comparisons. These also
exercise coordinate and strain finite differences, all cutoffs/activations,
Chebyshev versions 0/1/10, species-specific networks, periodic images, subsets,
empty/isolated systems and additive outputs. Array comparisons allow
`2e-10 + 2e-10*abs(reference)`; results agree numerically, not bit for bit.
Blackwell maximum absolute E/F/W differences were 1.024e-15 (synthetic) and
3.425e-12 (real Ti/O).

New lifetime checks cover stable upload/allocation counts, repeated buffer growth,
model copies, workspace assignment, same-shape parameter reloads, scope
finalization and explicit release. Every timing run fails if steady-state calls
change model-upload or buffer-growth counters; all final runs report `RESIDENCY 1 1`.
This assertion concerns application buffers, not runtime-internal allocation pools.

The CPU performance gate uses pre-change `c663146` and candidate GNU builds,
12 cases (Chebyshev/LJ/G4/G5 × 8/64/512 atoms), seven alternating-order samples,
0.2 seconds minimum per sample, pinned to CPU 0. Paired median candidate/baseline
time ratios were **0.9937–1.0089** (10% regression gate); all printed E/F/W values
were identical. This confirms preservation on the tested fixtures.

## Device-memory diagnostics

Compute Sanitizer memcheck with standard runtime settings passed the complete
59-case forced-moment synthetic suite, including finite differences and lifetime
checks, with **0 access errors** on [H100](h100-access-sanitizer.log) and
[Blackwell](blackwell-access-sanitizer.log). No API/access errors were suppressed.

An additional diagnostic build queried `omp_target_is_present` for every model
and scratch array immediately after each `target exit data map(delete:)`, before
host deallocation. All mappings were absent through the complete 59-case moment
suite, including growth, reloads, explicit release and scope finalization:
[log](mapping-release-check.log), [instrumented API](diagnostic-mapping-check.f90).
This diagnostic is separate from the performance binary and verifies OpenMP
mapping ownership; it does not claim the runtime returns its cached allocations
to the driver immediately.

Strict `--leak-check full` is **not clean** with NVHPC 25.3: it reports 88 retained
allocations (10,560,056 bytes) at process exit in the full suite, including runtime
and pooled data allocations ([original log](h100-sanitizer.log)). A standalone
OpenMP program without AccelNet, with `NV_ACC_MEM_MANAGE=0`, still reports four
runtime allocations (2,124,856 bytes): [source](runtime_array_probe.f90),
[log](runtime-array-sanitizer.log). Turning off memory management for AccelNet's
test instead produces a shutdown failure inside `libnvomp`'s
`finalizeDeviceMemoryPool`; this occurs with and without the sanitizer:
[normal run](blackwell-nopool-numerical.log), [debugger stack](nopool-gdb.log),
[sanitizer](h100-nopool-sanitizer.log). Therefore **a clean process-exit leak
check is not claimed**, and disabling this runtime setting is not recommended
as a workaround. Standard settings passed all numerical/access checks.

## Blackwell timings

Milliseconds per call; medians of five samples, each at least 0.2 seconds.
Both versions include transfers, use an already-built neighbor list for batch
timing, and check E/F/W before/after measurements. "Before" is the saved stages
2–3 executable ([saved sources and checksums](before-sources/checksums.json));
"after" uses current auto mode. Initial device setup and packing
are excluded; before still allocates/maps each call, after reuses buffers.
"Full" includes rebuilding the neighbor list on the CPU every call.
Speedups below divide median times; greater than one is faster.

| Atoms | Order | CPU batch ms | Before GPU ms | After GPU ms | Before/after | CPU/GPU batch | CPU full ms | GPU full ms |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 3 | 0.064 | 9.621 | 0.532 | 18.07× | 0.12× | 0.302 | 0.774 |
| 64 | 3 | 0.491 | 7.418 | 0.535 | 13.87× | 0.92× | 2.337 | 2.378 |
| 512 | 3 | 3.951 | 10.092 | 0.666 | 15.14× | 5.93× | 125.601 | 122.591 |
| 1024 | 3 | 7.142 | 9.494 | 0.766 | 12.40× | 9.32× | 558.036 | 551.866 |
| 512 | 8 | 16.091 | 12.631 | 1.393 | 9.07× | 11.55× | 137.666 | 123.161 |

[Raw medians](final-performance.json) and `final-before-*`/`final-after-*` logs in
this directory retain all samples. The main 512-atom/order-8 improvement is
**9.07× over the previous GPU backend**, or **11.55× versus the single-thread CPU
batch**. Full evaluation improves only to 123.16 ms versus CPU 137.67 ms
(1.12×): CPU neighbor construction remains the dominant cost. Small systems
still favor the CPU; existing CPU callers are not automatically redirected.

For 512 atoms/order 8/spacing 1.7, median phase times were:

| Phase | ms |
| --- | --- |
| Host validation/cache comparison | 0.0195 |
| Input packing/upload | 0.0888 |
| Geometry, radial and moment descriptors | 0.6670 |
| NN and input gradient | 0.1205 |
| Force contraction and scatter/virial | 0.4682 |
| Download and additive accumulation | 0.0212 |
| Total GPU API | 1.3925 |

These are synchronous wall-clock phases including launch/synchronization, not
CUDA-event kernel-only measurements. Independently computed phase medians need
not sum exactly to the median total. [Full profile samples](final-after-512-8.log).

At denser spacing 1.1 (512 atoms/order 8), current direct mode took **19.351 ms**,
forced moments **4.074 ms**, and auto **4.087 ms** including transfers: moments
were **4.75× faster than current pair enumeration**. This comparison isolates
the method choice within the optimized backend; the overall before/after gain
also includes data residency, cached geometry and increased parallelism.
[Direct](final-dense-mode1.log), [moment](final-dense-mode2.log),
[auto](final-dense-mode0.log).

Auto selects moments when `angular_neighbors * (order+1) >= number_of_monomials`.
This initial work estimate was informed by low/high-order and dense/sparse
fixtures. Forced modes remain available; optimal thresholds on other devices
have not been established.

## Reproduction and remaining scope

```sh
export CUDA_VISIBLE_DEVICES=GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9
export OMP_NUM_THREADS=1
ctest --test-dir build-target -R predictor_target --output-on-failure
taskset -c 6 build-target/bin/accelnet-target-benchmark 512 8 0.2 0 1.7
taskset -c 6 build-target/bin/accelnet-target-benchmark 512 8 0.2 1 1.1
taskset -c 6 build-target/bin/accelnet-target-benchmark 512 8 0.2 2 1.1
```

Only Chebyshev models with one component per element are supported on the GPU.
CPU neighbor construction, CLI/C/LAMMPS GPU integration, GPU n2p2/LJ/composite
models and AMD/Intel hardware validation remain outside this change. Scratch
storage grows with batch rows and neighbor edges; callers can partition the
CSR batch to bound memory. GPU performance is measured here for 8–1024 atoms,
not claimed for every size/model/device.
