# Unified potential-inference APIs — methods revision 1.12

Measured 2026-09-27 (JST), AccelNet 1.0.1. Baseline: `gpu` checkpoint
**543b180** (numerical checkpoint ea65290). Final source checksums are in
[source-sha256.json](source-sha256.json); compiler flags and executable hashes
are in [build-audit.json](build-audit.json).

## Scope

Structure/file, ordinary CSR batch, Fortran/C atomic energy/force/virial,
ænet-compatible SFB, and LAMMPS potential evaluation now reach the same numerical
kernels. The ordinary CPU compilation removes OpenMP directives. The target
compilation uses OpenMP host threads or GPU offload. Multiple Chebyshev
components and Chebyshev/LJ/Behler mixtures preserve each element's original
input offsets, scaling and NN, including different component lists by element.
There is no production fallback to the former CPU potential evaluator.

The former numerical routines live in [legacy/cpu-reference](../../../legacy/cpu-reference/README.md).
Five numerical include bodies match the baseline byte-for-byte except for one
new provenance comment; three complete pre-migration API files are exact
snapshots. See [reference-audit.json](reference-audit.json). Low-level standalone
descriptor/NN compatibility utilities and explicit reference tests still compile
these routines. They are not routed into the production potential-inference
path. Types, model loading, neighbor generation and shared scalar formulas stay
in the normal source directories.

## Correctness and runtime validation

| Final check | Result |
|---|---|
| GNU 11.4, OpenMP compiled out | 51/51 functional CTests passed |
| GNU array bounds/runtime checks | 14/14 focused CTests passed |
| H100 NVL, mandatory offload | 34/34 GPU CTests passed |
| RTX PRO 6000 Blackwell Max-Q, mandatory offload | 5/5 focused GPU CTests passed |
| GNU OpenMP host, 2 and 8 physical cores | 9 mixed-component cases each passed; maximum E/F/W absolute difference 1.11e-16 |
| H100 compute-sanitizer 13.1 memcheck | All 81 composite cases; zero errors |
| LAMMPS GPU 29 Aug 2024 Update 4 | 21 CPU/GPU comparisons passed, 1/2 MPI ranks, orthogonal/triclinic/empty rank and GPU neighbor no/yes/hybrid |
| Existing performance CTests | 3/3 passed, including the unchanged 1.10 force-inclusive CPU regression threshold |
| Identical public-API driver, baseline/current | 19 cases × 4 reversed rounds; energy and every force/virial component checked |

The full composite test spans descriptor versions 0/1/10, three component
families, auto/direct/moment and three geometries, plus coordinate finite
differences where applicable. Lifecycle tests cover reloads, copied models,
workspace reuse and explicit release. SFB is compared against the independent
retained descriptor evaluator. GPU suites include real ænet Ti/O networks and
independent n2p2 checks, not just comparisons of two common-kernel wrappers.

The H100 JUnit and memcheck transcripts are gzip-compressed without changing
their contents.

One earlier host run with `OMP_WAIT_POLICY=PASSIVE GOMP_SPINCOUNT=0` timed out
after 300 seconds; its empty output is retained as
`openmp-composite-passive.log`. The default-wait full 81-case run passed; final
2/8-thread checks use the bounded 9-case composite subset. These are correctness
checks, not a claim of efficient thread scaling or support for every wait policy.

## CPU public-API timings

GNU Fortran 11.4, `-O3`, both libraries compiled without OpenMP; CPU 6 on Xeon
Gold 6526Y. The exact same driver and fixtures are linked separately against the
immutable baseline and the final library. Four alternating/reversed rounds,
minimum 0.15 s per timing; medians below. Per-call E/F/W output logs are
retained in `public-performance-final/raw-outputs.tar.gz`. No benchmark overlapped our builds or
correctness runs. Clock frequency was not locked, so small differences are not
claimed as improvements. [Raw timings and commands](public-performance-final/report.json).

`structure` includes neighbor construction and E/F/W; `energy` includes neighbor
construction and energy only. `batch` uses prepared CSR and evaluates E/F/W;
`atomic`/`atomic-energy` use prepared environments with one call per center.
All times cover an entire structure. “n2p2” below means an imported n2p2 model
inside AccelNet, **not** execution time of the upstream n2p2 executable.

| Model | Atoms | API | Before (ms) | Common (ms) | Common / before |
|---|---:|---|---:|---:|---:|
| aenet | 192 | structure | 39.169 | 34.132 | 0.871 |
| aenet | 192 | energy | 24.758 | 28.569 | 1.154 |
| aenet | 192 | batch | 12.465 | 11.099 | 0.890 |
| aenet | 192 | atomic | 18.777 | 10.831 | 0.577 |
| aenet | 192 | atomic-energy | 4.271 | 5.363 | 1.256 |
| n2p2 | 512 | structure | 141.562 | 133.031 | 0.940 |
| n2p2 | 512 | energy | 141.072 | 133.645 | 0.947 |
| n2p2 | 512 | batch | 6.663 | 6.232 | 0.935 |
| n2p2 | 512 | atomic | 13.991 | 7.118 | 0.509 |
| n2p2 | 512 | atomic-energy | 10.135 | 5.770 | 0.569 |
| combined | 128 | structure | 10.692 | 9.048 | 0.846 |
| combined | 128 | energy | 8.175 | 8.783 | 1.074 |
| combined | 128 | batch | 2.643 | 1.005 | 0.380 |
| multi-chebyshev | 128 | structure | 9.662 | 9.538 | 0.987 |
| multi-chebyshev | 128 | energy | 8.913 | 8.799 | 0.987 |
| multi-chebyshev | 128 | batch | 1.739 | 1.139 | 0.655 |
| mixed-components | 128 | structure | 11.776 | 9.926 | 0.843 |
| mixed-components | 128 | energy | 9.345 | 10.013 | 1.072 |
| mixed-components | 128 | batch | 3.765 | 1.740 | 0.462 |

Forces-inclusive paths improved or remained close in this sweep. The remaining
regression is concentrated in energy-only paths: ænet per-atom energy is about
26% slower, and structure energy is about 15% slower. Combined/mixed synthetic
energy-only calls are about 7% slower. They are not hidden by an aggregate pass
flag: the JSON `passed` field denotes numerical equivalence, not a performance
threshold. The large original ~91% atomic-energy regression was reduced by
four-moment tiling, skipping unused derivatives and explicit private-cache
invalidation. Eight-moment tiling was rejected after paired CPU/GPU probes.

## GPU before/after timings

NVHPC 25.3, `-fast -O3 -mp=gpu -gpu=cc90,cc120`, H100 NVL, mandatory offload.
Host driver pinned to CPU 6. Two reversed rounds, five warmed timing samples
per executable per round, at least 0.05 s each. Same compiler/options and driver
in both revisions. [Raw timings, numerical errors and commands](gpu-performance-final/report.json).
These are target-batch E/F/W times without neighbor construction, not MD loop
times. Direct and moment are forced here; ordinary defaults remain auto.

| Family | Atoms | Order | Method | Before (ms) | Common (ms) | Common / before |
|---|---:|---:|---|---:|---:|---:|
| chebyshev | 2048 | 8 | direct | 2.629 | 2.504 | 0.952 |
| chebyshev | 2048 | 8 | moment | 2.310 | 2.233 | 0.967 |
| g5-series | 2048 | 4 | direct | 4.769 | 4.751 | 0.996 |
| g5-series | 2048 | 4 | moment | 2.977 | 2.920 | 0.981 |
| g5-high | 256 | 16 | direct | 2.486 | 2.507 | 1.009 |
| g5-high | 256 | 16 | moment | 1.971 | 1.963 | 0.996 |

The migration adds composite-model support; this timing sweep measures existing
single/grouped models to detect regressions. It does not establish composite
GPU speedups at every size, or performance on AMD/Intel GPUs.

## Reproduction

`AccelNetPredictor/benchmark/compare_public_apis.py` accepts `--before`, `--after`,
`--aenet`, `--n2p2` and `--output`. Build `accelnet-public-api-benchmark` for the
current tree. For the baseline, build libraries from `git archive 543b180`, then
compile the **current** `batch_test_support.f90` and `public_api_benchmark.f90`
against those libraries using identical compiler/options. This preserves the
same composite fixtures and driver in both executables. The JSON records exact
paths, hashes, commands, samples and errors.

The archived `gpu_compare.py` records this machine's exact GPU sweep; adapt
its executable/output paths for another checkout. `tile_probe.py` and the probe
JSON files preserve the tuning comparison. Build flags are archived separately.

The ordinary CPU performance CTest now selects the retained original structure
evaluator explicitly (`--baseline-mode reference`), so migration of the public
structure API cannot make the regression test compare the new implementation
against itself. The historical n2p2 type-21 1.101346 ratio versus a strict 1.10
parity limit remains documented in the earlier atomic-removal report; this
migration does not claim to resolve that separate gate.
