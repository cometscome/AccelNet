# Common-kernel energy-only recovery — methods revision 1.13

Measured 2026-09-27 (JST), AccelNet 1.0.1. CPU baseline: **543b180**, before
unifying the production potential APIs. GPU baseline: **e6a96d3**, the previous
common implementation (the older GPU API did not expose energy-only execution).
Both before/after pairs use identical benchmark sources and compiler options.
Source/build hashes and raw measurements accompany this report.

## Scope and method

The production CPU/GPU kernels remain common. Species-partitioned Chebyshev
moments accumulate each neighbor once, combine weighted channels algebraically,
and retain the existing force contraction. Explicit resident scratch supports
NVHPC device calls. The host uses a team-size limit of 32 and the GPU 4; serial
CPU compilation removes OpenMP directives. Energy-only skips unused force
clears/transfers, and unchanged validated workspace capacities skip shape scans.
Odd row capacities avoid unfavorable column strides in moment/descriptor/NN
loops. See [Section 24](../../../speedupmethods.md#24-restore-energy-only-cpu-performance-in-the-common-kernels-revision-113)
for the equations and experiments.

One-Chebyshev-component-per-element is **not** a public restriction since
`e6a96d3`. Complete basis blocks can repeat and mix with other descriptors;
the private component initializer accepts one already-separated block.

n2p2 types **13/21/24 are implemented** with exact direct values, analytic
forces and virial. Their general forms do not have the finite separable
single-neighbor moment representation used by G5. This is a moment-method
boundary, not a missing potential type. General third-distance/angle windows
cannot be replaced by a truncated expansion under the exactness requirement.

## Measurement conditions

* CPU: GNU Fortran 11.4, `-O3`, OpenMP compiled out in both libraries, Xeon
  Gold 6526Y physical CPU 6. The script inspects undefined symbols to reject
  linked OpenMP runtimes. Energies and all relevant force/virial components
  are compared to the independent archived evaluator.
* GPU: NVHPC 25.3, `-fast -O3 -mp=gpu -gpu=cc90,cc120`, mandatory offload.
  Timings use the idle RTX PRO 6000 Blackwell Max-Q, UUID
  `GPU-ebe7b5e0-03cd-5b5f-55e3-d6a16ed996f5`, with the host on CPU 6.
  Both H100s and the other Blackwell were occupied by unrelated jobs; no
  H100 speed claim is made. A contaminated H100 probe was excluded.
* Executables alternate in reversed rounds. No timed comparison overlaps our
  builds, correctness tests or another benchmark. Clocks are not locked;
  small differences should not be interpreted as certain speedups.
* `energy` includes neighbor construction; `atomic-energy` uses prepared
  environments, one public API call per center. Times cover the entire
  structure. GPU times exclude neighbor construction and model setup and
  include the ordinary per-call transfer/validation path.

## Reproduction

The archived JSON contains executable hashes, exact commands and samples.
Use an immutable 543b180 public-API executable and the current serial executable:

```sh
python3 AccelNetPredictor/benchmark/compare_public_apis.py \
  --before /path/to/543b180/accelnet-public-api-benchmark \
  --after /path/to/current/accelnet-public-api-benchmark \
  --aenet /path/to/Ti-O-models \
  --n2p2 AccelNetPredictor/test/data/n2p2-virial-angular \
  --aenet-energy-only --max-energy-ratio 1.10 \
  --rounds 6 --seconds .3 --output /tmp/energy-comparison
```

For CTest registration, set `ACCELNET_PUBLIC_API_BASELINE_EXECUTABLE` to that
archived executable and run `predictor_aenet_energy_performance`. It uses the
existing `ACCELNET_CPU_MAX_SLOWDOWN` limit (1.10), not a relaxed threshold.
`--aenet-energy-only` exercises 64, 192 and 512 atoms. Omitting it runs the
broader 19-case public API comparison. Both modes check numerical agreement.

`accelnet-target-energy-benchmark N ORDER SECONDS MODE MODEL_DIR_OR_DASH`
checks independent-reference E/F/W, then times energy-only while checking that
nonzero force/virial sentinels remain untouched. Compile this same source with
both e6a96d3 and current target libraries. The retained GPU comparison script
runs auto/direct/moment cases; force timings use the existing target driver.

## Rejected or superseded probes

Global-neighbor moment layouts, two-/eight-lane tiling, and plane work items did
not restore CPU energy performance. A whole-center typed-moment kernel restored
CPU speed but cost about 13--14% on Blackwell even after team-size tuning.
Partitioning by species recovered GPU parallelism while retaining the same
CPU arithmetic. An early assumed-shape device helper faulted; explicit-shape
arguments and resident scratch fixed the invalid device read. Failed prototype
runs are not counted as final validation.

Before row padding, the 512-atom structure-energy probe remained about 10%
slower than the original CPU evaluator. A paired experiment changed only row
capacity and reduced this to about 3%; force batch time also decreased. This
supports a cache-layout explanation, but hardware cache counters were not
collected. Final timings below supersede those probes.


## Final correctness and performance gates

| Check | Result |
|---|---|
| GNU serial functional/reference CTests | 52/52 passed |
| GNU bounds/runtime-check focused CTests | 22/22 passed |
| Blackwell mandatory-offload CTests | 35/35 passed, including upstream n2p2 multi-element validation |
| H100 new species-moment cases | 18 passed; maximum absolute E/F/W error 2.22e-16 |
| H100 compute-sanitizer 13.1 memcheck | Same 18 cases; zero errors |
| GNU host OpenMP, 2 and 8 physical cores | 18 species and 9 composite cases per thread count passed |
| LAMMPS 29 Aug 2024 Update 4 | 21 comparisons passed, 1/2 ranks, GPU neighbor no/yes/hybrid |
| Performance CTests | 4/4 passed, including the new energy-only 1.10 limit |
| Same public-API driver against archived/current libraries | All 19 cases × 4 reversed rounds numerically passed |

The H100 test and memcheck transcripts are gzip-compressed without editing
their original whitespace.

The new species cases span three/four elements, descriptor versions 0/1/10,
auto/direct/moment, rotated global/local species mappings, zero/negative/
fractional species weights, and finite differences. The existing broad tests
now also check energy-only with nonzero force/virial sentinels after a force
call, then reuse the workspace for further force calculations. This covers
composite blocks, n2p2 descriptors, reloads and shape changes as well.

## Final CPU timings

The [energy regression report](cpu-energy-final/report.json) uses four reversed
rounds, at least 0.15 s per timing, on CPU 6. Ratios are current/original; smaller
is faster. Structure energy includes neighbor construction. All six cases meet
the unchanged 1.10 regression limit. Clock variability prevents interpreting
percent-level differences as guaranteed improvements.

| Atoms | API | Original CPU (ms) | Common CPU (ms) | Ratio |
|---:|---|---:|---:|---:|
| 64 | energy | 4.696 | 4.726 | 1.007 |
| 64 | atomic-energy | 1.666 | 1.718 | 1.031 |
| 192 | energy | 24.880 | 25.081 | 1.008 |
| 192 | atomic-energy | 4.298 | 4.327 | 1.007 |
| 512 | energy | 146.466 | 150.547 | 1.028 |
| 512 | atomic-energy | 13.756 | 14.038 | 1.021 |

The earlier 15--26% energy-only loss at 192 atoms is removed in this comparison.
The maximum remaining measured difference in this size sweep is 3.1%. This does
not establish a bound for every model, order, density or species count.

The [broader public API report](public-final/report.json) separately covers
forces and mixed descriptor models. “n2p2” here means an imported model in
AccelNet, not the upstream n2p2 executable's runtime.

| Model | Atoms | API | Original CPU (ms) | Common CPU (ms) | Ratio |
|---|---:|---|---:|---:|---:|
| aenet | 192 | structure | 39.170 | 31.320 | 0.800 |
| aenet | 192 | energy | 24.814 | 24.947 | 1.005 |
| aenet | 192 | batch | 12.565 | 9.872 | 0.786 |
| aenet | 192 | atomic | 18.548 | 10.249 | 0.553 |
| aenet | 192 | atomic-energy | 4.305 | 4.282 | 0.995 |
| n2p2 | 512 | structure | 135.493 | 132.306 | 0.976 |
| n2p2 | 512 | energy | 137.368 | 132.484 | 0.964 |
| n2p2 | 512 | batch | 6.549 | 6.519 | 0.995 |
| n2p2 | 512 | atomic | 14.133 | 7.130 | 0.504 |
| n2p2 | 512 | atomic-energy | 10.221 | 5.764 | 0.564 |
| combined | 128 | structure | 11.045 | 8.864 | 0.802 |
| combined | 128 | energy | 8.249 | 8.524 | 1.033 |
| combined | 128 | batch | 2.599 | 0.948 | 0.365 |
| multi-chebyshev | 128 | structure | 9.707 | 9.153 | 0.943 |
| multi-chebyshev | 128 | energy | 8.310 | 8.581 | 1.033 |
| multi-chebyshev | 128 | batch | 1.722 | 1.267 | 0.736 |
| mixed-components | 128 | structure | 11.556 | 10.539 | 0.912 |
| mixed-components | 128 | energy | 9.173 | 9.498 | 1.035 |
| mixed-components | 128 | batch | 3.648 | 1.648 | 0.452 |

The retained G5/common-CPU performance tests also passed; their JSON reports
are included. Unlike the paired public-API table above, those tests use the
retained explicit reference inside the current benchmark executable.


## Final GPU timings

[Full samples, commands, binary hashes and errors](gpu-final/report.json).
Four reversed rounds; energy uses at least 0.15 s per invocation and force
uses five warmed samples of at least 0.10 s per invocation. Synthetic cases
have 64 neighbors per center and order 8. The Ti/O model uses its own orders,
cutoffs and constructed neighbors. Ratios are current/e6a96d3, on the same
Blackwell GPU; they are not H100-to-Blackwell or GPU-to-CPU speedup ratios.

| API | Model | Atoms | Mode | Previous common (ms) | Current (ms) | Ratio |
|---|---|---:|---|---:|---:|---:|
| energy | chebyshev | 2048 | direct | 2.861 | 2.821 | 0.986 |
| energy | chebyshev | 2048 | moment | 1.881 | 1.937 | 1.030 |
| energy | aenet | 192 | auto | 3.468 | 2.980 | 0.859 |
| energy | aenet | 192 | direct | 45.691 | 45.588 | 0.998 |
| energy | aenet | 192 | moment | 3.454 | 2.995 | 0.867 |
| force | chebyshev | 2048 | direct | 3.857 | 3.848 | 0.998 |
| force | chebyshev | 2048 | moment | 2.475 | 2.586 | 1.045 |
| force | g5-series | 1024 | direct | 4.850 | 4.854 | 1.001 |
| force | g5-series | 1024 | moment | 2.965 | 2.972 | 1.003 |

The real Ti/O energy-only auto/moment path takes 13--14% less time. Direct is
close to the previous version. The synthetic Chebyshev moment case still costs
about 3% more for energy-only and 4.5% more with forces. G5 is essentially
unchanged in this sweep. Thus the change fixes the measured CPU energy-only
regression and improves the real-model GPU energy path, but is not a universal
GPU speedup. No speed gate was weakened to hide the remaining synthetic GPU
cost. There is no separate legacy CPU production path.

`probes.tar.gz` records earlier candidate measurements; these must not be
mixed into the final tables. In particular, the unpadded `gpu-final` and
`cpu-energy-final` reports inside that probe archive are superseded by the
standalone final reports linked above.
