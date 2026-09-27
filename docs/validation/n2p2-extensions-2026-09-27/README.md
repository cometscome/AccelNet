# n2p2 extensions and multi-element validation — 2026-09-27

AccelNet **1.0.1**, methods document **1.8**, based on `gpu` commit `f8c19a0`.
Added types **12, 13, 20, 21, 22, 23, 24, 25** use common Fortran CPU/GPU
numerical kernels. Existing type 9/G5 retains both direct and exact integer
moment implementations. See [the formulas and scope](../../../speedupmethods.md#19-n2p2-weightedcompact-descriptors-and-multi-element-validation-revision-18).

## Numerical coverage

* **72 synthetic cases** pass with GNU serial (OpenMP compiled out), GNU bounds
  checking, H100, and Blackwell. The independent oracle is n2p2 v2.3.0 libnnp,
  GNU C++ 11.4.0 `-O3`, Eigen 3.4.0, with OpenMP disabled. No n2p2 numerical
  sources were modified. `n2p2-source-sha256.json` identifies the local source
  snapshot, which has no Git metadata.
* All nine compact subtypes and weighted cutoff types 0--8 are covered.
  Two-element mixed models have 66 descriptors per element. Four-element
  H/C/N/O models have 110 per element, all eleven SF types and all ten unordered
  chemical pairs, tested in auto, direct and forced G5-moment modes.
* Independent force finite differences have maximum error **1.05e-11**; all
  nine strain/virial components have maximum error **6.44e-12** (fixture units).
  Isolated atoms, nonperiodic molecules and axial collinearity are included.
* Forward conversion to native AccelNet and reverse conversion to n2p2 agree
  with the unconverted model. Global energy/length normalization, atomic energy
  offsets, node normalization, descriptor sorting and subtype sorting are used.
* The ordinary serial CTest suite passes **49/49 enabled correctness tests**;
  GPU-only tests are disabled in that build. Existing common-CPU and G5-moment
  performance gates pass separately with one CPU core and no OpenMP runtime.
* H100 existing GPU regression: **78 checks**, maximum E/F/virial error
  **7.98e-16**. GPU C API initialization/reload/validation accepts the mixed
  extended native model. H100 memcheck of the four-element mixed model reports
  **0 errors**. G5 moment regression adds **185 checks per GPU**, maximum
  E/F/virial error **2.23e-16**.

## Official trained models

The unmodified files are from n2p2's official `examples/nnp-predict` tree.
The implementation authors identify the SCAN examples as polynomial-SF tests
in [the E-CAM module documentation](https://e-cam.readthedocs.io/en/latest/Classical-MD-Modules/modules/n2p2/n2p2_polynomial_symfuncs/readme.html).
Model weights and structures remain external; their SHA256 hashes are in each
`official-*.json` report. Inference outputs, not supplied DFT target values,
are compared against the independently built n2p2 library.

| Model | Elements | Atoms | SF types | Descriptors per central element |
|---|---|---:|---|---|
| [Ethylbenzene_SCAN](https://github.com/CompPhysVienna/n2p2/tree/master/examples/nnp-predict/Ethylbenzene_SCAN) | H/C | 288 | 20,22 | H=184, C=184 |
| [Anisole_SCAN](https://github.com/CompPhysVienna/n2p2/tree/master/examples/nnp-predict/Anisole_SCAN) | H/C/O | 256 | 20,22 | H=351, C=354, O=331 |
| [DMABN_SCAN](https://github.com/CompPhysVienna/n2p2/tree/master/examples/nnp-predict/DMABN_SCAN) | H/C/N | 21 | 20,22 | H=334, C=333, N=219 |
| [H2O_RPBE-D3](https://github.com/CompPhysVienna/n2p2/tree/master/examples/nnp-predict/H2O_RPBE-D3) | H/O | 1080 | 2,3 | H=27, O=30 |

All four pass on GNU serial, H100, and Blackwell, both imported and native
converted; reverse-converted files also pass n2p2. Maximum absolute total-energy
error is **1.82e-11**, maximum force-component error **1.22e-13**, in each model's
physical units. The molecular-model energy unit is not assumed to be eV.
These tests caught and fixed the redundant per-element topology overrides and
the overly large neighbor-buffer bound described in the methods document.

## Known compatibility boundary

`endpoint-diagnostic/` intentionally preserves a failing upstream comparison:
type 22 with an exponential compact core and a window centered on 0/180 degrees,
in a small periodic cell containing collinear self-image pairs. Maximum energy
difference is **5.061e-6**, force difference **6.532e-9** for this synthetic NN.
Both implementations skip exactly endpoint cosine values, but floating-point
operation order can move a pair across that branch; the window need not vanish
there. Consequently exact compatibility is **not claimed for this degenerate
case**. It is not counted among the passing cases and tolerances were not
relaxed. The ordinary periodic fixture has cell vectors longer than Rc, while
separate axial collinear and nonperiodic checks cover the defined exclusion.

## Reproduce

Build the independent serial oracle from external n2p2 and Eigen sources:

```sh
python3 AccelNetPredictor/benchmark/compile_n2p2_reference.py \
  --source /path/to/n2p2 --eigen /path/to/eigen-3.4.0 --output /tmp/n2p2-oracle
cmake -S . -B /tmp/accelnet-serial -DCMAKE_BUILD_TYPE=Release \
  -DACCELNET_BUILD_OPENMP_TARGET=ON -DACCELNET_TARGET_SERIAL=ON \
  -DN2P2_REFERENCE_BENCHMARK=/tmp/n2p2-oracle/n2p2-reference \
  -DN2P2_PREDICT_EXAMPLES=/path/to/n2p2/examples/nnp-predict \
  -DN2P2_SCALING_EXECUTABLE=
cmake --build /tmp/accelnet-serial -j 6
ctest --test-dir /tmp/accelnet-serial --output-on-failure -R 'n2p2_(extended|official)_multi_element'
```

The local snapshot's precompiled `nnp-scaling` is a macOS executable and cannot
run here; its optional older tests were disabled with the explicit empty path.
The freshly compiled libnnp oracle above checks every new SF and cutoff 0--8.
Two pre-existing benchmark build issues were also fixed: GNU 11 crashes on
`!GCC$ ivdep` preceding `do concurrent`, and a G5 output string exceeded the
standard Fortran line length. Removing the redundant hint and splitting the
string lets the full build succeed without changing numerical kernels.

For GPU tests, pass the offload build's candidate and `--backend gpu` to
`check_n2p2_extended.py` or `check_n2p2_official.py`; both require `--output`,
`--reference` and `--converter` as above. The official script additionally
requires `--examples`. Set `CUDA_VISIBLE_DEVICES` and
`OMP_TARGET_OFFLOAD=MANDATORY`. CMake's opt-in `ACCELNET_TEST_N2P2_GPU=ON`
registers the extended GPU test. Without an external oracle, the default
synthetic CTest still checks the retained CPU reference, conversion and finite
differences, but its report explicitly records `independent_n2p2=false`.

The performance script `AccelNetPredictor/benchmark/compare_n2p2.py` accepts
`--reference`, `--cpu`, `--gpu`, and `--output`. It rejects CPU/reference
executables importing known OpenMP runtime symbols, pins one CPU core, checks
results, and records all five wall-clock samples and commands. An optional
`--baseline REPORT --max-slowdown 1.10` checks both CPU and GPU regressions on
the same GPU and measurement scope. Default timing scope is fixed neighbors;
`--scope full` includes neighbor construction.

## Measured speed

These are warmed **fixed-neighbor** energy + all-force evaluations. AccelNet
also produces virial. Host/device input/output transfers are included; initial
model loading/allocation is excluded. CPU uses the public batch API, including
per-call model packing; GPU uses a prepared resident model. One Xeon Gold 6526Y
core (CPU affinity 6), GNU `-O3` without OpenMP, is compared with an H100.
Five samples of at least 0.08 s are recorded; the table shows their median.
Models are synthetic H/O fixtures with six descriptors per type per element
(type 9 has three integer-power descriptors). This is not a training benchmark
or an end-to-end LAMMPS trajectory benchmark.

| Type | Atoms | n2p2 CPU (ms) | Common CPU (ms) | H100 (ms) | H100 speedup vs n2p2 |
|---:|---:|---:|---:|---:|---:|
| 9 | 512 | 28.823 | 18.367 | 0.973 | 29.63× |
| 12 | 512 | 8.828 | 6.762 | 1.078 | 8.19× |
| 13 | 512 | 83.051 | 192.521 | 8.206 | 10.12× |
| 20 | 512 | 3.922 | 3.189 | 0.897 | 4.37× |
| 21 | 512 | 27.715 | 60.035 | 7.645 | 3.63× |
| 22 | 512 | 42.276 | 70.384 | 7.867 | 5.37× |
| 23 | 512 | 19.156 | 3.309 | 0.894 | 21.43× |
| 24 | 512 | 96.695 | 141.041 | 8.506 | 11.37× |
| 25 | 512 | 172.509 | 164.796 | 8.506 | 20.28× |

For the added types, H100 gives roughly 3.6--21× speedup over n2p2 at 512 atoms.
At 64 atoms, launch/transfer overhead can make the GPU slower than n2p2 for
compact radial/narrow functions. CPU types 13/21/22/24 are about 1.5--2.3× slower
than n2p2 here. Source inspection shows that n2p2's groups compute the pair cosine,
angle and angular-gradient geometry once before looping over members; the new
flat descriptor-owned kernel repeats them. That is a concrete optimization
opportunity, not a claim that this first extension already matches all n2p2 CPU
performance. Existing G5 moment performance remains separately gated.

The forced type-9/G5 timings use the prepared common code on both CPU and GPU:

| Atoms | Common CPU direct (ms) | Common CPU moment (ms) | H100 direct (ms) | H100 moment (ms) |
|---:|---:|---:|---:|---:|
| 64 | 3.684 | 2.161 | 1.988 | 0.785 |
| 512 | 26.177 | 18.338 | 2.307 | 0.968 |

Blackwell was measured independently on the same 512-atom fixtures:

| Type | n2p2 CPU (ms) | Common CPU (ms) | Blackwell (ms) | Blackwell speedup vs n2p2 |
|---:|---:|---:|---:|---:|
| 9 | 39.960 | 25.339 | 1.465 | 27.27× |
| 12 | 8.890 | 6.795 | 1.869 | 4.76× |
| 13 | 82.866 | 204.605 | 23.707 | 3.50× |
| 20 | 3.923 | 3.178 | 1.265 | 3.10× |
| 21 | 27.722 | 60.102 | 18.964 | 1.46× |
| 22 | 42.301 | 70.326 | 19.973 | 2.12× |
| 23 | 19.081 | 3.320 | 1.267 | 15.06× |
| 24 | 97.048 | 140.933 | 20.956 | 4.63× |
| 25 | 172.829 | 164.932 | 21.446 | 8.06× |
