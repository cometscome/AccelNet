# CPU/GPU implementation status

Assessment: 2026-09-27; AccelNet 1.0.1, methods revision 1.12; migration baseline `gpu`
checkpoint `543b180`. This describes evaluation of models **inside AccelNet**,
not GPU support in the upstream ænet or n2p2 programs.

## Descriptor coverage

“Common” means the production CPU batch and OpenMP-target GPU paths compile
one maintained numerical source. The CPU batch instance strips OpenMP
directives. A dash in the moment column means no separate moment algorithm is
implemented; radial descriptors already use a single neighbor sum.

| Descriptor / model family | Common CPU/GPU direct | Common moment | Auto in the common batch path |
|---|---|---|---|
| ænet Chebyshev, versions 0/1/10 | Yes | Yes | Operation-count heuristic using angular neighbors, angular order and moment count |
| ænet Behler G1/G2/G3; n2p2 type 2 maps to G2 | Yes | — (radial sum) | Direct radial evaluation |
| ænet Behler G4 / n2p2 type 3 | Yes | — | Direct, grouped pair values and Jacobians |
| ænet Behler G5 / n2p2 type 9 | Yes | Exact integer powers 1--16 | Integer powers 1--10 and at least 16 angular neighbors; otherwise direct |
| n2p2 type 12 (weighted radial) | Yes | — (radial sum) | Direct |
| n2p2 type 13 (weighted narrow angular) | Yes | — | Direct |
| n2p2 types 20/23 (compact / weighted compact radial) | Yes | — (radial sum) | Direct |
| n2p2 types 21/24 (compact / weighted compact narrow angular) | Yes | — | Direct |
| n2p2 types 22/25 (compact / weighted compact wide angular) | Yes | — | Direct; no approximate angular moment expansion |
| AccelNet LJ extension | Yes | — (radial sum) | Direct |

Chebyshev auto selects moment when `N_angular * (angular_order + 1) >= M`,
where `M` is the number of angular moments. This is a heuristic, not a timed
per-device optimizer. Chebyshev modes are 0=auto, 1=direct, 2=moment.

G5 modes are 0=auto, 1=direct, 2=moment with the 16-neighbor threshold, and
3=forced eligible moments. Explicit modes 2/3 permit orders 11--16; auto keeps
them direct. Fractional, near-integer and orders above 16 remain direct in all
modes. Neighbors are counted within each component's maximum angular cutoff.
LAMMPS `g5 moment` selects mode 3. Chebyshev and G5 selectors are independent;
LAMMPS defaults both to auto. Neither auto rule chooses according to measured
CPU/GPU timing or thread count.

## What “ænet supported” includes

* Loading ænet/AccelNet ASCII networks and compatible native binary networks;
  native binary compatibility depends on Fortran record representation.
* Chebyshev and Behler2011 descriptors, affine input scaling, energy
  normalization/reference energies and native activation codes 0--4.
* Energy, analytic forces and virial in the common batch and LAMMPS paths.
* Multi-element networks, including independent element-specific topologies.
* Chebyshev versions 0/1/10. Standard ænet NN metadata does not encode this
  version: the default is 0; select another convention explicitly when needed.
* CPU reference tests against ænet descriptor/activation implementations, real
  Ti/O golden-model comparisons, and historical H100/Blackwell and LAMMPS
  comparisons of the same model. This does not establish every ænet model or
  binary compiler format as tested.

Multiple Chebyshev components and Chebyshev/LJ/Behler mixtures inside an element
now use the common component pipeline, preserving original descriptor offsets
and one NN per element. Different element-specific component lists are supported.
No production inference path silently falls back to the former CPU evaluator.

## Execution paths and remaining separate code

| Entry point | Numerical implementation | CPU threading / GPU |
|---|---|---|
| Default CPU CSR batch | Serial compilation of common kernels for supported configurations | OpenMP compiled out |
| LAMMPS `pair_style accelnet` | Common serial batch for supported models | Ordinary CPU path remains serial within each MPI rank |
| Explicit target API with `use_host=.true.` | Common kernels | OpenMP host threads, measured at 1/2/4/8 |
| Target API / LAMMPS `accelnet/gpu` | Common kernels | GPU offload; H100 and Blackwell validated |
| Structure/file, per-atom Fortran/C, and ænet-compatible SFB APIs | Common serial kernels; atomic environments become one CSR row | CPU, OpenMP compiled out |
| Explicit `evaluate_batch_reference` and standalone low-level descriptor/NN utilities | Former CPU evaluator in `legacy/cpu-reference/` | Independent reference / lower-level compatibility |

The production potential-inference paths now share their kernels. The repository
retains an explicit independent reference in [legacy/cpu-reference](../legacy/cpu-reference/README.md). That older G5
moment implementation retains its order-10 limit. “Atomic API” means per-atom
API and is unrelated to an OpenMP atomic instruction.

Current host force accumulation avoids fine-grained atomic updates. GPU
physical-atom force scatter still uses atomic additions; virial uses reduction.
AMD/Intel GPU execution has not been validated by the archived measurements.

## Validation and performance scope

The [revision 1.12 migration report](validation/unified-api-2026-09-27/README.md)
records public API, composite descriptor, lifecycle and performance checks.

* [Latest atomic-removal validation](validation/atomic-reduction-2026-09-27/README.md):
  CPU functional/reference tests, host 2/8-thread checks, H100/Blackwell
  independent n2p2 checks, LAMMPS comparisons and paired CPU/GPU timings.
* [ænet Ti/O GPU validation](validation/openmp-target-2026-09-25/README.md)
  and [later Chebyshev/LAMMPS tuning](validation/moment-horner-2026-09-26/README.md).
  These are measurements at their documented revisions, not new reruns at
  checkpoint ea65290. The latest CPU functional suite includes the Ti/O batch
  golden test and ænet activation reference test.
* [All eleven n2p2 types and trained multi-element models](validation/n2p2-extensions-2026-09-27/README.md).
* [Exact high-order G5 scope and direct/moment comparison](validation/high-g5-moments-2026-09-27/README.md).

The latest direct/moment timing table uses **forced methods**, not auto.
High-order G5 moment can lose on CPU even when it wins on GPU. No universal
moment speed advantage is claimed. The strict n2p2 type-21 CPU parity test
currently misses its 10% limit narrowly (10.13% slower); numerical checks pass.
For compact angular windows, the exact-collinearity/upstream endpoint limitation
is documented in [model compatibility](model-compatibility.md#10-weightedcompact-validation-methods-revision-18).

Not implemented: potential training, n2p2 4G/Q charge/electrostatic models, or
approximate moments for compact angular functions. AccelNet-specific extended
ASCII metadata is not guaranteed to be interpretable by upstream ænet.
