# Retained CPU reference code

Source checkpoint: `gpu` **543b180** (numerical implementation **ea65290**),
AccelNet 1.0.1, methods revision 1.11. This directory is the home of the former
CPU evaluators. The production potential-inference APIs no longer fall back to
them. Do not add an automatic dispatch from the common backend into this code.

| File | Purpose |
|---|---|
| `accelnet_cpu_reference.f90` | Retained CSR evaluator and its independent workspace; called explicitly through `evaluate_batch_reference` by validation/benchmarks |
| `legacy_chebyshev_evaluation.inc` | Former Chebyshev values, derivatives, moments and force contraction |
| `legacy_behler_evaluation.inc` | Former Behler/n2p2 descriptor and force loops, including the legacy G5 moment limit |
| `legacy_lj_evaluation.inc` | Former LJ values and derivatives |
| `legacy_model_evaluation.inc` | Former composition and contraction of descriptor components |
| `legacy_network_evaluation.inc` | Standalone NN evaluator used by the retained reference and NN compatibility tests |
| `*.f90.snapshot` | Exact pre-migration structure, per-atom, C and batch API sources; archival only, not compiled |

The descriptor/NN include files remain compiled in their original modules so
standalone low-level descriptor/NN compatibility APIs and independent reference
tests remain available. Their type definitions, initialization, model loading,
neighbor construction and shared scalar formula includes stay in the normal
source directories. This avoids duplicating the data model or adding a runtime
selection between old and new potential evaluators.

All production predictor structure/file calls, public atomic energy/force/
virial calls, CSR batch calls, ænet-compatible structural-fingerprint calls,
and LAMMPS calls use the common numerical backend. The standalone descriptor
library/descriptor-generation utilities are reference-compatible lower-level
interfaces; they are not the production potential-inference dispatch.

`evaluate_batch_reference` is retained as an explicit compatibility spelling in
`accelnet_batch`; its implementation lives here. The preserved formulas remain
independent so tests can detect changes in the shared kernels. Correctness
checks also use upstream ænet/n2p2 and coordinate/strain finite differences.
