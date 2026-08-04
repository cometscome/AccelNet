# AccelNet Compatibility with ænet and n2p2 Models

Assessment date: 2026-08-04<br>
AccelNet version assessed: `1.0.0`<br>
Reference implementations: ænet `2.0.4`, n2p2 `v2.3.0`

## 1. Summary

AccelNet can load existing models and calculate energies and analytical forces
within the following limits.

| Input | Loading | Energy | Forces | Main conditions |
|---|---:|---:|---:|---|
| ænet/AccelNet ASCII NN | Supported | Supported | Supported | Chebyshev, Behler2011, and the AccelNet LJ extension; ænet 2.0.4 activation codes 0--4 |
| AccelNet native binary NN | Supported | Supported | Supported | Fortran sequential-unformatted format; the reader must use a compatible record representation |
| ænet NN + explicit `.fingerprint.stp` | ASCII only | Supported | Supported | The number of NN inputs must match the number of setup descriptors |
| n2p2 2G-HDNNP directory | Conditional | Supported | Supported | SF types 2/3/9, global or per-element topology, with optional `normalize_nodes` |
| n2p2 4G-HDNNP / Q-HDNNP | Not supported | Not supported | Not supported | Charge NNs, charge equilibration, and electrostatics are not implemented |
| n2p2 weighted / compact SFs | Not supported | Not supported | Not supported | Types 12/13/20--25 are not implemented |

Here, “ænet/AccelNet ASCII” means that the files use the same basic record
layout. However, an extended ASCII file written by AccelNet may not have the
same semantics in upstream ænet 2.0.4 when it contains extended activation
functions or cutoff metadata. Being able to parse a file and being able to
evaluate it with numerically identical semantics are separate properties.

AccelNet is an inference engine. It does not reimplement the training,
dataset-processing, optimization, or MPI training functionality of ænet or
n2p2.

## 2. Scope and Assessment Method

The assessment cross-checked the following sources:

- AccelNet model loaders and inference paths:
  - [`aenet_network.f90`](../AccelNetPredictor/src/aenet_network.f90)
  - [`n2p2_network.f90`](../AccelNetPredictor/src/n2p2_network.f90)
  - [`accelnet_predictor.f90`](../AccelNetPredictor/src/accelnet_predictor.f90)
- AccelNet model converter:
  - [`model_conversion.f90`](../AccelNetModelConverter/fortran/model_conversion.f90)
- Upstream source trees available in the same workspace:
  - [ænet 2.0.4 `feedforward.f90`](../../aenet-master/src/ext/feedforward.f90)
  - [ænet 2.0.4 `sfsetup.f90`](../../aenet-master/src/sfsetup.f90)
  - [n2p2 v2.3.0 `Element.cpp`](../../n2p2-master/src/libnnp/Element.cpp)
  - [n2p2 v2.3.0 `Mode.cpp`](../../n2p2-master/src/libnnp/Mode.cpp)
  - [n2p2 v2.3.0 `NeuralNetwork.cpp`](../../n2p2-master/src/libnnp/NeuralNetwork.cpp)
- AccelNet tests and loading tests with real models distributed with n2p2
  v2.3.0.

The conclusions are based on the current implementation, not inferred only
from file-format names.

## 3. Support for ænet / AccelNet Neural Networks

### 3.1 File formats

Two formats are supported.

| Format | Detection | Notes |
|---|---|---|
| ASCII | The filename must end exactly in `.ascii` | Automatic detection by `read_aenet_network()` is case-sensitive |
| Binary | Any filename not ending in `.ascii` | Fortran sequential-unformatted records |

`load_predictor(setup_files, network_files, model)`, which accepts explicit
setup files, currently always invokes the ASCII reader. For a binary NN, use
`load_predictor_from_networks()`, which reconstructs the setup from embedded
descriptor metadata.

The binary format is not a self-describing, portable binary standard. Compiler
record markers, integer and real kinds, endianness, and related representation
details must be compatible between the writer and reader. ASCII is more
portable.

### 3.2 Network architecture

AccelNet supports the following properties of fully connected feed-forward
networks:

- Any number of hidden layers and any number of nodes per layer.
- A bias for every layer.
- A scalar atomic energy from the first output node, with atomic-environment
  descriptors as inputs.
- Analytical gradients with respect to the inputs and forces derived from
  those gradients.
- Affine descriptor transformation `(G - shift) * scale`.
- Inverse cohesive-energy scaling, atomic reference energies, and a per-atom
  energy shift.

There is no facility for interpreting multiple outputs as distinct physical
quantities. The evaluator uses only the first node in the final layer as the
atomic energy.

A multi-element model requires:

- One NN file per element.
- The number of NN files to equal the number of embedded global species.
- Files to be ordered in the same order as the embedded species.
- Each NN central-element name to match its position.
- Every central element's environment-species list to include all global
  species.

### 3.3 Embedded descriptors

The path that reconstructs a setup from `descriptor_name` in an NN file
supports the following descriptors.

| `descriptor_name` | Support | Restrictions |
|---|---|---|
| `Chebyshev` | Radial and angular Chebyshev terms with analytical derivatives | Version 0/1/10 must be specified outside the file |
| `Behler2011` | G1, G2, G3, G4, and G5 with analytical derivatives | Reordered into ænet canonical order |
| `LJ` | Per-element `sum(r^-6)` / `sum(r^-12)` with analytical derivatives | AccelNet extension; not a standard ænet 2.0.4 basis |
| Any other name | Not supported | Reports `unsupported embedded descriptor type` |

The Chebyshev version is not stored in standard NN metadata. Select version 0,
1, or 10 through `VERSION` / `CHEBYSHEV_VERSION` in `predict.in`, a Fortran API
argument, or the setter in the ænet-style API. The default is version 0.

Extended NNs can store `cutoff_type` and `cutoff_alpha` in rows 5 and 6 of
`descriptor_parameters`. Legacy NNs with fewer than six rows use a cosine
cutoff for Chebyshev and Behler descriptors and a hard cutoff for LJ
descriptors.

The explicit `.fingerprint.stp` path supports Chebyshev, LJ, Behler2011 G1--G5,
and `BASIS type=multi` combining Chebyshev and LJ. On this path, the setup file
takes precedence over the embedded `descriptor_name`, and the loader checks
that the number of setup descriptors matches the number of NN inputs.

### 3.4 ænet 2.0.4 activation functions

All native ænet 2.0.4 activation codes are implemented with the same semantics.

| Code | ænet 2.0.4 | Current AccelNet | Status as an ænet model |
|---:|---|---|---|
| 0 | linear | linear | Supported |
| 1 | tanh | tanh | Supported |
| 2 | logistic/sigmoid | logistic/sigmoid | Supported |
| 3 | modified/scaled tanh (`mtanh`) | modified/scaled tanh | Supported |
| 4 | `twist` | `mtanh + 0.1*x` | Supported |
| 5--10 | Undefined in ænet 2.0.4 | ReLU, Gaussian, cosine, reverse logistic, exponential, harmonic | AccelNet n2p2 extensions |
| 11 | Undefined in ænet 2.0.4 | softplus | AccelNet n2p2 extension |

An earlier AccelNet 0.1.0 implementation assigned n2p2 softplus to code 3,
which conflicted with ænet `mtanh`. Softplus now uses code 11. Files produced by
the earlier Fortran converter are identified as n2p2-derived from their
description and continue to interpret code 3 as softplus.

An AccelNet ASCII NN converted from n2p2 and containing codes 5--11 is not a
general-purpose NN for upstream ænet 2.0.4.

### 3.5 `predict.in` compatibility

`accelnet-predict predict.in` uses only these sections:

- `TYPES`
- `NETWORKS`
- `FILES`
- Optional `VERSION` / `CHEBYSHEV_VERSION`

It calculates energies and forces for each input file but is not a complete
replacement for ænet `predict.x`. For example, it does not reproduce structure
optimization, detailed output controls, or ænet-specific result files.

## 4. Support for n2p2 Models

### 4.1 Required files

A model directory must contain:

- `input.nn`
- One `weights.%03d.data` file per element:
  - `%03d` is the atomic number; for example, H uses `weights.001.data` and O
    uses `weights.008.data`.
- `scaling.data` when symmetry-function scaling is enabled.

Structures can be provided as AccelNet XSF files, as an `atomic_structure`
created through the Fortran API, or through
`accelnet-predict --n2p2-data [MODEL_DIR] input.data`. The `input.data` reader
supports multiple `begin`/`end` structures, nonperiodic molecules without
`lattice` records, and periodic structures with three `lattice` records.
Reference energies, charges, stored forces, and comments are not used as
inference inputs.

Element symbols must use the correct letter case and may be any of the 118
elements from H through Og. Internally, elements are reordered by atomic
number.

### 4.2 NNP generation

| n2p2 model | Status |
|---|---|
| Short-range 2G-HDNNP | Conditionally supported |
| 4G-HDNNP | Not supported |
| Q-HDNNP | Not supported |
| Charge NNs, electronegativity/hardness, charge equilibration, Ewald/electrostatic terms | Not supported |

In addition to conventional 2G files without an `nnp_type` record, the loader
accepts `2G`, `2g`, the official n2p2 v2.3.0 name `2G-HDNNP`
(case-insensitive), and the numeric form `2`. Names and numeric identifiers for
4G and Q models are rejected.

### 4.3 Symmetry functions

Upstream n2p2 v2.3.0 implements types 2, 3, 9, 12, 13, and 20--25. AccelNet can
directly load the following three types.

| n2p2 type | Meaning | Direct inference | Forces | Notes |
|---:|---|---:|---:|---|
| 2 | exponential radial | Supported | Supported | `eta`, `rshift`, `rcutoff` |
| 3 | narrow angular | Supported | Supported | An optional trailing radial shift is supported in direct inference |
| 9 | wide angular | Supported | Supported | An optional trailing radial shift is supported in direct inference |
| 12 | weighted radial | Not supported | Not supported |  |
| 13 | weighted angular | Not supported | Not supported |  |
| 20 | compact radial | Not supported | Not supported |  |
| 21 | compact narrow angular | Not supported | Not supported |  |
| 22 | compact wide angular | Not supported | Not supported |  |
| 23--25 | weighted compact | Not supported | Not supported |  |

AccelNet maps type 3 to its internal Behler G4 implementation and type 9 to G5.
Functions are sorted into n2p2 canonical order before the scaling arrays and
first-layer weights are associated with them.

### 4.4 Cutoff functions

All upstream n2p2 v2.3.0 cutoff types 0--8 are supported.

| Type | Function | Status |
|---:|---|---:|
| 0 | hard | Supported |
| 1 | cosine | Supported |
| 2 | unnormalized `tanh^3` | Supported |
| 3 | normalized `tanh^3` | Supported |
| 4 | exponential | Supported |
| 5 | polynomial 1 | Supported |
| 6 | polynomial 2 | Supported |
| 7 | polynomial 3 | Supported |
| 8 | polynomial 4 | Supported |
| 9 | fractional cutoff | Cutoff from Mori *et al.*, Appendix A; implemented as AccelNet extension type 9 and not present in upstream n2p2 v2.3.0 |

`0 <= cutoff_alpha < 1` is required. AccelNet extension type 9 additionally
requires `cutoff_alpha > 0`.

#### 4.4.1 The two meanings of “type 9”

The type number in `symfunction_short` and the number given to `cutoff_type`
belong to independent numbering systems.

**n2p2 symmetry-function type 9 = Behler G5 (wide angular)**

For central atom $i$ and neighbors $j$ and $k$, n2p2 type 9, which AccelNet
maps to G5, is

$$
G^{(9)}_i = 2^{1-\zeta}
\sum_{j<k}
\left(1+\lambda\cos\theta_{ijk}\right)^{\zeta}
\exp\!\left[-\eta\left\{(r_{ij}-r_s)^2+(r_{ik}-r_s)^2\right\}\right]
f_c(r_{ij})f_c(r_{ik}).
$$

Here, $r_{ij}=|\mathbf r_j-\mathbf r_i|$, $\theta_{ijk}$ is the angle at atom
$i$, $\lambda$ is normally $+1$ or $-1$, $\zeta$ controls angular resolution,
$\eta$ controls radial width, and $r_s$ is an optional angular radial shift.
Because this is the wide angular form, it does not include an exponential or
cutoff factor for $r_{jk}$, unlike the narrow angular form (n2p2 type 3 /
Behler G4).

**AccelNet cutoff type 9 = the fractional cutoff of Mori *et al.***

This is not a symmetry function. It is a selectable $f_c(r)$ in the equation
above. The function was introduced by H. Mori *et al.*, *Physical Review
Materials* **7**, 063605 (2023), Appendix A, Eqs. (A4)--(A6)
([doi:10.1103/PhysRevMaterials.7.063605](https://doi.org/10.1103/PhysRevMaterials.7.063605)).
The AccelNet-specific part is assigning the function extension number 9 in an
n2p2-compatible cutoff numbering scheme, not the cutoff function itself. With
$h=\alpha R_c$ and $X=(r-R_c)/h$,

$$
f_c^{(9)}(r)=
\begin{cases}
\displaystyle\frac{X^2}{1+X^2}, & 0\le r<R_c,\\[6pt]
0, & r\ge R_c,
\end{cases}
\qquad
X=\frac{r-R_c}{\alpha R_c},\quad 0<\alpha<1.
$$

Its radial derivative inside the cutoff is

$$
\frac{d f_c^{(9)}}{dr}
=\frac{2X}{\alpha R_c(1+X^2)^2},
\qquad 0\le r<R_c,
$$

and it is zero for $r\ge R_c$. Symmetry-function type 9 and cutoff type 9 may
therefore be selected simultaneously. In that case, the fractional cutoff is
substituted for the two G5 factors $f_c(r_{ij})f_c(r_{ik})$. The former type
number selects the overall descriptor form; the latter selects the truncation
function applied to each distance.

### 4.5 Network topology and activation functions

The following short-range topology settings are supported:

- `global_hidden_layers_short`
- `global_nodes_short`
- `global_activation_short`
- `element_hidden_layers_short`
- `element_nodes_short`
- `element_activation_short`
- `normalize_nodes`

The three `element_*` settings override the global hidden-layer count, node
counts, and activation functions for the named element. Zero hidden layers are
allowed, parser storage is available for up to 64 layers, and every network has
one atomic-energy output node.

With `normalize_nodes`, n2p2 evaluates each node as

$$
y_i^{(k)} = f_k\!\left(
\frac{b_i^{(k)} + \sum_{j=1}^{n_{k-1}}w_{ji}^{(k)}y_j^{(k-1)}}{n_{k-1}}
\right).
$$

AccelNet preserves these semantics by dividing every weight and bias feeding
layer $k$ by $n_{k-1}$ while loading the model. Its ordinary propagation and
analytical differentiation then produce the same energies and forces.

AccelNet also accepts an element-specific hidden-layer count. In the examined
n2p2 v2.3.0 source, a setup-loop boundary leaves one hidden-layer node count
unset when an element's depth differs from the global depth. The numerical
upstream comparison therefore exercises different per-element widths and
activations at a common depth; a separate AccelNet fixture covers the intended
hidden-layer-count override syntax.

All n2p2 v2.3.0 activation characters are implemented.

| Character | Function | Status |
|---|---|---:|
| `l` | identity/linear | Supported |
| `t` | tanh | Supported |
| `s` | logistic | Supported |
| `p` | softplus | Supported |
| `r` | ReLU | Supported |
| `g` | Gaussian | Supported |
| `c` | cosine | Supported |
| `S` | reverse logistic | Supported |
| `e` | exponential | Supported |
| `h` | harmonic (`x^2`) | Supported |

### 4.6 Symmetry-function scaling

The following modes are supported:

- No scaling.
- `scale_symmetry_functions`.
- `center_symmetry_functions`.
- Scaling plus centering.
- `scale_symmetry_functions_sigma`.
- `scale_min_short` / `scale_max_short`.

When scaling is enabled, `scaling.data` must contain a record for every
symmetry function of each relevant element. Missing or duplicate records,
out-of-range indices, zero sigma, and zero range produce errors.

### 4.7 Energy normalization and units

The following settings are supported:

- Per-element `atom_energy`.
- `mean_energy`.
- `conv_energy`.
- Corresponding inverse scaling of energies and forces.

As in upstream n2p2 v2.3.0, `mean_energy`, `conv_energy`, and `conv_length` must
be provided together. A partial set, `conv_energy = 0`, or `conv_length <= 0`
is rejected.

n2p2 multiplies coordinates, cutoffs, and shifts by `conv_length` and divides
`eta` by `conv_length^2` for evaluation in its internal units. These factors
cancel in the descriptor values of 2G types 2/3/9. AccelNet evaluates input
coordinates and model parameters directly in their original physical length
unit, producing the same descriptor values; derivatives with respect to the
physical coordinates already contain the required `conv_length` factor.
AccelNet therefore does not perform an additional coordinate conversion.

XSF coordinates and the length-dependent parameters in `input.nn` must use the
same physical length unit as was used to construct the model.

### 4.8 Ignored and rejected settings

Many training-only keywords are ignored because they are irrelevant to
inference. In contrast, an unknown key beginning with `global_`, `element_`,
`symfunction_`, or `nnp_type` is rejected as an
`unsupported n2p2 model setting`, because it could change model semantics.

Consequently, an `input.nn` that retains training settings such as epochs or
optimizer options will normally load, while new network- or descriptor-related
features are not silently ignored.

## 5. Inference Interfaces

| Feature | ænet/AccelNet NN | n2p2 directory |
|---|---:|---:|
| `accelnet-predict` CLI | Supported | Supported with `--n2p2` |
| Fortran `predictor_model` | Supported | Supported |
| Energy from XSF | Supported | Supported |
| Energy and forces from XSF | Supported | Supported |
| In-memory `atomic_structure` | Supported | Supported |
| ænet-style Fortran/C atomic API | Supported | Directly supported with `accelnet_init_n2p2` or `init` + `load_n2p2` |
| LAMMPS interface | Through NN files | Directly supported with `pair_style accelnet n2p2 DIR ELEMENT...` |

Direct loading of an n2p2 directory uses `input.nn`, the
`weights.%03d.data` files, and `scaling.data` when required. LAMMPS element
arguments are specified in LAMMPS atom-type order and automatically mapped to
the model's internal atomic-number order. Every MPI rank reads the same model
directory.

In addition to XSF, the structure-level CLI can directly read n2p2 `input.data`
with `accelnet-predict --n2p2-data [MODEL_DIR] input.data`. It supports multiple
`begin`/`end` structures, nonperiodic molecules without lattice records, and
periodic structures with three lattice records. Coordinates, elements, and the
cell are used for inference; reference energies, charges, existing forces, and
comments are not retained. General ASE formats and ænet training sets cannot be
read directly.

## 6. Model Conversion

### 6.1 n2p2 to AccelNet ASCII

For a 2G model accepted by the n2p2 loader, the Fortran converter writes:

- `<symbol>.nn.ascii` for each element.
- `networks.list`.
- Descriptor reordering from n2p2 order to ænet/AccelNet order.
- Corresponding permutations of the scaling arrays and first-layer weights.

The following restrictions apply:

1. An optional angular radial shift for n2p2 type 3/9 is stored in row 7 of
   AccelNet's extended metadata and restored when the file is reloaded.
2. n2p2 activations `p/r/g/c/S/e/h` use AccelNet extension integer codes, so
   upstream ænet 2.0.4 cannot evaluate the output NN with identical semantics.
3. Cutoff types 2--9 are also not generally equivalent to standard ænet 2.0.4
   Behler/Chebyshev evaluation.

The output is therefore extended ASCII intended for reloading by AccelNet. It
is not a general converter from arbitrary n2p2 models to native ænet 2.0.4
models.

### 6.2 AccelNet/ænet NN to n2p2

Conversion is possible only when all the following conditions hold:

- An NN is present for every element.
- `descriptor_name = Behler2011`.
- The descriptors contain only G2, G4, and G5:
  - Internal kinds 2, 4, and 5, corresponding to n2p2 types 2, 3, and 9.
- Species, atomic references, and energy scale/shift are consistent across all
  NNs.
- Every element uses the same cutoff type and alpha.
- Every element symbol appears in the H--Og table.

The following features cannot be converted:

- Chebyshev.
- LJ.
- Behler G1 and G3.
- 4G/Q and charge/electrostatic models.
- Weighted and compact symmetry functions.

The converter writes `input.nn`, `scaling.data`, and one
`weights.%03d.data` file per element. Descriptor order, scaling data, and
first-layer weights are permuted into n2p2 order. When element networks use
different architectures or activation functions, it writes global topology
defaults plus the corresponding `element_*_short` overrides. A model with an
element-specific depth may encounter the n2p2 v2.3.0 setup issue described in
Section 4.5 when loaded by that upstream version.

Upstream n2p2 has no activation equivalent to ænet `mtanh` or `twist`, so a
native ænet model containing code 3 or 4 is explicitly rejected during n2p2
conversion. Linear, tanh, and sigmoid can be converted. AccelNet extension
softplus code 11 and n2p2-derived codes 5--10 can also be converted back to
n2p2.

## 7. Validated Coverage

As of 2026-08-04, all 40 CTest cases in `build-safeopt` passed. The principal
compatibility tests cover the following behavior.

| Test | Coverage |
|---|---|
| `n2p2_model_loader` | Official `nnp_type 2G-HDNNP`, all three normalization entries, single-element H, type 2, cutoff type 7 + alpha, linear activation, atomic offset, and energy/forces |
| `n2p2_per_element_topology_and_normalize_nodes` | Per-element node counts and activations plus `normalize_nodes`, with energy and every force component checked against upstream n2p2 v2.3.0 reference values; a separate fixture covers a hidden-layer-count override |
| `n2p2_atomic_api_and_input_data` | Two-stage and one-call n2p2 loading through the Fortran atomic API; multiple `input.data` structures; molecular and periodic cells; agreement with the structure API |
| `n2p2_input_data_cli` | Energy/force output for multiple structures through `--n2p2-data` |
| `network_activation_compatibility` | Values and input derivatives for ænet `mtanh`/`twist` and n2p2 softplus |
| `network_ascii_roundtrip` | Loading of a Ti ænet/AccelNet ASCII model and gradients with respect to NN inputs |
| `network_binary_compatibility` | Energy agreement and reload for Ti/O Chebyshev ASCII and native binary models |
| `predictor_fortran_golden` | TiO2 energy and analytical forces against stored reference values |
| `predictor_original_input` | The `predict.in` path |
| `predictor_chebyshev_version*` | Chebyshev versions 0/1/10 |
| `fortran_converter_roundtrip` | Ti/O, types 2/3/9, angular shift, cutoff type 9, softplus, and a no-scaling round trip |
| `fortran_converter_per_element_roundtrip` | n2p2 → AccelNet → n2p2 → AccelNet round trip preserving per-element topology and the semantics of imported `normalize_nodes` |
| `descriptor_aenet_generate_chebyshev_match` | Atom-by-atom, coefficient-by-coefficient comparison of all radial, angular, and chemical Chebyshev channels between upstream `generate.x` and AccelNet |
| `descriptor_aenet_generate_behler_radial_match` | Behler G1/G2/G3 for every environment species and multiple parameter combinations against upstream ænet |
| `descriptor_aenet_generate_behler_angular_match` | Behler G4/G5 for every species pair, `lambda=+1/-1`, and multiple `zeta`/`eta`/cutoff values against upstream ænet |
| `network_aenet_activation_reference_match` | Values and first derivatives of ænet codes 0--4 against upstream `ff_activate` |
| `network_n2p2_activation_reference_match` | Values and first derivatives of all supported n2p2 activations against upstream `NeuralNetwork` |
| `descriptor_n2p2_cutoff_reference_match` | Values and derivatives of cutoff types 0--8 against upstream n2p2 and type 9 against its defining formula |
| `n2p2_descriptor_cutoff_{0..8}_match` | Upstream n2p2 v2.3.0 `nnp-scaling` compared with AccelNet directly reading the same `input.nn`; atom-by-atom and coefficient-by-coefficient checks for types 2 (G2), 3 (G4/narrow), and 9 (G5/wide), all species pairs, `lambda=+1/-1`, multiple `eta`/`zeta`/`rshift`/`rcutoff` values, and each cutoff type 0--8 on an asymmetric eight-atom H/O structure |

Because `function.data` from `nnp-scaling` contains ten digits after the
decimal point, the integrated comparison uses a tolerance of `7e-10`,
including relative scaling. The maximum observed scaled error over the nine
cases was at most `4.99e-11`, consistent with output rounding. This comparison
covers unscaled descriptor **values**. Cutoff values and radial derivatives are
compared separately with the upstream implementation by
`descriptor_n2p2_cutoff_reference_match`; there is not yet a
coefficient-by-coefficient comparison of the full coordinate derivatives of
types 3/9 against upstream n2p2.

With the integrated LAMMPS 29Aug2024 build, direct n2p2-directory loading was
validated with one and two MPI ranks. A four-atom fixture with LAMMPS type order
`Ti,O` and reversed internal model order `O,Ti` produced a potential energy
matching the standalone predictor value of `35.320559162555654 eV` within the
display precision at both rank counts. Existing n2p2, AccelNet, and ænet LAMMPS
smoke tests also passed with zero error in the energy and every force
component.

Loader/converter smoke tests were also run with real models distributed with
n2p2 v2.3.0.

| Bundled n2p2 model | Result | Notes |
|---|---|---|
| `H2O_RPBE-D3` | Loaded and wrote AccelNet ASCII successfully | Types 2/3 and all three normalization entries |
| `Cu2S_PBE` | Loaded and wrote AccelNet ASCII successfully | Types 2/3/9 and all three normalization entries |
| `Anisole_SCAN` | Rejected as expected | Compact types 20/22 |
| `H2O_RPBE-D3_4G` | Rejected as expected | 4G-related settings |

“Successfully” here means that parsing and conversion completed. A separate
corpus regression suite would be desirable to continuously guarantee
atom-by-atom energy and force agreement with upstream n2p2 for official
models.

## 8. Conservative Model-Selection Guidance

### ænet 2.0.4 models

A model is within the conservative compatibility subset when all of the
following conditions hold:

- It uses activation codes 0--4.
- Its descriptor is Chebyshev or Behler2011.
- It has one scalar atomic-energy NN per element.
- The correct Chebyshev version is specified.
- It is ASCII, or a native binary produced in a compatible environment.

### n2p2 v2.3.0 models

A model is within the conservative compatibility subset when all of the
following conditions hold:

- It is short-range 2G.
- It has no `nnp_type`, or uses `2G`, `2G-HDNNP`, or `2`.
- It uses only SF types 2, 3, and 9.
- It may use global topology defaults, per-element topology overrides, and
  `normalize_nodes`.
- It uses cutoff type 0--8.
- Type 3/9 with an angular shift is supported for direct inference and Fortran
  conversion.
- `mean_energy`, `conv_energy`, and `conv_length` are provided together.
- XSF coordinates and model parameters use the same physical length unit.

## 9. Recommended Next Improvements

The following order is reasonable for expanding and clarifying compatibility:

1. Add regression tests comparing energies and every force component of
   official n2p2 models with upstream `nnp-predict`.
2. Add real native ænet fixtures containing activation codes 3 and 4.
3. Add extended angular-shift metadata to the Julia converter to match the
   Fortran converter.
4. Add weighted types 12/13 to the Behler evaluator.
5. Add descriptor kernels for compact types 20--25.
