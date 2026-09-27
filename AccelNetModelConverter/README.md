# AccelNetModelConverter

Convert supported 2G-HDNNP models between n2p2 directories and AccelNet ASCII
networks. The **Fortran converter** uses the predictor's current loaders and
writers. The **Julia converter** has a narrower compatibility scope; see below.
Neither converter approximates unsupported descriptors.

For inference, AccelNet can read n2p2 directories directly through its CPU
CLI/APIs and LAMMPS CPU pair style. Conversion is needed for the current
LAMMPS GPU adapter, which accepts embedded network files. See the
[root README](../README.md) and [GPU guide](../docs/lammps-gpu.md).

## Fortran converter

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel --target accelnet-model-converter-fortran

# n2p2 directory -> one AccelNet ASCII network per element
build/bin/accelnet-model-converter-fortran n2p2-to-accelnet \
  /path/to/n2p2-model /path/to/accelnet-output

# AccelNet ASCII or compatible native binary -> n2p2 directory
build/bin/accelnet-model-converter-fortran accelnet-to-n2p2 \
  /path/to/n2p2-output Ti.nn.ascii O.nn.ascii
```

The first command also writes `networks.list`. Supply one network for every
embedded element when converting back. Species metadata, descriptor ordering,
first-layer weights and scaling arrays are transformed together to preserve
inference. Existing files with matching output names are replaced; use a new
output directory to retain the originals.

To build this component alone, run `cmake -S . -B build` from
`AccelNetModelConverter/`, followed by `cmake --build build --parallel`.
That build places the executable at `build/accelnet-model-converter-fortran`
and builds the sibling predictor with the same compiler. Override
`ACCELNET_PREDICTOR_SOURCE_DIR` if the predictor is located elsewhere.

## Conversion coverage

| Feature | Fortran | Julia |
|---|---|---|
| n2p2 types 2/3/9 (G2/G4/G5) | Supported | Supported |
| n2p2 weighted/compact types 12/13/20--25 | Supported through AccelNet extended metadata | Not supported |
| Nonzero angular radial shift for types 3/9 | Supported | Rejected |
| Per-element network topology | Supported in both directions | Not supported |
| n2p2 `normalize_nodes` | Folded into weights/biases on import | Not supported |
| Affine descriptor scaling, energy normalization, atomic references | Supported | Supported |
| Compatible native binary input | Supported | Supported |
| Chebyshev/LJ -> n2p2 | Not supported | Not supported |
| n2p2 4G/Q charge/electrostatic models | Not supported | Not supported |

Both converters handle n2p2 cutoff types 0--8 and the AccelNet fractional-cutoff
extension numbered 9. The normal cutoff-alpha interval is `0 <= alpha < 1`;
type 9 requires positive alpha (the Julia converter additionally requires
`alpha < 1`). Type 9 is an AccelNet extension, not an upstream n2p2 cutoff type.
Activation functions must have an equivalent in the destination format.

Fortran conversion preserves per-element topology and folds `normalize_nodes`
exactly into each layer's weights and biases. The output therefore does not
need that keyword. Descriptor data for weighted/compact functions are retained
in AccelNet's `n2p2_extended` representation; this is not a promise that upstream
ænet can read those extended files. Native binary input remains dependent on
Fortran record representation; use ASCII for portable interchange.

Reverse conversion requires consistent species, energy/reference and cutoff
metadata across element networks. See [model compatibility](../docs/model-compatibility.md)
for accepted parameters, ordering and format restrictions.

## Julia converter

Requires Julia 1.10 or compatible. Run from `AccelNetModelConverter/`:

```sh
julia --project=. bin/accelnet-model-converter.jl n2p2-to-accelnet \
  /path/to/n2p2-model /path/to/accelnet-output

julia --project=. bin/accelnet-model-converter.jl accelnet-to-n2p2 \
  /path/to/n2p2-output Ti.nn O.nn

julia --project=. test/runtests.jl
```

Use the Fortran converter for weighted/compact models, nonzero angular shifts,
per-element topology or `normalize_nodes`.

## Validation

From a complete repository-root build with tests enabled:

```sh
cmake --build build --parallel
ctest --test-dir build -R 'fortran_converter|n2p2_extended_multi_element' \
  -LE gpu --output-on-failure
```

The suite checks round trips, per-element networks, and mixed multi-element
weighted/compact models. An independent n2p2 executable can be enabled with
`N2P2_REFERENCE_BENCHMARK`; it is optional, not a converter dependency. The
[extension validation report](../docs/validation/n2p2-extensions-2026-09-27/README.md)
records comparisons against upstream n2p2 and trained model examples.
