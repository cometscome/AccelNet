# AccelNetModelConverter

Julia版とFortran版を備えた、AccelNet/ænetとn2p2 2G-HDNNPモデルの
相互変換ツールです。

## Fortran版

Fortran版は `AccelNetPredictor` を直接リンクし、Predictorに実装された
n2p2ローダー、AccelNet ASCII/binaryローダー、およびASCIIライターを利用します。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure

# n2p2 -> AccelNet ASCII
build/accelnet-model-converter-fortran n2p2-to-accelnet \
  /path/to/n2p2-model /path/to/accelnet-output

# AccelNet ASCIIまたはnative binary -> n2p2
build/accelnet-model-converter-fortran accelnet-to-n2p2 \
  /path/to/n2p2-output Ti.nn O.nn
```

標準では隣接する `../AccelNetPredictor` を同一コンパイラでビルドします。
別の場所にある場合は
`-DACCELNET_PREDICTOR_SOURCE_DIR=/path/to/AccelNetPredictor` を指定します。

## Julia版

n2p2 → AccelNet ASCII networks:

```sh
julia --project=. bin/accelnet-model-converter.jl n2p2-to-accelnet \
  /path/to/n2p2-model /path/to/accelnet-output
```

AccelNet native binary or ASCII networks → n2p2:

```sh
julia --project=. bin/accelnet-model-converter.jl accelnet-to-n2p2 \
  /path/to/n2p2-output Ti.nn O.nn
```

両実装ともBehler–Parrinello 2G types 2, 3, 9、n2p2 cutoff types 0--8、
および[Mori *et al.*のfractional cutoff](https://doi.org/10.1103/PhysRevMaterials.7.063605)を
AccelNet拡張番号type 9として扱う形式に対応します。
`0 <= cutoff_alpha < 1`（type 9では`0 < cutoff_alpha < 1`）、affine scaling、
活性化関数、energy normalization、atomic reference energiesに対応します。
n2p2とænetの対称関数順序に合わせ、第一層の重みとscaling配列も並べ替えます。

Fortran版はtype 3/9のangular radial shiftをAccelNet拡張メタデータの
7行目に保存して往復できます。Julia版は現時点では非ゼロangular shiftを
明示的に拒否します。

Chebyshev、LJ、4G/charge、weighted/compact symmetry functionsは、近似せず
明示的に拒否します。Fortran版はper-element n2p2 topologyを双方向で保持します。
n2p2からの読込み時には、`normalize_nodes`を各層のweightとbiasへ厳密に
取り込むため、変換後のモデルにこのキーワードは不要です。
Julia版は現時点ではper-element topologyと`normalize_nodes`に未対応です。

Run tests:

```sh
julia --project=. test/runtests.jl
```
