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

両実装ともBehler–Parrinello 2G types 2, 3, 9、cosine cutoff、affine scaling、
活性化関数、energy normalization、atomic reference energiesに対応します。
n2p2とænetの対称関数順序に合わせ、第一層の重みとscaling配列も並べ替えます。

Chebyshev, LJ, 4G/charge, weighted/compact or shifted-angular symmetry
functions, non-cosine cutoffs, `normalize_nodes`, and per-element n2p2
topologies are rejected rather than approximated.

Run tests:

```sh
julia --project=. test/runtests.jl
```
