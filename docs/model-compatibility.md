# AccelNet の ænet / n2p2 モデル互換性

調査日: 2026-08-03  
対象 AccelNet: `0.1.0`, commit `02e7cb9b65bb945454a566b70a059a7e92aeac51`  
比較対象: ænet `2.0.4`, n2p2 `v2.3.0`

## 1. 結論

AccelNet は、次の範囲で既存モデルを読み込み、エネルギーと解析的な力を計算できる。

| 入力 | 読み込み | エネルギー | 力 | 主な条件 |
|---|---:|---:|---:|---|
| ænet/AccelNet ASCII NN | 対応 | 対応 | 対応 | Chebyshev、Behler2011、AccelNet拡張LJ。ænet 2.0.4の活性化コード0--4に対応 |
| AccelNet native binary NN | 対応 | 対応 | 対応 | Fortran sequential-unformatted形式。同じレコード表現を読める環境が必要 |
| ænet NN + 明示的 `.fingerprint.stp` | ASCIIのみ対応 | 対応 | 対応 | NN入力数とsetupの記述子数が一致すること |
| n2p2 2G-HDNNPディレクトリ | 条件付き対応 | 対応 | 対応 | SF type 2/3/9、共通トポロジー、`normalize_nodes`なし |
| n2p2 4G-HDNNP / Q-HDNNP | 非対応 | 非対応 | 非対応 | 電荷NN、電荷平衡、静電相互作用を実装していない |
| n2p2 weighted / compact SF | 非対応 | 非対応 | 非対応 | type 12/13/20--25は未実装 |

ここでいう「ænet/AccelNet ASCII」は同じ基本レコード配置を指す。ただし、AccelNetが書き出す拡張ASCIIは、活性化関数やcutoffメタデータによっては本家ænet 2.0.4では同じ意味で実行できない。ファイルを読めることと、数値的に同じモデルを実行できることは区別する必要がある。

AccelNetは推論器であり、ænetやn2p2の学習機能、データセット処理、最適化器、MPI学習などを再実装しているわけではない。

## 2. 調査範囲と判定方法

以下を相互に照合した。

- AccelNetのモデルローダーと推論経路
  - [`aenet_network.f90`](../AccelNetPredictor/src/aenet_network.f90)
  - [`n2p2_network.f90`](../AccelNetPredictor/src/n2p2_network.f90)
  - [`accelnet_predictor.f90`](../AccelNetPredictor/src/accelnet_predictor.f90)
- AccelNetの変換器
  - [`model_conversion.f90`](../AccelNetModelConverter/fortran/model_conversion.f90)
- 同じワークスペースにある本家側ソース
  - [ænet 2.0.4 `feedforward.f90`](../../aenet-master/src/ext/feedforward.f90)
  - [ænet 2.0.4 `sfsetup.f90`](../../aenet-master/src/sfsetup.f90)
  - [n2p2 v2.3.0 `Element.cpp`](../../n2p2-master/src/libnnp/Element.cpp)
  - [n2p2 v2.3.0 `Mode.cpp`](../../n2p2-master/src/libnnp/Mode.cpp)
  - [n2p2 v2.3.0 `NeuralNetwork.cpp`](../../n2p2-master/src/libnnp/NeuralNetwork.cpp)
- AccelNetのテストと、n2p2 v2.3.0に同梱された実モデルによるロード試験

判定は現在の実装に対するものであり、ファイル形式名だけから推測したものではない。

## 3. ænet / AccelNet NNへの対応

### 3.1 ファイル形式

対応する形式は次の2種類である。

| 形式 | 判定方法 | 備考 |
|---|---|---|
| ASCII | ファイル名末尾が厳密に `.ascii` | `read_aenet_network()`による自動判定は大文字小文字を区別する |
| binary | `.ascii`以外 | Fortran sequential-unformattedレコード |

明示的setupファイルを同時に渡す `load_predictor(setup_files, network_files, model)` は、現状では常にASCIIリーダーを呼ぶ。binary NNを使う場合は、埋め込み記述子メタデータからsetupを再構成する `load_predictor_from_networks()` を使う。

binary形式は自己記述的・可搬な標準バイナリ形式ではない。コンパイラのレコードマーカー、整数・実数kind、エンディアンなどが作成側と読み込み側で互換である必要がある。ASCIIの方が移植性は高い。

### 3.2 ネットワーク構造

AccelNetは、全結合フィードフォワードネットワークについて以下を扱う。

- 任意の隠れ層数と各層ノード数
- 各層のバイアス
- 原子環境記述子を入力とする、出力先頭ノードのスカラー原子エネルギー
- 入力に対する解析的勾配と、それを用いた力
- 記述子のaffine変換 `(G - shift) * scale`
- cohesive energyの逆スケーリング、原子参照エネルギー、原子当たりenergy shift

一方、複数出力を意味のある別物理量として扱う機能はない。評価器は最終層の第1ノードだけを原子エネルギーとして使用する。

複数元素モデルでは、次が必要である。

- 元素ごとに1個のNNファイル
- NNファイル数と埋め込みglobal species数が同じ
- ファイルの並びが埋め込みspecies順と同じ
- 各NNの中心元素名とその位置が一致
- 各中心元素のenvironment speciesにglobal speciesがすべて含まれる

### 3.3 埋め込み記述子

NNファイル内の `descriptor_name` からsetupを再構成する経路は次に対応する。

| `descriptor_name` | 対応内容 | 制約 |
|---|---|---|
| `Chebyshev` | radial/angular Chebyshev、解析微分 | version 0/1/10をファイル外から指定 |
| `Behler2011` | G1、G2、G3、G4、G5、解析微分 | ænet内部のcanonical orderingに並べ直す |
| `LJ` | 元素別 `sum(r^-6)` / `sum(r^-12)`、解析微分 | AccelNet拡張。本家ænet 2.0.4の標準basisではない |
| その他 | 非対応 | `unsupported embedded descriptor type` |

Chebyshev versionは標準NNメタデータに保存されていないため、`predict.in`の `VERSION` / `CHEBYSHEV_VERSION`、Fortran API引数、またはaenet-style APIのsetterで0、1、10のいずれかを指定する。省略値は0である。

拡張NNでは `descriptor_parameters` の5、6行目に `cutoff_type` と `cutoff_alpha` を保存できる。行数が6未満の従来NNでは、Chebyshev/Behlerはcosine、LJはhard cutoffとして扱う。

明示的 `.fingerprint.stp` を使う経路では、Chebyshev、LJ、Behler2011 G1--G5、およびChebyshev+LJの `BASIS type=multi` を利用できる。この経路では埋め込み `descriptor_name` よりsetupファイルが優先され、setupの記述子数とNN入力数の一致を検査する。

### 3.4 ænet 2.0.4の活性化関数

ænet 2.0.4のネイティブ活性化コードをすべて同じ意味で実装する。

| コード | ænet 2.0.4 | 現在のAccelNet | ænetモデルとしての判定 |
|---:|---|---|---|
| 0 | linear | linear | 対応 |
| 1 | tanh | tanh | 対応 |
| 2 | logistic/sigmoid | logistic/sigmoid | 対応 |
| 3 | modified/scaled tanh (`mtanh`) | modified/scaled tanh | 対応 |
| 4 | `twist` | `mtanh + 0.1*x` | 対応 |
| 5--10 | ænet 2.0.4では未定義 | ReLU、Gaussian、cos、reverse logistic、exp、harmonic | AccelNetのn2p2拡張 |
| 11 | ænet 2.0.4では未定義 | softplus | AccelNetのn2p2拡張 |

以前のAccelNet 0.1.0はn2p2 softplusをコード3へ割り当ててænet `mtanh`と衝突していた。現在はsoftplusをコード11へ分離した。以前のFortran変換器が生成したファイルはdescriptionからn2p2由来と判定してコード3をsoftplusへ読み替える。

n2p2から変換してコード5--11を含むAccelNet ASCII NNは、本家ænet 2.0.4用の汎用NNとはみなせない。

### 3.5 `predict.in`互換性

AccelNetの `accelnet-predict predict.in` は、次のセクションだけを実際に使用する。

- `TYPES`
- `NETWORKS`
- `FILES`
- 任意の `VERSION` / `CHEBYSHEV_VERSION`

各ファイルについてエネルギーと力を計算するが、ænet `predict.x` の全機能互換ではない。例えば、構造最適化、詳細な出力制御、ænet独自の結果ファイル生成などは再現しない。

## 4. n2p2モデルへの対応

### 4.1 必要なファイル

モデルディレクトリに次が必要である。

- `input.nn`
- 元素ごとの `weights.%03d.data`
  - `%03d` は原子番号。例: Hは `weights.001.data`、Oは `weights.008.data`
- scalingを有効にした場合は `scaling.data`

`input.data`は読まない。構造入力はAccelNet側のXSF、またはFortran APIで構築した `atomic_structure` を使う。

元素記号はHからOgまでの118元素を、正しい大文字小文字で指定する必要がある。内部では原子番号順に並べ替える。

### 4.2 NNP世代

| n2p2モデル | 判定 |
|---|---|
| 2G-HDNNP short-range model | 条件付き対応 |
| 4G-HDNNP | 非対応 |
| Q-HDNNP | 非対応 |
| 電荷NN、電気陰性度・硬度、電荷平衡、Ewald/静電項 | 非対応 |

`nnp_type`行がない通常の2Gファイルに加え、`2G`、`2g`、本家n2p2 v2.3.0の正式名 `2G-HDNNP`（大文字小文字を区別しない）、数値表記 `2` を受理する。4G/Qの名称と数値は拒否する。

### 4.3 対称関数

n2p2 v2.3.0本体はtype 2、3、9、12、13、20--25を実装しているが、AccelNetが直接読み込めるのは次の3種類である。

| n2p2 type | 内容 | 直接実行 | 力 | 備考 |
|---:|---|---:|---:|---|
| 2 | exponential radial | 対応 | 対応 | `eta`, `rshift`, `rcutoff` |
| 3 | narrow angular | 対応 | 対応 | 任意の末尾radial shiftも直接実行では対応 |
| 9 | wide angular | 対応 | 対応 | 任意の末尾radial shiftも直接実行では対応 |
| 12 | weighted radial | 非対応 | 非対応 |  |
| 13 | weighted angular | 非対応 | 非対応 |  |
| 20 | compact radial | 非対応 | 非対応 |  |
| 21 | compact narrow angular | 非対応 | 非対応 |  |
| 22 | compact wide angular | 非対応 | 非対応 |  |
| 23--25 | weighted compact | 非対応 | 非対応 |  |

type 3はAccelNet内部のBehler G4、type 9はG5へ対応付ける。関数はn2p2のcanonical orderへソートしてから、scaling配列と第1層重みに対応させる。

### 4.4 cutoff関数

本家n2p2 v2.3.0の全cutoff type 0--8に対応する。

| type | 関数 | 対応 |
|---:|---|---:|
| 0 | hard | 対応 |
| 1 | cosine | 対応 |
| 2 | unnormalized `tanh^3` | 対応 |
| 3 | normalized `tanh^3` | 対応 |
| 4 | exponential | 対応 |
| 5 | polynomial 1 | 対応 |
| 6 | polynomial 2 | 対応 |
| 7 | polynomial 3 | 対応 |
| 8 | polynomial 4 | 対応 |
| 9 | fractional cutoff | AccelNet独自拡張。n2p2 v2.3.0本体にはない |

`0 <= cutoff_alpha < 1`を要求し、AccelNet拡張type 9ではさらに `cutoff_alpha > 0`を要求する。

#### 4.4.1 「type 9」の二つの意味と定義式

`symfunction_short`のtype番号と`cutoff_type`の番号は、互いに独立した番号体系である。

**n2p2 symmetry-function type 9 = Behler G5（wide angular）**

中心原子を$i$、二つの隣接原子を$j,k$とすると、AccelNetがG5へ対応付ける
n2p2 type 9は次式である。

$$
G^{(9)}_i = 2^{1-\zeta}
\sum_{j<k}
\left(1+\lambda\cos\theta_{ijk}\right)^{\zeta}
\exp\!\left[-\eta\left\{(r_{ij}-r_s)^2+(r_{ik}-r_s)^2\right\}\right]
f_c(r_{ij})f_c(r_{ik}).
$$

ここで$r_{ij}=|\mathbf r_j-\mathbf r_i|$、$\theta_{ijk}$は$i$を頂点とする角、
$\lambda$は通常$+1$または$-1$、$\zeta$は角度分解能、$\eta$は動径方向の幅、
$r_s$は任意のangular radial shiftである。wide angularなので、narrow angular
（n2p2 type 3 / Behler G4）と異なり$r_{jk}$に対する指数因子とcutoff因子を含まない。

**AccelNet cutoff type 9 = fractional cutoff**

こちらは対称関数ではなく、上式中の$f_c(r)$として選べるAccelNet独自のcutoffである。
$h=\alpha R_c$、$X=(r-R_c)/h$と置くと、

$$
f_c^{(9)}(r)=
\begin{cases}
\displaystyle\frac{X^2}{1+X^2}, & 0\le r<R_c,\\[6pt]
0, & r\ge R_c,
\end{cases}
\qquad
X=\frac{r-R_c}{\alpha R_c},\quad 0<\alpha<1.
$$

cutoff内側での距離微分は、

$$
\frac{d f_c^{(9)}}{dr}
=\frac{2X}{\alpha R_c(1+X^2)^2},
\qquad 0\le r<R_c,
$$

であり、$r\ge R_c$では0とする。したがって、symmetry-function type 9と
cutoff type 9は同時に指定でき、その場合はG5の二つの因子
$f_c(r_{ij})f_c(r_{ik})$へfractional cutoffを代入する。前者が記述子全体の形、
後者が距離ごとの打ち切り関数を決める。

### 4.5 ネットワークトポロジーと活性化関数

対応する設定はglobal short-range topologyだけである。

- `global_hidden_layers_short`
- `global_nodes_short`
- `global_activation_short`
- 隠れ層0個も可
- 最大64層分のパーサー領域
- 出力は1ノードの原子エネルギー

次は非対応で、明示的にエラーにする。

- `element_hidden_layers_short`
- `element_nodes_short`
- `element_activation_short`
- `normalize_nodes`

n2p2 v2.3.0の全活性化文字を実装している。

| 文字 | 関数 | 対応 |
|---|---|---:|
| `l` | identity/linear | 対応 |
| `t` | tanh | 対応 |
| `s` | logistic | 対応 |
| `p` | softplus | 対応 |
| `r` | ReLU | 対応 |
| `g` | Gaussian | 対応 |
| `c` | cosine | 対応 |
| `S` | reverse logistic | 対応 |
| `e` | exponential | 対応 |
| `h` | harmonic (`x^2`) | 対応 |

### 4.6 symmetry-function scaling

次のモードに対応する。

- scalingなし
- `scale_symmetry_functions`
- `center_symmetry_functions`
- scale + center
- `scale_symmetry_functions_sigma`
- `scale_min_short` / `scale_max_short`

scalingを有効にした場合、対象元素・全対称関数の `scaling.data` 行が必要で、欠落、重複、範囲外index、ゼロsigma、ゼロrangeはエラーになる。

### 4.7 energy normalizationと単位

次に対応する。

- 元素別 `atom_energy`
- `mean_energy`
- `conv_energy`
- 上記に応じたエネルギーと力の逆スケーリング

本家n2p2 v2.3.0と同様に、`mean_energy`、`conv_energy`、`conv_length`を一組として要求する。部分的にしか存在しない場合、`conv_energy = 0`、または `conv_length <= 0` は拒否する。

n2p2本体は座標、cutoff、shiftを `conv_length` 倍し、`eta`を `conv_length^2` で割って内部単位で評価する。2G type 2/3/9ではこれらが記述子値で相殺される。AccelNetは入力座標とモデルパラメータを元の物理長さ単位のまま評価するため、同じ記述子値になり、物理座標に対する微分には必要な `conv_length` 因子が既に含まれる。このため、追加の座標変換は行わない。

XSF座標と `input.nn` の長さ依存パラメータは、モデル作成時と同じ物理長さ単位で与える必要がある。

### 4.8 読み飛ばす設定と拒否する設定

学習専用の多くのキーワードは推論に不要なため読み飛ばす。一方、未知のキーが `global_`、`element_`、`symfunction_`、`nnp_type` で始まる場合は、モデルの意味を変え得るため `unsupported n2p2 model setting` として拒否する。

この方針により、epochsやoptimizerなどの学習設定が残った `input.nn` は通常読めるが、新しいネットワーク・記述子関連機能は黙って無視しない。

## 5. 実行インターフェース

| 機能 | ænet/AccelNet NN | n2p2 directory |
|---|---:|---:|
| `accelnet-predict` CLI | 対応 | `--n2p2`で対応 |
| Fortran `predictor_model` | 対応 | 対応 |
| XSFからenergy | 対応 | 対応 |
| XSFからenergy + forces | 対応 | 対応 |
| メモリ上の `atomic_structure` | 対応 | 対応 |
| aenet-style Fortran/C atomic API | 対応 | 直接ロードAPIなし |
| LAMMPS interface | NNファイル経由 | n2p2 directoryの直接ロードなし |

n2p2モデルをaenet-style C APIや現在のLAMMPS interfaceから使う場合は、対応範囲内のモデルをAccelNet ASCIIへ変換してからロードする方法がある。ただし、次節の変換制約を必ず確認すること。

構造レベルCLIはXSFの `PRIMVEC` と `PRIMCOORD` を読む。n2p2の `input.data`、ænetの学習セット、一般的なASE形式を直接読む機能はない。

## 6. モデル変換

### 6.1 n2p2からAccelNet ASCII

Fortran変換器は、n2p2ローダーで受理できる2Gモデルについて次を出力する。

- 元素ごとの `<symbol>.nn.ascii`
- `networks.list`
- n2p2順からænet/AccelNet順への記述子並べ替え
- 並べ替えに対応したscaling配列と第1層重みの置換

ただし、次の制約がある。

1. n2p2 type 3/9の任意angular radial shiftは、AccelNet拡張メタデータの7行目へ保存して再ロードする。
2. n2p2活性化 `p/r/g/c/S/e/h` はAccelNet拡張整数コードになるため、出力NNを本家ænet 2.0.4で同じ意味で実行できない。
3. cutoff type 2--9も本家ænet 2.0.4の標準Behler/Chebyshev評価とは一般に同義ではない。

したがって出力は「AccelNetが再ロードするための拡張ASCII」であり、任意のn2p2モデルを本家ænet 2.0.4モデルへ変換する一般変換器ではない。

### 6.2 AccelNet/ænet NNからn2p2

変換できるのは以下をすべて満たす場合である。

- 全元素分のNNがある
- `descriptor_name = Behler2011`
- 記述子がG2、G4、G5だけ
  - 内部kind 2、4、5、n2p2 type 2、3、9に対応
- 全元素で同じ層構造と活性化関数
- species、atomic reference、energy scale/shiftが全NNで整合
- 全元素で同じcutoff type/alpha
- 元素記号がH--Ogの表にある

以下は変換できない。

- Chebyshev
- LJ
- Behler G1、G3
- 元素別トポロジー
- 4G/Q、charge/electrostatics
- weighted/compact symmetry functions

出力は `input.nn`、`scaling.data`、元素ごとの `weights.%03d.data` である。記述子順、scaling、第1層重みをn2p2順へ置換する。

本家n2p2にはænet `mtanh`と`twist`の同等活性化がないため、コード3/4を含むネイティブænetモデルのn2p2変換は明示的に拒否する。linear、tanh、sigmoidは変換できる。AccelNet拡張softplusコード11とn2p2由来のコード5--10もn2p2へ戻せる。

## 7. 検証済み範囲

2026-08-03時点で `build-safeopt` のCTest 36件はすべて成功した。互換性に直接関係する主な検証は次の通り。

| テスト | 実際に確認する範囲 |
|---|---|
| `n2p2_model_loader` | 正式な`nnp_type 2G-HDNNP`、正規化3項目、1元素H、type 2、cutoff type 7 + alpha、linear、atomic offset、energy/force |
| `network_activation_compatibility` | ænet `mtanh`/`twist`とn2p2 softplusの値・入力微分 |
| `network_ascii_roundtrip` | Tiのænet/AccelNet ASCIIロードとNN入力勾配 |
| `network_binary_compatibility` | Ti/O Chebyshev ASCIIとnative binaryのenergy一致、reload |
| `predictor_fortran_golden` | TiO2のenergyと解析力を保存済み参照値と比較 |
| `predictor_original_input` | `predict.in`経路 |
| `predictor_chebyshev_version*` | Chebyshev version 0/1/10 |
| `fortran_converter_roundtrip` | Ti/O、type 2/3/9、angular shift、cutoff type 9、softplus、scalingなしの往復 |
| `descriptor_aenet_generate_chebyshev_match` | 本家`generate.x`とAccelNetでChebyshevのradial/angular/chemical全チャネルを原子・係数単位で比較 |
| `descriptor_aenet_generate_behler_radial_match` | Behler G1/G2/G3、全environment species、複数パラメータを本家と比較 |
| `descriptor_aenet_generate_behler_angular_match` | Behler G4/G5、全species pair、`lambda=+1/-1`、複数`zeta`/`eta`/cutoffを本家と比較 |
| `network_aenet_activation_reference_match` | ænetコード0--4の値と1階微分を本家`ff_activate`と比較 |
| `network_n2p2_activation_reference_match` | n2p2の全対応活性化の値と1階微分を本家`NeuralNetwork`と比較 |
| `descriptor_n2p2_cutoff_reference_match` | cutoff type 0--8の値・微分を本家n2p2と比較し、type 9を定義式と比較 |
| `n2p2_descriptor_cutoff_{0..8}_match` | 本家n2p2 v2.3.0の`nnp-scaling`と、同じ`input.nn`を直接読むAccelNetを比較。H/Oの非対称8原子構造に対し、type 2（G2）、type 3（G4/narrow）、type 9（G5/wide）、全species pair、`lambda=+1/-1`、複数`eta`/`zeta`/`rshift`/`rcutoff`を、cutoff type 0--8のそれぞれで全原子・全係数比較 |

`nnp-scaling`の`function.data`は小数点以下10桁であるため、統合比較の許容差は相対スケール込みで`7e-10`とした。9ケースで観測された最大scaled errorは`4.99e-11`以下であり、出力丸め誤差の範囲内だった。この比較は未スケーリングの記述子**値**を対象とする。cutoffの値と距離微分は上記の`descriptor_n2p2_cutoff_reference_match`で本家実装と別途比較しているが、type 3/9の全座標微分を本家n2p2と係数単位で直接比較するテストはまだない。

さらに、同梱n2p2 v2.3.0実モデルに対してローダー/変換器のsmoke testを行った。

| n2p2同梱モデル | 結果 | 備考 |
|---|---|---|
| `H2O_RPBE-D3` | ロード・AccelNet ASCII出力成功 | type 2/3、正規化3項目 |
| `Cu2S_PBE` | ロード・AccelNet ASCII出力成功 | type 2/3/9、正規化3項目 |
| `Anisole_SCAN` | 期待通り拒否 | compact type 20/22 |
| `H2O_RPBE-D3_4G` | 期待通り拒否 | 4G関連設定 |

「成功」はパースと変換完了を意味する。公式モデルについて本家n2p2との全原子energy/force一致を継続的に保証するには、別途corpus回帰テストが望ましい。

## 8. 安全に利用できるモデルの目安

### ænet 2.0.4モデル

次をすべて満たすモデルが安全側である。

- 活性化コード0--4
- descriptorはChebyshevまたはBehler2011
- 元素ごとに1個のスカラー原子エネルギーNN
- Chebyshev versionを正しく指定
- ASCII、または作成環境と互換なnative binary

### n2p2 v2.3.0モデル

次をすべて満たすモデルが安全側である。

- short-range 2G
- `nnp_type`なし、`2G`、`2G-HDNNP`、または `2`
- SF type 2、3、9のみ
- global topologyのみ
- `normalize_nodes`なし
- cutoff type 0--8
- angular shift付きtype 3/9も直接実行・Fortran変換に対応
- `mean_energy`、`conv_energy`、`conv_length`は3項目を揃える
- XSFとモデルパラメータの物理長さ単位を揃える

## 9. 優先的に改善すべき点

互換性を拡大・明確化するなら、優先度は次の順が妥当である。

1. 公式n2p2モデルについて、本家 `nnp-predict` とenergy/全force成分を比較する回帰テストを追加する。
2. native ænetのコード3/4を含む実ファイルfixtureを追加する。
3. Julia版変換器にもangular shift用の拡張メタデータを実装し、Fortran版と揃える。
4. weighted type 12/13をBehler評価器へ追加する。
5. compact type 20--25の新しい記述子kernelを追加する。
