# LAMMPS・n2p2・CabanaMDのGPU実装調査

調査日: 2026-09-26。公式文書・著者論文・公開ソースを確認した。今回の調査では
CabanaMDのビルドや速度測定、AccelNetの実装変更は行っていない。

## LAMMPSがGPUを使う仕組み

LAMMPSはpair style等に用意されたGPU用実装を選択して実行する。
GPU実行オプションがCPU専用pair styleを自動変換するわけではない。

| 方式 | 計算とデータの置き方 | 選択方法の例 |
| --- | --- | --- |
| GPUパッケージ | 力計算等をGPUへ移し、座標・力は毎ステップCPU/GPU間で転送。近傍リスト構築はCPU/GPUを選べる | 対応pairの `/gpu`、`-sf gpu -pk gpu 1` |
| KOKKOSパッケージ | pair・近傍処理・fix等をKokkosで実装。対応する処理で揃えればデータを複数ステップGPUに保持できる | GPUバックエンドでビルドし、`-k on g 1 -sf kk` |
| 外部GPUライブラリを呼ぶpair | LAMMPSのpairからライブラリ側のGPUカーネルを呼ぶ。転送や同期の責任分担を個別に設計する | AccelNetに提案した `/target` はこの方式、未実装 |

GPUパッケージはCUDA/HIP/OpenCL等のGPUライブラリ経由、KOKKOSはC++の並列処理と
メモリ抽象化を用いる。LAMMPSのOPENMPパッケージの `/omp` はCPUスレッド並列で、
Fortran OpenMP `target` によるGPU offloadとは異なる。
KOKKOSでもCPU専用fix/computeや出力を使えば転送が必要で、「常に転送ゼロ」ではない。

出典: [GPU](https://docs.lammps.org/Speed_gpu.html)、
[KOKKOS](https://docs.lammps.org/Speed_kokkos.html)、
[OPENMP](https://docs.lammps.org/Speed_omp.html)。

## n2p2の通常のLAMMPS連携とCabanaMD連携

通常の `pair_style hdnnp` はML-HDNNPパッケージからn2p2を呼ぶ。
調査時点のLAMMPS developツリー `ec02ed0f8b0b347d3693d96892a87499b24f40f9` には
`src/ML-HDNNP/pair_hdnnp.cpp/.h` があるが、`hdnnp/kk`・`hdnnp/gpu` の実装は
見つからなかった。公式pair文書にもaccelerated variantは示されていない。
したがって標準のn2p2連携に `-sf kk` 等を付けるだけでNNP全体がGPU化するとは扱えない。
外部forkや別のNNP実装の存在まで否定するものではない。
[hdnnp公式説明](https://docs.lammps.org/pair_hdnnp.html)、
[確認したソース](https://github.com/lammps/lammps/tree/ec02ed0f8b0b347d3693d96892a87499b24f40f9/src/ML-HDNNP)。

n2p2公式のGPU実装例はCabanaMD向けインターフェース。
CabanaMDはLAMMPSとは別のMD研究用アプリケーションで、CabanaとKokkosを使う。
n2p2のモデル読込等を再利用しつつ、計算部分をGPU対応カーネルに置き換えている。
LAMMPS風の入力を使うことと、LAMMPS内部で動くpair styleであることは区別する。
[n2p2 CabanaMD公式説明](https://compphysvienna.github.io/n2p2/interfaces/if_cabanamd.html)、
[CabanaMD](https://github.com/ECP-copa/CabanaMD)。

## 実装で確認できたGPU計算の流れ

確認した版:

- n2p2: `29b9c9f10b1ac9ea45631b3d556a129ca3b701c2`
- CabanaMD: `b79c780c4edcbba5ea3a2ac1dff4ffdf1ea21956`

CabanaMDの `ForceNNP::compute()` は次を連続して呼ぶ。

```text
calculateSymmetryFunctionGroups  記述子G
calculateAtomicNeuralNetworks   原子エネルギーEとdE/dG
calculateForces                 dE/dGと記述子の座標微分から力を集約
```

出典: [ForceNNPの呼出元](https://github.com/ECP-copa/CabanaMD/blob/b79c780c4edcbba5ea3a2ac1dff4ffdf1ea21956/src/force_types/force_nnp_cabana_neigh_impl.h)、
[n2p2のGPU計算本体](https://github.com/CompPhysVienna/n2p2/blob/29b9c9f10b1ac9ea45631b3d556a129ca3b701c2/src/libnnpif/CabanaMD/ModeCabana_impl.h)。

| 対象 | 実装 |
| --- | --- |
| 動径記述子 | `Cabana::neighbor_parallel_for`、中心iと近傍jを処理 |
| 角度記述子 | 同じAPIの `SecondNeighborsTag`、i-j-kの近傍三つ組を処理 |
| 正規化・NN | `Kokkos::parallel_for`。NNの順伝播と入力に対する微分を計算 |
| 力 | 動径・角度用カーネルで幾何と座標微分を計算し、dE/dGと縮約。atomicなforce sliceへ加算 |
| データ | 粒子・G・dE/dG等はCabana AoSoA、モデル係数等はKokkos View |
| 近傍リスト | 対象deviceのmemory spaceを持つ `Cabana::VerletList` |
| NVE積分 | 座標・速度の更新も `Kokkos::parallel_for` |

[近傍リスト実装](https://github.com/ECP-copa/CabanaMD/blob/b79c780c4edcbba5ea3a2ac1dff4ffdf1ea21956/src/neighbor_types/neighbor_verlet.h)、
[NVE実装](https://github.com/ECP-copa/CabanaMD/blob/b79c780c4edcbba5ea3a2ac1dff4ffdf1ea21956/src/integrator_nve_impl.h)、
[G・dE/dG・Eのデータ配置](https://github.com/ECP-copa/CabanaMD/blob/b79c780c4edcbba5ea3a2ac1dff4ffdf1ea21956/src/system_types/system_nnp_3aosoa.h)。

GPUバックエンドでビルドすれば、近傍・記述子・NN・力・積分でdevice上のデータを
使える構造になっている。I/O・初期化・通信に関する全転送がなくなるという意味ではない。

## メモリと並列化で参考になる点

1. 記述子GとdE/dGは保持するが、全原子×近傍×記述子の座標微分を巨大配列として
   保存しない。力カーネルで必要な微分を再計算し、その場で縮約する。
2. 動径と角度を分け、原子だけに並列化する場合と近傍にも並列化する場合を選べる。
   並列度を増やすとatomic競合も増えるので、最大並列度が常に最速とは限らない。
3. AoSoAは粒子データを小さなまとまりに分け、その中で同じ成分を連続配置するもの。
   CPU/GPUに応じた配置を試せる。単に同じループをGPUで実行するだけではない。

これらは論文の節3.2–3.3、4.3にも説明されている。2020年arXiv版では、GPU上の大きい系で
近傍方向の追加並列化がatomic競合により遅くなる例も報告されている。
CPUの少スレッド・小規模条件では既存n2p2の方が速い例もある。
「移植したコードに全面的に置き換えればCPUも必ず速くなる」という結果ではない。
[著者論文・2020年arXiv版](https://arxiv.org/pdf/2002.00054)。

同研究の出版版はCPC 270 (2022), 108156。
arXiv初稿と出版版では報告する規模や性能値が異なるため数値を混ぜない。
古いV100/POWER9等の測定をH100や現在のAccelNetの速度予測に使わない。
[出版版の書誌・概要](https://doi.org/10.1016/j.cpc.2021.108156)。

## そのまま採用しない方がよい部分

確認した `ModeCabana_impl.h` の主要な記述子分岐はtype 2/3、NN活性化分岐は
線形/tanhであり、n2p2の全モデル形式への互換性を仮定しない。
NN作業用のViewは評価関数内で作られるので、このコード全体が「割当も毎ステップ完全ゼロ」
というわけではない。またCabanaMDのenergy集約には原子energy offsetや単位変換について
一般化が必要というTODOが残っている。性能設計の例として学び、既存モデル・単位・元素数の
仕様を独立に検証する。

## AccelNetへの判断

既存のFortran + OpenMP target路線は維持できる。CabanaMDから取り入れるべき点は
GPUへ載せる計算範囲、巨大Jacobianを避けること、データ配置と並列度の実測。
現在のAccelNetは記述子→NN→力までGPU化済みで、巨大Jacobianを作らず、モデルと
作業領域を常駐化している。CabanaMD例とは角度記述子の数式も異なり、GPU向けの設計原則を
参考にするのであって、そのカーネルを直接置き換えるわけではない。

最初のLAMMPS連携はCPUの近傍リストとMPI・積分を利用し、Fortran GPUバッチAPIを呼ぶ。
この構成は転送の点でLAMMPS GPUパッケージに近い。
次に座標と近傍indexからGPUで変位を作り、近傍indexの転送を更新時だけにする。
MD全体のGPU常駐化が必要になった場合は、KOKKOSとのdeviceメモリ共有・同期・MPIを
含む追加設計を行う。Fortran GPUカーネルをC ABIで呼べることだけではこの共有は成立しない。

既存CPU経路とGPU経路は維持し、CPU少数原子・CPU複数rank・GPUのそれぞれで速度回帰を
測る。この方針はCabanaMD論文のCPU/GPU間の性能差とも整合する。
具体的な接続手順は [LAMMPS実装案](lammps-openmp-target-plan.md) を参照。
