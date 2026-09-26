# LAMMPSからOpenMP target GPUを使うための設計検討

2026-09-25。コード調査に基づく提案。以下のpair style、C API、オプションは未実装。

2026-09-26追記: GPUパッケージを使う方針については、後続の
[GPUパッケージ接続案](lammps-gpu-package-plan.md) を優先する。

## 結論

既存Fortran GPUバックエンドを呼ぶ `pair_style accelnet/target` を追加する。
最初は29Aug2024 Update 4向けCMake連携を対象とし、LAMMPSのCPU近傍リストと
MPI通信を利用する。記述子、NN、力・ビリアルは既存のOpenMP targetカーネルを使う。
既存 `pair_style accelnet` とGPU無効ビルドにはGPUランタイム依存を加えない。

今回確認した既存コード:

- `interfaces/lammps/29Aug2024/ACCELNET/pair_accelnet.cpp` は `REQ_FULL` を要求し、
  LAMMPSの `ilist/numneigh/firstneigh` を使って原子ごとに
  `accelnet_atomic_energy_and_forces` を呼んでいる。
- `AccelNetPredictor/src/accelnet_batch_target.f90` のGPU入口はFortran専用。
  CSR形式の近傍、中心原子の部分集合、ghostへの力、加算型の出力を扱える。
  CからGPUバッチを呼ぶ入口はない。
- `AccelNetPredictor/src/accelnet.f90` の現行C APIは非公開のプロセス共通
  `global_model` を持つ。GPUコンテキストからこの変数を直接参照する構成は避ける。
- LAMMPS用 `cmake/ACCELNET.cmake` は現状CPUの静的ライブラリ2個とFortran
  ランタイムをリンクするだけで、GPU用ライブラリやdevice linkは扱っていない。

## 実行の流れ

```text
LAMMPS: 座標更新・ghost座標通信・必要なときに近傍リスト再構築
  ↓ local原子を中心とするfull listをCSRへ変換
C ABI → Fortran: MPI rankごとのモデル・常駐GPUワークスペース
  ↓ バッチ評価（記述子 → NN → 力・ビリアル）
LAMMPS: local + ghostの力へ加算、エネルギー・圧力へ反映
  ↓ LAMMPS既存のreverse communication
所有rankに力を集約、積分を継続
```

AccelNet独自の `build_neighbor_list` は呼ばない。
LAMMPSはskin付きVerletリストを複数ステップで再利用し、通常は空間binを使って
構築する。そのため、前回の「AccelNet近傍探索を毎回含めた全体123 ms」はLAMMPSの
予測時間ではない。LAMMPSのpacking、通信、積分、近傍再構築頻度を含めて再測定する。
[LAMMPS neighbor-list設計](https://docs.lammps.org/Developer_par_neigh.html)。

## 入口とビルド

GPUライブラリ側に `bind(C)` のcreate/evaluate/destroyを設ける。
外部にはopaque handleを公開し、内部で `predictor_model`、`target_model`、
`target_workspace` を保持する。モデル読込は既存の読込関数を再利用し、
ロード直後にGPUスナップショットを作る。pairインスタンス間で状態を共有しない。
Fortranのallocatableを持つ型をC++側に公開・コピーさせない。

C APIで固定する契約:

- atom数は有効な `nlocal + nghost`、中心は `ilist` のlocal原子。
- indexはrank内配列index。グローバルatom IDではない。Cの0始まりからFortranの
  1始まりへの変換は境界で一度だけ行う。`NEIGHMASK` を適用する。
- 変位は `x[j] - x[i]`。LAMMPSのghost座標に含まれる周期イメージを使い、
  minimum-image変換を重ねない。同じatom IDの複数周期イメージを統合しない。
- speciesはLAMMPS typeからモデルspeciesへ変換する。skin領域の近傍は保持し、
  物理cutoffで評価を切る。近傍順序を無断で変更しない（version 10の互換動作にも注意）。
- 出力energyは各中心原子、forceはlocal+ghost、virialは当該rankの中心原子分。
  Forceは他の相互作用を消さずに加算する。MPIでghostのenergyを重複計上しない。
- 空rank、zero edge、NULLと長さ0、整数幅、容量拡大を明示的に扱う。
- 不正入力・未対応モデル・GPU不在はstatus/messageを返し、LAMMPS側で終了処理する。
  現行GPU APIの `error stop` をそのままMPI境界に露出させないため、エラー処理の
  整理も実装範囲に含める。

GPUオプション有効時だけGPU pair styleとC ABIをビルドする。
現在の `AccelNet::Target` はNVHPCのoffloadリンクフラグを伝播するため、
GNU C++の最終リンクにそのまま渡しても使えるとは限らない。
初期構成はNVHPCのC++/Fortranドライバでdevice linkとMPIを確認する。
GNU C++から呼ぶ構成はNVHPCでリンクした共有ライブラリ等を別途検証する。
異なるFortranコンパイラの `.mod` や内部派生型は混在させない。

利用形の案（現時点では実行できない）:

```lammps
newton on
pair_style accelnet/target auto Ti.nn O.nn
pair_coeff * *
```

初版は起動環境でrankごとのGPU可視性を設定し、各rankで可視device 0を使う方式を
基本とする。明示device指定も用意する。1 GPUあたり1 MPI rankから検証し、rank数を
増やした場合のGPU共有は別測定する。global rankをそのままGPU番号にはしない。

LAMMPSの `-sf gpu` は既存GPUパッケージ対応styleを選ぶ機構であり、現行AccelNetを
自動的にGPU化するものではない。GPUパッケージのNewton制約はバージョン依存であり、
後続のGPUパッケージ接続案に詳細仕様とソースの確認結果を記す。
またLAMMPSの `/omp` はCPUスレッド並列なので、今回のOpenMP targetには `/target`
という明示名を提案する。
[GPUパッケージ](https://docs.lammps.org/Speed_gpu.html)、
[OPENMPパッケージ](https://docs.lammps.org/Speed_omp.html)。

## MPI・energy・virialで確認すべき点

初版は `newton pair on` を要求する。各rankがlocal原子の原子エネルギーを一度評価し、
ghostに生じた力もLAMMPSの `atom->f` に加算する。その後のLAMMPS標準の力のreverse
communicationに任せる。pair内から同じ力を独自送信して二重集約しない。
GPU評価・ホストへのコピーを完了してからLAMMPSに戻る。
[LAMMPS communication設計](https://docs.lammps.org/Developer_par_comm.html)。

GPU APIの全体virialをLAMMPSの6成分へ正しい符号・順序で反映し、`fdotr` と
二重計上しない。全体圧力はひずみ差分でも検証する。現在のGPU APIには原子別virialの
出力がないので、初版の `compute stress/atom` は明示的に未対応として扱い、
必要なら局所応力の定義と分配方法を決めたうえで別途実装する。

既存CPU pairには `ev_tally(0,0,...)` が各原子のenergyに対して呼ばれている。
29Aug2024の `Pair::ev_tally` 実装と照合すると、原子別energy要求時にはatom 0へ
集約される構造になっている。GPU版は正しい `eatom[ilist[row]]` に入れる。
原子別energyの検証はこの既存CPU pair出力だけを正解とせず、Fortran参照値でも行う。
既存CPU pair側の修正は原因を分離したテスト付き変更にする。

## 転送削減は二段階に分ける

1. 最初は現在のCSR＋変位APIをそのままCから使う。LAMMPSの近傍を毎回安全にpackし、
   GPUモデルと作業領域を再利用する。独自の近傍探索は不要だが、edge数に比例する
   packingと変位転送は残る。
2. 実測で必要なら「座標＋近傍index」の入口を追加する。GPUで変位を作り、近傍indexは
   リスト再構築・atom sorting・migration等の無効化時に更新する。座標は毎ステップ更新。
   これにより座標転送をedge数比例からlocal+ghost atom数比例に減らせる。

初版からGPU近傍探索やKokkosとのdevice pointer共有を追加する必要はない。
Kokkos側のメモリ所有権、同期、device対応付けには別の検証が必要。
大きい系は中心原子のバッチ分割でGPUメモリを制限できるが、現APIで分割ごとに
全ghost/local出力を転送するとコストが増えるため、バッチ容量と転送量を併せて測る。

## 実装順と合格条件

| 段階 | 内容 | 主な検証 |
| --- | --- | --- |
| 1 | GPU C ABIとコンテキスト寿命、エラー処理、C++リンク | Fortran直接呼出とのE/F/W一致、作成/破棄/再読込、空入力、GPU不在 |
| 2 | 29Aug2024の `accelnet/target`、1 rank/1 GPU | `run 0`、原子別energy、全force、圧力、直交/三斜晶、有限差分、短いNVE |
| 3 | 複数MPI rank・GPU、粒子移動と近傍更新 | 1/2/4 rankの同一snapshot、ghost力、sorting、空rank、周期境界越え |
| 4 | LAMMPS内の速度測定と必要な転送削減 | Pair/Neigh/Comm/Loop、packing/転送/kernel内訳、GPU常駐領域の再確保回数 |

対象はまず既存GPUバックエンドが扱える単一Chebyshev成分/元素、FP64。
LJ・n2p2/Behler・複合記述子をGPUで使うには追加カーネルが必要。
4Feb2020対応、restart、hybrid、原子別応力、r-RESPA等は実装・検証が済むまで
初版の対応機能として宣言しない。

数値比較はまず同じsnapshotでE/F/Wを比較し、現在のGPU試験の絶対/相対許容誤差を
起点に原子数・総和の規模を考慮する。長時間MDの軌跡がbit一致することを要求せず、
短時間NVEのエネルギードリフトを比較する。
性能はCPU 1 rankとの比較に加えてCPU複数rankの実用構成とも比較し、初期化と
ウォームアップを除いて繰り返し測る。8〜数百原子に加え数千〜数万原子、複数次数を
測定する。既存CPU pairの結果と性能の回帰チェックも実施する。

前回のBlackwell 512原子・次数8の計算部分約11.6倍という数字は、LAMMPS全体の
速度倍率ではない。今回の設計なら近傍探索の負担を軽くできると見込むが、改善幅は
LAMMPS上の測定前には確定しない。AMD/Intel GPUは別途実機検証が必要。
