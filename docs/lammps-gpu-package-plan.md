# LAMMPS GPUパッケージへのAccelNet接続案

2026-09-26。ユーザーのGPUパッケージ利用方針を受けた設計検討。
`accelnet/gpu`、bridge、device入力APIは実装・検証済み。操作方法と対応範囲は
[実装ドキュメント](lammps-gpu.md)を参照。以下は実装前の設計記録。

## 推奨する構成

`PKG_GPU` と `package gpu` を使う `pair_style accelnet/gpu` を追加する。
GPUの割当・設定と近傍リストはLAMMPS GPUパッケージ側、記述子・NN・力計算は
既存Fortran + OpenMP target側に分担する。まずCUDA/FP64で検証する。
これはソース確認に基づく設計案であり、2つのランタイムの相互運用は最初の検証項目。

```text
package gpu / FixGPU / LAMMPS_AL::Device
  ├─ GPUの選択・MPI rankへの割当・近傍設定
  └─ PairAccelNetGPU → C++ bridge → Fortran bind(C)
                                  └─ 常駐モデル・作業領域
                                     記述子 → NN → 力・ビリアル
```

既存CPUの `pair_style accelnet` とCPU専用ビルドは維持する。
前案の独立した `/target` より、今回はこちらの `/gpu` 接続を優先する。
GPUパッケージのGPU_APIはCUDA/HIP/OpenCLで、Fortran OpenMP targetを選ぶ既存スイッチは
ない。従ってC++側の接続層とビルド変更が必要。
[公式ビルド説明](https://docs.lammps.org/Build_extras.html#gpu-package)。

## 実装方式の比較

| 方式 | 再利用できるもの | 課題 | 判断 |
| --- | --- | --- | --- |
| GPUパッケージとFortranのbridge | 既存GPUカーネル、数値テスト、モーメント法、CPU経路 | device選択、CUDA context、同期、device配列の受渡し | Fortranを維持する今回の第一候補 |
| `lib/gpu/lal_accelnet.*` に計算全体を移植 | モデル形式、アルゴリズム、参照テスト | C++/CUDA/HIP/OpenCL向けカーネルの追加実装・二重保守 | 将来GPUパッケージ各バックエンドへの深い統合が必要な場合に再検討 |

bridgeの初期版でHIP/OpenCLも自動対応するとは宣言しない。GPUパッケージ側の
対応GPUとFortran offloadコンパイラの対応GPUは別条件であり、メモリ相互運用も別に検証する。

## コード調査で分かった制約

手元の `lammps-29Aug2024` は `29 Aug 2024 Update 4` と表示するソースツリー。
加えて公開developの `ec02ed0f8b0b347d3693d96892a87499b24f40f9` を確認した。
実装時は対象版を固定し、その版のAPIで試験する。

### GPUの設定・メモリ所有権

- `/gpu` pairは `Suffix::GPU` と `GPU_EXTRA::gpu_ready` 等でpackageとの接続を行う。
- GPUパッケージは `LAMMPS_AL::Device` にdevice割当やneigh設定を保持する。
  bridgeはそこから実際の割当を取得し、OpenMP側と同じ物理GPUかを確認する。
  GPU番号を別々に自動選択して偶然の一致に依存しない。
- CUDA側のUCLにはprimary contextを使う経路がある。OpenMP側とのcontext同一性、
  device切替、stream同期、初期化と破棄の順序を小さな相互運用テストで確認する。
  primary contextが使えることだけで互換性が証明されたとは扱わない。
- lib/gpuの所有メモリを借りるdevice入力APIではFortran側から解放しない。
  library側の再確保や近傍再構築でポインタが変わったら更新する。

### ghost力を返すタイミング

GPUパッケージの通常のanswer queueは `FixGPU::post_force()` で結果を回収する。
Verletループでは通常のforce reverse communicationがその前に実行される。
現在のAccelNet方式はlocal中心原子からghost原子にも力を加えるので、通常のqueueに
そのまま積むとghost力の通信に間に合わない。

初版は `PairAccelNetGPU::compute()` 内で同期・download・力への加算を完了し、
LAMMPSの通常reverse communicationへ渡す。AccelNetの同じ結果をanswer queueにも
登録して二重に加算しない。GPUとCPUの非同期計算の重ね合わせは初版の目標に含めない。
将来非同期化するなら、pre_reverseで完了を待つ仕組みを設計する。

出典: [FixGPU](https://github.com/lammps/lammps/blob/ec02ed0f8b0b347d3693d96892a87499b24f40f9/src/GPU/fix_gpu.cpp)、
[Verlet実行順](https://github.com/lammps/lammps/blob/ec02ed0f8b0b347d3693d96892a87499b24f40f9/src/verlet.cpp)、
[GPU Device管理](https://github.com/lammps/lammps/blob/ec02ed0f8b0b347d3693d96892a87499b24f40f9/lib/gpu/lal_device.cpp)。

### newtonの扱い

「GPUパッケージは必ずnewton off」とはしない。2026-09-02のpackage文書では
newton onの制約解除とCPU/GPU particle splitの廃止が明記されている。
一方Speed_gpuの説明にはoff必須という記述が残るため、対象版の詳細仕様とソースで判断する。
手元の2024版 `FixGPU` はnewton onかつsplit < 1の場合を拒否する実装で、split 1なら
このチェックは通る。初版はnewton on、旧版ではsplit 1を明示する。
これは独自AccelNet pairをその条件に実装する方針であり、任意の既存GPU pairについて
newton onで正しいことを保証するものではない。
[package gpuの詳細・変更履歴](https://docs.lammps.org/package.html#gpu-package-settings)。

`tersoff/gpu`・`sw/gpu` は多体相互作用の参考になるが、ghostの近傍まで評価する構成で
通信cutoffを概ね2Rc+skinへ拡大する部分がある。AccelNetのlocal中心＋ghostへのscatter
方式にそれを無条件にコピーしない。必要な通信範囲は計算方式とcutoffから決める。

## 推奨する実装順

### 0. ビルドとランタイムの相互運用

`PKG_GPU=ON`, `GPU_API=cuda`, `GPU_PREC=double` を基本に、既存GPU package初期化と
Fortran OpenMP targetの小さな呼出を同一processで動かす。GPU割当、context、終了処理を
確認する。GPU packageの既定はmixedなのでdoubleを明示する。
NVHPC FortranとC++/CUDAの最終リンクを確認し、GPUフラグをCPU用ターゲットへ伝播させない。
H100用にはsm_90、手元Blackwell用にはsm_120等、対象とtoolchainに合うコードを用意する。

### 1. GPUパッケージに接続し、近傍はCPUで構築

`PairAccelNetGPU`、C++ bridge、Fortran C ABIを追加する。
初版は `neigh no` のfull listをCSRへpackし、既存GPU APIを呼ぶ。
モデル・workspaceは常駐化し、原子ごとのGPU呼出はしない。
`neigh yes/hybrid` など未実装設定は明示的に拒否し、勝手に無視しない。

入力の完成形案（まだ実行不可、手元の2024版を想定）:

```lammps
package gpu 1 neigh no newton on split 1
pair_style accelnet/gpu auto Ti.nn O.nn
pair_coeff * *
```

2026-09-02以降ではsplit 1は不要。1 rank/1 GPUから始め、同一snapshotのE/F/W、
原子別energy、短いNVE、GPU不在、モデル非対応、生成・破棄を確認する。
Fortranの `error stop` をC ABI/MPI境界へそのまま露出させず、statusへ整理する。

### 2. GPUパッケージの近傍リストを直接使う

`neigh yes` を追加する。lib/gpuのnear-neighbor配列は現APIのCSRと同じ形ではなく、
packed/pitched形式なので、GPU上の変換カーネルまたはその形式を読む入力adapterを用意する。
OpenMPの外部deviceポインタ利用に対応する入口も追加し、ホストの `map` と区別する。

```text
ホスト座標 → lib/gpu座標配列
           → GPU近傍構築（更新時）
           → GPU上で入力形式変換・変位計算
           → Fortranの記述子・NN・力
           → ホストのlocal+ghost力 → MPI集約
```

近傍リストをGPU→CPUへ戻してからFortranへ再uploadする経路は、切り分け用には使えても
最終の高速経路にはしない。座標・力の毎ステップ転送はGPUパッケージの設計上残る。
近傍のみGPUで作ることが常に速いとは限らず、neigh noと同じ系で比較する。
三斜晶・skin・再構築・atom sorting・migration・周期自己像を個別に試験する。

### 3. MPI・性能・対応範囲の拡大

1/2/4 rank、複数GPU、空rank、モデル再読込と容量増加を試す。
1 GPUを複数rankで共有する構成は正しさ確認後に測定する。
CPU/GPUのE/F/W、有限差分、CPU性能回帰を継続する。
Pair/Neigh/Comm/Loopに加え、CSR変換・転送・記述子・NN・力を分けて計測する。
8/64/512原子だけでなく数千〜数万原子を含め、CPU複数rankとの比較も行う。

最初はFP64・単一Chebyshev成分/元素・通常Verlet・全体圧力を対象とする。
原子別応力、hybrid、r-RESPA、restart、LJ/n2p2等は個別に対応・検証する。
既存Fortran GPU APIの詳細と既知のNVHPCメモリ検査上の制約は
[OpenMP target説明](openmp-target.md) と検証報告を参照。

## 判断

第一段階を「package gpuの設定に従うGPU力計算」、第二段階を「GPU packageの近傍探索も
利用」に分けるのが良い。両方を完成対象として明示し、neigh noの段階でGPU近傍探索まで
実装済みと報告しない。既存Fortran資産を維持しながら、LAMMPS GPUパッケージの利用を
段階的に広げられる。速度倍率は接続後のLAMMPS実測まで確定しない。
