# G1–G5/LJ GPU対応と、同じカーネルのCPU比較（2026-09-26）

G1–G5/LJの記述子・NN・力・virialをFP64 OpenMP targetに追加した。
既存CPUの数式（cutoff、angular power、activation）は共通includeから生成する。
既存CPUの最適化ループは維持し、targetカーネルをCPUで実行する選択肢も追加した。

「GPU向けだからCPUで遅い」と仮定せず、**OpenMPなし**の同一カーネルと
OpenMP有効・1スレッドの両方を既存CPUと比較した。遅い箇所については、
数値を変えずに重複計算を取り除く比較実装で原因を検証した。

## 測定条件

- Intel Xeon Gold 6526Y、CPU 6に固定、`OMP_NUM_THREADS=1`、他の本作業のベンチと同時実行しない。
- GNU Fortran 11 Release `-O3`、NVHPC 25.3 `-fast -O3`。コンパイラごとに別集計。
- 既存CPUとtargetは同じビルド内のライブラリを使う。OpenMP版のみ
  GNU `-fopenmp`、NVHPC `-mp=gpu -gpu=cc90,cc120`を追加。
- OpenMPなしは`ACCELNET_TARGET_SERIAL=ON`：同一数値ソースの指示文を無視し、
  runtime facadeもCPU定数と`system_clock`に置換する。GNUにはOpenMPライブラリ依存がない。
  NVHPCは標準Fortran runtime経由でlibnvompもリンクするが、カーネルにOpenMP起動の呼出しはない。
- 64 / 512 / 4096原子、2元素、間隔1.7の同じ合成構造、NNは入力–8–4–1。
  Chebyshevは次数5、LJは4成分、G4/G5は各4成分、BehlerはG1–G5の14成分。
  各ケースでCPUとtargetを交互に5サンプル測定し、**対になった時間比の中央値**を使う。
- 近傍CSRは作成済み。モデル初期化・初回確保・近傍構築を除外する。
  ステップごとの領域起動、同期、バッファ確認、入出力コピーは含める。
- 全E/F/virial成分を毎サンプル検査。許容差は絶対・相対とも`2e-10`。
- Chebyshevは同じdirect/momentを双方に指定。G5も性能比較では双方direct。
  この合成モデルの倍率を、任意の実用ポテンシャルへ一般化しない。

## OpenMPを無効にした比較

時間比 = 共通target実装 / 既存CPU。1未満が高速、1より大きい値は低速。

| 記述子 | GNU 64 | GNU 512 | GNU 4096 | NVHPC 64 | NVHPC 512 | NVHPC 4096 |
|---|---:|---:|---:|---:|---:|---:|
| chebyshev direct | 0.872 | 1.069 | 1.032 | 0.595 | 0.708 | 0.707 |
| chebyshev moment | 0.754 | 0.832 | 0.829 | 0.759 | 0.830 | 0.820 |
| lj | 1.511 | 1.618 | 1.586 | 1.451 | 1.499 | 1.486 |
| g4 | 4.892 | 4.801 | 4.875 | 3.657 | 3.590 | 3.611 |
| g5 | 3.289 | 3.281 | 3.237 | 2.110 | 2.092 | 2.077 |
| behler | 2.599 | 2.561 | 2.577 | 1.757 | 1.726 | 1.696 |

Chebyshev momentは今回の条件で両コンパイラとも高速化した。
directはNVHPCで速く、GNUの512/4096原子では数%遅い。
LJ/G4/G5の遅さはOpenMPを除去しても残るため、領域起動だけでは説明できない。

## OpenMP有効・1スレッドとの差（512原子）

| 記述子 | GNU OpenMPあり | GNU なし | NVHPC OpenMPあり | NVHPC なし |
|---|---:|---:|---:|---:|
| chebyshev direct | 1.190 | 1.069 | 0.795 | 0.708 |
| chebyshev moment | 0.896 | 0.832 | 0.962 | 0.830 |
| lj | 2.206 | 1.618 | 1.824 | 1.499 |
| g4 | 5.653 | 4.801 | 5.461 | 3.590 |
| g5 | 3.802 | 3.281 | 3.004 | 2.092 |
| behler | 3.001 | 2.561 | 2.638 | 1.726 |

NVHPCでは`device(omp_get_initial_device())`だけでは実GPUへ実行されたため、
初期の誤ったhost測定は採用していない。最終実装は`target if(...)`を使い、
初期化時にも`omp_is_initial_device()`で実行先を検査する。
Nsight Systemsでhost経路にCUDAカーネル起動がないことも確認した。

原子1個・近傍なし（NNとバッファ処理は実行）の追加時間は次の通り。
これはOpenMP起動・同期を含む固定コストの目安であり、厳密にfork単体の時間ではない。

| ケース | GNU OpenMP追加時間 [µs/call] | NVHPC OpenMP追加時間 [µs/call] |
|---|---:|---:|
| chebyshev direct | 2.06 | 2.66 |
| chebyshev moment | 3.28 | 4.47 |
| g4 | 2.99 | 4.04 |

実計算の差はミリ秒単位。OpenMP有効/無効で生成コードやatomic処理も変わるので、
その差のすべてを「起動コスト」とは呼ばない。

## 遅さの原因と、数値を保った比較実装

### G4/G5: 同じradial項と原子対の再計算

- GPU向け共通実装は、各記述子・各原子対で`exp(-eta*(r-Rs)^2)*cutoff(r)`と
  その微分を繰り返し計算する。既存CPUはcutoff/exponentialのグループを作り、近傍ごとに再利用する。
- 値はunordered pairを1回、力は独立した各edgeから両方向に2回評価する。
  GPUでedge単位に並列化でき、完全なdescriptor Jacobianを保存せずに済む一方、計算量が増える。
  既存CPUのG4は1回で値と両側のJacobianを計算し、G5 directは値・力の各段階で
  unordered pairを1回ずつ評価する。Chebyshev directの力も既存CPUは両側を一度に計算する。
- 初期のgeneric値カーネルは値/微分共通helperを呼び、不要な微分も計算する。
  CPUの合成G4/G5では、targetの大半の時間が記述子と力にあり、NNは1%未満。
- `cached`実験はparallelな仕事の分け方・NN・scatterを変えず、q(r),q'(r)だけを近傍ごとに保存する。
  キャッシュ作成・確保の時間も含める。さらにCPUと全E/F/virial・有限差分の一致を確認した。

| 記述子・512原子 | GNU 元の共通コード | GNU radial再利用 | NVHPC 元の共通コード | NVHPC radial再利用 |
|---|---:|---:|---:|---:|
| g4 | 4.801 | 2.558 | 3.590 | 2.336 |
| g5 | 3.281 | 1.509 | 2.092 | 1.285 |
| behler | 2.561 | 1.299 | 1.726 | 1.112 |

再利用だけで差の大きな部分が消えるため、これは演算の重複が原因であるという直接の証拠。
残る原子対の重複評価やパラメータグループ共有の差を、キャッシュミスと推測して置き換えない。
ハードウェアカウンタによるcache miss/bandwidth測定は行っていない。

### LJ: 2成分を別々に計算するコスト

既存CPUはr^-6/r^-12、cutoff、微分を一度に計算する。
generic targetは2成分×値/力で同じcutoff等を繰り返す。
`lj-powers`は累乗の式だけを変え、`lj-fused`は2成分をまとめて計算する。

| 512原子 | 元の共通コード | 累乗だけ置換 | 2成分まとめ計算 |
|---|---:|---:|---:|
| gnu | 1.618 | 1.639 | 1.056 |
| nvhpc | 1.499 | 1.424 | 0.999 |

GNUでは累乗だけの変更は改善にならず、両コンパイラともまとめ計算で既存CPUに近づいた。
したがってLJの主因は重複した成分計算。累乗関数だけのせいとは結論しない。

これらの比較実装はCPU専用の診断コードとして保存した。**改善版のGPU性能は未測定**。
GPUでも演算削減は有望だが、中間配列の容量・帯域、カーネル追加コスト、並列度を含めて測る必要がある。
GPU用キャッシュは単純なfeature×edge全保存より、同じradialパラメータのグループ共有が候補になる。

## GPU機能・回帰検証

- H100: 29 GPUテスト + 2 hostテスト、全31件成功。Blackwell: 全29 GPUテスト成功。
- GNU host、GNU OpenMPなし、NVHPC OpenMPなし: 各2スイート成功。
- 新記述子スイートは72ケース。全cutoff、混合species/family、local species mapとパラメータのreload、
  周期画像、孤立原子、力とvirialの有限差分を含む。
- H100 Compute Sanitizer: 新記述子72ケース、`ERROR SUMMARY: 0 errors`。
  NVHPC runtimeの終了時プールについてleak-freeという主張はしない（`--leak-check no`）。
- LAMMPS: LJ / Behler G1–G5 / 変換n2p2の計63比較。
  CPU / GPU近傍no・yes・hybrid、1/2 MPIランク、直交・三斜晶・空rank、短時間MDのE/F/virial/atom Eを確認。
- 既存CPU: 数値42テスト + 性能2テスト成功。12ケースの変更前との時間比は0.999–1.015。

G4/G5のGPU計算はdirect。G5 moment強制指定は拒否する。
Fortran APIではn2p2を直接ロードできる。LAMMPSではn2p2→AccelNet変換後の埋め込みモデルを使う。
Chebyshevと他成分を同じ元素内で組み合わせるモデル、複数Chebyshev成分は未対応。

## 既存ChebyshevのGPU速度

H100、実Ti/Oモデル、24,000原子、20 warmup + 100 measured steps、3サンプルの中央値。
各実行の全E/F/atom E/pressureをdirect基準と比較済み。

| mode | 変更前 [ms/step] | 変更後 [ms/step] |
|---|---:|---:|
| direct | 27.767 | 27.681 |
| moment | 10.791 | 10.528 |
| auto | 10.788 | 10.532 |

## 再現

ビルド・比較・診断のコマンドは[OpenMP target手順](../../openmp-target.md)を参照。
GPUベンチは`benchmark_gpu_modes.py --replicates 10`、新モデルのLAMMPS検証は
`check_gpu_descriptors.py`を使った。各ディレクトリに生ログ・JSON・入力を保存している。
診断用Fortranの実ソースは`experiments/`に保存した。

既存CPUの既定経路は維持し、共通カーネルのCPU実行とOpenMPなしビルドを選択可能にした。
今回の結果だけを根拠に、すべての記述子・コンパイラ・モデルで速い経路を固定しない。
