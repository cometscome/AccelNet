# LAMMPS GPU package + Fortran OpenMP target

`interfaces/lammps/29Aug2024` に `pair_style accelnet/gpu` を追加した。
LAMMPS GPU package がGPUを割り当て、既存Fortran OpenMP targetが記述子・NN・力を計算する。
`neigh no`、`neigh yes`、`neigh hybrid` に対応する。

moment経路は原子ごとに成分をまとめて構築し、力用係数の事前計算と多項式勾配の縮約を行う。
H100・24,000原子のTi/Oでは25.07→10.81 ms/step、Blackwellでは54.25→17.19 ms/step。
`auto`でmomentを選んだ場合にも適用される。[実装・数値検証・速度測定](validation/gpu-moment-force-2026-09-26/README.md)。

## 対応範囲

- 対象: LAMMPS 29 Aug 2024 Update 4、NVHPC 25.3、CUDA、`GPU_PREC=double`。
- GPU: H100 / Blackwellを検証済み。[実測結果と検証報告](validation/lammps-gpu-2026-09-26/README.md)を参照。
- モデル: `.nn` / `.nn.ascii` の埋込みChebyshev・LJ・Behler G1〜G5モデル。LAMMPSのtype順は
  ネットワークに埋め込まれたglobal species順と一致させる。バージョンは既存LAMMPS CPU版と同じ0。
- Chebyshevの`auto` / `direct` / `moment`と、独立した`g5 auto` / `g5 direct` / `g5 moment`。G4はdirect計算。
- 局所原子・ghost原子、周期境界、直交セル・制限三斜晶。
- 全体エネルギー・力・全体virial・原子別エネルギー。`newton on`、旧GPU packageでは`split 1`が必要。
- `pair_style accelnet` とGPUを無効にしたCPUビルドは引き続き利用可能。
- 原子別stress、pair hybrid、r-RESPA、分子トポロジー、近傍除外は未対応で拒否する。
  n2p2ディレクトリの直接指定は未対応だが、下記の変換後モデルを利用できる。
  restartへのモデル保存は非対応。AMD/HIP/OpenCLの相互運用は未実装。

## ビルド

既存ビルドと分けたディレクトリを使う。以下では`$ACCELNET_SOURCE`がこのリポジトリ、
`$LAMMPS_SOURCE`がLAMMPS 29Aug2024 Update 4ソース、`$NVHPC_ROOT`が
`/opt/nvidia/hpc_sdk/Linux_x86_64/25.3` を指す。

```bash
cmake -S "$ACCELNET_SOURCE" -B "$ACCELNET_SOURCE/build-gpu" \
  -DCMAKE_Fortran_COMPILER="$NVHPC_ROOT/compilers/bin/nvfortran" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DACCELNET_BUILD_OPENMP_TARGET=ON \
  -DACCELNET_OPENMP_TARGET_FLAGS="-mp=gpu -gpu=cc90,cc120"
cmake --build "$ACCELNET_SOURCE/build-gpu" --target AccelNetTarget -j 8

python3 "$ACCELNET_SOURCE/interfaces/lammps/29Aug2024/install.py" "$LAMMPS_SOURCE"
cmake -S "$LAMMPS_SOURCE/cmake" -B "$LAMMPS_SOURCE/build-accelnet-gpu" \
  -DCMAKE_CXX_COMPILER="$NVHPC_ROOT/compilers/bin/nvc++" \
  -DCMAKE_Fortran_COMPILER="$NVHPC_ROOT/compilers/bin/nvfortran" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG" \
  -DBUILD_MPI=ON -DBUILD_OMP=OFF \
  -DPKG_ACCELNET=ON -DPKG_GPU=ON -DGPU_API=cuda -DGPU_PREC=double \
  -DGPU_ARCH=sm_90 -DCUDA_BUILD_MULTIARCH=OFF \
  -DCUDA_NVCC_FLAGS="-gencode=arch=compute_120,code=sm_120" \
  -DCUDA_TOOLKIT_ROOT_DIR="$NVHPC_ROOT/cuda/12.8" \
  -DBIN2C="$NVHPC_ROOT/cuda/12.8/bin/bin2c" \
  -DACCELNET_DIR="$ACCELNET_SOURCE/build-gpu" \
  -DACCELNET_TARGET_FLAGS="-mp=gpu -gpu=cc90,cc120"
cmake --build "$LAMMPS_SOURCE/build-accelnet-gpu" -j 16
```

MPIが自動検出されない場合は`-DMPI_CXX_COMPILER=/path/to/mpicxx`を指定する。
H100のみならBlackwell向け`cc120` / `compute_120`指定は省ける。
`install.py` はACCELNET/GPU pair、lib/gpu adapter、CMake moduleをコピーし、
CMakeのpackage登録と三斜晶ソートの修正を適用する。
既存のCPUインターフェースファイルもこのリポジトリの版へ更新する。

Fortranランタイムを手動で`-lnvf -lnvomp ...`と並べず、NVHPCの`-fortranlibs`を使う。
手動のライブラリ順ではGPUが存在しても`omp_get_num_devices()`が0を返した。
付属CMakeはLAMMPS本体によるランタイムの重複追加も抑制する。

## 入力例

```lammps
package gpu 1 neigh yes newton on split 1
units metal
atom_style atomic
read_data TiO2.data
pair_style accelnet/gpu auto Ti.nn.ascii O.nn.ascii
pair_coeff * *
neighbor 0.6 bin
neigh_modify every 10 delay 0 check yes
fix integrate all nve
run 100
```

`package gpu`はboxを作る前に置く。`CUDA_VISIBLE_DEVICES`で使用GPUを限定できる。
LAMMPSがその範囲から各MPI rankへ割り当てたdevice番号をFortranにも渡す。
CPU近傍との比較は`neigh no`へ変更する。小規模系ではCPUの方が速い場合がある。

## データと通信

`neigh no`ではfull neighbor listをCPUでCSRに詰め、一括でGPUへ送る。
`neigh yes/hybrid`ではlib/gpuの座標と近傍配列を借り、GPU上でCSRと変位に変換する。
近傍本体をCPUへdownloadして再uploadする処理はない。CSRの辺数を表す整数1個は回収する。
モデル・NN作業領域・descriptor/force用配列は常駐し、容量不足時だけ再確保する。

このLAMMPS版の標準CUDA設定では、`neigh yes`でもlib/gpu内部のbinningはCPUで、
近傍候補の探索はGPUで行う。`neigh hybrid`と内部実装が共通になる場合がある。
「全binningもGPU」とは扱わない。

各`compute()`内でGPU計算とlocal+ghost力のdownloadを完了する。
その後にLAMMPSがghost力のreverse communicationを行う。
FixGPUのpost_force answer queueにはAccelNetの結果を重複登録しない。
GPUライブラリとOpenMPの間は同期してから外部device pointerを渡す。

LAMMPS 2024のGPU用原子ソートは三斜晶で未初期化の直交セル用sublo/subhiを読む箇所がある。
付属patchは三斜晶では通常のbboxに基づく原子ソートを使わせる。sorting自体は有効のまま。
NVHPC 25.3でelemental finalizerの誤った呼出を確認したため、workspaceの終了処理は
scalar/rank-1用の明示的なfinalizerにした。C handleごとにモデルとworkspaceを解放する。
OpenMPランタイムのキャッシュはLAMMPSの`clear`より長く生存するため、利用GPUのprimary
CUDA contextはプロセス終了まで保持する。`clear`後のモデル再読込も検証済み。

## 検証・速度チェック

`interfaces/lammps/29Aug2024/tests/check_gpu.py` はCPU 1 rankを基準に、
エネルギー・原子別エネルギー・全原子の力・virial圧力6成分・短いNVE軌道を比較する。
既定は1/2/4 rank、直交・三斜晶・空rank、毎stepのneighbor rebuildとsorting。
`--gpus 2 --ranks 2 4`で複数GPUも検証できる。

```bash
export OMP_NUM_THREADS=1
export OMP_TARGET_OFFLOAD=MANDATORY
python3 interfaces/lammps/29Aug2024/tests/check_gpu.py \
  --lammps /path/to/lmp --golden /path/to/fortran_predict --output /tmp/gpu-check
python3 interfaces/lammps/29Aug2024/tests/check_gpu_errors.py \
  --lammps /path/to/lmp --input /tmp/gpu-check/orthogonal-auto-no-1rank/in.test \
  --output /tmp/gpu-errors
python3 interfaces/lammps/29Aug2024/tests/benchmark_gpu.py \
  --lammps /path/to/lmp --input /tmp/gpu-check/orthogonal-auto-cpu-1rank/in.test \
  --output /tmp/gpu-benchmark
```

比較スクリプトはPython + NumPy、MPI launcherを必要とする。
benchmarkはwarmup後のLAMMPS Loop time、Pair timeを記録し、複数回の中央値を出す。
CPU 1/4 rank、GPU 1 rank、24/192/5184/24000原子が既定。ログにはNeigh/Comm等も残る。
単一rankの速度測定はCPU 6へ固定するため、別環境ではスクリプトのCPU番号を変更する。

`direct` / `moment` / `auto` の比較には次を使う。
GPU近傍、1 rank、20 step warmup + 100 stepの計測を各3回行い、モードの順序を交代する。
計測後に全原子の力・原子別エネルギー・座標を比較し、thermoのエネルギー・virialも確認する。
dump出力は計測区間に含めない。CPU固定先は`--cpu`で変更できる。

```bash
python3 interfaces/lammps/29Aug2024/tests/benchmark_gpu_modes.py \
  --lammps /path/to/lmp --input /tmp/gpu-check/orthogonal-auto-cpu-1rank/in.test \
  --output /tmp/gpu-modes
```

[H100 / Blackwellでのdirect・moment比較結果](validation/lammps-gpu-modes-2026-09-26/README.md)。
その後のmoment最適化を含む最新結果は[こちら](validation/gpu-moment-force-2026-09-26/README.md)。
今回のTi/Oモデルではmomentが速く、autoもほぼ同じ速度となった。
一般には角度次数と近傍数、GPU、原子数に依存する。autoは近傍数と角度次数に基づく
演算数の見積りであり、実測による自動チューニングは行わない。

`check_gpu_lifecycle.py`は原子数増加と`clear`後の再読込を検証する。
`check_gpu.py --cases migration --rebuild-every 5 --modes direct moment`で
rank間移動と近傍再利用も比較できる。

Fortranの既存数値テストとCPU性能回帰に加え、`predictor_target_c_api`は複数handleの
独立性、不正CSR、不正モデルファイル、繰り返し生成・破棄を確認する。
`read_aenet_network*`にはoptional statusを追加し、C APIからのファイル読込エラーを
LAMMPS側でMPI-safeに処理する。既存Fortran呼出のstatus省略時の失敗動作は維持する。

## LJ・Behler・n2p2モデル

埋め込み記述子付きのLJ・BehlerネットワークはChebyshevと同じ入力方法で使える。
n2p2 2G-HDNNPは既存コンバーターで記述子順序、重み、正規化を保持して変換する。

```sh
build/bin/accelnet-model-converter-fortran n2p2-to-accelnet n2p2-model converted
```

```lammps
package gpu 1 neigh yes newton on split 1
pair_style accelnet/gpu auto converted/H.nn.ascii converted/O.nn.ascii
pair_coeff * *
```

type順は変換後モデルのglobal species順に合わせる。G5のmoment法は次のように指定する。

```lammps
pair_style accelnet/gpu auto converted/H.nn.ascii converted/O.nn.ascii g5 moment
```

末尾の`g5 moment`がG5を選択し、先頭の`auto`はChebyshev用である。
G5の`auto`は従来と同じく、各成分の角度カットオフ内に16以上の近傍があると
moment法を選ぶ。整数次数1〜10が対象で、非整数・高次数はdirect法を併用する。
`g5 moment`は近傍数の閾値を解除するが、次数の上限は解除しない。
`tests/check_gpu_descriptors.py` はLJ、G1〜G5、変換n2p2について、CPU・GPU近傍各方式、
1/2 MPIランク、直交・三斜晶・空ランクの短時間MDを比較する。
