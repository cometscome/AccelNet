# direct・moment・aenet 比較ベンチマーク

この手順は、同じ学習済みポテンシャルを使って次の3方式を比較します。

| ラベル | LAMMPS pair style | AccelNet内部アルゴリズム |
|---|---|---|
| `aenet` | `aenet` | original aenet |
| `accelnet-direct` | `accelnet` | 角度近傍ペアを直接列挙 |
| `accelnet-moment` | `accelnet` | 厳密なCartesian moment展開 |

1個のAccelNet/LAMMPS実行ファイルに対して、AccelNetの実行時APIで`direct`または
`moment`を明示的に指定します。LAMMPSログに実際のモード名を出力し、測定スクリプトは
期待したモードと一致しなければ停止します。`auto`は通常利用できますが、この3方式比較
では使いません。`pair_style aenet`実行時にAccelNetライブラリは計算に使われません。

LAMMPS入力では次のように選択できます。モードを省略した場合は`auto`です。

```lammps
pair_style accelnet auto   H.ann O.ann
pair_style accelnet direct H.ann O.ann
pair_style accelnet moment H.ann O.ann
```

## 1. Linuxでビルド

各ソースの場所を環境変数で明示できます。

```sh
cd "$HOME/AccelNet-lammps"

ACCELNET_SOURCE="$HOME/AccelNet" \
AENET_SOURCE="$HOME/aenet-master" \
AENET_LAMMPS_SOURCE="$HOME/aenet-lammps" \
LAMMPS_SOURCE="$HOME/lammps-4Feb20" \
BENCH_BUILD_JOBS=4 \
./linux/build_variants_linux.sh
```

生成物は次の1個です。

```text
lammps-4Feb20/src/lmp_mpi_accelnet_modes
```

同一バイナリ内でモードを切り替えるため、異なる名前の同一バイナリを比較する問題は
生じません。GNU Fortran 11で並列コンパイル時にinternal compiler errorが出る場合は、
`BENCH_BUILD_JOBS=4`または`1`に下げてください。

## 2. energy・force一致の確認

本測定の前に、H2OとTiO2について3方式のenergy・forceを比較します。

```sh
AENET_SOURCE="$HOME/aenet-master" \
LAMMPS_SOURCE="$HOME/lammps-4Feb20" \
./linux/smoke_test_variants.sh
```

すべて`PASS`になった場合だけ速度測定へ進んでください。

## 3. バックグラウンド実行

次のコマンドはSSH端末から切り離して実行し、結果ディレクトリを表示して直ちに戻ります。

```sh
AENET_SOURCE="$HOME/aenet-master" \
LAMMPS_SOURCE="$HOME/lammps-4Feb20" \
./linux/run_variants_background.sh
```

既定条件は、H2OとTiO2、1,152原子、固定座標と動的NVT、1/2 MPI、3反復、
warm-up 100 step、測定500 stepです。各trial内でaenet、direct、momentの順に
実行して、マシン負荷や温度ドリフトの影響を抑えます。

5反復で測定する例:

```sh
AENET_SOURCE="$HOME/aenet-master" \
LAMMPS_SOURCE="$HOME/lammps-4Feb20" \
BENCH_MPI_RANKS='1 2 4 8' \
BENCH_TRIALS=5 \
BENCH_WARMUP_STEPS=100 \
BENCH_MEASURE_STEPS=500 \
./linux/run_variants_background.sh
```

実行時に表示された`RESULT_DIR`を使って監視します。

```sh
RESULT_DIR="$HOME/AccelNet-lammps/linux/results/variant-run-YYYYMMDD-HHMMSS"
cat "$RESULT_DIR/status.txt"
tail -f "$RESULT_DIR/driver.log"
ps -p "$(cat "$RESULT_DIR/background.pid")"
```

`status.txt`は`STARTING`、`RUNNING`、`COMPLETED`の順に変化します。エラー時は
`FAILED exit_code=...`になります。停止する場合はPIDと結果ディレクトリを確認してから
`kill "$(cat "$RESULT_DIR/background.pid")"`を実行します。

## 4. 条件の選択

```sh
BENCH_SYSTEMS='water' \
BENCH_MODES='dynamic' \
BENCH_IMPLEMENTATIONS='aenet accelnet-direct accelnet-moment' \
BENCH_MPI_RANKS='1 2' \
BENCH_TRIALS=3 \
BENCH_MEASURE_STEPS=500 \
./linux/run_variants_background.sh
```

主な変数は`BENCH_SYSTEMS` (`water tio2`)、`BENCH_SIZES`
(`small medium large`)、`BENCH_MODES` (`fixed dynamic`)、
`BENCH_IMPLEMENTATIONS`、`BENCH_MPI_RANKS`、`BENCH_TRIALS`、
`BENCH_WARMUP_STEPS`、`BENCH_MEASURE_STEPS`です。

## 5. 結果

結果ディレクトリには次が保存されます。

```text
variant-run-YYYYMMDD-HHMMSS/
├── background.pid
├── driver.log
├── runner.pid
├── status.txt
├── system-info.txt
├── summary.tsv
└── logs/
```

`summary.tsv`には各条件のloop time中央値・IQR、steps/s、atom-step/s、LAMMPS timing
breakdown、`speedup_vs_aenet`、`speedup_moment_vs_direct`が入ります。
`system-info.txt`にはCPU、コンパイラ、MPI、Git revision、バイナリ・モデル・方式選択に
関わるソースのSHA-256を記録します。

## SLURMを使う場合

`nohup`はログインノードや、終了する対話allocationを越えて計算資源を保証しません。
管理クラスタではバッチスクリプト内で次を実行してください。

```sh
export BENCH_LAUNCHER=srun
export BENCH_MPI_RANKS='1 2 4 8'
export BENCH_RESULT_DIR="$PWD/linux/results/variant-${SLURM_JOB_ID}"
./linux/run_variant_benchmarks.sh
```

この部分をサイト固有の`#SBATCH`指定とともにジョブスクリプトへ入れ、`sbatch`で投入します。
