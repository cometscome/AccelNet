# LAMMPS GPU direct / moment比較 — 2026-09-26

GPU計算コードは変更せず、公開済みの実行ファイルでモードを切り替えて測定した。
このTi/Oモデルでは両GPUともmomentが速いが、24,000原子では差が小さくなる。
autoはmomentとほぼ同じ時間だった。これは今回のモデル・構造・実装についての結果。

## 条件

- LAMMPS 29 Aug 2024 Update 4、NVHPC 25.3、CUDA倍精度。
- executable SHA256: `6fb237a8c912d811423d1d23d5e74a949f2e8fc66aebe307173c22f1688d28b3`
- H100 NVL: `GPU-2644154d-7268-af42-6631-59e1f3c6e7f3`。
- RTX PRO 6000 Blackwell Max-Q: `GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9`。
- 各GPUで順次実行。開始時に選択GPUが空いていることを確認した。
  別のH100は他のプロセスが使用中であり、ホスト全体は専有していない。
- 1 MPI rank、CPU 6に固定、OMP_NUM_THREADS=1。
- Ti/O実モデル、TiO2 24原子構造を各軸に1/2/6/10倍複製。
  Blackwellは192/5184/24000原子を測定。
- radial cutoff 6.5 Å / order 20、angular cutoff 5 Å / order 6。
  周期像を含めた初期角度近傍数は48〜54、平均50.17。
- `package gpu 1 neigh yes newton on split 1`。近傍skin 0.6 Å、
  `neigh_modify every 10 delay 0 check yes`、`atom_modify sort 100 0.0`。
- 初期温度100 K、dt=0.0001 ps、NVE 20 step warmup + 100 step計測。
- 各3回、実行順はdirect→moment→auto、moment→auto→direct、auto→direct→moment。
  各sampleは新規プロセスで同じ初期条件から開始。
- 表はLAMMPS Loop time / 100の中央値。近傍更新、転送、記述子、NN、力、積分を含む。
  初期化、終了後の力・原子別エネルギーdumpは含まない。
  記述子カーネル単独の速度ではない。

## 結果

単位ms/step、倍率はdirect / moment（1より大きいほどmomentが速い）。

| GPU | 原子数 | direct | moment | auto | momentの倍率 |
| --- | ---: | ---: | ---: | ---: | ---: |
| h100 | 24 | 8.327 | 1.573 | 1.572 | 5.29× |
| h100 | 192 | 9.073 | 1.607 | 1.605 | 5.65× |
| h100 | 5,184 | 13.703 | 5.153 | 5.158 | 2.66× |
| h100 | 24,000 | 29.093 | 25.069 | 25.071 | 1.16× |
| blackwell | 192 | 17.246 | 2.689 | 2.690 | 6.41× |
| blackwell | 5,184 | 25.193 | 9.109 | 9.116 | 2.77× |
| blackwell | 24,000 | 57.780 | 54.249 | 54.254 | 1.07× |

## 数値一致

各サイズ・GPUで最初のdirect実行を基準に、すべての実行について
最終全原子の力・原子別エネルギー・座標と、thermoに出力されたエネルギー・virial圧力6成分を比較。
原子配列は`atol=2e-8, rtol=2e-9`、thermoは`atol=2e-5, rtol=2e-9`。
今回の比較はGPUモード間。CPUとの比較は[先の検証報告](../lammps-gpu-2026-09-26/README.md)を参照。

全63実行が完了し、比較はすべて合格。基準自身の7実行も含む。最大絶対差：

- 力: 1.662e-13 eV/Å
- 原子別エネルギー: 2.27374e-13 eV
- 出力された全体エネルギー: 0 eV
- virial圧力: 8.14907e-10 bar

## 解釈と使い分け

`direct`は各中心原子の近傍ペアについて角度項を直接加算する。
`moment`は近傍方向の単項式の和を先に計算し、その積から同じ角度項を得る。
近傍数をK、角度次数をLとすると、主要項はdirectのK²(L+1)に対して
momentではK×M（M=(L+1)(L+2)(L+3)/6）となる。momentは近似による省略ではないが、
演算順序が異なるため浮動小数点の結果はbit単位では一致しない。

今回L=6、M=84。GPUのautoは`K*(L+1) >= M`でmomentを選ぶので、初期構造では
すべての中心原子がmoment側になる。autoは実時間を測る方式ではなく、
GPUの機種や総原子数を使った性能モデルもまだ持っていない。

[追加のカーネル計測](profile/README.md)で、moment構築が5,184→24,000原子で
1.23→12.36 msに増加することを確認した。原子数4.63倍に対して約10倍であり、
moment構築の原子あたり効率が低下している。一方、momentの近傍ごとの力計算は
ほぼ原子数に比例した。倍率縮小をdirectの並列度向上だけで説明するのは不十分だった。
近傍データの再読込とキャッシュ容量の影響が原因候補だが、ハードウェアカウンタは
権限制限で取得できず、メモリ帯域・キャッシュミスへの原因特定は未完了。

その後の[処理順序だけを変えた比較実験](locality/README.md)では、24,000原子の
構築カーネルが12.36→1.04 ms、LAMMPS全体が25.21→13.86 msに短縮した。
同じ原子の複数成分を近くに配置して近傍データを再利用する効果が大きく、
現行構築カーネルのデータ再利用・メモリアクセス配置が主要因と判断する。
診断版での結果であり、既存の実行ファイルや製品コードは変更していない。

このTi/Oモデルではautoを使ってよい。高次数でMが増える場合や近傍が少ない場合まで
moment優位とは限らず、別モデルでは同様の実測で選ぶ。

## 再現と記録

[測定スクリプト](../../../interfaces/lammps/29Aug2024/tests/benchmark_gpu_modes.py)を追加した。
`check_gpu.py`でCPU用in.testを生成した後、次を実行する。

```bash
CUDA_VISIBLE_DEVICES=GPU-2644154d-7268-af42-6631-59e1f3c6e7f3 \
python3 interfaces/lammps/29Aug2024/tests/benchmark_gpu_modes.py \
  --lammps /home/nagai/AccelNetGPU/lammps-accelnet-gpu/lmp \
  --input /tmp/accelnet-lammps-gpu/validation-h100/orthogonal-auto-cpu-1rank/in.test \
  --output /tmp/accelnet-lammps-modes-h100
```

BlackwellはUUIDを変更し、`--replicates 2 6 10`で実行。
JSONは全sample、Loop/Pair時間、最大誤差、command、実行ファイルhashを含む。

- [H100](h100.json) / [Blackwell](blackwell.json)
- `h100-logs/` / `blackwell-logs/`: 全sampleの入力と標準出力。
- [モデル・構造・ソースのSHA256](sources.json)、[初期構造](orthogonal.data)。
- 元のログにある絶対パスはこのマシンの測定時の値。再現時はモデルとdataのパスを合わせる。

GPU/CPU計算コードやautoの選択条件はこの測定で変更していない。
