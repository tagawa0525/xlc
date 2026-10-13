---
status: open
labels: [lattice, performance]
opened: 2026-10-13
closed:
milestone:
blocked-by: []
resolution:
---

# 軸方向の平板の格子計算が 32 スレッドでも 1 コアしか使わない

## 問題

examples/reflector/dvc.yaml の段階 `axial-fine`（入力 examples/reflector/axial-3A-fine.toml。軸方向の平板を 1080 の列に
刻んだ 1 次元の格子）の `lattice` のプロセスが、排他の 32 スレッドのノードで 1 コアしか使っていない。

観察（2026-10-13、Slurm のジョブ 261）:

- `ps -o pcpu` で CPU 使用率 101%（走り始めて 33 分の時点）。`pcpu` は走り始めからの累積の CPU 時間を経過時間で割った
  平均なので、今この瞬間の使用率ではない。33 分のうち大半が 1 スレッドの処理だったことは示すが、今の段が何かは
  示さない（`top -H` などで確認していない）
- ログ（axial-3A-fine.log 相当）は `moc{regions=1080 tracks=211956}`。対称操作は「mirror across the horizontal axis (1/2)」
  だけで、representative_tracks=105978、representative_regions=1080
- examples/reflector/verification.md の「軸方向の反射体の平板」の表（2026-10-13 時点）で、この段階は 3316 秒。列 540 の
  `16` と `32` は 465 秒
- 入力は `azimuthal = 32`、`spacing = 0.05`（examples/reflector/axial-3A-fine.toml）
- 比較: 集合体の MOC（2A、領域は約 3 万）は 32 スレッドで 7 秒。ただし main.rs の `use_physical_cores` は
  `RAYON_NUM_THREADS` がなければ rayon を物理コアの数（16）に抑える。ジョブ 261 と 2A の測定で `RAYON_NUM_THREADS` を
  設定したかは確かめていない

領域の数は平板のほうが集合体より少ないのに、並列が効いていないように見える。原因は確かめていない。

### 並列の構造（コードを読んだ範囲）

crates/lattice/src/moc.rs の `sweep` は rayon の `par_iter_mut` で、代表のトラックの組（`representatives.chunks`）ごとに
スイープする。組の数は crates/lattice/src/moc/symmetric.rs の `balanced_chunks(tracking, &tracks, tracking.azimuths.len())` で
決まり、方位角の数に等しい（対称性を使わないときも方位角ごと）。組の中は逐次で、組の部分和は組の順に足す（結果が
スレッドの数に依らないため）。部分和の足し合わせ（`par_chunks_mut`）は領域・面ごとに並列。自己遮蔽（shielding.rs）は
ピンごとの `par_iter` で、平板の燃料のセルは少ない。main.rs は `RAYON_NUM_THREADS` がなければ物理コアの数（16）で
抑える。組の数、スレッド数は、ジョブ 261 で確かめていない。

### 候補（いずれも確かめていない）

1. 並列の単位が組（方位角の数）で、方位角は 32 なので組の数は足りるはずだが、1 次元の平板では方位角ごとのトラックの
   区間の数や長さが大きく偏り、少数の組が全体の時間を決める（`balanced_chunks` で揃えた結果は確かめていない）
2. 時間の大半が MOC のスイープでなく、並列にしていない別の段（自己遮蔽、トラッキング、CMFD、平衡値の計算、対称性の
   確認など）にある
3. 領域・セルの列の数（1080）に依って時間が伸びる部分（CMFD の行列の解き、列数の高次に比例する処理など）が逐次である。
   540 列から 1080 列で 465 秒から 3316 秒（約 7 倍）と、列数の 2 倍に対して見かけ上 2 乗を超えて伸びている。ただし
   トラックの数や反復の回数、全体の時間に含まれる段が同じかは確かめていないので、見かけの指数にすぎない

## 完了の条件

- 平板の 1 状態（`axial-fine` の入力を 1 状態に絞ったもの）を `perf` または `tracing` の span で計測し、時間がどの段にあるか、
  並列の単位（組の数、スレッド数、組ごとの時間）を確認する。計測は 1 本ずつ Slurm で流す
- 同じスレッド数で集合体の MOC 並みにスケールするようにするか、しないと決める。しない場合は、理由を ADR に書く

関係する ADR: なし

## 経緯

- 近い決定は ADR-0112（対称操作を確かめて代表のトラックだけをスイープする。ADR-0113 を置き換えた）と ADR-0128
  （軸方向の平板の格子計算）。ADR は確定済みで、前付けの `issues` に足せないので、「関係する ADR:」には挙げない

- 2026-10-13: ジョブ 261 の `axial-fine` を `ps` で見て起こした
