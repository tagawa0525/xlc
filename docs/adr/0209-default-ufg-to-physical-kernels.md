---
status: accepted
date: 2026-10-28
requires: [ADR-0070, ADR-0076, ADR-0185, ADR-0203, ADR-0204, ADR-0205, ADR-0207, ADR-0208]
supersedes: []
superseded-by:
amends: [ADR-0203, ADR-0204, ADR-0205, ADR-0207]
issues: ["#0152"]
---

# ADR-0209: UFG の減速計算の既定を、自由ガスと S(α,β) の散乱核とリングの中性子束の移行にする

## 背景

自己遮蔽の UFG の減速計算には、散乱核の標的と群の移行の形について、入力で選べる方法が並んでいる。

- 散乱核: 静止した標的（[ADR-0070](0070-treat-targets-at-rest.md)）、分離共鳴を持たない核種の自由ガス（[ADR-0203](0203-scatter-light-nuclides-by-free-gas-in-ufg.md)）、
  分離共鳴を持つ核種の 400 kT より下の自由ガス（[ADR-0204](0204-scatter-resonant-nuclides-by-free-gas-in-ufg.md)）、H in H2O の H の B(4) より下の
  S(α,β) の移行（[ADR-0205](0205-scatter-h-in-h2o-by-sab-in-ufg.md)）
- UFG の範囲の群の散乱行列の移行の形: 無限希釈の形（[ADR-0076](0076-store-legendre-scattering-moments.md)）、等価なセルのリングの UFG の中性子束で
  数えた移行（[ADR-0207](0207-collapse-ufg-transfers-by-ring-flux.md)。間引いた格子と組めることは [ADR-0208](0208-count-ring-flux-transfers-on-thinned-grids.md)）

これらの ADR は方法を入力で選べる形にし、既定（入力の省略）を静止した標的と無限希釈の形に残して、既定の切り替えをユーザーに残した
（各 ADR の「既定の案」）。参照解の OpenMC 0.15.2 は、H in H2O の H を B(4) より下で S(α,β)、ほかの核種を 400 kT より下で断面積一定の
自由ガスとして扱い（ADR-0204 の背景）、群の移行を連続エネルギーの散乱で数える。静止した標的と無限希釈の形は、その扱いからの近似である。

ユーザーの決定（2026-10-28。「なぜ精度が悪いものがデフォルトになるのですか？」の後に、既定を物理的に正しい方法にする案に「それでいいです」）:
既定は物理的に正しい方法にする。速度のための近似は、入力でそれを選んだときだけ使う。参照そのもの（軽い核種の自由ガスの上限なし、
スイープを収束まで）は時間がかかり既定にできないので、参照からの差を測った速い設定を既定にする。

判断の軸: 正確さ > 速度 > 保守。後方互換を作らない（古い既定を保つための黙ったフォールバックや、入力への自動の追記をしない）。

## 検討した案

- 案 D1: 既定を静止した標的と無限希釈の形に残す（今まで）。既定の結果の記録はそのまま使えるが、入力を省いた計算が、参照解と異なる、
  精度の低い近似になる。静止した標的の k は、散乱核の参照（S(α,β)、重い核種は 400 kT まで、軽い核種は上限なし）から 1C の新しい燃料で
  −94.9 pcm、2A-40 で −46.7 pcm 離れる（examples/depletion/verification.md「H in H2O の S(α,β) の UFG の散乱核（ADR-0205）」の表の
  静止した標的と参照の行の差。どちらも無限希釈の形）。近似を選んだことが入力に現れない
- 案 D2: 既定を参照そのもの（`{ upper = inf, cells = 4, sweeps = inf, heavy = 400, sab = true }`）にする。参照は全範囲の格子点でスイープを
  6 回ほど繰り返し、上限のある速い方法より重い（ADR-0203 の帰結。`--timed` では測っていない）。速い方法でも 2A の集合体の燃焼の自己遮蔽は
  1 時間 42 分（「確認」）なので、参照の集合体の燃焼はそれ以上かかる。速い方法との差は −4.3、−3.0 pcm（ADR-0205）で、1C と 2A の燃焼の
  OpenMC の標準偏差（9〜22 pcm）より小さい
- 案 D3: 既定を、散乱核は ADR-0205 の速い方法 `{ upper = 400, cells = 2, sweeps = 3, heavy = 400, sab = true }`、移行はリングの中性子束
  （ADR-0207）にする（採用）

### 参照の条件の書き方

条件の項目ごとの既定を速い方法の値にすると、参照の値（上限なし、収束まで、重い核種は静止した標的）を書く方法が要る（今までは省略が参照の
値だった）。

- `upper = inf`、`sweeps = inf` を「上限なし」「収束まで」とし、重い核種を静止した標的にするのは `heavy = false`（採用）。`heavy = inf` は
  ADR-0204 のとおり上限なしの自由ガス。`inf` は TOML の浮動小数点の無限大で、`heavy` がすでに受け付けている書き方と揃う
- 「上限なし」を 0 などの特別な数で書く: 0 は上限 0（自由ガスにしない）と読めるので、上限なしと逆の意味になる
- 項目が省略されたら参照、`ufg-free-gas` 全体が省略されたら速い方法: 一部だけ書いた条件と全く書かない条件で、書いていない項目の値が
  変わる

## 決定と理由

案 D3 にする。入力の既定（省略したときの値）を次にする。

- `ufg-scattering`: `"free-gas"`（今まで `"static"`）
- `ufg-free-gas` の各項目: `upper = 400`、`cells = 2`、`sweeps = 3`、`heavy = 400`、`sab = true`（今まで上限なし、4、収束まで、重い核種は静止した
  標的、false）。参照の値は `upper = inf`、`sweeps = inf`、`heavy = false` で書く（上の「参照の条件の書き方」）。`ufg-free-gas` は
  `ufg-scattering` を省略しても書ける（`"static"` と組めば止める）
- `ufg-transfer`: `"ring-flux"`（今まで `"infinite-dilution"`）

これは ADR-0203、ADR-0204、ADR-0205、ADR-0207 の「既定」の部分（ADR-0203 の「既定は静止した標的のまま」と条件の省略の値、ADR-0204 の
「省略すると重い核種は静止した標的」、ADR-0205 の「省略すると false」、ADR-0207 の「既定は `infinite-dilution`」）を置き換える。方法そのもの
（核、区間の機構、スイープ、S(α,β) の移行の表、移行の数え方、間引いた格子での数え方）は変わらない。ADR-0203、ADR-0204、ADR-0205 は、
この ADR で既定に採るので accepted にする。

理由: 速い方法は参照（OpenMC と同じ扱い）からの k の差が −4.3、−3.0 pcm（1C の新しい燃料、2A-40。ADR-0205）で、静止した標的と無限希釈の形
より参照に近い。移行は、状態で変わる群の中のスペクトルを行に入れ、1C の新しい燃料から 2A-40 への k の差の増えを 42.0 pcm 縮めた（ADR-0207）。
間引いた格子との組は ADR-0208 で 1 pcm 以内。静止した標的と無限希釈の形は入力で選べる近似として残す。

時間と、切り替えの前後の k の差は「確認」に書く。

## 帰結

- 入力に `ufg-scattering`、`ufg-free-gas`、`ufg-transfer` を書かない計算は、新しい既定で解く。既定を保つための入力の書き換えや読み替えは
  しない。前の既定の結果は、入力に `ufg-scattering = "static"` と `ufg-transfer = "infinite-dilution"` を書けば同じになる（自己遮蔽の条件
  `ShieldingOptions` が前の既定と同じになり、解く側のコードは変えていない）
- 自己遮蔽の時間が増える（`--timed` のジョブで測った。下の「確認」）: 1 回の自己遮蔽は、1C のピンセルの 1 リングで 8.5 → 76〜77 秒（約 9 倍）、
  2A の集合体の燃焼で 9.7 → 74〜77 秒、2A の燃焼の自己遮蔽の和は 787.3 → 6,131.5 秒（7.8 倍）、ジョブの経過時間は 18 分 58 秒 → 1 時間 47 分
  59 秒。2A（解くピン 39 本、格子点の和 370 万）と 1C（1 本、88 万）の 1 回の自己遮蔽がほぼ同じで、大きさに釣り合わない。候補は S(α,β) の
  散乱源の直列の評価と表のコピー（#0180）で、確かめていない
- 作り直した記録: 1C の燃焼（examples/depletion の `lattice-1C`、`compare-1C`）と 2A の集合体の燃焼（examples/depletion/p2 の `lattice@2A`、
  `compare@2A`。入力は `ufg-tolerance = 1e-4` のまま）。p2 のライブラリは ADR-0205 の形式で作り直した（90 核種の出力がすべて変わる）ので、
  p2 のほかの集合体（`lattice@2F` など）と、2A の履歴（`p2/lattice-2A.msgpack`）を deps にする examples/depletion/branches のステージは changed になる
- 作り直していない記録（入力が `ufg-scattering` と `ufg-transfer` を書かず、自己遮蔽のある lattice の結果）は、古い既定（静止した標的、
  無限希釈の形）の値で、記録の先頭にそう書いた: examples/depletion（1C と 2A のほか。2A 以外の集合体、MOX-1 など）、branches、build-alignment、
  colorset、core、mox、split、vera（verification.md、verification-endfb70.md）。m0、m1、reflector は自己遮蔽をしない入力だけ、c5g7 は巨視的
  断面積で、UFG を使わない
- 前の測定の試行の入力（examples/depletion/issue-0133 の、`ufg-free-gas` の一部の項目だけを書いたもの）は、省略した項目が今は速い方法の値に
  なる。それらの記録の値は、記録に書いたコミットで再現する。集合体の自由ガスの試行の表は、静止した標的の本番の要約を試行の隣にコピーして比べる
  （`issue-0133/ufg-free-gas-assemblies/compare-<問題>-static.json`）
- 参照（`{ upper = inf, cells = 4, sweeps = inf, heavy = 400, sab = true }`）は重いまま（#0178、#0179）

## 確認

- 入力の試験（`crates/lattice/src/input/tests.rs`）: `ufg_scattering_defaults_to_the_fast_free_gas_with_sab`（省略した項目が速い方法の値、
  `inf` と `heavy = false` で参照の値、`"static"` と `ufg-free-gas` を組めば止める）、`ufg_transfer_defaults_to_the_ring_flux_rows`、
  `static_targets_and_infinite_dilution_reproduce_the_former_defaults`（静止した標的と無限希釈の形を書いた入力の自己遮蔽の条件が、前の既定と
  同じ）。`crates/lattice/src/shielding/tests.rs` の `ring_flux_transfers_are_the_default_and_infinite_dilution_carries_no_rows`。ADR-0207 と
  ADR-0208 の「確認」の `infinite_dilution_transfers_are_the_default_and_carry_no_rows` と `ufg_transfer_defaults_to_the_infinite_dilution_rows` は、
  この 2 つに名前を変えた
- 1C の燃焼（examples/depletion/verification.md「格子計算の燃焼」。xlc f434af311、2026-10-28、Slurm の timed のジョブ 563）: k の差（xlc − OpenMC）は
  1 リングで −61〜+73 → +29〜+136 pcm、4 リングで −10〜+50 → +78〜+126 pcm。xlc の k は 0 の点で +97.7（1 リング）、+97.2 pcm（4 リング）、
  5 MWd/kg で +63.7、+63.2 pcm 動き、0 の点から 5 MWd/kg への変化は +117.6 → +83.6 pcm（1 リング）、+44.9 → +11.0 pcm（4 リング）に縮む。
  自己遮蔽の和は 1 リングで 110.5 → 994.5 秒、4 リングで 284.0 → 1,232.2 秒（前の既定は Slurm 457 の排他のジョブ）、ジョブの経過時間は
  37 分 41 秒
- 2A の集合体の燃焼（examples/depletion/verification.md「集合体の燃焼」の「2A」と「時間」。xlc f434af311、2026-10-28、Slurm の timed のジョブ 565。
  間引いた格子と移行の組で `UfgError::Transfer` で止まらなかった）: k の差は +5〜+203 → +56〜+164 pcm。xlc の k は 0 の点で +54.0 pcm、
  18 MWd/kg で +1.2 pcm、48 MWd/kg で −39.2 pcm 動き、0 の点からの変化は 40 MWd/kg で +158.4 → +67.3 pcm、60 MWd/kg で +169.0 → +81.2 pcm に
  縮む。60 MWd/kg の数密度の差は U-235 +1.27 → +0.92%、Pu-239 +0.01 → −0.01%、Xe-135 +0.84 → +0.74%、Sm-149 +1.92 → +1.82%、ピン出力の
  rms の変化は 0.0006 ポイント以下
