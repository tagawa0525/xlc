---
status: proposed
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
  精度の低い近似になる（ADR-0205 の表で、1C の新しい燃料の k を S(α,β) の参照から −12.5〜−16.8 pcm（自由ガスだけ）、ADR-0203 の前の静止した
  標的はさらに離れる）。近似を選んだことが入力に現れない
- 案 D2: 既定を参照そのもの（`{ upper = inf, cells = 4, sweeps = inf, heavy = 400, sab = true }`）にする。参照は全範囲の格子点でスイープを
  6 回ほど繰り返し、ピンセルの 1 回の自己遮蔽で数分、集合体の燃焼で数時間かかる（ADR-0203 の帰結。`--timed` では測っていない）。
  速い方法との差は −4.3、−3.0 pcm（ADR-0205）で、集合体の燃焼の他の誤差より小さい
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
  しない。前の既定の結果は、入力に `ufg-scattering = "static"` と `ufg-transfer = "infinite-dilution"` を書けばビットで同じになる（`ShieldingOptions`
  が前の既定と同じになり、解く側のコードは変わらない）
- 作り直した記録: 1C の燃焼（examples/depletion の `lattice-1C`、`compare-1C`）と 2A の集合体の燃焼（examples/depletion/p2 の `lattice@2A`、
  `compare@2A`。入力は `ufg-tolerance = 1e-4` のまま）。p2 のライブラリは ADR-0205 の形式で作り直した
- 作り直していない記録（入力が `ufg-scattering` と `ufg-transfer` を書かず、UFG の自己遮蔽を使う lattice の結果）は、古い既定（静止した標的、
  無限希釈の形）の値で、記録の先頭か該当の節にそう書いた。examples/depletion の集合体の 2A 以外、MOX-1、分岐、core、examples/vera、m1、mox など
- 前の測定の入力（examples/depletion/issue-0133 の `ufg-free-gas` の一部の項目だけを書いた試行の入力）は、省略した項目が今は速い方法の値に
  なる。それらの記録の値は、記録に書いたコミットで再現する
- 参照（`{ upper = inf, cells = 4, sweeps = inf, heavy = 400, sab = true }`）は重いまま（#0178、#0179）

## 確認

- 入力の試験: `crates/lattice/src/input/tests.rs` の `ufg_scattering_defaults_to_the_fast_free_gas_with_sab`、
  `ufg_transfer_defaults_to_the_ring_flux_rows`、`static_targets_and_infinite_dilution_reproduce_the_former_defaults`
