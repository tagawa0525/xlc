---
status: accepted
date: 2026-11-03
requires: [ADR-0187, ADR-0195, ADR-0200, ADR-0216, ADR-0221, ADR-0224, ADR-0226]
supersedes: []
superseded-by:
issues: ["#0213", "#0214"]
---

# ADR-0225: 燃焼の履歴からの分岐と制御棒の分岐を 1 つの格子計算の入力で組み、値ごとの形状で燃焼の領域を履歴と照合し、値ごとの燃焼点を [history] に書く

## 背景

[ADR-0224](0224-write-rod-branches-as-cell-replacements.md) は、制御棒の分岐を、値ごとに基準の形状のセルを置き換えたセルで書き（`[branches.rod]`）、値ごとに形状、
トラッキング、対称操作、自己遮蔽を作り直して 1 つの格子計算で解くと決めた。燃焼の履歴からの分岐（`[history]`。[ADR-0187](0187-tabulate-burnup-branches-from-history.md)）とは
組まず（入力のエラー）、「燃焼の領域（燃料のセル）は置き換えないので、履歴の領域と値ごとの形状の燃焼の領域の対応は変わらない
見込み（確かめていない。#0213）」と残した。そのため、燃焼度の軸を持つ制御棒の組（[ADR-0200](0200-tabulate-burnup-branches-by-cut-hdmr.md) の「lattice の入力」）は、制御棒を
入れた格子の入力で同じ履歴を読んで別に解き、examples/core/rod_table.py でまとめている（examples/depletion の
`branches/branch-2A-micro-rods.toml` と `branch-2A-micro-boron-rods.toml`）。

決まっていないのは次の 3 つである。

- 値ごとの形状の燃焼の領域と、履歴の燃焼の領域の対応。燃焼の領域（`DepletionRegion`）は、物質、組の代表のセル、組のセルの
  位置（ADR-0226）、体積と、MOC の領域の番号を持つ。案内管のセルを置き換えると、そのセルの帯の数が変わり、燃料の MOC の領域の番号がずれる（#0196。
  `micro_cross_sections_take_the_regions_of_the_branch_lattice_when_the_guide_tubes_differ`）
- 値ごとに違う燃焼点の書き方。ADR-0200 は、公称でない状態（制御棒を入れた組の公称の状態を含む）を一部の燃焼点だけで解けるように
  し、「組ごとに燃焼点を入力で選ぶ」とした。制御棒を入れた組は、測るまですべての燃焼点で解く（#0176）。今の入力の
  `off-nominal-burnup` は 1 つの列で、棒のない組の公称でない状態の燃焼点を選ぶ
- 値ごとの格子で求める量（ミクロ燃焼の均質化した数密度。[ADR-0216](0216-tabulate-micro-depletion-cross-sections-per-state.md)）と、燃焼点を間引いた値の状態の反復を始める解
  （[ADR-0195](0195-start-moc-from-the-previous-transport.md)）を、別の入力の格子計算とどう揃えるか

判断の軸: 制御棒の値の状態が、セルを置き換えた別の入力で同じ履歴を読んだ格子計算と同じ問題を解くこと（ADR-0224 の「値 `aic` の
状態は、`rows` の gt を gt-aic に書き替えた別の入力と同じ問題を解く」）、ADR-0200 の組と燃焼点の選び方を 1 つの入力で書けること、
ADR-0224 の `[branches.rod]` の書き方を変えないこと。

## 検討した案

燃焼の領域の対応:

- 案 R1: 基準の形状の燃焼の領域（MOC の領域の番号）を、すべての値に使う: 案内管のセルの帯の数が値で違うと、燃料の MOC の領域の
  番号がずれ、制御棒の値の状態で燃料でない領域に燃焼の組成を置く
- 案 R2: 値ごとの形状から燃焼の領域を求め直し、[ADR-0226](0226-match-depletion-regions-by-cell-orbits.md) の比べ方（物質、組の代表のセル、組のセルの位置、体積、MOC の
  領域の数。`region_mismatch`）で履歴と同じことを確かめ、MOC の領域の番号はその形状のものを使う（採用）。燃焼の輸送と違う格子の
  分岐で、MOC の領域の番号を入力の格子から求めている（ADR-0187、ADR-0226、#0196）のと同じ扱い。燃料のセルを置き換える値や、
  燃料のセルの組（対称性）を変える値はエラーになる

値ごとの燃焼点の書き方:

- 案 B1: 制御棒の値のすべての状態を `off-nominal-burnup` で解き、値ごとの列を持たない: ADR-0200 の公称でない状態の定義（制御棒を
  入れた組の公称の状態を含む）とは合うが、棒のない組の公称でない状態を間引き、制御棒を入れた組はすべての燃焼点で解く形（ADR-0200 の
  「測るまで」の形）が書けない
- 案 B2: `[branches.rod]` の値に燃焼点を書く（`aic = { cells = { gt = "gt-aic" }, burnup = [...] }`）: `[branches.rod]` は履歴のない
  分岐（M3 のテーブル）にも使い、燃焼点は履歴からの分岐にしかない。ADR-0224 の書き方（値 = 置き換えるセル）を変える
- 案 B3: `[history]` に値ごとの燃焼点 `rod-burnup = { aic = [...] }` を足し、書かない値は `off-nominal-burnup`（採用）。燃焼点の選び方は
  すべて `[history]` に置き、`[branches.rod]` は変えない。既定は ADR-0200 の公称でない状態の定義（制御棒を入れた組の公称の状態は
  公称でない）に従う
- 案 B4: 値ごとの別の入力で解き、rod_table.py でまとめる（今の形）: ADR-0224 の検討した案の最初の案と同じく、[ADR-0097](0097-tabulate-branches-on-a-full-grid.md) が採らないと
  決めている

燃焼点を間引いた値の状態の反復を始める解:

- 案 S1: 今の規則のまま（ADR-0195 の、添字が 0 でない最後の軸を 1 つ戻した状態を解いていなければ、直前に解いた状態の解。ADR-0200）:
  制御棒の値が燃焼度 0.5 を間引くと、燃焼度 2.0、ホウ素 600 ppm の状態は、1 つ戻した燃焼度 0.5 の状態を解かないので、直前に解いた
  燃焼度 0、ホウ素 1300 ppm の状態の解から始める。別の入力（燃焼点 0、2.0）は燃焼度 0、ホウ素 600 ppm の状態から始めるので、収束の
  誤差の程度で別の入力と違う値になる（下の結合テストの格子で、値の燃焼度 2.0 の状態の少数群定数がビットで一致しなかった。
  2026-11-03、手元の実行）
- 案 S2: 戻す点を、前の軸が同じ解いた状態がその軸に持つ点のうち 1 つ前のものにする（採用）。制御棒の値の燃焼度の軸では、その値の
  1 つ前の燃焼点になり、別の入力と同じ始める解になる。前の軸が同じ状態がその軸のどの点も持つ形（直積、ADR-0200 の部分格子、棒の
  ない値の燃焼点。棒のない値の公称の状態はすべての燃焼点で解く）では 1 つ前の点で、今の規則と同じ

## 決定と理由

`[history]` と `[branches.rod]` を同じ格子計算の入力で組む。値ごとに、基準の形状のセルを置き換えた形状（ADR-0224）で、同じ燃焼度の
点の履歴の組成を解き、テーブルは燃焼度の軸を最初に、制御棒の軸を最後に持つ。

```toml
[history]
file = "../p2/lattice-2A.msgpack"
burnup = [0, 0.1, 0.5, 10, 20, 40, 60]
off-nominal-burnup = [0, 0.1, 0.5, 20, 60]            # 棒のない値の公称でない状態（ADR-0200）
rod-burnup = { aic = [0, 0.1, 0.5, 10, 20, 40, 60] }  # 値 aic のすべての状態を解く燃焼点（書かない値は off-nominal-burnup）

[branches.rod]
aic = { gt = "gt-aic" }
```

- 燃焼の領域（案 R2）: 値ごとの形状とその形状で選んだ対称操作から燃焼の領域を求め、ADR-0226 の比べ方（物質、組の代表のセル、
  組のセルの位置、体積（ビット）、MOC の領域の数）で履歴の燃焼の領域と同じことを確かめる（違えば値の名前を添えたエラー）。
  一部の案内管にだけ棒を入れて燃料のセルの組を分ける値は、組のセルの位置が違って止まる（値ごとに対称性を揃えて履歴の組成を
  分けて当てることはしない。履歴の組の数密度は組の平均で、組を分けた位置ごとの組成を持たない）。燃焼の組成を置く MOC の領域、自己遮蔽と MOC に渡す物質の名前
  （`<物質>@<セル>`）、微視的断面積を重み付ける領域は、その値の形状のもの
- 値ごとの燃焼点（案 B3）: `rod-burnup` は `[branches.rod]` の値（`out` でない）ごとに、`burnup` から選んだ増える列で、最初と最後の点を
  含む（core は公称の状態との差を外挿しない。ADR-0200）。棒のない値は今のとおり、公称の状態をすべての燃焼点で、ほかの状態を
  `off-nominal-burnup` で解く。制御棒の値は、公称の状態を含むすべての状態を、その値の `rod-burnup`（書かなければ `off-nominal-burnup`）で
  解く。どれかの値が燃焼点を間引けば、公称の状態（入力の物質が軸の格子点にあること）が要る（ADR-0200 の検査を値ごとの燃焼点に
  広げる）
- 均質化した数密度（ADR-0216）: 値の微視的断面積は、その値の格子の体積で均質化した数密度で割る。格子の体積は値の形状の領域の体積の
  和で、帯の数が違うと丸めの程度で違う。テーブルの燃焼度の点ごとの数密度は棒のない値のもの（rod_table.py と同じ）。値の状態の
  数密度の微分（[ADR-0221](0221-correct-micro-linearization-by-number-density-derivatives.md)）も、その値の数密度で割り、摂動の組成の自己遮蔽は値の形状ごとに持つ
- 反復を始める解（案 S2）: ADR-0195 の、添字が 0 でない最後の軸を戻す点を、前の軸（制御棒の値を含む）が同じ解いた状態がその軸に持つ
  点のうち 1 つ前のものにする。その状態を解いていなければ、今のとおり直前に解いた状態の解から始める
- ADR-0224 の、値ごとの解く順（`out` から）、値の最初の状態を平らな初期値から始めること、チェックポイントの番号、標準出力の状態の
  番号（テーブルの軸の直積の順）は、履歴からの分岐でも同じ

これで、制御棒の値の状態は、セルを置き換え、燃焼点をその値の燃焼点にした別の入力で同じ履歴を読んだ格子計算と同じ問題を、同じ
順と同じ始める解で解き、ADR-0200 の「組ごとに燃焼点を入力で選ぶ」を 1 つの入力で書ける。

### 帰結

- 燃焼度の軸を持つ制御棒の組は 1 つの入力で作れる。examples/depletion の制御棒の組（branch-2A-micro-rods.toml と
  branch-2A-micro-boron-rods.toml を rod_table.py でまとめたテーブル）は、1 つの入力（`branch-2A-micro-rod-axis.toml`、
  `branch-2A-micro-boron-rod-axis.toml`）に置き換えた。置き換えた入力のテーブルはまだ作り直しておらず、記録の数値は前の形のもの（#0213）
- 置き換えた `branch-2A-micro-rod-axis.toml` は、棒のない組の入力（branch-2A-micro.toml）の数密度の微分（ADR-0221 の `derivative`）を
  制御棒の値にも持つ。前の制御棒の組の入力（branch-2A-micro-rods.toml）は `derivative` を持たず、rod_table.py は持つ値と持たない
  値をそのまま結合していた。branch-2A-micro-rods.toml は、部分挿入のノードの確認の一様に棒の入った層の比べる先として残し、値 aic と
  同じ問題になるよう `derivative` を足した
- rod_table.py は、examples/vera/p5 の HZP のテーブルを 1 つの入力の格子に置き換え（ノードごとに制御棒の軸を付け分ける仕組みは
  [ADR-0227](0227-give-the-rod-axis-to-nodes-with-replaced-cells.md)）、確認の手順（ADR-0224 の examples/core/rod-branch-check、この ADR の examples/depletion/rod-axis-check）とともに
  消した（#0213、2026-11-04）。この ADR の確認の計算の入力と手順は、タグ `record/adr-0225-rod-axis-check` にある
- 制御棒の履歴の補正（[ADR-0215](0215-correct-untracked-history-by-history-branches.md) の制御棒の変数。#0190）は、この組み合わせの上に作る別のトピック

### 確認

- 入力: `history_branches_take_the_rod_axis`、`rod_burnups_are_history_burnups_of_rod_values`（crates/lattice/src/branches/tests.rs）
- 燃焼の領域と燃焼点: `rod_values_take_the_history_densities_in_the_regions_of_their_own_geometry`、
  `rod_values_are_solved_at_their_own_burnups`、`rod_values_must_keep_the_depletion_regions_of_the_history`、
  `rod_values_that_change_the_fuel_orbits_are_errors`（crates/lattice/src/history/tests.rs）
- 反復を始める解: `states_step_back_to_the_previous_point_their_value_solves`、`states_outside_the_subgrids_start_from_the_last_solved_state`
  （crates/lattice/src/branches/tests.rs）
- 別の入力との一致: `history_rod_branches_match_separate_lattices`（crates/lattice/tests/history_rod_branch.rs）。合成のライブラリの
  2 × 2 の格子で燃焼させた履歴から作った、制御棒の値の燃焼点を間引いた 1 つの入力のテーブルの各状態（少数群定数、追う核種の微視的
  断面積と数密度の微分）が、棒のない入力と、案内管を置き換えてその値の燃焼点だけを読む入力の格子計算の状態とビットで一致する
- VERA 2A（自己遮蔽あり、32 スレッド）: 燃焼度 0、10、40 MWd/kg（値 aic は 0 と 40）の 1 つの入力のテーブルと、値ごとの 2 つの
  入力を rod_table.py でまとめたテーブルは、軸、点、少数群定数、微視的断面積と数密度の微分のすべての項目でビットまで等しい
  （examples/depletion/verification.md の「制御棒の分岐と値ごとの格子計算の結合の一致（燃焼の履歴からの分岐）」、2026-11-03）。
  計算は ADR-0226 の前のコード（b588ab8d。タグ `record/adr-0225-rod-axis-check`、分岐点 3dda2078）と前の形式の燃焼の履歴
  （md5 7efea664）で流し、ADR-0226 の値ごとの `region_mismatch` の経路は通っていない（単体テストが通す）。今のコードは使った履歴を
  読めない（再現の手順は記録の「計算」）
