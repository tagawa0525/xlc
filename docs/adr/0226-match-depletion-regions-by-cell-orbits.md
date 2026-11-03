---
status: accepted
date: 2026-11-03
requires: [ADR-0173, ADR-0187, ADR-0215, ADR-0224]
supersedes: []
superseded-by:
issues: ["#0190", "#0196"]
---

# ADR-0226: 続ける燃焼と履歴からの分岐は、燃焼の領域を物質、組のセルの位置、体積、MOC の領域の数で比べ、MOC の領域の番号は比べない

## 背景

格子計算の燃焼の領域（crates/lattice/src/depletion.rs の `DepletionRegion`）は、燃焼させる物質と、対称操作で移り合うセルの組
（ADR-0173）ごとに 1 つで、組の代表のセル（組の最小の番号）、その物質の MOC の領域、体積の和を持つ。燃焼の履歴
（`xlc-depletion-history`）はこれを書き、次の 2 つの経路が、入力の格子から求め直した燃焼の領域を履歴のものと照らし、履歴の
数密度（燃焼の領域と核種の 2 次元の列）を領域の順に当てる。

- 履歴からの分岐（ADR-0187。history.rs の `history_branches`）: 制御棒を入れた組（ADR-0200、ADR-0224）などで、履歴の格子と
  案内管の中身が違う格子を解く
- 途中で条件を変えた燃焼（ADR-0215 の「確認」。depletion.rs の `continued_start`、入力の `[depletion]` の `continue`）: Pyrex を
  抜いた燃焼の案 Q2 と制御棒の履歴（#0190）で、案内管の中身を変えて続ける

`continued_start` は燃焼の領域を MOC の領域の番号まで比べていたので、案内管のセルの領域の数が違う格子は、燃料のセルが同じでも
止まった（#0190 の 2026-10-30 の項）。`history_branches` は（物質、組の代表のセル、体積）で比べ、番号は入力の格子のものを
使っていた（#0196）。ADR-0215 は「燃焼の領域…が履歴と違うならエラー」と書くが、何が同じなら同じ燃焼の領域かを決めていない。

組の代表のセルと体積で比べる形には反例がある（2026-11-03、#0190 の 2026-11-03 の項。今はテスト
`continuing_needs_the_cells_of_each_depletion_region_of_the_history` の格子）。2 × 2 の格子で、履歴を `rows = ["pin gt", "pin gt"]`（燃料の組 {0, 2}、MOC の領域 [0, 4]）で燃焼させ、
`rows = ["pin gt", "gt pin"]`（組 {0, 3}、[0, 6]）で続けると、組の代表 0、体積、MOC の領域の数が同じなのでエラーにならず、
セル 2 で燃焼させた組成をセル 3 に置いて続けた。燃料の組成を別の位置に当てる誤りを、黙って出す経路である。

決まっていないのは、2 つの格子の燃焼の領域を同じとみなす条件と、そのために履歴が持つものである。判断の軸は、怪しい入力を
推測せずにエラーで止めること（CLAUDE.md の「方針」）、案内管の中身を変えた続き（Q2、制御棒）を止めないこと、2 つの経路で
同じ条件を使うこと。

## 検討した案

- MOC の領域の番号まで比べる（前の `continued_start`）: 組のセルの違いも止めるが、案内管のセルの領域の数が違う格子は、燃料の
  セルが同じでも番号がずれて止まる。Q2 と制御棒の続きが流せない。番号は燃料のセルの外の離散化（案内管のリングの数）で変わり、
  燃焼の領域の対応とは関係がない
- 物質、組の代表のセル、体積で比べる（前の `history_branches`、#0196）: 上の反例（組の代表と大きさが同じで、組のセルの違う
  格子）を通す。燃料のリングの数を変えた格子は、体積が MOC の領域の体積の和なので、組に複数のセルがあるとリングの分け方で
  和の丸めが変わり、通るか止まるかが丸め次第になる（このブランチの最初の版のレビュー、2026-11-03）
- 代表のセルで比べ、組のセルの違う格子は使う側の運用で避ける（Q2 と制御棒の入力で対称性が変わらないことを、流す前に確かめる）:
  形式を変えずに済むが、確かめ忘れた入力は誤った組成で黙って流れる。強制する対称性（`symmetry`）や、Pyrex の配置で、対称性の
  位数が同じで別の操作になる格子は、入力を読むだけでは気づきにくい
- 履歴に組のセルの番号を書き、物質、組のセルの番号の集合、体積（ビット）、MOC の領域の数で比べる（このブランチの 2 つ目の版）:
  形の違う格子では同じ番号のセルが別の位置にあり、対称性を使わない（`symmetry = "none"`）2 × 2 の格子の履歴を 1 × 4 の格子で
  続けると、番号 0〜3、体積、MOC の領域の数が同じで通り、セル 2 の組成を 2 行目の 1 列目から 1 行目の 3 列目に移した（3 周目の
  レビュー、2026-11-03。テスト `continuing_needs_the_positions_of_the_cells_of_the_history`）
- 履歴に組のセルの位置（[行, 列]）を書き、物質、組の代表のセル、組のセルの位置、体積（ビット）、MOC の領域の数で比べる（採用）

## 決定と理由

燃焼の領域は、物質、組の代表のセル（番号）に加えて組のセルの位置（[行, 列]、上から、左から。行の順）を持ち、燃焼の履歴に
`cells` として書く。続ける燃焼（`continued_start`）と履歴からの分岐（`history_branches` と、履歴の組成で自己遮蔽する
`History::shield_with`）は、入力の格子から求めた燃焼の領域と履歴の燃焼の領域を、順に（物質、組の代表のセル、組のセルの位置、
体積（ビット）、MOC の領域の数）で比べ、1 つの関数
（depletion.rs の `region_mismatch`）を使う。MOC の領域の番号は比べず、入力の格子のものを使う。どれかが違うか、燃焼の領域の
数が違えばエラーで止め、最初に違う燃焼の領域を両方の値で示す。

- 組のセルの位置を比べるので、組の代表と大きさが同じで組のセルの違う格子（上の反例）と、番号が同じで位置の違う形の格子を
  止める。燃料のセルの組が同じなら、案内管の中身を変えた格子（Q2、制御棒）は通る
- MOC の領域の数を比べるので、燃料のリングの数を変えた格子は、体積の和の丸めによらず止まる。燃焼の領域の中の離散化が違う格子は、
  同じ燃焼の領域とみなさない
- 体積はビットで比べる。同じセルの同じ帯なら、体積は解析的な面積（`Pin::volumes`）の同じ順の和で、ビットで同じになる
- 読むときに、組のセルの位置が空でなく、行の順に増えて重ならないことを確かめる。組の代表の番号は、格子の列の数を履歴が持たない
  ので、位置と照らさない
- 後方互換を作らない（CLAUDE.md の「方針」）。`cells` を持たない前の履歴は読めない（`DepletionRegion` は欠けたキーをエラーに
  する）。燃焼の履歴を埋め込む lattice のチェックポイントも同じ。Python のツールは履歴を map のキーで読むので読めるが、`cells` の
  有無で結果が変わるものがある（帰結）

### 帰結

- 燃焼の履歴の形式が変わり、lattice の実行ファイルが変わる。前の履歴は lattice（`[history]`、`continue`、`shielding-file`）では
  読めない。DVC で追跡する燃焼の履歴は 23 本（2026-11-03 に数えた）: examples/depletion/dvc.yaml の 4 本（lattice-1C-rings1、
  lattice-1C-rings4、lattice-2A、lattice-mox1-rings1）、examples/depletion/p2/dvc.yaml の 8 本（2A、2F、2M、2N、2O、2P、2A-density、
  2A-rods）、`.dvc` の 7 本（branches/lattice-2A-6272e937、issue-0133/ufg-free-gas-assemblies/lattice-2A-fg と lattice-2O-fg、
  history/lattice-residual-density、lattice-switch-density、lattice-weight-0661、lattice-weight-0702）、make_hybrid.py が履歴から
  作る history/hybrid-*.msgpack の 4 本。作り直しは、それを lattice で読む作業が来るまで後回しにする（CLAUDE.md の「計算」）
- 作り直す時間の見込み（main の作業ツリーのステージのログの始めと終わりの時刻の差。割り込まれたかは確かめていない。2026-10-23 の
  ログは、その後の高速化の前の版）: dvc.yaml の 4 本の和 1335 秒（1C-rings1 449、1C-rings4 648、2A 147、mox1 91）、p2 の 8 本の和
  29021 秒（2A 2221（2026-10-28）、2F 2195、2M 10661、2N 5194、2O 2360、2P 2945、2A-density 1106、2A-rods 2339）、history の
  residual-density-history 2489 秒、switch-density-history 972 秒。weight-history-0661 と 0702、lattice-2A-6272e937、lattice-2A-fg と
  2O-fg のログは見つけておらず、2A の燃焼と同じ程度（各 2221〜2489 秒）と置いた（確かめていない）。hybrid の 4 本は make_hybrid.py
  で数秒。合計はおよそ 9.4 時間（ログのあるもの）と、見込みの 5 本で 3.1〜3.5 時間
- 分岐のテーブル（`xlc-table`）と core の履歴は燃焼の領域を持たないので、形式は変わらない
- `cells` の有無で結果が変わる Python のツール（2026-11-03 に見つけたもの）: examples/depletion/p2/compare_checkpoint.py は燃焼の
  領域の map を丸ごと比べるので、前の形式と今の形式の履歴を同じとしない。examples/depletion/micro/make_shifted.py は、前の履歴からは
  前の形式の履歴を書き、lattice が読むときに止まる
- examples/depletion/history/make_hybrid.py の `identity` も同じ項目（物質、組の代表のセル、組のセルの位置、体積、MOC の領域の数）で
  2 つの履歴を比べ、組のセルを持たない前の履歴を止める（ハイブリッドを作る経路で、組成を別の位置に置かないため）

### 確認

- crates/lattice/src/depletion/tests.rs: `continuing_with_other_guide_tube_contents_keeps_the_fuel_regions_of_the_history`（案内管の
  中身を変えて続く）、`continuing_needs_the_fuel_materials_cells_volumes_and_rings_of_the_history`、
  `continuing_needs_the_cells_of_each_depletion_region_of_the_history`（上の反例）、`histories_write_the_cells_of_each_depletion_region`、
  `continuing_needs_the_positions_of_the_cells_of_the_history`、`histories_whose_region_cells_are_not_increasing_positions_are_errors`
- examples/depletion/history/test_make_hybrid.py: `test_histories_with_other_depletion_regions_are_errors`、
  `test_histories_without_the_cells_of_the_regions_are_errors`
- crates/lattice/src/history/tests.rs: `the_branch_lattice_must_hold_the_cells_and_rings_of_the_history_regions`、
  `micro_cross_sections_take_the_regions_of_the_branch_lattice_when_the_guide_tubes_differ`（#0196）
