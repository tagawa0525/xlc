---
status: open
labels: [lattice, core]
opened: 2026-11-03
closed:
milestone: M4
blocked-by: []
resolution:
---

# rod_table.py を、燃焼の履歴からの分岐の制御棒の軸に置き換えて消す

## 問題

格子計算は制御棒の分岐を 1 つの入力で解けるようになった（ADR-0224、#0093）が、燃焼の履歴からの分岐（`[history]`。ADR-0187）と
制御棒の軸は組めない（入力をエラーにする）。次がまだ examples/core/rod_table.py で値ごとのテーブルをまとめている。

- examples/depletion の制御棒の組（`branches/branch-2A-micro-rods.toml`、`branch-2A-micro-boron-rods.toml` とそれをまとめた
  テーブル。ADR-0200 の「lattice の入力」の rod_table.py の項）。燃焼度の軸、ミクロ燃焼の微視的断面積（ADR-0216）、履歴の補正
  （ADR-0215 の制御棒の変数。#0190）と絡む
- examples/vera/p5 の HZP のテーブル（`p5_inputs.py` が書く dvc.yaml の 8 本。燃料と支持格子の層、I 型の反射体の格子の燃料）。
  燃焼のない M3 のテーブル。燃料と支持格子の層は今の格子計算の制御棒の分岐で書けるが、I 型の反射体の格子はノードごとの出力で、
  制御棒の分岐と組む入力はエラーになる（反射体のノードのテーブルは制御棒の軸を持てない。ADR-0127、ADR-0224）。ノードごとに
  制御棒の軸を持たせる扱いが要る

関係する ADR: ADR-0224、ADR-0225

## 完了の条件

`[history]` と `[branches.rod]` を組め（値ごとの形状の燃焼の領域と履歴の領域の対応、値ごとに違う燃焼点（ADR-0200 の制御棒を入れた組の
間引き）の書き方を含む）、径方向の反射体の格子のノードごとの出力で、反射体のノードのテーブルに制御棒の軸を入れない仕組みを作り（ADR-0224）、
examples/depletion と examples/vera/p5 の入力を制御棒の分岐に置き換え、rod_table.py とそのテストを消す。

## 経緯

- 2026-11-03: #0093 で格子計算の制御棒の分岐（燃焼なし）と examples/core の置き換えを行い、履歴からの分岐と残りの例を分けた
- 2026-11-03: #0093 を分岐の部分で閉じた。ADR-0200（213 行目）と ADR-0215（423 行目）の「#0093 まで」の rod_table.py の
  つなぎは、凍結する部分なので直さず、この issue が引き継ぐ
- 2026-11-03: 前半（`[history]` と `[branches.rod]` の組み合わせ、値ごとの燃焼点 `rod-burnup`、値ごとの形状の燃焼の領域の照合）を
  ADR-0225 で行い、examples/depletion の制御棒の組の入力を 1 つの入力（`branches/branch-2A-micro-rod-axis.toml`、
  `branch-2A-micro-boron-rod-axis.toml`）に置き換え、`branch-2A-micro-boron-rods.toml` を消した。`branch-2A-micro-rods.toml` は、
  core-2A-micro-rodded の一様に棒の入った層の比べる先（制御棒の軸を持たないテーブル）として残し、値 aic と同じ問題になるよう
  Pu-240 の数密度の微分（`derivative`）を足した
- 2026-11-04: p2 の 2A の履歴を今の形式で作り直し（chore/rerun-after-adr-0226、Slurm 664）、置き換えた入力のテーブルと棒のない組、
  制御棒を入れた組を今の bin で作り直して `dvc add` した（Slurm 673〜676。md5 は `table-2A-micro.msgpack` a0776f62（前の 646d0ec2）、
  `table-2A-micro-rods.msgpack` dd5ed4fe（前の 70b204f3）、`table-2A-micro-rod-axis.msgpack` a33f1edf（前のまとめた形 9689f898）で、
  `table-2A-micro-boron-rod-axis.msgpack` 283e440b だけが前のまとめた形とバイトまで同じ）。
  記録の「部分挿入のノードのミクロ燃焼（2A の集合体 1 体）」と「ステップの大きさの誤差と中点の臨界ホウ素（小さな炉心）」の手順と
  数値を、作り直したテーブルの core の燃焼（Slurm 677）で書き直した。部分挿入のノードの k は、前の記録から 52 MWd/kg で棒の入った燃焼が
  +100.9 pcm（0.73414963 → 0.73515853）、棒のない燃焼が +0.1 pcm（0.79047490 → 0.79047634）動いた。候補は、値 aic が数密度の微分を
  持つようになったこと（棒の入った燃焼）と、前のテーブルを作った bin（87f5dd7b）からの bin の変更（両方）で、分けていない

## 残り

- examples/depletion/verification.md の「サイクルのステップの core の燃焼（2A の集合体 1 体）」は前の形のテーブル（646d0ec2、
  70b204f3、まとめたもの 9689f898）、「core のミクロ燃焼（2A の集合体 1 体）」は前の棒のないテーブル（be25725e）の数値のまま。
  作り直したテーブル（経緯の 2026-11-04）で core を再実行して書き直す。「制御棒の分岐と値ごとの格子計算の結合の一致」も、今の形式の
  2A の履歴（`p2/lattice-2A.msgpack` a2cc87df）で再実行して書き直す（下の項）。記録のほかの節にも、2A の履歴を前の md5 7efea664 で
  引く文が残る
- 径方向の反射体の格子のノードごとの出力で、反射体のノードのテーブルに制御棒の軸を入れない仕組み（ADR-0224。#0214 の
  ノードごとの付け分けの項）と、examples/vera/p5 の HZP のテーブルの置き換え
- rod_table.py とそのテスト（examples/core/test_rod_table.py）の削除と、直下の `dvc.yaml` の `python-ast` の foreach からの削除
- 前半のレビューの 3 周目（2026-11-03）で、直さずに残した指摘: examples/depletion/verification.md の「制御棒の分岐と値ごとの
  格子計算の結合の一致（燃焼の履歴からの分岐）」と ADR-0225 の確認の計算は、分岐点 dcc9036a の上の ADR-0226 の前のコード
  （a0de3041。タグ `record/adr-0225-rod-axis-check`）で流した。今のコードでは使った燃焼の履歴（前の形式）を読めず、2A の比較は
  ADR-0226 の `region_mismatch` を値ごとに通す経路を通っていない（単体テスト `rod_values_that_change_the_fuel_orbits_are_errors` が
  通す）。テーブルを作り直すときに、今の形式の履歴で比較を再実行して記録を書き直す。core-2A-micro-rods-assemblies.toml と
  core-small-assemblies.toml の説明は新しい入力を出どころに書くが、今あるテーブルは rod_table.py でまとめた前のもの（テーブルを
  作り直せば合う）。examples/core/rod_table.py の説明は「#0213 まで」のつなぎと書くが、今は p5 の 8 本のテーブル（燃料と支持格子の層の 6 本と、
  反射体の格子の 2 本）、確認の手順（examples/core/rod-branch-check、examples/depletion/rod-axis-check）、examples/depletion の記録の
  前の形のテーブルの手順（「部分挿入のノードのミクロ燃焼」と小さな炉心の節）、examples/core/verification.md の「方式の誤差と
  微分価値」の前の形の説明、examples/core/dvc.lock に残る前の形の `table` のステージの記載（#0214）が使うか指す（rod_table.py を
  消すときに、これらの手順と記載も書き直す）
