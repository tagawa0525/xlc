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

関係する ADR: ADR-0224

## 完了の条件

`[history]` と `[branches.rod]` を組め（値ごとの形状の燃焼の領域と履歴の領域の対応、値ごとに違う燃焼点（ADR-0200 の制御棒を入れた組の
間引き）の書き方を含む）、径方向の反射体の格子のノードごとの出力で、反射体のノードのテーブルに制御棒の軸を入れない仕組みを作り（ADR-0224）、
examples/depletion と examples/vera/p5 の入力を制御棒の分岐に置き換え、rod_table.py とそのテストを消す。

## 経緯

- 2026-11-03: #0093 で格子計算の制御棒の分岐（燃焼なし）と examples/core の置き換えを行い、履歴からの分岐と残りの例を分けた
- 2026-11-03: #0093 を分岐の部分で閉じた。ADR-0200（213 行目）と ADR-0215（423 行目）の「#0093 まで」の rod_table.py の
  つなぎは、凍結する部分なので直さず、この issue が引き継ぐ
