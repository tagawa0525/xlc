---
status: open
labels: [lattice, shielding]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# セルをまたぐ構造材の塊（バッフル）の自己遮蔽

## 問題

構造材の塊の自己遮蔽（ADR-0086）は、セルの中で接する、核種の組が同じ領域を 1 つの塊にする
（`crates/lattice/src/shielding.rs` の `shield_lumps`）。VERA 5-2D の反射体の格子（examples/vera/p5）で、SS304 の
バッフル（厚さ 2.85 cm）を examples/reflector と同じく 0.95 cm の 3 列のセルに分けて自己遮蔽をすると、中の列の
セルは上下左右を SS304 のセルに囲まれ、`lump 1 has no face open to non-black regions` で止まる（b89a69e8 の
`reflector-I.toml`、2026-10-12）。

5-2D では、バッフルを厚さ 2.85 cm の 1 列のセルにし、セルの中を `grid = 3` で分けて解いた（examples/vera/p5/p5_inputs.py）。
このとき塊は「ピンの行の高さ（1.26 cm など）× 2.85 cm」のバッフルの切れ端で、上下の隣のセルのバッフルは黒体に
せず、無限希釈の全断面積の影として扱う（ADR-0086 の、隣のセルとつながる枠と同じ近似）。物理的な塊はバッフル
全体で、この近似が実効断面積と反射体の定数にどれだけ効くかは確かめていない。

関係する ADR: なし

## 完了の条件

セルをまたいでつながる構造材の領域を 1 つの塊にする扱いの要否を、5-2D の反射体の格子で、今の近似との差（バッフルの
実効断面積、反射体の定数、炉心の k と外周の出力）を測って決める。

## 経緯

- 2026-10-12: 5-2D の反射体の格子を自己遮蔽ありで解こうとして分かった（feat/p5-2d）
