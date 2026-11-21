---
type: Issue
status: open
labels: [lattice, sensitivity]
opened: 2026-11-22
closed:
milestone:
blocked-by: []
resolution:
---

# 一般化随伴の少数群定数の感度が、同じビルドの 2 回の実行で 1e-14 の桁で揺れる

## 問題

examples/depletion/issue-0279/lattice-3x3-few-group.toml を核種 3 つにし、`direct` と `shielded` を除いた入力（3 × 3 の格子、
`[sensitivity.few-group]`）を、同じ lattice のビルド（main 630459cba）で続けて 2 回流すと、感度の 60 個の値の相対差の最大が 1.3×10⁻¹⁴ になる
（2026-11-22、手元、32 スレッド。ADR-0275 の確認の「一般化随伴」。入力はコミット 56e4978b9 の examples/depletion/trial-0241-gen/）。MOC の
固有値計算の出力（VERA 2A などの標準出力、ピン出力、少数群定数）は、ビルドが同じなら 2 回の実行でビットで同じだった（同じ確認）。

原因は確かめていない。候補は、一般化随伴の反復か GMRES の中の、スレッドの数や仕事の分け方で順が変わる並列の和。

読み返す場面: 出力が変わらないことで再実行を省く仕組み（#0150 の、結果の変わらないステージの付け直し）と、感度の記録を再実行で
確かめるとき。揺れがあると、コードを変えていなくても出力が「変わった」と見える。

関係する ADR: [ADR-0275](../adr/0275-solve-lattice-coarse-mesh-eigenvalue-with-xlc-cmfd.md)

## 完了の条件

- [ ] 揺れの出どころを最小の入力で特定し、スレッドの数によらずビットで同じ結果にするか、揺れを許す理由と大きさを記録する

## 経緯

- 2026-11-22: #0241（ADR-0275）の前後の比較で、同じビルドの 2 回の差として見つけた
