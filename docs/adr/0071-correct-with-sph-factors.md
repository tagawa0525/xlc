---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: []
---

# ADR-0071: 等価なセルの減速計算の反応率に合わせる SPH 因子で、MOC の角度の扱いの違いを直す

## 背景

ピンごとの減速計算（[ADR-0068](0068-solve-pinwise-slowing-down.md)）はリングごとのスカラー中性子束の重みで群定数を作る。この群定数を MOC で使うと、角度の
扱いの違いによる誤差が出る（文献では約 −200 pcm。docs/plans/m2-assembly.md、2f127db、2026-10-07）。

## 検討した案

計画には、SPH 因子以外の補正を比べた記録はない。

## 決定と理由

同じセルを多群で解いた反応率が UFG と一致するよう求めた SPH 因子で直す。SPH 因子は UFG の中性子束を参照にして求める
ので、参照のない範囲の外（[ADR-0074](0074-shield-outside-range-by-equivalence.md)）には掛けない。

## 確認

`factors_reproduce_the_reference_shape`、`consistent_reference_needs_no_correction`、`rejects_references_without_flux`、
`rejects_rings_without_flux`。VERA のピンセルの結果（examples/vera/README.md）。

## 結果

- 等価なセルの被覆管と減速材には掛けない（[ADR-0074](0074-shield-outside-range-by-equivalence.md)）
- 枠のリングの SPH 因子は 1 のままにする（[ADR-0088](0088-add-frame-ring-to-equivalent-cell.md)）
