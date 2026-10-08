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

SPH 因子（Hébert の SPH 等価）で直す。リング r・群 g の反応に掛ける因子 μ を、同じセルを μ を掛けた多群の群定数で
解いた中性子束が、UFG のリングごとの群の中性子束 φ_ref に比例するよう反復で求める。再現できるのは群ごとのリング間の
形で、大きさは、因子を掛けるリングの平均の中性子束を保つよう規格化して決める。群の釣り合い（反応率の和 = 上の群からの
線源）は μ では変えられず、多群の散乱行列は UFG と同じ群間の移行を与えないためである（`crates/lattice/src/sph.rs` の
冒頭の説明）。M2 の計画は「同じセルを多群で解いた反応率が UFG と一致するよう求めた SPH 因子」と書いていたが、一致する
のはリング間の形と平均の中性子束までである。

因子は、UFG の範囲の群で、等価なセルの真空でないすべてのリング（燃料、被覆管、減速材）の物質に掛ける
（`crates/lattice/src/shielding.rs` の `shield_cell`）。SPH 因子は UFG の中性子束を参照にして求めるので、参照のない範囲の外（[ADR-0074](0074-shield-outside-range-by-equivalence.md)）には掛けない。真空のリングと、因子を
掛ける先のないリング（格子に使わない等価なセルの枠のリング。[ADR-0088](0088-add-frame-ring-to-equivalent-cell.md)）は μ = 1 とする。

## 確認

`factors_reproduce_the_reference_shape`、`consistent_reference_needs_no_correction`、`rejects_references_without_flux`、
`rejects_rings_without_flux`。VERA のピンセルの結果（examples/vera/README.md）。

## 結果

- 範囲の外の群では、どのリングにも掛けない（[ADR-0074](0074-shield-outside-range-by-equivalence.md)）
- 枠のリングの SPH 因子は 1 のままにする（[ADR-0088](0088-add-frame-ring-to-equivalent-cell.md)）
