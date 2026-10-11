---
status: accepted
date: 2026-10-08
supersedes: []
superseded-by:
issues: ["#0006", "#0007"]
requires: []
---

# ADR-0086: 燃料のピンの外の構造材の塊は、黒体の Dancoff 係数と平均弦長で、等価原理により遮蔽する

## 背景

燃料のピンの UFG の外にも、分離共鳴を持つ核種の領域（支持格子の枠、案内管と計装管の被覆管など）がある。無限希釈の
ままでは、枠の Zircaloy の吸収を k で約 380 pcm 多く見積もった（1B、OpenMC と比較。docs/plans/m3-core-hzp.md、0e8af68、
2026-10-08）。

- 隣のセルとつながる、露出した面が一部だけの塊にも使えること
- 計算の量（黒体ごとに MOC を解く）

## 検討した案

- 無限希釈のまま: 上のとおり吸収を多く見積もる
- 核種ごとに黒体を分ける: 正確だが、組の数だけ MOC を解く
- 分離共鳴を持つ核種の組が同じ物質の領域をまとめて黒体にする（採用）

## 決定と理由

全体計画のとおり Bondarenko の表を等価原理の σ0 で引く（[ADR-0074](0074-shield-outside-range-by-equivalence.md)）。領域を黒体にした MOC で、黒体（セルの中で接する、
分離共鳴を持つ核種の組が同じ物質の領域）の Dancoff 係数と平均弦長（等方な中性子束 1 での入射率 R_iso から
l̄ = V / R_iso）を求め、黒体の中の物質ごとに燃料と同じ Carlvik の2項有理近似を使う。隣のセルとつながる枠のように、
露出した面が一部だけの塊にも使えるため。

同時に黒体にするのは、核種の組が同じ物質の領域（共鳴が同じエネルギーにあり影を作る）で、燃料は SPH 因子を掛けない
実効断面積、ほかは無限希釈の全断面積とする（燃料のピンの Dancoff 係数と同じ。[ADR-0087](0087-use-dilute-totals-in-shadows.md)）。組が一部だけ重なる物質
（Zircaloy とステンレスの Fe など）の影は入れない。黒体の中で密度の違う物質（囲まれた薄い層など）も、黒体の Dancoff
係数と平均弦長で遮蔽する近似とする。

## 確認

`structural_frames_are_shielded_as_lumps`、`touching_materials_with_the_same_resonant_nuclides_are_one_body`、
`materials_with_the_same_resonant_nuclides_shadow_a_lump`、`other_resonant_materials_do_not_shadow_a_lump`、
`connected_frames_see_only_their_open_faces`、`isolated_annulus_has_its_own_chord`。遮蔽した 1B の枠の吸収断面積は、OpenMC の
枠の中性子束で重みづけした平均で 0.00195 /cm（OpenMC 0.00194 /cm。examples/vera/verification.md、612913c、2026-10-08）。

## 帰結

- 核種ごとに黒体を分ける方式は、必要になったところで見直す
- 円筒でない塊への Carlvik の係数の適用は #0006、AIC は #0007
