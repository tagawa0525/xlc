---
status: proposed
date: 2026-10-23
requires: [ADR-0173, ADR-0194, ADR-0195]
amends: [ADR-0173]
supersedes: []
superseded-by:
issues: ["#0134"]
---

# ADR-0197: 燃焼の自己遮蔽を、入力で選べばステップの始まりだけで解く

## 背景

格子計算の燃焼（[ADR-0173](0173-deplete-lattice-pins-by-ce-cm.md)。CE/CM で 1 ステップに輸送 2 回）は、輸送のたびに自己遮蔽を解き直す
（案 F1）。ADR-0173 は、中点で始まりの微視的断面積を使う案 F2 の誤差と時間を集合体で測って決めるとした（#0134）。燃料のピンの
Dancoff 係数の持ち越し（[ADR-0194](0194-reuse-pin-dancoff-while-inputs-match.md)）と MOC の前の解からの開始
（[ADR-0195](0195-start-moc-from-the-previous-transport.md)）の後の 2A では、1 輸送の自己遮蔽が 15.6 秒、MOC が 3.8 秒で、燃焼の時間の
約 8 割が自己遮蔽である（examples/depletion/verification.md の「集合体の燃焼」の「時間」）。F2 なら自己遮蔽の回数は、40 点（輸送 81 回）で
81 回から 41 回になる。

判断の軸:

- 近似は入力で選べる項目にし、誤差を測って記録する。既定は原理どおり（F1）のままにし、既定の切り替えは誤差を見てから別に決める
- 中点の輸送の物質は、始まりの自己遮蔽と同じ作り方で、数密度だけを中点のものにする（作り方の違いを誤差に混ぜない）
- shielding.rs の自己遮蔽の手順は変えない（並行する高速化の作業が触る）

## 検討した案

頻度の選び方:

- 案 S1: 入力の `[depletion]` に `shielding = "every-transport" | "step-start"` を足す。省略すると `every-transport`（採用）
- 案 S2: F2 を既定にする: 誤差を測る前に既定を変えることになる
- 案 S3: 自己遮蔽を N ステップに 1 回にする: 中点の F2 の誤差（組成の変化が半ステップ分）を測る前に、より大きな近似を足すことになる

中点の物質の作り方（F2）:

- 案 M1: 始まりの自己遮蔽の結果（`FuelShielding`）から、燃焼の領域の MOC の領域ごとに、核種の実効断面積
  （`region_cross_sections`。UFG の範囲はリングの減速計算の値、範囲の外は等価原理の値）から群定数（`shielded_nuclide`）を作り、中点の
  数密度で混ぜ（`Material::from_nuclides`）、その領域の SPH 因子を掛ける。自己遮蔽が燃料の物質を作る手順と同じ式で、燃焼させない
  物質（被覆管、減速材、構造材の塊）は始まりのまま（採用）
- 案 M2: 始まりの巨視的断面積をそのまま使う: 数密度の変化（Xe、Sm、Gd の吸収）も消え、F2 より粗い近似になる
- 案 M3: 自己遮蔽の関数に「微視的断面積を固定して混ぜ直す」口を足す: shielding.rs を変える。M1 は公開された結果（実効断面積、
  SPH 因子を戻す `physical_flux`）だけで同じことができる

## 決定と理由

S1、M1 にする。

- `[depletion]` の `shielding`: `"every-transport"`（ADR-0173 の案 F1。既定）か `"step-start"`（案 F2）。ほかの値はエラー
- `step-start` では、ステップの始まりと最後の点の輸送で自己遮蔽を解き、その結果を次の中点の輸送で使う（`depletion::remix_depleted`）。
  SPH 因子は、1 の中性子束を `FuelShielding::physical_flux` に通して領域ごとに引く。1 つの物質を別の燃焼の領域か燃焼させない領域と
  共有していれば、作り直せないのでエラーにする
- 作り直しの前提（自己遮蔽の燃料の物質が、領域の実効断面積と数密度から作った群定数に SPH 因子を掛けたもの）は、最初の輸送で、
  始まりの組成で作り直した物質が自己遮蔽の物質とビットで同じことを確かめ、違えば止める
- 履歴に自己遮蔽を解いた回数（`shieldings`）を書く
- これは ADR-0173 の「決定と理由」の「自己遮蔽は輸送のたびに解き直す（案 F1）」を、「自己遮蔽は既定で輸送のたびに解き直し（案 F1）、
  入力の `shielding = "step-start"` で中点の輸送の自己遮蔽を省く（案 F2）」に置き換える。ほかの部分は変わらない

## 確認

- 入力: `input::tests::depletion_shields_at_every_transport_unless_told_to_shield_at_the_step_start`
- 燃焼: `depletion::tests::every_transport_solves_the_shielding_by_default`（自己遮蔽 7 回）、
  `step_start_shielding_skips_the_midpoint_transports`（輸送 7 回、自己遮蔽 4 回。群定数が σ0 によらないテストのライブラリでは、F1 と
  相対 1e-12 で同じ）、`remixing_the_shielded_composition_reproduces_the_shielded_materials`、
  `remixed_materials_take_the_new_densities_and_keep_the_microscopic_cross_sections`、`remixing_needs_the_shielded_nuclides_of_each_region`
- 集合体（examples/depletion/verification.md の「集合体の燃焼」）: 測定の後に書く

## 帰結

- 測定の後に書く
