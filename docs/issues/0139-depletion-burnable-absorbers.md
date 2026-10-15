---
status: open
labels: [lattice, depletion]
opened: 2026-10-15
closed:
milestone: M4
blocked-by: []
resolution:
---

# 格子計算の燃焼は燃料のリングの物質しか燃やせず、Pyrex、IFBA、WABA の吸収体を燃やせない

## 問題

集合体の燃焼の検証（ADR-0178）の xlc の入力のうち、`[depletion] materials` に Pyrex（2F）、IFBA の被覆（2M、2N）、WABA の吸収体（2N）を
挙げたものは、lattice が燃料のリングの物質だけを燃焼の領域にできないので止まる（1 ステップの試行、2026-10-15、examples/depletion/p2/lattice-<ジョブ>.toml
から燃焼度の点を 0.1 だけにしたもの）:

- 2F: `cell 63: the depleted material is not shielded fuel`（crates/lattice/src/depletion.rs の `depletion_regions` の後。Pyrex は自己遮蔽する燃料でない）
- 2M、2N: `cell 20: region 60 is outside the fuel rings`（IFBA の被覆は燃料のリングの外の領域）

OpenMC の参照解はこれらの吸収体も燃やす（openmc_depletion_assembly.py）ので、xlc の領域に吸収体を足すか、決めて比べる範囲を絞るまで、2F、2M、2N の
xlc と OpenMC の比較ができない。2A、2A-density、2A-rods（燃料だけ）と、2O、2P（Gd 燃料のリング。`track-table = false`。#0129）は走る。

関係する ADR: ADR-0178

## 完了の条件

燃料でない物質（B-10 などの吸収体、Pyrex の Si と O）を燃焼の領域にでき（反応率は自己遮蔽しない微視的断面積、または吸収体の自己遮蔽を決める）、
2F、2M、2N が走る。IFBA の 10 μm の被覆の中の flux の取り方（薄い領域の MOC の平均）を、OpenMC の燃焼の B-10 の数密度で確かめる。

## 経緯

- 2026-10-15: 集合体の燃焼の入力（ADR-0178）の xlc の試行で分かった
- 2026-10-15: 吸収体を燃焼の領域にできるようにした（ブランチ fix/0139-deplete-absorbers）。燃焼の領域の分け方（対称操作の組 × 物質）は
  ADR-0173 のままで、OpenMC の openmc_depletion_assembly.py の `regions` と同じ。領域の微視的断面積は、その物質を遮蔽した値
  （`FuelShielding::region_cross_sections`）: 燃料のピンの等価なセルのリングは、UFG の範囲の群がそのリングの減速計算の値で、範囲の外は
  燃料のリングが等価原理、ほかのリング（IFBA の被覆など）が無限希釈。遮蔽した塊（Pyrex など）は塊の等価原理の値（中心のリングを
  埋める棒は UFG の値）。遮蔽しない物質は MOC と同じ無限希釈。燃やす物質の温度が違うと止める（輸送する核種の集合を決めていない）。
  古い p2 のライブラリ（ADR-0177 の前）で 2F と 2M を 1 ステップ（0.1 MWd/kg）流し、走ることと、Pyrex（5 領域）の B-10 が −1.3%、
  IFBA の被覆（19 領域）の B-10 が −2.0〜−2.2% 減ることを確かめた（記録に書かない動作確認。Slurm のジョブ 405）
- 残り: IFBA の被覆の中の中性子束の取り方を、作り直したライブラリでの OpenMC との比較の B-10 の数密度で確かめる。燃やす吸収体は、
  OpenMC と同じく輸送できる核種すべてを数密度 0 から持つので、IFBA の被覆のリングの UFG の縮約と Pyrex の塊の反復が核種の数だけ
  重い（#0137）。中心のリングを埋める吸収体の棒を燃やすと、UFG の範囲（`material_range`）が数密度 0 の核種の分離共鳴で広がる
  （p2 にはない。確かめていない）
