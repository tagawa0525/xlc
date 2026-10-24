---
status: closed
labels: [lattice, branches, shielding, performance]
opened: 2026-10-24
closed: 2026-10-24
milestone: M4
blocked-by: []
resolution: ADR-0201
---

# 分岐の自己遮蔽をホウ素の軸だけで持ち越す案の誤差と時間を測る

## 問題

ADR-0201 の `shielding = "nominal-moderator"`（燃焼度の点と燃料の温度の点ごとに公称のホウ素と減速材で自己遮蔽を 1 度だけ解く）は、
状態ごとに解き直すテーブルとの k∞ の差が、ホウ素だけを変えた状態で 0.2〜0.4 pcm、減速材の密度を 3.6〜4.2% 変えた状態で 69〜118 pcm
だった（examples/depletion/verification.md の「分岐の自己遮蔽（ADR-0201）」。xlc d3e0b64e、2026-10-24、Slurm 531、532）。差の大半は減速材の
密度による（2O の温度を変えない密度の変化で −94〜−102 pcm）。

参照を燃焼度の点 × 減速材の点 × 燃料の温度の点ごとに持ち、ホウ素の軸の物質だけを作り直せば、誤差はホウ素の分（上の 0.2〜0.4 pcm）に
留まる見込みで、自己遮蔽の回数は ADR-0097 の M3 の 27 状態で燃焼度の点ごとに 9 回（今の every-state の 1/3）になる（見込みで、測っていない）。
ADR-0201 の作り直し（`history::remix_varied`）は、参照の取り方（`Branches::references` と `HistoryBranches::reference`）と作り直す物質（`Branches::varied`）を変えれば使える。

関係する ADR: ADR-0201

## 完了の条件

ホウ素だけを持ち越す方式を入力で選べる項目にするかを決め、作るなら ADR-0201 と同じ入力（examples/depletion/branches/shielding-*.toml）で
誤差と自己遮蔽の回数を測って記録する。

## 経緯

- 2026-10-24: ADR-0201 の測定で、ホウ素と減速材の密度の誤差が 2 桁以上違うと分かった
- 2026-10-24: ユーザーが「ホウ素の軸は再利用しましょう」と決め、ADR-0201 で既定（`shielding = "nominal-boron"`）にした。2A と 2O で、ホウ素だけを
  変えた状態の差は −0.38〜+0.30 pcm、自己遮蔽の回数は 27 状態で 9 回の比（examples/depletion/verification.md の「分岐の自己遮蔽（ADR-0201）」。
  xlc 667c3a51、Slurm 533）
