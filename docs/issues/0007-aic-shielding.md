---
status: open
labels: [lattice, shielding]
opened: 2026-10-08
closed:
milestone: M3
blocked-by: []
resolution:
---

# AIC の自己遮蔽（C4）

## 問題

制御棒の AIC（Ag-In-Cd）は、今は構造材の塊としてしか扱えない。Ag と In の共鳴は強く、Ag-109 の 5.2 eV は U-238 の
6.67 eV に近い。影の計算の全断面積を無限希釈にする近似（docs/plans/m3-core-hzp.md）の前提も、ここで崩れる
おそれがある。In-115 の 1.46 eV は UFG の範囲の下にある（#0001）。

## 完了の条件

制御棒のセルを入れる段階で、OpenMC と比べて AIC の吸収と制御棒の反応度を確かめる。影の計算の全断面積は、
燃料の吸収率で重みづけした平均（UFG の範囲の群。examples/vera/README.md、
docs/plans/m3-core-hzp.md。実装はブランチ feat/dancoff-potential に残してある）で測り直す。

## 経緯

- 2026-10-08: 自己遮蔽の前提の洗い出しで挙がった。制御棒の段階で扱う（ユーザーと合意）
