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
6.67 eV に近い。影の計算の全断面積を無限希釈にする近似（[ADR-0087](../adr/0087-use-dilute-totals-in-shadows.md)）の前提も、ここで崩れる
おそれがある。In-115 の 1.46 eV は UFG の範囲の下にある（#0001）。

関係する ADR: [ADR-0086](../adr/0086-shield-structural-lumps.md)、[ADR-0087](../adr/0087-use-dilute-totals-in-shadows.md)

## 完了の条件

制御棒のセルを入れる段階で、OpenMC と比べて AIC の吸収と制御棒の反応度を確かめる。影の計算の全断面積は、
燃料の吸収率で重みづけした平均（UFG の範囲の群。examples/vera/verification.md、
[ADR-0087](../adr/0087-use-dilute-totals-in-shadows.md)。実装はブランチ feat/dancoff-potential に残してある）で測り直す。

## 経緯

- 2026-10-08: 自己遮蔽の前提の洗い出しで挙がった。制御棒の段階で扱う（ユーザーと合意）
- 2026-10-11: VERA 2G（AIC の制御棒 24 本）で、xlc と OpenMC の k の差は +27 pcm、吸収体のない 2E-empty（+247 pcm）から −220 pcm
  動く（examples/vera/verification.md の「OpenMC との比較」、c0662d26）。B4C は逆に +100 pcm と正の向きに動かす（Pyrex 24 本の +41 pcm は、統計誤差とトラック間隔の振れと区別できない）。
  xlc が AIC の吸収を多く見積もるかは、吸収体の中の反応率を OpenMC と比べていないので確かめていない
- 2026-10-12: ブランチ feat/dancoff-potential には、自動レビューの指摘 10 件が未対応のまま残っている。指摘の中身はどこにも
  残していないので、測り直しに使うときにレビューを回し直して対応する
