---
status: open
labels: [verification, xs, shielding]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# 燃料の温度が高いほど OpenMC との k の差が小さくなる

## 問題

[ADR-0133](../adr/0133-scatter-h-in-h2o-up-to-the-evaluated-sab-limit.md) のライブラリで、OpenMC との k の差（xlc − OpenMC）は燃料の温度が高いほど小さい
（examples/vera/verification.md の「OpenMC との比較」。xlc は 1f4fcaaf、2026-10-12）:

| 問題 | 燃料の温度 | 差      |
| ---- | ---------- | ------- |
| 1B   | 600 K      | +60 pcm |
| 1C   | 900 K      | −1 pcm  |
| 1D   | 1200 K     | −69 pcm |
| 2B   | 600 K      | +45 pcm |
| 2C   | 900 K      | −8 pcm  |
| 2D   | 1200 K     | −61 pcm |

1B から 1D で −129 pcm、2B から 2D で −106 pcm。xlc は燃料の温度による k の変化（ドップラーの反応度）を OpenMC より
大きく見積もる。OpenMC の k の統計誤差は 1 本あたり約 5〜6 pcm で、2 本の差の統計誤差は約 8 pcm。

候補（どれも確かめていない）:

- UFG の散乱核で標的を静止として扱うこと（[ADR-0070](../adr/0070-treat-targets-at-rest.md)）。ADR-0070 は U-238 の 36.7 eV と 20.9 eV の群の吸収が
  600 K で +0.3〜+0.8%、1200 K で +1.6〜+2.6% 多いことを記録している（前の版のライブラリでの値）
- 共鳴の上方散乱（U-238 の共鳴の近くの標的の熱運動。OpenMC の側で何を使ったかも確かめていない）

関係する ADR: なし

## 完了の条件

1B と 1D の差の変化（−129 pcm）を、燃料の核種ごと・群ごとの吸収の差（`compare_fuel_xs.py`、`compare_balance.py`）に
分け、原因を特定して直すか、定量化して記録する。

## 経緯

- 2026-10-12: ADR-0133 の測り直し（#0070）で、温度による差の傾きが残ることが分かった。2G の AIC の調査（#0082）で起こした
