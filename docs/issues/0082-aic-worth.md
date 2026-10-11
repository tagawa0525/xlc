---
status: open
labels: [verification, lattice, shielding]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# VERA 2G の AIC の制御棒の価値を OpenMC より大きく見積もる

## 問題

[ADR-0133](../adr/0133-scatter-h-in-h2o-up-to-the-evaluated-sab-limit.md) のライブラリで、2G（AIC の制御棒 24 本）の OpenMC との k の差（xlc − OpenMC）は −212 pcm で、
吸収体のない 2E-empty（+45 pcm）から −257 pcm 動く（examples/vera/verification.md の「吸収体による差の変化」。xlc は e4dbc0de、
2026-10-12）。xlc は AIC の制御棒の価値を OpenMC より 257 pcm 大きく見積もる。変化の統計誤差は大きめに 12 pcm、方位角 30・
間隔 0.03 cm の離散化の誤差は方位角 128・間隔 0.01 cm に対して −31 pcm で、どちらでも説明できない。B4C（2H）は逆に
+63 pcm 動く。

候補（記録と ADR にあるもの。どれも確かめていない）: Ag、In、Cd の共鳴の自己遮蔽（構造材の塊としての等価原理、
[ADR-0086](../adr/0086-shield-structural-lumps.md)、UFG の範囲の下の In-115 の 1.46 eV。#0007、#0001）、影の計算の全断面積（[ADR-0087](../adr/0087-use-dilute-totals-in-shadows.md)）、
トラック間隔（[ADR-0111](../adr/0111-use-30-azimuths-for-vera-assemblies.md)）、等価セルの枠（[ADR-0088](../adr/0088-add-frame-ring-to-equivalent-cell.md)）、2 項有理近似の重み（[ADR-0089](../adr/0089-weight-two-terms-with-potential-scattering.md)）、MF4 の角度分布（[ADR-0133](../adr/0133-scatter-h-in-h2o-up-to-the-evaluated-sab-limit.md)）。

関係する ADR: [ADR-0086](../adr/0086-shield-structural-lumps.md)

## 完了の条件

2G の k の差を物質ごと・群ごとの正味の吸収の差に分け（`openmc_balance.py`、`compare_balance.py`）、AIC の核種と群の
寄与を特定し、小さい体系で候補を 1 つずつ切り分けて原因を定量化する。直しが同じトピックに収まれば ADR と TDD で直し、
収まらなければ測った値と修正方法の案をこの issue（か別の issue）に書いて閉じる。

## 経緯

- 2026-10-12: ADR-0133 の測り直し（#0070）で、2G の差が +27 pcm から −212 pcm になり、AIC による変化が目立った。#0007
  （AIC の自己遮蔽）は原因の候補の 1 つで、この issue は差の原因を特定する
