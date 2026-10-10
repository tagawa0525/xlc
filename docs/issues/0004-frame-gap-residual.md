---
status: open
labels: [lattice, shielding, verification]
opened: 2026-10-08
closed:
milestone:
blocked-by: []
resolution:
---

# 支持格子による k の変化の OpenMC との差が +40 pcm 残る

## 問題

支持格子の支持格子による k の変化は、等価なセルに支持格子のリングを入れて OpenMC との差が 1B で −170 → +40 pcm、2B で
−126 → +38 pcm になった（examples/vera/verification.md、前の値は f54c62d3、後の値は 3964be55 の記録。2026-10-08）。OpenMC の2つの計算の差の統計誤差は 1σ で約 10 pcm なので、
約 4σ にあたるが、原因は分けていない。

関係する ADR: [ADR-0088](../adr/0088-add-frame-ring-to-equivalent-cell.md)

## 完了の条件

残りの差の原因を特定し、統計誤差の 2σ 以内まで直すか、近似として ADR に記録する。

## 経緯

- 2026-10-08: 等価なセルへの支持格子のリング（C1）で −170 → +40 pcm
