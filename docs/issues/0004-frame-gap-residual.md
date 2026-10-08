---
status: open
labels: [lattice, shielding, verification]
opened: 2026-10-08
closed:
---

# 支持格子による k の変化の OpenMC との差が +40 pcm 残る

## 問題

支持格子の支持格子による k の変化は、等価なセルに支持格子のリングを入れて OpenMC との差が 1B で −170 → +40 pcm、2B で
−126 → +38 pcm になった（examples/vera/README.md）。OpenMC の統計誤差（約 ±10 pcm）より大きいが、原因は分けて
いない。

## 完了の条件

残りの差の原因を特定し、直すか、近似として計画に記録する。

## 経緯

- 2026-10-08: 等価なセルへの支持格子のリング（C1）で −170 → +40 pcm
