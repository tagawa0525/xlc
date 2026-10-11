---
status: accepted
date: 2026-10-08
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0090: 炉心は対称性を使わず全体を解く

## 背景

VERA Problem 5 の制御棒バンク SC と SD は 1/8 対称でない。鏡面反射の 1/4 炉心では別の配置を解くことになる（仕様書。
KENO は 1/4 炉心で解いて −53/−39 pcm を補正している。docs/plans/m3-core-hzp.md、efdfc8f、2026-10-08）。

## 検討した案

- 鏡面反射の 1/4 炉心: 別の配置を解くことになり、補正が要る
- 全炉心（採用）

## 決定と理由

炉心は全体を解く。NEM なら全炉心でも計算は軽い。

## 確認

まだ確かめていない（core の装荷図と全炉心の入力は M3 の作業）。M3 の Problem 5 を、仕様書の実測値と KENO の値と比べて確かめる。

## 結果

なし。
