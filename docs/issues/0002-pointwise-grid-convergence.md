---
status: open
labels: [lattice, shielding]
opened: 2026-10-08
closed:
milestone:
blocked-by: []
resolution:
---

# lattice のポイントワイズの減速計算の格子の収束

## 問題

lattice の減速計算（crates/lattice/src/ufg.rs）は、xs の flux calculator と同じ格子の作り方（散乱の窓に最低
WINDOW_POINTS = 16 点）を使っている。xs では格子を細かくしても中性子束が変わらないことを確かめているが、lattice の
セル（多リングの CP、全核種）では確かめていない（コードのコメントのとおり）。

## 完了の条件

VERA 1B などで WINDOW_POINTS を増やし、燃料の実効断面積と k が収束していることを確かめる。収束していなければ
値を見直す。#0001 で範囲を下に延ばすときも同じ確認をする。

## 経緯

- 2026-10-08: ポイントワイズの意味を確かめる中で、未確認であることを確かめた
