---
status: open
labels: [core, depletion, branches]
opened: 2026-10-23
closed:
milestone: M4
blocked-by: []
resolution:
---

# 燃焼度の点の指針を Pyrex の 2F で確かめる

## 問題

ADR-0192 の燃焼度の点の指針（26 点）は、Gd の 2O で決め、2A で確かめた（examples/depletion/verification.md の「燃焼度の内挿の方式と
点の選び方」、daaff322e、2026-10-23）。Problem 9 のサイクル 1 の装荷に入る Pyrex の集合体（2F。B-10 の燃え尽きの折れ曲がりがある）では
確かめていない。全点のテーブルの入力（examples/depletion/branches/branch-2F-all.toml、公称の 1 状態 × 履歴の 41 点）はあり、Slurm の
排他のジョブで 25〜30 分の見込み（2F の 1 輸送は自己遮蔽 28 秒、MOC 8 秒。ADR-0195 の前の値）。優先度の低い計算として、親の指示まで
投げないことにした（2026-10-23）。

関係する ADR: ADR-0192

## 完了の条件

2F の全点のテーブルを作り、examples/depletion/branches/burnup-interpolation.toml に 2F（「ADR-0190 の点」と「指針」）を足して比べ、
指針の点で点の間の誤差が ±10 pcm に入るかを記録に書く。入らなければ、入るのに要る点だけを指針に足す（ADR を改める）。

## 経緯

- 2026-10-23: ADR-0192 のブランチで入力を用意した。ジョブ（491、507）は、高速化の測定を先に実行するため取り消した
