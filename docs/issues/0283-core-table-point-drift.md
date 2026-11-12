---
status: open
labels: [lattice, depletion, core]
opened: 2026-11-13
closed:
milestone:
blocked-by: []
resolution:
---

# ミクロ燃焼なしの core の k が、テーブルの点で格子計算の燃焼の履歴からずれる

## 問題

- #0279 の「段 5 の結果」の core（examples/depletion/issue-0279/core-micro/）で、ミクロ燃焼なしの core の k が、テーブルの燃焼度の点で
  格子計算の燃焼の履歴（examples/depletion/p2/lattice-2A.msgpack、md5 c028d773）の k∞ から −13.30〜−2.31 pcm ずれる。ADR-0221 の時点の
  同じ比べ方（examples/depletion/verification.md の「core のミクロ燃焼（2A の集合体 1 体）」）は −0.071〜−0.026 pcm だった
- ずれは core の前の、テーブルを作る格子計算で既にある: table-d0.msgpack の状態の k∞ と履歴の k∞ の差は燃焼度 0、20、60 MWd/kg で
  −2.855、−7.509、−13.299 pcm（core の k との差と同じ）。燃焼度 0 は組成が入力のままなので、組成の履歴によらず、格子計算の版か
  ライブラリの違いで生じている
- 履歴の版（examples/depletion/p2/dvc.lock の `lattice@2A` の deps）は bin/lattice md5 598a45ce とライブラリ md5 a96f4f23 で、
  テーブルは今の bin とライブラリ（md5 975be970）。その間に、自己遮蔽の自由ガスの上限の既定が変わり（ADR-0243、ADR-0246）、ライブラリに
  温度の点が足された（ADR-0239）

関係する ADR: なし

## 完了の条件

- どの違いがずれのどれだけを説明するかを数値で示し、正しい比べ方を記録に書く

## 経緯

- 2026-11-13: #0279 の段 5 の比較で、参照とのずれが持ち方の間の差（最大 2.7 pcm）より大きく、持ち方の良し悪しを決められなかったので起こした
