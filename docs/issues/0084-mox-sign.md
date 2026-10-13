---
status: open
labels: [verification, xs, shielding]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# MOX のピンセルで OpenMC との k の差の符号が UO2 と逆になる

## 問題

[ADR-0133](../adr/0133-scatter-h-in-h2o-up-to-the-evaluated-sab-limit.md) のライブラリで、MOX のピンセルの k の差（xlc − OpenMC）は MOX-1 で −225 pcm、MOX-2 で −113 pcm
（examples/mox/verification.md の「結果」。xlc は 1f4fcaaf、2026-10-12）。ADR-0133 の前は +79、+178 pcm で、UO2 のピンセル
（1B、+60 pcm）とは符号が逆になった。

分かっていること（同じ記録）:

- 燃料の核種の 14〜69 群の実効断面積の差を 1 次の近似で k に換算すると MOX-1 −162 pcm、MOX-2 −2 pcm。MOX-1 の大半は
  Pu-242 の 2.67 eV の共鳴（29 群、捕獲 +11%、k で −170 pcm）で、UFG の範囲（4 eV）の下にある（#0001）
- 残り（k の差 − 燃料の核種の計）は MOX-1 で −63 pcm、MOX-2 で −111 pcm で、燃料の断面積の差では説明できない

候補（どれも確かめていない）:

- 4 eV より下の共鳴の自己遮蔽（#0001。MOX-1 の Pu-242）
- 熱群の水の移行（ADR-0133 で変わったもの）と、Pu の熱群の大きな共鳴（0.3 eV）のあるスペクトルでの散乱行列の群の中の重み
- MOX-2 の残り −111 pcm は、水と被覆管の吸収と移行の差（MOX では balance を取っていない）

関係する ADR: なし

## 完了の条件

MOX-1 と MOX-2 で `openmc_balance.py` と `compare_balance.py` の内訳を取り、残りの −63、−111 pcm を物質と範囲に分け、
原因を特定して直すか、定量化して記録する。

## 経緯

- 2026-10-12: ADR-0133 の測り直し（#0070）で符号が反転した。2G の AIC の調査（#0082）で起こした
