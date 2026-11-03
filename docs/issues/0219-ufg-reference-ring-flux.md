---
status: open
labels: [lattice, shielding]
opened: 2026-11-03
closed:
milestone: M4
blocked-by: []
resolution:
---

# リングの中性子束の移行で、UFG の散乱核の既定が参照から 28〜46 pcm 離れる

## 問題

ADR-0209 は、UFG の散乱核の既定（速い方法 `{ upper = 400, cells = 2, sweeps = 3, heavy = 400, sab = true }`）を、参照
（`{ upper = inf, cells = 4, sweeps = inf, heavy = 400, sab = true }`）からの k の差 −4.3、−3.0 pcm（1C の新しい燃料、2A-40。ADR-0205）を
理由に採った。この差は無限希釈の形の移行（`ufg-transfer = "infinite-dilution"`）で測ったもので、既定のリングの中性子束の移行
（`"ring-flux"`、ADR-0207）では測っていなかった。

リングの中性子束の移行のまま散乱核だけを参照にすると、0 MWd/kg の k は 1C の 1 リングで 1.172387 → 1.172111（−27.6 pcm）、MOX-1 の
1 リングで 1.196631 → 1.196167（−46.4 pcm）動く（examples/depletion/issue-0152/ の `1c-r1.toml` と `1c-r1-ref.toml`、`mox1-r1.toml` と
`mox1-r1-ref.toml`。xlc は e5969df2a（ブランチ feat/0152-mox1-gap、分岐点は main の 54ab434ff）、2026-11-03、Slurm 661、662。#0152 の
2026-11-03 の項）。どちらも OpenMC との差を負の側へ動かす（1C の 1 リングの 0 の点 +52.7 ± 22.0 pcm、MOX-1 −253.5 ± 16.4 pcm。
examples/depletion/compare-1C-rings1.json、compare-mox1-rings1.json）。

候補（確かめていない）: 縮約に使うリングの UFG の中性子束が、スイープの回数（3 回と収束まで）や区間の細かさ（2 と 4）で変わり、その差が
移行の行に入る。無限希釈の形の移行では行が中性子束によらないので、この分は ADR-0205 の差に入らない。

関係する ADR: なし

## 完了の条件

参照と既定の差を、散乱核の項目（`upper`、`cells`、`sweeps`）ごとに 1 つずつ戻してリングの中性子束の移行で測り、どの項目の差かを
示す。差が ADR-0209 の理由（参照に近い）を崩すなら、既定の条件を ADR で決め直す。

## 経緯

- 2026-11-03: #0152 の MOX-1 の 0 の点の切り分け（ブランチ feat/0152-mox1-gap）で、自己遮蔽の設定の候補として参照の条件を解いて見つけた
