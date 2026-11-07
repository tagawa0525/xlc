---
status: open
labels: [depletion, verification]
opened: 2026-10-24
closed:
milestone: M4
blocked-by: []
resolution:
---

# 集合体をまたぐ燃焼の比較の要約の表が検証の記録に無い

## 問題

examples/depletion/verification.md の「集合体の燃焼」には、集合体ごとの節（2A、2F、2N、2O、2P、2M、2A-density、2A-rods）しか
なく、種類をまたぐ傾向（k∞ の差の燃焼度 0 の点、最大とその燃焼度、60 MWd/kg の値、ピン出力の rms の範囲）を 1 か所で見る表が無い。
傾向を読むには 8 つの節を行き来する必要がある。

2026-11-03〜04 の値（examples/depletion/p2/compare-*.json の `k[].diff_pcm` と `pins[].rms`。xlc の計算は 06d6bd71（分岐点は main の a0fd48c0）、
2026-11-03〜04、JSON のコミットは 8488f5dfa。8 本とも ADR-0238 の前のライブラリの値で、2A は格子計算だけを
作り直し（examples/depletion/p2/dvc.lock の `compare@2A` の deps の lattice-2A.msgpack が `lattice@2A` の outs と違う）、ほかの 7 本の格子計算の燃焼は
ADR-0238 の「帰結」のとおり changed で、再実行していない。表を作るときは JSON から出し直す）:

- k∞ の差は 8 本の全体で +17〜+202 pcm（2A +56〜+164、2F +68〜+171）。2A-rods（+17〜+77）のほかの 7 本は ADR-0032 の目標の ±100 pcm を一部の点で
  超え、許容の ±150 pcm を 2〜25 点で超える
- ピン出力の rms の最大は 2P の 0.235% と 2M の 0.221% が目標の 0.2% を超え、ほかは 0.171% 以下

親のセッションが今は直さず issue に残す案を出し、ユーザーが 2026-10-24 09:17 にそれを認めた。

関係する ADR: なし

従う決定: ADR-0176（検証の記録の表は比較のステージの JSON からツールが描く）、ADR-0178（集合体の燃焼の検証）

## 完了の条件

`examples/verification_tables.py` に、`examples/depletion/p2/compare-*.json` から上の項目を描く表（印は
`<!-- table: depletion-assembly-summary -->` のような名前）を足し、「集合体の燃焼」の節に置く。
`python3 examples/verification_tables.py --check` が通り、表の数値が各節の表と JSON に一致する。

## 経緯

- 2026-10-24: p2 の再実行のあとの評価で、種類をまたぐ表が無いことが分かり、今は直さないと決めた
- 2026-11-08: M4 の計画の「進め方」の 5 から、2026-11-03〜04 の再実行の要約の値を移し、2026-10-24 の値と置き換えた
