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
なく、種類をまたぐ傾向（k∞ の差の 0 の点、最大とその燃焼度、60 MWd/kg の値、ピン出力の rms の範囲）を 1 か所で見る表が無い。
傾向を読むには 8 つの節を行き来する必要がある。

2026-10-24 の値（xlc 78a9a9b64 の examples/depletion/p2/compare-*.json から、親のセッションが集計したもの。表を作るときは
JSON から出し直す）:

- k∞ の差が 0 になる燃焼度の点（pcm）: 2A +11、2A-density +13、2A-rods +36、2F +63、2N +62、2O +66、2P +93、2M +10
- 60 MWd/kg の k∞ の差は +176〜+210 pcm。2A-rods だけ +66 pcm
- 差の最大は 44〜52 MWd/kg

ユーザーが 2026-10-24 に「今は直さない」と決めた。

関係する ADR: なし

従う決定: ADR-0176（検証の記録の表は比較のステージの JSON からツールが描く）、ADR-0178（集合体の燃焼の検証）

## 完了の条件

`examples/verification_tables.py` に、`examples/depletion/p2/compare-*.json` から上の項目を描く表（印は
`<!-- table: depletion-assembly-summary -->` のような名前）を足し、「集合体の燃焼」の節に置く。
`python3 examples/verification_tables.py --check` が通り、表の数値が各節の表と JSON に一致する。

## 経緯

- 2026-10-24: p2 の再実行のあとの評価で、種類をまたぐ表が無いことが分かり、今は直さないと決めた
