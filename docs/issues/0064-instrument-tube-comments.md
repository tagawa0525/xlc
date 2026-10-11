---
status: open
labels: [verification]
opened: 2026-10-11
closed:
milestone:
blocked-by: []
resolution:
---

# VERA 2F、2G、2H の入力で計装管のセルの注釈が「案内管」になっている

## 問題

examples/vera/p2/assembly-2F.toml、assembly-2G.toml、assembly-2H.toml の `[cell.it]` の注釈は「案内管（中は減速材）」だが、
このセルは中央の計装管の位置にだけ置かれ、計装管の半径（0.559 cm と 0.605 cm。案内管は 0.561 cm と 0.602 cm）を持つ。
2B〜2E などの入力では「計装管（中は減速材）」と書いている。計算には影響しない（2K のレビューで見つかった）。

関係する ADR: なし

## 完了の条件

3 つの入力の注釈を「計装管（中は減速材）」に直す。

## 経緯

- 2026-10-11: 2K の入力（2F から作った）のレビューで見つけ、2K の分だけ直した（ブランチ feat/vera-2k）
