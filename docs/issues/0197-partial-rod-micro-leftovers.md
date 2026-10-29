---
status: open
labels: [core, depletion]
opened: 2026-10-29
closed:
milestone: M4
blocked-by: []
resolution:
---

# 部分挿入のノードのミクロ燃焼の残り（数密度が割合について単調でないことと、#0192 のレビューの残り）

## 問題

部分挿入のノードのミクロ燃焼（ADR-0217）は、#0192 で 2A の集合体 1 体を確かめて閉じた。次が残っている。

- 2A の集合体 1 体の上のノードで、追う核種 31 のうち 17 の数密度が、挿入の割合について単調でない点を持つ
  （examples/depletion/verification.md の「部分挿入のノードのミクロ燃焼（2A の集合体 1 体）」。k とノードの燃焼度は単調）。
  候補は、ノードの燃焼度（割合について下がる）と、棒の入った区間の重みで硬くなるスペクトルが逆の向きに効くこと。燃焼度あたりの変化は
  U-238 で単調だが、数密度が単調な核種を含む 26 核種で単調でなく、候補を区別しない（確かめていない）。割合ごとの上のノードの履歴を同じノードの燃焼度に内挿して比べれば、燃焼度の
  差を除いて確かめられる（`compare_core_partial_rods.py` は燃焼度で割った値だけを見る）。区間の重みの不具合なら、Problem 9 の
  サイクルの燃焼に入る
- #0192 のレビューの 3 周目の残り（#0192 の「経緯」）:
  - 部分挿入のノードの区間が別のテーブルを使うときのエラーと、`mix_micro` の核分裂率の和が正でないときのエラーに試験がない
  - 区間の収率の境界より下の割合は、`interpolate_micro` が群の和の σf で内挿してから、`mix_micro` が中性子束で重みを付けた率で
    混ぜるので、1 つの量を 2 つの規則で混ぜる（大きさは測っていない）
  - examples/core/rod_table.py の `DENSITY_ROUNDING`（1e-12）のコメントは「数個の ulp」と書くが約 4500 ulp。2A の制御棒の組と
    棒のない組の差は相対 2.0e-14（#0192 の「経緯」）
  - `rods::weighted` は `average` を呼ぶだけで、`mix_micro` は群ごとの重みの和を核種ごとに作り直す
- examples/depletion の比較のツールは、燃焼度の一致の許容 `BURNUP_TOLERANCE`（1e-9）と `same_burnup`、`rho_pcm` を
  compare_core_depletion.py、compare_depletion.py、compare_depletion_assembly.py、compare_core_micro.py、compare_core_partial_rods.py の
  5 つに複製している

関係する ADR: ADR-0217

## 完了の条件

- 数密度が割合について単調でない原因を、同じノードの燃焼度で比べて確かめるか、区間の重みの不具合として直す
- ほかの各項目を直すか、直さない理由を書いて閉じる
