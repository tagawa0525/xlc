---
status: open
labels: [lattice, shielding, depletion]
opened: 2026-10-30
closed:
milestone: M4
blocked-by: []
resolution:
---

# 自己遮蔽の履歴（shielding-file）と分けるツールのレビューの残り

## 問題

[#0191](0191-micro-correction-first-order-error.md) の候補を分ける輸送のために足した `[history]` の `shielding-file`（crates/lattice/src/history.rs の `HistoryBranches::shield_with`）と
examples/depletion/micro/make_shifted.py の `--kind` の、ブランチ feat/0191-linearization-error-split の /code-review（1 周目）で、
直さなかった指摘（再現する誤動作がなく、効率、重複、使い勝手のもの）:

- `shielding = "every-state"`（またはホウ素が 1 点の分岐）と `shielding-file` を組むと、状態ごとに参照の番号（`references.len() + b`）が
  違い、`BranchShielder` の `held` が燃焼度の点の中のすべての状態の `FuelShielding` を持つ。#0191 の計算は燃焼度の点ごとに状態 1 つなので
  増えない。状態の多い分岐で使うなら、1 つの状態だけが使う参照を使った後に捨てる
- `HistoryBranches::reference_state` は 1 つの番号で 2 つの並び（ホウ素の軸の参照と、分岐の状態）を引き、範囲の外の番号は
  panic する。呼ぶのは `reference` の返した番号だけ。列挙型か、`shield_with` で参照の並びを伸ばす形にすれば番号の約束が要らない
- `BranchShielder` に渡す作り直す物質は、呼び出し側（main.rs と history/tests.rs）がホウ素の軸の物質に `remixed` を足して作る。
  `HistoryBranches` に全体を返すメソッドを足せば重複がなくなる
- `shield_with` は `history_branches` の燃焼の領域の照合と燃焼度の点の引き方を繰り返し、食い違う領域を示さない。共通の関数にする
- `a_shielding_history_must_match_the_branch_history` の燃焼度の点のない場合は、エラーの文が "2" を含むことだけを見る。核種の並びの
  検査は、それだけを外しても輸送する核種の検査で止まるので、単独で確かめていない
- make_shifted.py の `--kind previous` は追う核種を使わないが、`micro` のある入力を求める
- examples/verification_tables.py の `micro_linearization_single` は、`--single` なしの要約で TableError でなく KeyError で止まる
- compare_split.py の `report` は、1 核種の組成の内訳の交差の項（最大 0.84 pcm）を出さない
- compare_split.py の `NUCLIDE_THRESHOLD_PCM`（5 pcm）に理由を書いていない。Am-243 は最大 5.04 pcm で閾値の近くにあり、再計算で
  0.05 pcm 変われば表の列が消える
- compare_split.py は各分の値（行の値と `split` の全体）を同じ式で 2 度書く。1 度求めて両方に渡せば食い違わない
- test_compare_split.py の `test_the_single_nuclide_error_splits_into_its_term_and_the_rest` は `compare` を 2 度呼び、1 核種の内訳の
  χ と交差の項を確かめない
- `shielding-file` の形（ADR-0201 のホウ素の軸の参照の仕組みで燃焼の領域の物質を作り直す）を ADR にしていない。#0191 の検査のための
  入力の項目で、既定の挙動は変えない

関係する issue: #0191

関係する ADR: なし

## 完了の条件

- 上の各項を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-10-30: ブランチ feat/0191-linearization-error-split の /code-review の 1 周目と 3 周目で挙がった。「自己遮蔽なしで `shielding-file` を
  黙って無視する」は、`[history]` が自己遮蔽を求める（crates/lattice/src/input.rs）ので前提が誤り、「compare_split.py がコミット
  されていない」は、レビューの後にコミットした
