---
status: open
labels: [depletion, examples, maintenance]
opened: 2026-10-22
closed:
milestone:
blocked-by: []
resolution:
---

# ピン出力の比較のツールとテストに残る重複

## 問題

燃焼点ごとのピン出力の比較（ADR-0181。ブランチ feat/0135-depletion-pin-powers）のレビューで、構成と重複だけの指摘が 4 件残った。
精度と速度に関わらないので、マージの前に直さず、ここにまとめる（/git-merge の打ち切りの規則）:

1. examples/verification_tables.py の `ASSEMBLY_DEPLETION_JOBS` が examples/depletion/p2/dvc.yaml の `jobs` のコピーで、一致は
   examples/test_verification_tables.py が正規表現で dvc.yaml を読んで確かめている。問題を足すと 2 か所を直す
2. examples/depletion/test_compare_depletion_assembly.py の `pins(edge, corner)` と `std(edge, corner)` の本文が同じ
3. crates/lattice/src/depletion/tests.rs の `written_pin_powers` が、`history_is_written_as_a_map_with_the_format_name` の
   map の引き方（`get`）を書き直している
4. crates/lattice/src/depletion.rs の `State.pin_powers` が `Option` で、`run` の終わりで `anyhow` のエラーに開く。残す輸送
   （`kept`）では必ず `Some` なので、このエラーは出ない。`few_group` も同じ条件なので、残す値を 1 つの `Option<Kept>` に
   まとめれば実行時の検査がなくなる

関係する ADR: ADR-0181

## 方針の案

- 1: 比較のステージの JSON に問題の名前が入っている（`problem`）ので、表の名前を JSON から取れば一覧が要らなくなる。
  表の定義が記録のファイルの `<!-- table: -->` の名前で引かれるので、名前の規則と JSON の対応を決めてから
- 2、3: 共通の helper にする。小さいので、同じファイルを次に触るときに合わせて直す
- 4: `ce_cm` が `kept` を型で返す形にする

## 完了の条件

4 件を直すか、直さない理由を書いて閉じる。

## 経緯

- 2026-10-22: ADR-0181 のレビュー（/code-review high）の指摘から
