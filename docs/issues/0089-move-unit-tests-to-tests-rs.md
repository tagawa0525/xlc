---
status: open
labels: [refactor]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# 単体テストを tests.rs に移し、本番のファイルのテストのコードを検査で止める

## 問題

main（14793ea7、2026-10-12）で、`crates/*/src` の 64 個のファイルがテストをファイルの中に持つ。テストだけの項目が
本番のファイルに 4 個ある（crates/format/src/fewgroup.rs の `sample`、`sample_pins`、crates/numeric/src/pins/node.rs の
`modes_at`、`flux`）。置き方と確かめ方は ADR-0143 で決めた。

関係する ADR: [ADR-0143](../adr/0143-put-unit-tests-in-tests-rs.md)

## 完了の条件

- crates/docs の検査が、本番のファイルの中のテストのコードを見つけて失敗し、移したあとの main で通る
- 64 個のファイルのテストと 4 個の項目が `tests.rs` にあり、ADR-0143 の「確認」に実際の結果（実行ファイルのビット一致、
  node.rs の差分が `#[cfg(test)]` の項目だけであること、`passed` の件数）を書いて `status: accepted` にしてある
  （docs/adr/README.md の「運用」）

## 経緯

- 2026-10-12: コードのレビューで、大きなファイルの半分以上がテストと分かった。ADR-0143 で置き方を決めた
