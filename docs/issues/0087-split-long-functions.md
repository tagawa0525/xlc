---
status: open
labels: [refactor, lattice, core, xs]
opened: 2026-10-12
closed:
milestone:
blocked-by: ["#0089"]
resolution:
---

# 長い関数を段ごとに分ける

## 問題

clippy の `too_many_lines`（100 行）を超える関数が、main（14793ea7、2026-10-12）で 36 個ある（本番のコード 26 個、
テスト 6 個、examples 4 個）。長いものは反復の段（初期化、スイープ、収束の判定、出力）を 1 つの関数に持ち、段ごとに
テストやレビューができない。

| 関数                                               | 行数 |
| -------------------------------------------------- | ---- |
| crates/lattice/src/moc.rs `solve_with`             | 349  |
| crates/lattice/src/main.rs `run`                   | 316  |
| crates/lattice/examples/swap_transfer.rs `main`    | 273  |
| crates/lattice/src/shielding.rs `shield_fuel_with` | 265  |
| crates/lattice/src/input.rs `parse`                | 244  |

（行数は clippy が数える、空行とコメントを除いたもの。残りは `cargo clippy --workspace --all-targets -- -W clippy::too_many_lines` で出る）

関係する ADR: [ADR-0144](../adr/0144-refactor-and-tune-with-bit-identical-outputs.md)

## 完了の条件

- `too_many_lines` を拒否にする lint があり、36 個の関数の `#[expect(clippy::too_many_lines, reason = "#0087")]` がすべて
  外れている
- 分けたブランチごとに、ADR-0144 の判定（DVC のステージの出力のビット一致と、時間の比較）の結果がここにある

## 経緯

- 2026-10-12: コードのレビューで、長い関数が保守性の最大の弱点と指摘された。ADR-0144 で進め方を決めた
