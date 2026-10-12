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

clippy の `too_many_lines`（100 行）を超える関数が、main（1842033d、2026-10-12）で 36 個ある（本番のコード 26 個、
テスト 6 個、examples 4 個）。長い関数は、いくつもの段（`solve_with` では初期化、スイープ、収束の判定、出力。`run` や `parse` では入力の
部分ごとの処理）を 1 つの関数に持ち、段ごとにテストやレビューができない。

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

- `cargo clippy --workspace --all-targets -- -W clippy::too_many_lines` が本番のコードに 100 行を超える関数を出さないか、
  出る関数ごとに分けない理由がここに書いてある（テストと examples の関数は、分けるかをここで判断する）
- 長い関数が再び育つのを止めるか（clippy の `too_many_lines` を拒否にし、残る関数に `#[expect]` を付ける案など）を、
  ADR で決めてある
- 本番のコードを分けたブランチ（ADR-0144 の「確認」に書く最初のブランチを除く）ごとに、ADR-0144 の判定（全桁の msgpack のビット一致と、時間の比較）の結果がここにある

## 経緯

- 2026-10-12: コードのレビューで、長い関数が保守性の最大の弱点と指摘された。ADR-0144 で進め方を決めた
