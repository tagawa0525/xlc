---
status: open
labels: [lattice, test]
opened: 2026-11-14
closed:
milestone:
blocked-by: []
resolution:
---

# 自由ガスの核の計数のテストが、テストの実行の順で落ちる

## 問題

- `ufg::free_gas::tests::the_census_counts_the_table_bytes_of_its_points_by_nuclide`（crates/lattice/src/ufg/free_gas/tests.rs）が、
  ほかのテストと並べて実行すると `assertion failed: counted > 0` で落ちることがある。単独では通る。main の 1765e5381 で
  `cargo test -p xlc-lattice --lib ufg::free_gas` を 2 回実行して 1 回目は通り、2 回目は落ちた（2026-11-14）。#0285 のブランチの
  run-checks.sh も同じテストで落ちた（ブランチは crates/lattice/src/ufg/ を変えていない）
- テストは `tracing::subscriber::with_default` でスレッドの購読者を置き、DEBUG の計数の行から格子点の大きさを数える。候補
  （確かめていない）: 別のスレッドのテストが購読者の無いまま同じ呼び出し箇所を先に通ると、tracing が呼び出し箇所の関心を
  「無い」と持ち、このテストの購読者に行が届かない。計数が DEBUG の行を出すかどうかに頼るので、行が届かないと 0 になる
- run-checks.sh の全体の `cargo test` が順によって落ち、マージの前のチェックが通ったり落ちたりする

関係する ADR: なし

## 完了の条件

原因を確かめ、順によらず通るようにする（呼び出し箇所の関心を作り直す、計数を tracing の行でなく戻り値から取る、など。どちらに
するかは原因を確かめてから決める）。`cargo test -p xlc-lattice --lib ufg::free_gas` を繰り返して落ちないこと。

## 経緯

- 2026-11-14: #0285 のブランチのチェックで落ち、main でも再現した
