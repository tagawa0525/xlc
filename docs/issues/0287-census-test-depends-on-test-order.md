---
status: closed
labels: [lattice, test]
opened: 2026-11-14
closed: 2026-11-14
milestone:
blocked-by: []
resolution: コミット c97a8fb4d
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

## 原因

tracing は呼び出し箇所の関心（interest）をプロセスで 1 つ持ち、購読者（dispatcher）が 1 つしかない間は、関心を登録する
スレッドの既定の購読者だけから求める（tracing-core 0.1.36 の `callsite.rs` の `Dispatchers::rebuilder` の `JustOne`）。
このテストの購読者が作られると全体の最大のレベルが DEBUG になり、その後に購読者の無いほかのテストのスレッドが
`Thermal::new` の `tracing::enabled!(DEBUG)` を先に通ると、関心は「無い」で登録され、このテストの購読者の下でも false を
返す。計数が切れたまま減速計算が走るので、`counted` が 0 になる。2 つのスレッドの最小の再現（tracing 0.1.44）で、別の
スレッドが先に通ると購読者の下で false、通らなければ true を確かめた。main の 1765e5381 で
`ufg::free_gas` のテストを 20 回実行して 14 回落ちた。

## 修正方法

`Thermal::new` が計数するかを引数で受け、本番の減速計算は今までどおり全体の購読者の下の `tracing::enabled!(DEBUG)` を渡す
（DEBUG でないときは数えない）。テストは tracing を介さずに計数を入れる。テストの側で関心を作り直す
（`rebuild_interest_cache`）案は、作り直した後に別のスレッドが登録すると同じく「無い」になるので採らなかった。

crates/ でスレッドの購読者（`with_default`、`set_default`）を置くテストはこの 1 件だけだった。xs の
`parallel::tests` は全体の購読者（`set_global_default`）で、関心は全体の購読者から求められるので同じ原因では落ちない。

## 経緯

- 2026-11-14: #0285 のブランチのチェックで落ち、main でも再現した
- 2026-11-14: 原因を確かめ、計数の切り替えを `Thermal::new` の引数にして閉じた。直した後に 30 回実行して落ちなかった
