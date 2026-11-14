---
status: closed
labels: [lattice, test]
opened: 2026-11-14
closed: 2026-11-14
milestone:
blocked-by: []
resolution: コミット ff5ddb7b3
---

# 自由ガスの核の計数のテストが、テストの実行の順で落ちる

## 問題

- `ufg::free_gas::tests::the_census_counts_the_table_bytes_of_its_points_by_nuclide`（crates/lattice/src/ufg/free_gas/tests.rs）が、
  ほかのテストと並べて実行すると `assertion failed: counted > 0` で落ちることがある。単独では通る。main の 1765e5381 で
  `ufg::free_gas` のテストを 20 回実行して 14 回落ちた（2026-11-14）。#0285 のブランチの run-checks.sh も同じテストで落ちた
  （ブランチは crates/lattice/src/ufg/ を変えていない）
- テストは `tracing::subscriber::with_default` でスレッドの購読者を置き、`Thermal::new` の `tracing::enabled!(DEBUG)` で
  計数を入れていた。原因: tracing は呼び出し箇所の関心（interest）をプロセスで 1 つ持ち、購読者（dispatcher）が 1 つしかない
  間は、関心を登録するスレッドの既定の購読者だけから求める（tracing-core 0.1.36 の `callsite.rs` の `Dispatchers::rebuilder`
  の `JustOne`）。このテストの購読者が作られると全体の最大のレベルが DEBUG になり、その後に購読者の無いほかのテストのスレッドが
  同じ呼び出し箇所を先に通ると、関心は「無い」で登録され、このテストの購読者の下でも false を返す。計数が切れたまま減速計算が
  走るので、`counted` が 0 になる。2 つのスレッドの最小の再現（tracing 0.1.44）で、別のスレッドが先に通ると購読者の下で
  false、通らなければ true を確かめた。crates/ でスレッドの購読者（`with_default`、`set_default`）を置くテストはこの 1 件だけで、
  xs の `parallel::tests` は全体の購読者（`set_global_default`）なので、関心は全体の購読者から求められ、同じ原因では落ちない
- run-checks.sh の全体の `cargo test` が順によって落ち、マージの前のチェックが通ったり落ちたりする

関係する ADR: なし

## 完了の条件

原因を確かめ、順によらず通るようにする。`cargo test -p xlc-lattice --lib ufg::free_gas` を繰り返して落ちないこと。

## 検討して捨てた形

- テストの購読者を置いた後に関心を作り直す（`tracing::callsite::rebuild_interest_cache`）: 作り直した時点で呼び出し箇所が
  まだ登録されていなければ、その後に購読者の無い別のスレッドが登録して同じく「無い」になる

## 経緯

- 2026-11-14: #0285 のブランチのチェックで落ち、main でも再現した
- 2026-11-14: 計数するかを `Thermal::new` の引数にし（減速計算は共有の表があり DEBUG のとき。ログも同じ値で決める）、テストは
  tracing を介さずに計数を入れて閉じた。直した後に 25 回実行して落ちなかった
