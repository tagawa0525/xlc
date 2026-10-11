---
status: open
labels: [xs, lattice]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# 計算の条件を環境変数でなく入力で受ける

## 問題

[ADR-0122](../adr/0122-run-trials-from-committed-inputs.md) は、試行を本番のコードに入力ファイルを与えて行い、環境変数で
切り替えないと決めた。examples にはまだ環境変数で計算や出力を変えるものが残っている（2026-10-12 に `std::env::var` を grep）。

- crates/xs/examples/compare_njoy_shielding.rs: `XLC_SLOWING_DOWN`、`XLC_URR_ONLY`
- crates/xs/examples/compare_njoy_gendf.rs、compare_njoy_thermr.rs: `XLC_ROWS`

環境変数は入力ファイルにも DVC の deps（ADR-0123）にも残らないので、同じ入力から違う結果が出る。clippy.toml の
`disallowed-methods` で `std::env` の読み取り（`var`、`var_os`、`vars`、`vars_os`）を禁じれば、新しく足すことも止められる。
`disallowed-methods` は変数の名前で許せないので、残す読み取りには場所ごとに `#[allow]` と理由を付ける。残す候補は
`RUST_LOG`（crates/{xs,lattice,core}/src/main.rs のログの設定）、`RAYON_NUM_THREADS`（crates/lattice/src/main.rs）、
crates/docs/tests/consistency.rs の `vars_os`。examples が引数で選ぶ計算の条件（crates/lattice/examples/homogeneous.rs の
`--swap` など）も入力ファイルで受けるかを一緒に決める。

関係する ADR: なし

## 完了の条件

方式を ADR にし、上の環境変数を入力で受ける形に直し、clippy の禁止で新しい環境変数を止める。

## 経緯

- 2026-10-11: ADR-0122 の後の作業として決めた（メモリに置いていたものを 2026-10-12 に issue にした）
