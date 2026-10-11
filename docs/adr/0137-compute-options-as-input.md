---
status: accepted
date: 2026-10-12
requires: [ADR-0067, ADR-0122, ADR-0123]
supersedes: []
superseded-by:
issues: ["#0074"]
---

# ADR-0137: 計算の条件を環境変数でなく入力で受け、環境変数の読み取りを clippy で禁じる

## 背景

ADR-0122 は、試行を環境変数で切り替えず、本番のコードに入力ファイルを与えて行うと決めた。環境変数は入力ファイルにも
DVC の deps（ADR-0123）にも残らないので、同じ入力から違う結果が出る。2026-10-12 に `std::env::var` を grep すると、
examples に次の読み取りが残っていた。

- compare_njoy_shielding: `XLC_SLOWING_DOWN`（分離共鳴の重みを超微細群の減速計算で作る）、`XLC_URR_ONLY`（比べる群を
  非分離共鳴にかかる群に絞る）
- compare_njoy_gendf、compare_njoy_thermr: `XLC_ROWS`（全群の値を表示する）

ADR-0122 は禁止を文章にしただけで、新しい読み取りを止める仕組みがない。

## 検討した案

- 案 A: 環境変数を残し、文書に書く。ADR-0122 が退けた形のままで、新しい読み取りも止まらない
- 案 B: examples の選択肢を入力の TOML で受ける。比べる道具の入力は NJOY の出力のパスと核種で、位置の引数で受けている。
  選択肢だけのために TOML を足すと、入力が 2 つに分かれる。DVC の段階にする比べる道具は、`cmd` に選択肢が書かれ
  （`dvc.lock` にも残る）、再現できる
- 案 C: 選択肢を `--` で始める引数で受ける。引数は `cmd` と履歴に残り、環境のように黙って残らない。`--` で始まる引数で受ける点は
  `crates/lattice/examples/homogeneous.rs` の `--balance`、`--swap` と同じ（そちらは値を取る引数で、順も決まっている。
  比べる道具の選択肢は真偽だけなので、位置の引数から分けて読む）
- 案 D: 環境変数の読み取りを `clippy.toml` の `disallowed-methods` で禁じる。変数の名前では許せないので、読む場所ごとに
  `#[allow(clippy::disallowed_methods, reason = …)]` を付ける

## 決定と理由

案 C と案 D を採る。

- examples の選択肢は引数にする: `--slowing-down`、`--urr-only`（compare_njoy_shielding）、`--rows`（compare_njoy_gendf、
  compare_njoy_thermr）。知らない `--` の引数はエラーにする
- 本番の 3 つの実行ファイル（xs、lattice、core）の計算の条件は、これまでどおり入力の TOML に書く
- `std::env` の読み取り（`var`、`var_os`、`vars`、`vars_os`）を禁じる。残す読み取りは、計算の結果を変えないものに限り、
  読む関数に `#[allow]` と理由を付ける:
  - `RUST_LOG`（xs、lattice、core の `init_logging`）: 表示する量だけを選ぶ
  - `RAYON_NUM_THREADS`（lattice の `use_physical_cores`）: 並列の数だけを選ぶ（ADR-0067。結果が数によらないことは #0060）
  - crates/docs/tests/consistency.rs の `vars_os`: テスト用のリポジトリから、呼び出し元の `GIT_*` を外す

`RAYON_NUM_THREADS` は、結果がスレッド数によらないことを #0060 が確かめるまで、計算の条件でない前提で許す（確認前の
許しなので、`#[allow]` の理由にもそう書く）。`homogeneous.rs` の `--balance`、`--swap` は、すでに引数なので、そのままにする。

## 帰結

- 新しい環境変数の読み取りは、clippy（`run-checks.sh`）で止まる。許すには、理由を書いた `#[allow]` が要り、レビューで見える
- clippy が止めるのは `std::env` の読み取りで、`env!`、`option_env!` と、依存のクレートが自分で環境変数を読む場合は止めない。
  後者は、読むクレートを入れる時点のレビューで見る
- 引数の選択肢は `--` で始まる引数を位置の引数から分ける。引数の数や順は変わらない
- examples/m3-nuclides/verification.md の `XLC_SLOWING_DOWN` は `--slowing-down` に直した

### 確認

`cargo clippy --workspace --all-targets` が通ること（`run-checks.sh`）。禁止が効くことは、examples に `std::env::var` を
仮に足して clippy が `disallowed_methods` で警告することで確かめた（2026-10-12。足した行は戻した）。
