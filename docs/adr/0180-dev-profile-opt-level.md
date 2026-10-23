---
status: accepted
date: 2026-10-21
requires: [ADR-0109, ADR-0157, ADR-0159]
supersedes: []
superseded-by:
issues: []
---

# ADR-0180: dev プロファイルを opt-level 1 でビルドしてテストを速くする

## 背景

`cargo test` が最適化なし（opt-level 0）の dev ビルドで走り、数値計算のテストが遅い。nextest で測ると、`xlc-xs` の
333 テストは 35.5 秒（`processes_each_temperature_and_sigma0` 35.4 秒、`free_gas_matrix_satisfies_detailed_balance`
31.5 秒）、ワークスペース全体の 1152 テストは 46.8 秒かかった。チェックと /git-merge のたびにこの時間を待つ。

`.cargo/config.toml` はすでに Zen 5 の命令と `-align-loops=64` でコンパイルする（ADR-0109、ADR-0157）。これは opt-level
とは独立で、最適化をかけたテストにもそのまま効く。

- 判断の軸は、チェックの待ち時間（実行とビルド）と、テストの意味を変えないこと

## 検討した案

- 案 A: `[profile.dev] opt-level = 1`（採用）
- 案 B: `[profile.dev] opt-level = 2`。実行は 7.8〜8.1 秒で案 A より約 10% 速いが、フルビルドが 48 秒で案 A の 34 秒より
  14 秒長い。実行の差が小さく、ビルドとデバッグ（変数が最適化で消えにくい）で案 A が有利
- 案 C: `[profile.test]` だけを最適化する。`cargo test` が使う設定は同じで効果も同じだが、設定が 2 か所になり、
  `cargo build` との差が増える。利点がないので採らない
- 案 D: 遅いテストの入力を縮める。テストが調べる内容（詳細つり合い、温度ごとの処理など）の網羅が減る。入力を変えずに
  済む案 A で足りるので採らない
- 案 E: テストの並列をさらに増やす。測定はすでに nextest（テストごとに並列）で、最も遅い 1 テストが 35 秒なので、
  並列を増やしても全体はそれより縮まない
- 案 F: 依存だけを最適化する（`[profile.dev.package."*"]`）。遅いのは依存でなくワークスペースの数値計算なので、効かない
  見込み。測っていない
- opt-level 3 は測っていない。opt-level 1 で足りる以上、ビルドが長くなる側に進む理由がない

## 決定と理由

案 A。同じ機械（手元の Ryzen 9 9950X）で、nextest のワークスペース全体（実行した 1152 テストと、元から skip の 1 テスト）を測った。
測定は main の `a74381e90` に Cargo.toml の変更を足した worktree（2026-10-21）。

| opt-level | テストの実行 | フルビルド（`cargo test --no-run`） |
| --------- | ------------ | ----------------------------------- |
| 0（既定） | 46.8 秒      | 測っていない                        |
| 1         | 8.6 秒       | 34 秒                               |
| 2         | 7.8〜8.1 秒  | 48 秒                               |

`xlc-xs` だけ（333 テスト）では opt-level 0 が 35.5 秒、opt-level 1 が 6.9 秒（最長は `free_gas_matrix_satisfies_detailed_balance`
6.8 秒）、opt-level 2 が 6.8 秒。上の時間は 16 コアでの並列の実行時間で、1 コアに限った時間は測っていない。`debug-assertions` と
`overflow-checks` は dev の既定（有効）のまま残るので、オーバーフローと `debug_assert!` の検査は失われない。

## 確認

opt-level 1 でワークスペースの実行した 1152 テストがすべて通る（skip の 1 テストは opt-level 0 のときと同じ）。許容誤差は
変えていない。opt-level 0 と 1 で出力の数値そのものは比べていない。

## 帰結

- `cargo test`、`run-checks.sh`、/git-merge のチェックが、同じ設定で速くなる
- debug の `cargo build` も opt-level 1 になる。デバッガで変数が消えるときは `--config 'profile.dev.opt-level=0'` で戻す
- ステージの実行ファイルは release（`target/dvc`、ADR-0159）でビルドするので、この設定の影響を受けない見込み。
  `dvc status` で確かめてはいない（worktree には `bin/` がなく比べられない）
- 本番の release と dev で最適化の度合いが違う点は変わらない。release の結果を確かめるテストは、この ADR の対象外
