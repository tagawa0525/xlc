---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0021: 核データ処理は NJOY と同じ流れとし、共鳴の形式は対象核種で使うものから順に対応する

## 背景

ENDF から多群ライブラリを作る処理の流れと、対応する共鳴パラメータの形式の範囲を決める必要があった。

- NJOY とステージごとに比べて検証できること
- 必要になるまで作らないこと

## 検討した案

計画には、ほかの流れを比べた記録はない。

## 決定と理由

処理の流れは NJOY と同じく、共鳴の再構成 → ドップラー拡がり → 熱中性子散乱 → 群定数化 とする。ステージごとに NJOY の
同じ条件の出力と比べられる。

共鳴パラメータの形式（MLBW / Reich-Moore / R-matrix limited）は、対象核種で実際に使われているものから順に対応する。

- M0: 対象核種（H-1, O-16, U-235, U-238）の分離共鳴は Reich-Moore（LRF=3）だけなので、それだけを扱う
- M1: 対象核種の分離共鳴は Reich-Moore（U-235, U-238）と MLBW（U-234, U-236, Zr）なので、MLBW を追加する。それ以外の
  形式は、対象核種で必要になるまで追加しない
- M2: Fe-54 と Fe-57 の R-matrix limited 形式（LRF=7）を追加する（[ADR-0055](0055-reconstruct-r-matrix-limited-subset.md)）
- M3: 新しい核種（C、Si、P、Ti、Mn、Ni、Ag、Cd、In）の共鳴とエネルギー分布の形式は、見出しを調べた範囲
  （2026-10-08）では対応済みのものだけである。実際に処理して足りない形式が出たら、そのとき対応する

## 確認

ステージごとに NJOY の出力と比べる例（`crates/xs/examples/compare_njoy_pendf.rs`、`compare_njoy_broadr.rs`、
`compare_njoy_thermr.rs`、`compare_njoy_gendf.rs`、`compare_njoy_shielding.rs`）。結果は examples/m1/verification.md。

## 帰結

対応していない形式はエラーにする（`unsupported_formalism_is_error`）。
