---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0017: 本番のコードは panic しない

## 背景

数値計算の途中で panic すると、どの入力のどこで失敗したかの文脈が失われる。

## 検討した案

- `unwrap` や `assert!` で前提を確かめる
- 型で表し、表せなければエラーとして返す（採用）

## 決定と理由

本番のコードは panic しない（`unwrap`、`expect`、`panic!`、`unreachable!`、`assert!` を使わない）。起こりえない
場合は型で表せないかを先に考え（検査済みの値を enum や空でない型で持つなど）、表せなければエラーとして返す。

## 確認

clippy の lint でワークスペース全体を検査する（`Cargo.toml` の `[workspace.lints.clippy]` と `clippy.toml`。
テストでは許す）。

## 帰結

数値の失敗もエラーとして返す（`huge_count_is_error_not_abort`、`nan_pivot_is_singular` など）。
