---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: []
---

# ADR-0018: 計算結果は標準出力に、経過と診断は tracing のログとして標準エラーに出す

## 背景

計算結果と、経過・診断の出力を分け、性能の調査にも使えるようにしたい。

## 検討した案

- log と env_logger: 区間（span）を持たない
- tracing と tracing-subscriber（採用）

## 決定と理由

計算結果（k など）は標準出力に、経過と診断はログとして標準エラーに出す。ログは tracing とし、ライブラリは
マクロと区間（span）だけを使い、出力の設定は各コマンドで tracing-subscriber を使って行う。量は環境変数
`RUST_LOG` で選ぶ（既定は info: 段階ごとの時間と反復の回数。debug: 反復ごとの経過と内訳の時間。trace:
原因の調査用の細部）。log と env_logger より、区間で段階ごとの時間を測れるので性能の調査にも使える。

## 確認

約束事で、確かめるテストはない。

## 結果

なし。
