---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0053: xs の入力はコマンドライン引数をやめて入力ファイルにする

## 背景

M1 では、核種ごとに ENDF、温度、熱中性子散乱のファイルを指定する。

## 検討した案

- コマンドライン引数（M0 の形）: 核種ごとの指定が書きにくい
- 入力ファイル（採用）

## 決定と理由

核種ごとに ENDF、温度、熱中性子散乱のファイルを指定するため、コマンドライン引数をやめて入力ファイルにする。

## 確認

`parses_sample`、`nuclides_need_temperatures`、`parses_input_relative_to_its_directory`。

## 帰結

入力ファイルの形式は TOML とした（[ADR-0010](0010-write-inputs-in-toml.md)）。
