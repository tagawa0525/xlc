---
type: Issue
status: open
labels: [xs, depletion]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# xs が MF5 の LF=11（エネルギーによる Watt 分布）を読めない（ENDF/B-VII.0、VII.1 の U-233）

## 問題

- crates/xs/src/endf/mf5.rs は LF=1、5、7、9、12 を読み、LF=11（Watt の分布 exp(−E'/a) sinh(√(b E'))、a と b が入射エネルギーの表）で
  「LF=11 is not supported」で止まる
- 止まる評価は ENDF/B-VII.0 と VII.1 の U-233 の MF5 MT18（走査の数と一覧の出どころは examples/depletion/verification.md の「中性子のファイルの読み込みの走査」（`xs check`、ADR-0266））。VIII.0 と VIII.1 の U-233 は止まらない。#0292 の 9 核種には当たらない
- NJOY2016 の GROUPR は LF=11 を扱う（references/NJOY2016/src/groupr.f90 の `getsed` の LF=11 の分岐。式は読んでいない）
- 再現: `bin/xs check <報告> references/endf/ENDF-B-VII.0_neutrons/n-092_U_233.endf`

関係する ADR: なし

## 完了の条件

- [ ] LF=11 を ENDF-102 の式で読み、NJOY と比べる

## 経緯

- 2026-11-20: 4 つの版の中性子のファイルの全部を `xs check` に通し、起こした
