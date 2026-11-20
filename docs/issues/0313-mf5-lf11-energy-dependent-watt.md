---
type: Issue
status: closed
labels: [xs, depletion]
opened: 2026-11-20
closed: 2026-11-20
milestone:
blocked-by: []
resolution: コミット cfc8f6943
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

- [x] LF=11 を ENDF-102 の式で読み、NJOY と比べる（examples/depletion/verification.md の「版ごとのライブラリのエネルギーによる Watt の核分裂
  スペクトル」。χ の差は 4.7e-7 以下）

## 経緯

- 2026-11-20: 4 つの版の中性子のファイルの全部を `xs check` に通し、起こした
- 2026-11-20: 4 つの版の全部の核種が流れるようにすることが求められ（ユーザー 2026-11-20。別のセッションが伝えたもの）、読むことにした。a(E) を θ とする x = E'/a の
  累積の差で群に分け、x < 1 は桁落ちしない正の級数にした（cfc8f6943）。閉じた式が a b/4 の大きい評価で桁落ちするので、
  精度を確かめた a b/4 ≤ 5 の評価だけを読む（62814be3b。U-233 は 0.71 以下）。読み込みの走査で VII.0 と VII.1 の U-233 が「読めた」に移り、
  600 K の `xs nuclide` は終了コード 0 で終わる。ADR は無い（ENDF-102 の式のとおりで、累積の式は NJOY の anased と同じ。anased が √(E'/a) < 0.01 で級数の初項だけを使う
  近似は、相対 3.5e-5 以下の打ち切りの差で、誤りでないので真似ず warn も出さない。/adr の「ADR にするもの」に当たらない）
