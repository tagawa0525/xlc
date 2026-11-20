---
type: Issue
status: closed
labels: [xs, depletion]
opened: 2026-11-20
closed: 2026-11-20
milestone:
blocked-by: []
resolution: ADR-0272
---

# ENDF/B-VII.0 の Rb-87 の MF5 MT17 の TAB1 の補間の区切り（NBT）が点の数と食い違い、xs が読めない

## 問題

- ENDF/B-VII.0 の Rb-87 の MF5 MT17 の、最初の入射エネルギー（18.8049 MeV）の出射エネルギーの分布の TAB1 は NR=1、NP=3 で、区切りが
  NBT=1（ファイルの 1193〜1194 行）。ENDF-102 の TAB1 は最後の NBT が NP に等しいとするので、xs は「last NBT must equal NP」で止まる
  （`xs nuclide` も `xs check` も読み込みで止まる）
- 4 つの版でこの形で止まるのは VII.0 の Rb-87 だけ（走査の数と一覧の出どころは examples/depletion/verification.md の「中性子のファイルの読み込みの走査」（`xs check`、ADR-0266））。#0292 の 9 核種には当たらない
- NJOY2016 は止まらず、範囲の番号が NR を超えると、範囲の表の後ろの語（最初の点の y = 0）を補間則として読む（コードを読んだもの。
  ADR-0272 の背景）。リポジトリの方針（怪しい入力は止める）では、黙って直して読まない
- 正しい区切りは 3 だけ（NR = 1）で、ENDF/B-VII.1 の同じ評価はこの欄だけを 3 に修正した
- 再現: `bin/xs check <報告> references/endf/ENDF-B-VII.0_neutrons/n-037_Rb_087.endf`

関係する ADR: ADR-0272

## 完了の条件

- [x] 評価の誤りとして扱いを決める（VII.0 の Rb-87 を輸送しない、正しい区切りを入力で与える、NJOY の読み方に合わせて warn を出す、のどれか）。
  正しい値を入力の核種の項目 `corrections` に宣言して読む（ADR-0272）

## 経緯

- 2026-11-20: 4 つの版の中性子のファイルの全部を `xs check` に通し、起こした
- 2026-11-20: 4 つの版のライブラリの全部の核種が流れるようにすることをユーザーが求めた（ユーザー 2026-11-20。別のセッションが伝えたもの）ので、
  輸送しない案は採らず、入力に修正を宣言する形にして閉じた（ADR-0272。宣言した `xs nuclide` の結果は examples/depletion/verification.md の
  「走査が見逃す種類」）
