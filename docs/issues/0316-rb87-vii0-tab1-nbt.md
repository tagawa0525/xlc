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

# ENDF/B-VII.0 の Rb-87 の MF5 MT17 の TAB1 の補間の区切り（NBT）が点の数と食い違い、xs が読めない

## 問題

- ENDF/B-VII.0 の Rb-87 の MF5 MT17 の、最初の入射エネルギー（18.8049 MeV）の出射エネルギーの分布の TAB1 は NR=1、NP=3 で、区切りが
  NBT=1（ファイルの 1193〜1194 行）。ENDF-102 の TAB1 は最後の NBT が NP に等しいとするので、xs は「last NBT must equal NP」で止まる
  （`xs nuclide` も `xs check` も読み込みで止まる）
- 4 つの版でこの形で止まるのは VII.0 の Rb-87 だけ（走査の数と一覧の出どころは examples/depletion/verification.md の「中性子のファイルの読み込みの走査」（`xs check`、ADR-0266））。#0292 の 9 核種には当たらない
- NJOY がこの TAB1 をどう読むかは確かめていない。リポジトリの方針（怪しい入力は止める）では、黙って直して読まない
- 再現: `bin/xs check <報告> references/endf/ENDF-B-VII.0_neutrons/n-037_Rb_087.endf`

関係する ADR: ADR-0272

## 完了の条件

- [ ] 評価の誤りとして扱いを決める（VII.0 の Rb-87 を輸送しない、正しい区切りを入力で与える、NJOY の読み方に合わせて warn を出す、のどれか）

## 経緯

- 2026-11-20: 4 つの版の中性子のファイルの全部を `xs check` に通し、起こした
