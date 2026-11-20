---
type: Issue
status: closed
labels: [xs, depletion]
opened: 2026-11-20
closed: 2026-11-20
milestone:
blocked-by: []
resolution: コミット 967431abe
---

# xs が核分裂の中性子の分布を MF6 MT18 で与えた評価（Th-232、Pa-231、Pa-233）を輸送する核種として処理できない

## 問題

- Th-232、Pa-231、Pa-233 は 4 つの版のどれでも、核分裂の中性子の出射の分布を MF5 MT18 でなく MF6 MT18 で与える（ファイルの MF5 の
  節に MT18 が無い。2026-11-20 に references/endf の 4 つの版のファイルの制御欄で確かめた）
- xs は χ の材料を MF5 MT18 から取り（crates/xs/src/process.rs の `fission_parts`）、MF6 の核分裂を出射の分布に入れない（`read_distributions` は
  中性子を放出する反応のうち `emits_neutrons` のものだけを読む）。核分裂の断面積があるので、温度ごとの処理で「fission needs MF1 MT452 and
  MF5 MT18」で止まる（`xs check` も同じ条件で止める）
- 止まる評価は 走査の数と一覧の出どころは examples/depletion/verification.md の「中性子のファイルの読み込みの走査」（`xs check`、ADR-0266）。#0292 の 9 核種には当たらない。2A のライブラリでは放射化の核種で、放射化の経路は χ を作らないので止まらない
- 再現: `bin/xs check <報告> references/endf/ENDF-B-VIII.0_neutrons/n-090_Th_232.endf`

関係する ADR: ADR-0186

## 完了の条件

- [x] MF6 MT18 の核分裂の中性子の分布から χ を作る（NJOY と比べる）か、これらの核種を輸送しないことを決める（作った。比べた結果は
  examples/depletion/verification.md の「版ごとのライブラリの MF6 MT18 の核分裂スペクトル」）

## 経緯

- 2026-11-20: 4 つの版の中性子のファイルの全部を `xs check` に通し、起こした
- 2026-11-20: 4 つの版の全部の核種が流れるようにすることが求められ（ユーザー 2026-11-20。別のセッションが伝えたもの）、χ を作ることにした。
  12 の評価の MF6 MT18 は JP=0（ENDF-102 の 6.3.9.3 の平均の分布、y(E) = 1）なので、ADR-0186 の即発の分布の出どころを MF5 MT18 か JPN=0 の
  MF6 MT18 にした（ADR-0186 の「変更」。成分と重みは変えない）。読み込みの走査で 12 の評価が「読めた」に移り、
  600 K の `xs nuclide` は終了コード 0 で終わる。NJOY との χ の差は重心系から実験室系への変換に残り、xs の側の候補は #0325 に分けた
