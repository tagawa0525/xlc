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

# xs が MF6 の中性子の LAW=6（位相空間）と LAW=7（実験室系の角度とエネルギー）を読めない

## 問題

- crates/xs/src/endf/mf6.rs は中性子の出射の分布の LAW=1、2 だけを読み、LAW=6（n 体の位相空間の分布）と LAW=7（実験室系の角度と
  エネルギーの表）で「neutron LAW=6 is not supported」「neutron LAW=7 is not supported」で止まる
- 止まる評価（走査の数と一覧の出どころは examples/depletion/verification.md の「中性子のファイルの読み込みの走査」（`xs check`、ADR-0266））: LAW=6 は H-2 の MT16（4 つの版）と ENDF/B-VIII.1 の Li-6 の MT41、LAW=7 は Be-9 の MT16（4 つの版）。
  どれも燃焼チェーンの核種で、#0292 の 9 核種には当たらない
- 今の燃焼のライブラリ（65 核種）には入らない。2A のライブラリ（examples/depletion/p2/xs.toml）は H-2 と Be-9 を放射化の核種にしており、
  放射化の経路は MF6 を読まないので止まらない。輸送する核種にすると止まる
- NJOY2016 の GROUPR は LAW=6 と LAW=7 を扱う（references/NJOY2016/src/groupr.f90 の `getsix` の周り。式は読んでいない）
- 再現: `bin/xs check <報告> references/endf/ENDF-B-VIII.0_neutrons/n-001_H_002.endf`

関係する ADR: [ADR-0271](../adr/0271-group-mf6-law6-and-law7-by-endf-102.md)

## 完了の条件

- [ ] LAW=6 と LAW=7 を ENDF-102 の式で読み、散乱行列に入れる（NJOY と比べる）か、輸送しない核種に限ることを決める

## 経緯

- 2026-11-20: 4 つの版の中性子のファイルの全部を `xs check` に通し、起こした
