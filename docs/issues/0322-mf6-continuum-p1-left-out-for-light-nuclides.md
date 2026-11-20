---
type: Issue
status: open
labels: [xs]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# H-2、Li-6、Be-9 の MF6 の LAW=6、7 の反応の P1 以上を 0 にしている（ADR-0082）

## 問題

- MF6 の連続的な放出の P1 以上は 0 とする（ADR-0082）。ADR-0271 は LAW=6、7 にも同じ扱いを当てた
- ADR-0082 の調査（2026-10-08）は、重い核と O-16、Si-28 などで、重心の運動による P1 はわずかとした。#0311 で読めるようになった軽い核では、
  NJOY2016 GROUPR の P1 の行の和が P0 の 4〜9 割ある（H-2 の VIII.0 の MT16 で 0.89、Li-6 の VIII.1 の MT41 で 0.56、Be-9 の MT16 で 0.41〜0.51。
  examples/depletion/verification.md の「MF6 の LAW=6、7 の散乱行列（NJOY との比較）」）
- どれも 1.7 MeV 以上のしきい反応で、PWR の計算への影響は測っていない

関係する ADR: [ADR-0082](../adr/0082-leave-out-mf6-angular-data.md)、[ADR-0271](../adr/0271-group-mf6-law6-and-law7-by-endf-102.md)

## 完了の条件

- [ ] P1 以上を入れたときの影響（重水、Be、Li を含む体系、または全核種を輸送する燃焼）を測り、ADR-0082 を置き換えるかを決める

## 経緯

- 2026-11-20: #0311 で LAW=6、7 を読み、NJOY の P1 との差が大きいので、ADR-0082 を置き換えるかの材料として起こした
