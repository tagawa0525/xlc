---
status: superseded
date: 2026-10-08
supersedes: [ADR-0056]
superseded-by: ADR-0082
issues: []
---

# ADR-0091: MF6 の角度分布（連続的な放出の Kalbach-Mann など、LAW=2 の重心系の分布）を散乱行列に使う

## 背景

M2 で MOC を非等方散乱で解き（[ADR-0075](0075-solve-moc-with-anisotropic-scattering.md)）、ライブラリの散乱行列を Legendre のモーメントで持つことにした（[ADR-0076](0076-store-legendre-scattering-moments.md)、
2026-10-07）。それまで MF6 の LAW=2 の角度分布は、P0 の近似と同じく無視していた（[ADR-0056](0056-treat-mf6-law2-law3-as-two-body.md)）。M2 の計画では [ADR-0076](0076-store-legendre-scattering-moments.md) と
同じ項目に書いていたものを、別の決定として分けた（番号はこの分割のときに振った）。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

MF6 の角度のデータも散乱行列の P1 以上に使う（docs/plans/m2-assembly.md、579cde7、2026-10-07）。

- MF6 の LAW=2, 3 の離散準位は、[ADR-0056](0056-treat-mf6-law2-law3-as-two-body.md) のとおり Q 値の2体反応とする。LAW=3 は重心系で等方、LAW=2 は重心系の角度
  分布を使う
- 連続的な放出のうち、角度分布を持つ形式（MF6 の Kalbach-Mann など）はその角度分布を使う。角度分布が等方のときだけ
  P1 以上が 0

## 確認

なし。実装していない。

## 結果

コードには入らなかった。MF6 の読み込みは LAW=2 の角度のデータを読み飛ばし（`crates/xs/src/endf/mf6.rs`）、Kalbach-Mann
の角度を扱うコードはない（`git log -S Kalbach -- crates` で見つからない）。M3 で、対象の核種を調べて MF6 の角度のデータを
使わないことにした（[ADR-0082](0082-leave-out-mf6-angular-data.md)、2026-10-08）。
