---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0028: 炉心計算の固有値は非線形 CMFD で加速する

## 背景

炉心計算の固有値の反復を加速する方法を決める必要があった。M2 の NEM は CMFD を持たない。

## 検討した案

計画には、ほかの加速法を比べた記録はない。

## 決定と理由

固有値は非線形 CMFD で加速する。粗メッシュのソルバは lattice と共通のクレートのものを使う（[ADR-0003](0003-extract-shared-computation-on-second-use.md)）。

## 確認

[ADR-0232](0232-accelerate-nem-with-node-cmfd.md) で実装した。加速の有無で固有値と中性子束が許容の 2 倍の中で一致し、反復が減ることは
`cmfd_converges_to_the_unaccelerated_answer_in_fewer_iterations`（crates/core）が確かめる。

## 帰結

演算子と固有値の反復は分けておく（[ADR-0037](0037-separate-operators-from-iterations.md)）。粗メッシュの単位、補正、ADF と外周、粗メッシュの固有値問題と
許容、既定値、lattice との共有は [ADR-0232](0232-accelerate-nem-with-node-cmfd.md) で決めた。粗メッシュのソルバは共通のクレート `xlc-cmfd` に置き、lattice は
次に lattice の実行ファイルが変わる作業で移す（#0241）。
