---
status: closed
labels: [xs]
opened: 2026-10-08
closed: 2026-10-08
milestone: M3
blocked-by: []
resolution: ADR-0100
---

# 運動学のしきい値より下から始まる反応の扱いを決める

## 問題

ENDF の MF3 のしきい反応には、最初の点が運動学のしきい値 −Q (A + 1)/A より下にあるものがある（Cd-112 の MT51 は
623.0002 keV から始まるが、Q と AWR からは 623.5708 keV）。NJOY RECONR は最初の点を運動学のしきい値に上げ
（"changed threshold"）、xlc はファイルのとおりに扱う。xlc の2体の散乱行列は運動学のしきい値より下で中性子を出さず、
P0 の行の和を断面積に合わせ直す。

M3 の新しい核種の無限希釈の群定数では、しきい値を含む群でこの差が出る（Cd-112 の MT51 の 6 群で 3.0e-3、
ほかのしきい反応で最大 1.5e-2。examples/m3-nuclides/verification.md）。NJOY の PENDF から xlc で群平均すると差は消え、
xlc の値は ENDF の表の厳密な積分と 7 桁一致する。どちらが正しいかは評価済みデータの不整合の扱いの問題で、決めていない。

関係する ADR: [ADR-0100](../adr/0100-move-tables-to-kinematic-thresholds.md)

## 完了の条件

運動学のしきい値より下の断面積を、ファイルのとおり残すか、NJOY のように上げるか、エラーにするかを決めて ADR にする。

## 経緯

- 2026-10-08: M3 の新しい核種のライブラリを NJOY と比べて見つかった
- 2026-10-08: しきい値より下の断面積はありえないので、最初の点をしきい値に移して warn を出し、しきい値以下のほかの点は
  捨てることにした（NJOY RECONR と同じ。NJOY はほかの点をしきい値のすぐ上に詰める）。4cdc221 で実装し、Cd-112 の MT51 の
  6 群の差は 3.0e-3 から 1.1e-4 になった（examples/m3-nuclides/verification.md）
- 2026-10-09: 方針が認められ、[ADR-0100](../adr/0100-move-tables-to-kinematic-thresholds.md) にした。警告は核種ごとに1回にまとめた
