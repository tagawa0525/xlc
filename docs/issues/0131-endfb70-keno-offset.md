---
status: open
labels: [xs, vera, verification]
opened: 2026-10-15
closed:
milestone:
blocked-by: []
resolution:
---

# ENDF/B-VII.0 のライブラリでの KENO-VI との k の差が、VIII.0 での OpenMC との差より約 +130 pcm 大きい

## 問題

ADR-0175 の検証用ライブラリ（ENDF/B-VII.0。examples/vera/endfb70/xs.toml）で xlc が止まった核種はなく、1B〜1D、2A〜2K、2Q、
3A、3B を解けた。ただし xlc − KENO-VI の k の差は、ENDF/B-VIII.0 で xlc − OpenMC が示す手法の差より大きい
（examples/vera/verification-endfb70.md、数値は 43661306、2026-10-15）。

- 15 の格子（ピンセルと集合体）で、VII.0 の xlc − KENO-VI は平均 +165 pcm（+115〜+270）。VIII.0 の xlc − OpenMC は平均 +39 pcm
  （−69〜+108）。核データの版を揃えれば、差は VIII.0 の OpenMC との差と同程度になると見込んでいたが、そうならなかった
- 同じ xlc の k が、ライブラリを VIII.0 から VII.0 に替えると平均 +75 pcm 動く（−78〜+173）。OpenMC（VIII.0）と KENO-VI（VII.0）の
  k の差から見込む動き（KENO-VI − OpenMC）は平均 −52 pcm（−129〜+37）で、差し引くと平均 +127 pcm が説明できない。
  1B、1C、1D、2A〜2D、2I、2Q では xlc の動きが +98〜+173 pcm
- 1B〜1D は 600 K 以上の等温で、S(α,β) の温度の内挿も天然の C も効かない（C を含まない）。ほかの核種は燃料（U、O）と被覆管、水

## 候補（確かめていない）

- xlc の VII.0 のファイルの扱い（共鳴の再構成、ドップラー拡がり、非分離共鳴、MF4/MF5/MF6 の違い）。VII.0 でだけ起きる処理の誤り
- KENO-VI の連続エネルギーのライブラリ（AMPX）と NJOY の処理の差。同じ評価済みデータでも処理が違う
- KENO-VI の入力（公開されていない）とモデルの解釈の差
- 仕様書の KENO-VI が 565 K（2A、2Q）で S(α,β) の温度を補正していること

## 完了の条件

NJOY で処理した ENDF/B-VII.0 で OpenMC を解き（2D の小さな問題、1B と 2B など。ADR-0175 が許す診断）、OpenMC（VII.0）− KENO-VI が
小さければ原因は xlc の VII.0 の処理、大きければ KENO-VI との処理の差かモデルの解釈、と分けられる。分けた結果を ADR と記録に書く。

関係する ADR: ADR-0175
