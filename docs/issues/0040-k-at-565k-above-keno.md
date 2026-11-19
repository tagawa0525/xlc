---
type: Issue
status: closed
labels: [verification, lattice]
opened: 2026-10-09
closed: 2026-11-19
milestone: M3
blocked-by: []
resolution: "wontfix: examples/vera の同じライブラリと入力の形の 2A（565 K）と 2B（600 K）が KENO-VI と +26、+19 pcm で、OpenMC の 2A とも +18 pcm。170 pcm は ADR-0209 の前の古い設定の examples/branches の状態 13 だけに残り、作り直すまで原因を分ける対象が無い"
---

# 565 K の集合体の k∞ が KENO の Problem 2A より約 170 pcm 大きい

## 問題

examples/branches のブランチのテーブルの状態 13（1300 ppm、燃料・被覆管・減速材とも 565 K、減速材の密度 0.74271
g/cm³）は、VERA の Problem 2A（565 K、0.743 g/cm³、1300 ppm。仕様書の表 P2-5）と同じ形状と条件にあたる。xlc の
k∞ は 1.184552 で、KENO-VI の 1.182175 ± 0.000017 より反応度で +170 pcm 大きい。600 K の Problem 2B では、xlc
（examples/vera、1.183868）は KENO（1.183360）より +36 pcm だった。

等温で 565 K から 600 K に上げたときの KENO の反応度の変化は +85 pcm（約 +2.4 pcm/K）。xlc の 565 K と 600 K の
値は別のライブラリと別の減速材の組成の入力から求めたので、xlc の同じ変化はまだ求めていない。テーブルの
560〜570 K の等温温度係数は +4.9 pcm/K（examples/branches/verification.md）。

関係する ADR: なし

## 疑うところ

- H-1 の S(α,β) の温度の内挿（ADR-0080）。565 K は 550 K と 573.6 K の表の間にあり、600 K は表の上にある
- KENO は S(α,β) を温度で内挿できない（仕様書の付録 G）。Problem 2A の「565 K の結果」（仕様書の改訂 2 の変更点）で
  H-1 の S(α,β) をどの温度で扱ったかは書かれていない。付録 G は、565 K の炉心の KENO の固有値を、S(α,β) が温度に
  よらないことに対して約 −42 pcm 補正するとしている（炉心の問題の値。集合体の問題に同じ補正をしたかは書かれていない）
- 減速材の組成の作り方の違い（ADR-0106 のホウ酸水と、examples/vera の仕様書の表の組成）
- ライブラリの違い（examples/branches の 550〜580 K と、examples/vera の 600 K）

## 完了の条件

- 同じライブラリと同じ入力の形で 565 K と 600 K の k∞ を求め、Problem 2A と 2B の両方と比べる
- 565 K の OpenMC の参照（#0017）で、差が xlc の側にあるかを切り分け、原因と大きさを検証の記録に書く

## 経緯

- 2026-10-09: ブランチのテーブルの確認のレビューで、仕様書の参照値と比べていないことが指摘され、比べて分かった
- 2026-11-19: 閉じた。examples/vera/verification.md の「OpenMC との比較」（xlc 05714cea、2026-10-23）で、同じライブラリと入力の形の 2A は
  1.182543（KENO-VI 1.182175 より +26 pcm）、2B は 1.183627（KENO-VI 1.183360 より +19 pcm）。565 K の OpenMC の 2A は KENO-VI より +19 pcm で、
  xlc − OpenMC は +18 pcm。600 K から 565 K への変化の差は −4 pcm（同じ記録の「温度と支持格子による差の変化」）。完了の条件の 2 つは
  そこで満たした。170 pcm は examples/branches の状態 13（c6ffe4f、2026-10-09。方位角 32・間隔 0.05 cm、ADR-0133 と ADR-0209 の前の
  ライブラリと自己遮蔽）だけで、その記録は作り直していない（記録の冒頭）。作り直す作業が来たときに、2A と比べ直す
