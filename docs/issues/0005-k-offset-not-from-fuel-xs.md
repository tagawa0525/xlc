---
status: open
labels: [xs, lattice, verification]
opened: 2026-10-08
closed:
---

# 燃料の断面積で説明できない k の差（VERA +295 pcm、MOX +180〜+240 pcm）

## 問題

VERA 1B の k は OpenMC より +295 pcm 高い（examples/vera/README.md）。MOX のピンセルでは、9.1 keV より下の燃料の
核種の断面積の差から見積もった k の差（MOX-1 −162 pcm、MOX-2 −2 pcm）を除いても、+241 pcm（MOX-1）と
+180 pcm（MOX-2）が残る（examples/mox/README.md）。同じ原因の疑いがある。

## 完了の条件

差の内訳（燃料以外の物質、9.1 keV より上の群、輸送の近似など）を分け、原因を特定する。

## 経緯

- 2026-10-08: MOX のピンセルで、燃料の断面積で説明できない分を見積もった
