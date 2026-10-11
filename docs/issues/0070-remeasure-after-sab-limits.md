---
status: open
labels: [verification, xs, lattice]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# ADR-0133 のライブラリで、ピンセルの 1C・1D、集合体と他の例を測り直す

## 問題

[ADR-0133](../adr/0133-scatter-h-in-h2o-up-to-the-evaluated-sab-limit.md) で H in H2O の S(α,β) の上限、表の外の扱い、その上の自由ガスを変えた。
OpenMC と比べ直したのは、1B のピンセルと 1B を混ぜた均質な無限媒質だけ（examples/vera/verification.md「均質な無限媒質」
「k の差の内訳」）。次の記録は、ADR-0133 の前のライブラリ（S(α,β) の上限 4 eV）の値のまま。

- examples/vera/verification.md の「OpenMC との比較」の 1C、1D、2A〜2K、2Q と、その下の節（温度と支持格子、吸収体による
  差の変化、2K の釣り合い、離散化の誤差の基準など）
- `thermal` を使う他の例のライブラリ（examples/m1、mox、reflector、branches の xs.toml。既定値が変わる）と、その記録

関係する ADR: なし

## 完了の条件

上の記録を ADR-0133 のライブラリで測り直し、今の姿に書き直す。集合体は Slurm で 1 本ずつ流す。

## 経緯

- 2026-10-12: ADR-0133 のブランチで、1B と均質な媒質だけを測り直したので起こした（集合体は別の作業にするという指示）
