---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0031: 参照解は同じ核データから手元で作った OpenMC で求め、OpenMC は本体の依存に含めない

## 背景

格子計算と炉心計算の参照解が要る。公開されている参照解（VERA の KENO-VI）は ENDF/B-VII.0 で計算されており、
核データの違いを含む。

- 手法の誤差だけを取り出せること
- 本体の依存を増やさないこと

## 検討した案

- VERA の公開値（KENO-VI、ENDF/B-VII.0）だけと比べる: 核データの差を含む
- 同じ ENDF/B-VIII.0 から手元の NJOY で作った連続エネルギーのデータで、OpenMC を手元で解く（採用）

## 決定と理由

参照解を作るためのモンテカルロコード（OpenMC など）は検証にだけ使う外部ツールとし、本体の依存には含めない。
ソースから `references/` にビルドする。

参照解は OpenMC で手元に作る。ENDF/B-VIII.0 から手元の NJOY で連続エネルギーのデータを作り、同じ組成・同じ温度で
計算する。核データの違いを除き、手法の誤差だけを取り出すためである。これはコード同士の検証（verification）であり、
実測との比較による妥当性確認（validation）ではない。VERA の公開値（KENO-VI、ENDF/B-VII.0）とも比べるが、核データの
差を含むので参考とする。

## 確認

手順は docs/verification.md の「VERA の参照解」、結果は examples/vera/README.md。

## 結果

- 3D の全炉心の OpenMC は計算時間が大きいので作らない。Problem 5 は仕様書の実測値（妥当性確認）と KENO の値で比べる
- 参照解の組成は計算と同じにする（[ADR-0054](0054-model-vera-materials.md)）
