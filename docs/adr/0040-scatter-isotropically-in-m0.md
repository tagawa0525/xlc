---
status: superseded
date: 2026-10-07
supersedes: []
superseded-by: ADR-0075
issues: []
---

# ADR-0040: 散乱は等方（P0）とし、輸送補正はしない

## 背景

M0 の縦切り（2026-10-06）で、散乱の非等方性をどこまで扱うかを決める必要があった。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

散乱は等方（P0）とし、輸送補正はしない。格子計算の MOC も散乱を等方とする。

## 確認

M0 の例（examples/m0/README.md）。

## 結果

M2 で非等方散乱に変えた（[ADR-0075](0075-solve-moc-with-anisotropic-scattering.md)、2026-10-07）。VERA 2B のピンの核分裂率が OpenMC より、水の穴の上下左右の隣で平均
+0.45%、斜めの隣で −0.16%、それ以外で −0.42% ずれ（rms 0.43%。2C と 2D も 0.44%）、水素の前方に偏る散乱を表せない
ことによると見込んだため（examples/vera/README.md、docs/plans/m2-assembly.md、579cde7、2026-10-07）。
