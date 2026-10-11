---
status: accepted
date: 2026-10-07
supersedes: [ADR-0040]
superseded-by:
issues: []
requires: []
---

# ADR-0075: MOC を非等方散乱で解き、既定ではライブラリにある全ての次数（P3）を使う

## 背景

M0 から散乱を等方（P0）とし輸送補正もしていなかった（[ADR-0040](0040-scatter-isotropically-in-m0.md)）。VERA 2B のピンの核分裂率が OpenMC より、水の穴の
上下左右の隣で平均 +0.45%、斜めの隣で −0.16%、それ以外で −0.42% ずれた（rms 0.43%。2C と 2D も 0.44%。
examples/vera/verification.md、docs/plans/m2-assembly.md、579cde7 と aad5bc7、2026-10-07）。自己遮蔽とトラックが原因でないことを
確かめたうえで、水素の前方に偏る散乱を表せないことによると見込んだ。

- 精度（近似を持ち込まないこと）
- 時間とメモリ

## 検討した案

- 輸送補正（TCP0）: MOC の時間が変わらず CASMO などで広く使われるが、近似なので採らない
- 非等方散乱の散乱源を球面調和関数で展開する（採用）
- スイープで区間ごとにモーメントの和を取る: モーメントの数に比例して遅くなる

## 決定と理由

MOC を非等方散乱で解く。散乱源を、z 方向に対称な実の球面調和関数で L 次まで展開し、領域ごとに (L + 1)(L + 2)/2 個の
モーメントで持つ。次数は格子の入力で選べ、省略するとライブラリにある全ての次数（P3）を使う。VERA 2B で P1〜P3 の
k の幅は 10 pcm、ピンの核分裂率でほぼ 0 で、時間の差も小さい（MOC 60 → 71 秒。docs/plans/m2-assembly.md、3eb2d9a と
12b039e、2026-10-08）ため。

スイープでは、方位角ごとに方向別の線源を一度作り、減衰の和も方向別に貯めて方位角の終わりにモーメントへ畳み込む。区間ごとの
計算が P0 と同じになり、2B の MOC は P0 の 46 秒から P3 で 71 秒に収まる。方位角ごとの集計がモーメントの数だけ増え、
2B のメモリは 9.2 GB から P3 で 16.3 GB になる（12b039e、2026-10-08）。

## 確認

`anisotropy_does_not_change_an_infinite_medium`、`forward_scattering_acts_like_the_transport_correction`、
`addition_theorem_holds`、`anisotropic_result_does_not_depend_on_the_thread_count`。VERA 2B〜2D の水の穴のまわりのピンの
核分裂率を OpenMC と比べる（examples/vera/verification.md）。

## 結果

ライブラリの散乱行列のモーメントは [ADR-0076](0076-store-legendre-scattering-moments.md)。SMT のスレッドでのメモリは [ADR-0067](0067-default-threads-to-physical-cores.md)。
