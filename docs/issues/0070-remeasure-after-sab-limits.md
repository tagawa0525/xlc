---
status: open
labels: [verification, xs, lattice]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# ADR-0133 のライブラリで、ピンセルの 1C・1D、集合体と他の例を再測定する

## 問題

[ADR-0133](../adr/0133-scatter-h-in-h2o-up-to-the-evaluated-sab-limit.md) で H in H2O の S(α,β) の上限、表の外の扱い、その上の自由ガスを変えた。
再測定した記録は、examples/vera/verification.md の 1B、均質な無限媒質、k の差の内訳（1B と 2K）、1C〜2Q の「OpenMC との比較」、
温度と支持格子、吸収体による差の変化と、方位角 128・間隔 0.01 cm との差、examples/m0、m1、mox の verification.md（コミット e4dbc0de）。

まだ再測定していないもの（ADR-0133 の前のライブラリの値のまま、または未測定と書いたもの）:

- examples/vera/verification.md の、吸収体の集合体のトラック間隔のスイープ（2E-empty、2F、2G、2H、2K の方位角 30、間隔 0.02〜0.05 cm の
  6 点ずつ計 30 本。表は「未測定」と書いた）と、吸収体のリングの数、SS304 の被覆管と案内管の間の水のリングの数による k の動き
- examples/vera/verification.md の、「支持格子の枠と構造材の自己遮蔽」「散乱の次数」「同じ多群定数での OpenMC との比較」「離散化の誤差」
  （2B と 1B の方位角と間隔の表、メッシュ）「B1 の臨界スペクトルの拡散係数」「対称性」の表の xlc の値（それぞれ別の入力と設定で測った
  記録で、数値を出したコミットが書いてある）
- examples/branches/verification.md（分岐のテーブルの 27 状態で約 2 時間、直接の格子計算 6 本。共有モードでは約 10 時間）
- examples/reflector/verification.md（dvc.yaml の全段階、排他で約 3 時間。#0034 のマージ後に `dvc repro examples/reflector/dvc.yaml`
  を流し、lattice-I、lattice-corner、core-corner* の比較を今の姿に書き直す。main の `dvc checkout` は `--force` が要る）

関係する ADR: ADR-0133（フロントマターの `issues` にはこの issue を挙げていない。マージ済みで書き換えない）

## 完了の条件

上の記録をすべて ADR-0133 のライブラリで測り直し、今の姿に書き直す。集合体は Slurm で 1 本ずつ流す。

## 経緯

- 2026-10-12: ADR-0133 のブランチで、1B と均質な媒質だけを再測定したので起こした（集合体は別の作業にするという指示）
- 2026-10-12: feat/0070-remeasure で、vera の 1C〜2Q、2K の内訳、方位角 128 の参照、m0、m1、mox を再測定した（k の差は 2G を除く
  15 体系で +157〜+347 pcm から −69〜+108 pcm、2G は +27 から −212 pcm）。OpenMC の本番（共有モードの 24 コア）が走る間は
  排他のジョブが始まらないので、トラック間隔のスイープ（30 本）、branches、reflector を残した
