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
測り直した記録は、examples/vera/verification.md の 1B、均質な無限媒質、k の差の内訳（1B と 2K）、1C〜2Q の「OpenMC との比較」、
温度と支持格子、吸収体による差の変化と、方位角 128・間隔 0.01 cm との差、examples/m0、m1、mox の verification.md（コミット 1f4fcaaf）。

まだ測り直していないもの（ADR-0133 の前のライブラリの値のまま、または未測定と書いたもの）:

- examples/vera/verification.md の、Pyrex の環のリングの数と、SS304 の被覆管と案内管の間の水のリングの数による k の動き
  （吸収体の集合体のトラック間隔のスイープは 2026-10-14 に測り直した。下の経緯）
- examples/vera/verification.md の、「支持格子の枠と構造材の自己遮蔽」「散乱の次数」「同じ多群定数での OpenMC との比較」「離散化の誤差」
  （2B と 1B の方位角と間隔の表、メッシュ）「B1 の臨界スペクトルの拡散係数」「対称性」の表の xlc の値（それぞれ別の入力と設定で測った
  記録で、数値を出したコミットが書いてある）
- examples/branches/verification.md（分岐のテーブルの 27 状態で約 2 時間、直接の格子計算 6 本。共有モードでは約 10 時間）
- examples/reflector/verification.md（dvc.yaml の全段階、排他で約 3 時間。#0034 のマージ後に `dvc repro examples/reflector/dvc.yaml`
  を流し、lattice-I、lattice-corner、core-corner* の比較を今の姿に書き直す。main の `dvc checkout` は `--force` が要る）

関係する ADR: なし

## 完了の条件

上の記録をすべて ADR-0133 のライブラリで測り直し、今の姿に書き直す。集合体は Slurm で 1 本ずつ流す。

## 経緯

- 2026-10-12: ADR-0133 のブランチで、1B と均質な媒質だけを測り直したので起こした（集合体は別の作業にするという指示）
- 2026-10-12: feat/0070-remeasure で、vera の 1C〜2Q、2K の内訳、方位角 128 の参照、m0、m1、mox を測り直した（k の差は 2G を除く
  15 体系で +157〜+347 pcm から −69〜+108 pcm、2G は +27 から −212 pcm）。OpenMC の本番（共有モードの 24 コア）が走る間は
  排他のジョブが始まらないので、トラック間隔のスイープ（30 本）、branches、reflector を残した
- 2026-10-14: feat/0070-remainder で、トラック間隔のスイープと branches の測り直しを始めたが、格子の自己遮蔽の高速化（出力が変わる案を
  検討中）が先なので、調整役が待ち行列のジョブを取り消して中断した。測れたのは 2E-empty の方位角 128・間隔 0.01 cm の k = 1.181466
  だけで、記録の参照（1.181466）と一致した（8271ca18、ADR-0142 の後のコード）。ほかは未測定。
  再開の手順（高速化のマージ後。出力が変わるなら、上の測定をすべて測り直す）:
  - 入力は 8271ca18 の trials/disc/（`lattice_variant.py` で p2/assembly-{2E-empty,2F,2G,2H,2K}.toml から作った 40 本。方位角 30 の間隔
    0.02、0.025、0.03、0.035、0.04、0.045、0.05 cm と、方位角 128・間隔 0.01 cm）。再開時は main に rebase し、ライブラリと
    パスを `lattice_variant.py` で作り直す（入力にライブラリと出力の絶対パスが入る）
  - スイープ: 問題ごとに 1 ジョブ（8 本を逐次、`-t 45`、名前 `vera-disc-<問題>`）を 5 本、同時に投げて待ち行列に並べる
    （1 本ずつ投げ直すと、ほかのジョブの後ろに並び直して 1 本あたり数時間待つ）。k は標準出力の `k <値>`
  - branches: `branches-xs-table`（xs 13 分と assembly-2B.toml の 27 状態の約 2 時間、`-t 135`）のあとに、direct 6 本と core 7 本
    （`-t 70`、1 ジョブ）。dvc.yaml はないので手で流す（verification.md の手順）
  - 2G/2H の行はスイープの結果が出たあとに、ADR-0142（51eb3e73）の後の値に揃える
- 2026-10-14: feat/0070-remainder で、トラック間隔のスイープ（2E-empty、2F、2G、2H、2K の方位角 30 × 間隔 7 種と方位角 128・間隔 0.01 cm、
  計 40 状態）を測り直した（cedff2cf、分岐点 39e501e3、排他のジョブ 355〜359、各 4〜9 分）。幅は 2E-empty 7.0、2F 28.9、2G 65.8、
  2H 33.1、2K 29.3 pcm で、2G と 2H の行も ADR-0142 の後の値に揃えた（examples/vera/verification.md の「吸収体を入れた集合体のトラック間隔」）。
  branches の測り直しは後回しにした（M4-3 のライブラリ作成と待ち行列で競合する。燃焼の縦の一本を優先する指示）。
  手順と見積もり: `branches-xs-table`（xs 13 分と assembly-2B.toml の 27 状態。1 状態は 1B〜2K の今の値から 20〜50 秒の見込みで、
  ADR-0156 の前の 256 秒より短い。見積もりは測ってから決める）のあとに、direct 6 本と core 7 本を、verification.md の手順で流す。
  xs の排他で 1〜2 時間。dvc.yaml はなく、`bin/xs`、`bin/lattice`、`bin/core` を `examples/branches` で走らせる
