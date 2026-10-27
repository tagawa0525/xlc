---
status: proposed
date: 2026-10-27
requires: [ADR-0173]
supersedes: []
superseded-by:
issues: ["#0001", "#0152"]
---

# ADR-0206: MOX-1 のピンセルを 1C の燃焼の条件で燃焼させる

## 背景

格子計算の燃焼の参照解（[ADR-0173](0173-deplete-lattice-pins-by-ce-cm.md)）は VERA 1C のピンセル（UO2）だけで、燃焼で増える k の差の
切り分け（#0133、#0152）に Pu の多い燃料の題材がない。#0001 も、4 eV より下の共鳴の扱いを MOX のピンセルで検証するとしている。
そこで、examples/mox の LWR MOX ピンセル MOX-1 を、1C と同じ形の燃焼（燃料 1 リング、分岐なしの 1 つの状態）で OpenMC と xlc の
両方で燃焼させる。そのときの条件を決める。

前提:

- 1C の燃焼の条件: 燃料 900 K、被覆管と減速材 600 K、初めの重金属あたり 38.41 W/g、燃焼度 0.1、0.5、1、2、3、5 MWd/kg、OpenMC は
  1 世代 20 万粒子・130 世代（最初の 30 世代を捨てる）、乱数のシード 1（examples/depletion/openmc-1C-rings1.toml）
- examples/mox の MOX-1 は、手元のデータに合わせてすべての領域を 600 K にしている（examples/mox/verification.md の「仕様からの変更」。
  仕様は 300 K と、燃料だけ 560 K）
- 燃焼の OpenMC のライブラリ（examples/depletion/openmc-library.toml）は、燃料にある核種（U、Pu、Am、チェーンの核種）を 900 K だけで、
  被覆管と減速材の核種（Zr、H-1、O-16 など）を 600 K で持つ（2026-10-27 に `openmc-library/cross_sections.xml` の各核種の `kTs` で、
  Pu-238〜242、Am-241、U-235、U-238 は 900 K だけと確かめた）。xlc のライブラリ（examples/depletion/xs.toml）は燃料の核種を 565、600、900、
  1200 K で持つ
- 燃料では、燃料の温度そのもののデータのある核種だけを輸送する（ADR-0173）。燃料を 600 K にすると、OpenMC のライブラリの燃料の核種が
  すべて輸送から落ちる
- NEA の仕様（references/mox/pin_cell.pdf の付録 2）は、燃焼の条件（出力密度、燃焼度）を定めていない

判断の軸:

- 1C の参照と同じ条件にして、1C と MOX-1 の差を燃料の組成と形状の違いだけにする
- 重いステージ（1C の OpenMC の燃焼、1 回の輸送 約 410 秒 × 13 回）を古くしない

## 検討した案

- 案 T1: 燃料 900 K、被覆管と減速材 600 K（1C と同じ）。出力密度、燃焼度の点、粒子数、世代数、シードも 1C と同じにする（採用）。
  燃焼のライブラリ（xlc と OpenMC）をそのまま使える。燃焼の初めの k は、examples/mox の 600 K の MOX-1（k の差 −219 pcm）とは
  温度が違うので、そのまま比べられない
- 案 T2: examples/mox と同じくすべて 600 K。OpenMC のライブラリに燃料の核種（チェーンの核種を含む約 50）の 600 K が要る。
  openmc-library.toml に足すと、そのディレクトリを deps にする 1C の OpenMC の燃焼（`openmc-1C`）が古くなり再実行が要る。別の
  ライブラリにすると、同じ核種の組の入力とステージが 2 つになり、NJOY を約 50 核種ぶん流す。1C との比較では燃料の温度という違いが増える
- 出力密度を MOX の仕様に合わせる案は、仕様に燃焼の条件がないので立たない

## 決定と理由

案 T1。1C と同じ条件にすれば、1C と MOX-1 の燃焼の k の差の推移の違いは燃料の組成（Pu の多さ）と形状（ギャップのない被覆管、
ホウ素のない減速材）の違いだけになり、#0152 と #0001 の切り分けの題材として比べやすい。ライブラリを作り直さず、1C の参照解も
古くならない。燃焼の初めの k の差は、温度の違う examples/mox の値でなく、この燃焼の 0 の点の値を使う。

形状と組成は examples/mox の MOX-1（examples/mox/openmc_model.py と pincell-mox1.toml）のまま変えない。

### 帰結

- ステージは examples/depletion に置く。燃焼に要るもの（チェーン、燃焼の多群ライブラリ、OpenMC のライブラリとチェーン）が
  そこのステージの出力で、examples/mox の `mox.msgpack` はチェーンの核種を持たない
- OpenMC の燃焼は、別のモジュール examples/depletion/openmc_depletion_mox.py で組み立てる。`openmc_depletion.py` のコードを直すと
  `openmc-1C` が古くなり、再実行が要る（#0132 と同じ理由。`dvc commit` で付け直すのは、コメントだけの変更のとき（acfe57b0）に限る）。
  温度は入力の項目（`fuel-temperature`、`moderator-temperature`）にし、温度のモードはライブラリの表の温度だけを使う `nearest` にする。
  初めの燃料の核種のどれかを燃料の温度で輸送できなければ、輸送の前に止める
- 1C と MOX-1 の燃焼の k の差の推移を比べるときは、両方の xlc を同じライブラリと `bin/lattice` で流す（2026-10-27 の時点の 1C の xlc の
  結果は ADR-0205 の前のライブラリのもの）
- 参照解（OpenMC の燃焼）は、入力とステージを用意しただけで、まだ流していない

### 確認

- examples/depletion/test_openmc_depletion_mox.py: 入力の検査、等しい面積のリング、MOX-1 の組成と温度、ギャップのない形状、
  温度のモード、初めの燃料の核種が輸送されないときに止めること
- 粒子数を落とした OpenMC の試しの実行で、燃焼の 1 ステップが回ること（出力は残さない）
- 参照解を流したあと、`compare-mox1` の要約を examples/depletion/verification.md に書く
