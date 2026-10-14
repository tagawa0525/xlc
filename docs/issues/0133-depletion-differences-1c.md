---
status: open
labels: [lattice, depletion]
opened: 2026-10-15
closed:
milestone: M4
blocked-by: []
resolution:
---

# 格子計算の燃焼の 1C の OpenMC との差の切り分け

## 問題

VERA 1C のピンセルの燃焼（ADR-0173。examples/depletion/verification.md の「格子計算の燃焼」）で、k∞ の差の変化は 1 リングで
−58〜+70 pcm に収まるが、数密度に原因を確かめていない差が残る。

1. 最初のステップ（0.1 MWd/kg）の終わりで U-239、Np-239、Pu-239 が +4.7%、Xe-135 が +2.5%。2 ステップ目からは U-239 の差は
   0.2% 以下で、最初のステップの修正子（中点の組成。始まりに Xe-135 がなく、Xe-135 の燃え尽きなしに積もった組成）の U-238 の
   捕獲率の差が候補
2. I-135 −0.9%、Cs-133 −0.5%、Rh-103 −1.5%、Nd-149 と Pm-149 −0.5% が燃焼度によらずほぼ一定。核分裂の数は 0.011% 以内で
   合うので、収率の重みの違いが候補: xlc は ADR-0164 の区間（500 keV と 14 MeV の幾何平均 2.65 MeV で分ける）で U-238 などに
   14 MeV の収率を使い、OpenMC の `cutoff`（`fission_yield_mode="cutoff"`）は 112.5 eV より上をすべて 500 keV の収率にする
3. Sm-149 が +0.9%（先行核は −0.5% なので、xlc の Sm-149 の捕獲率が約 1.4% 小さい計算になる）。0.0973 eV の共鳴の熱群での
   縮約が候補
4. Pm-148m が +27%（5 MWd/kg）。Pm-148m の捕獲率か収率の扱いが候補
5. 4 リングの 2 MWd/kg の点の k の差の変化が +122 pcm で、前後から外れる。OpenMC の 4 リングと 1 リングの k の差が −35〜+62 pcm
   にばらつくので、OpenMC の燃焼で反応率の統計誤差が組成に入る揺らぎが候補

関係する ADR: ADR-0173

## 方針の案

- 1: OpenMC の `depletion_results.h5` の反応率（ステップごと）と、xlc の領域ごとの反応率を同じ組成で比べる（xlc の反応率を
  履歴に書き出す）
- 2: OpenMC の収率を xlc と同じ区間の重みにする（`AveragedFissionYieldHelper` か、自前の helper）か、xlc の区間を OpenMC の
  `cutoff` に合わせた試行で、I-135 と Cs-133 の差が消えるかを見る
- 3、4: 均質な無限媒質（examples/vera の `homogeneous_medium.py`）で Sm-149、Pm-148m の微視的断面積を比べる
- 5: OpenMC の燃焼を種を変えて数本流し、点ごとの k の揺らぎを測る

## 完了の条件

5 つの差の原因を、測った数値で説明するか直す。

## 経緯

- 2026-10-15: ADR-0173 の実装の検証（ブランチ feat/lattice-depletion-rates）で見つけた
