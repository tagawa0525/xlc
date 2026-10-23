---
status: open
labels: [depletion, verification, openmc]
opened: 2026-10-23
closed:
milestone: M4
blocked-by: []
resolution:
---

# OpenMC の抽出のピン出力の規格化が吸収体のセルを数える

## 問題

OpenMC の燃焼の参照解の抽出（examples/openmc_run.py の `fuel_mean`、examples/depletion/openmc_depletion_assembly.py の
`bos_points`）は、ピン出力を「値が正のセルの平均」で割る。OpenMC の燃焼はチェーンの全核種を CRAM で解くので、アクチニドのない
吸収体（Pyrex、WABA）にも丸めの誤差の大きさのアクチニドができ（2F の pyrex@6-3 は 60 MWd/kg で Pu-239 が 1.2e5 個。負の値も
出る）、核分裂の tally が 7e-24〜4.6e-16 の値を持つ。その結果、吸収体のセルが平均の分母に入り、2F のピン出力の相対差が一様に
−5.7〜−8.3% ずれていた（1 − 264/288 = 8.33% が測った rms と一致）。

本番の比較（chore/p2-depletion-production）は、比べる側（examples/depletion/compare_depletion_assembly.py）で核分裂する組を
物質の配置（fuel と gd の領域）で決め、OpenMC のピン出力を燃料のセルの平均で割り直し、吸収体の最大の核分裂を `absorber_max`
として JSON に残して 1e-6 を超えれば止める形にした（2F の rms は 0.062〜0.112% になった）。抽出の側は直していない。statepoint の
多くが残っていないので、再抽出はできない。

残っている判断:

- 抽出（`fuel_mean`、`bos_points`）の規格化を、値の正負でなく物質で決める形に直すか。直すなら、次に参照解を作るとき（2P、2M、
  2A-density、2A-rods の取り出しの後）から適用し、比べる側の割り直しは残すか消すかを決める
- ADR-0181 の「核分裂のある組」の定義（値が正の組）を、物質で決める定義に改める ADR（amends）が要るか
- xlc の側の `cell_fission_rates`（crates/lattice/src/power.rs）も r > 0 の平均で規格化している。今の xlc の吸収体のセルは厳密に 0
  （2F で燃料 264 セルの平均が 1 ± 2e-15）だが、すべての問題で確かめていない
- 採らなかった案（大きさのしきい値で 0 とみなす）は、しきい値が黙った既定になるので採らなかった

関係する ADR: ADR-0181

## 完了の条件

抽出の規格化を物質で決める形にするかを決め、直すなら次の参照解から適用して、比べる側との二重の割り直しを無くす。
ADR-0181 の定義の扱いを決める。

## 経緯

- 2026-10-23: 本番の比較が 2F の Pyrex の位置の 1e-24 の核分裂で止まり、原因を調べて見つけた。比べる側で割り直して本番を進めた
