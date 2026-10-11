---
status: open
labels: [verification]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# Problem 3 と 5-2D の OpenMC のライブラリに Ti と 565 K の構造材の核種を足す

## 問題

Problem 3 と 5-2D の参照解の入力（examples/vera/p3/openmc-3A.toml など）は、`examples/vera/openmc-library/cross_sections.xml`
を指す。いまの `examples/vera/openmc-library.toml` では、この入力は走らない。

- Inconel-718（3A、3B の端の支持格子。仕様書の表 P3-3）の Ti-46〜50 がない
- 565 K で解く問題（3B、5A〜5C-2D）は温度の内挿のモードで解くが、SS304、Inconel、ノズルと炉心板、炭素鋼、Pyrex、AIC、B4C
  の核種のうち C-12、C-13、P-31、Si-28〜30、Mn-55、Ni-58〜64、Ag-107、Ag-109、In-113、In-115、Cd-106〜116 は 600 K しか
  持たない。OpenMC 0.15.2 は、温度を 1 つしか持たない核種があると全体を最も近い温度のモードに戻し、H in H2O の S(α,β) に
  565 K から 1 K 以内の表がないので止まる（examples/vera/verification.md の「参照解」の 565 K の節）

短い実行（1 世代あたりの時間の測定）は、これらを足したライブラリ（`examples/vera/openmc-library-timing.toml`。
references/openmc/p3-p5-timing に置いた）で実行した。足した 29 核種は NJOY で 54 秒（Slurm で並列、2026-10-12）。

## 完了の条件

`examples/vera/openmc-library.toml` に、Ti-46〜50 を 565 K と 600 K で、上の核種に 565 K を足し、段階（examples/vera/dvc.yaml の
openmc-neutron、openmc-index）で作り、`dvc.lock` をコミットする。p3/ と p5/ の入力が走ることを確かめる。

## 経緯

- 2026-10-12: Problem 3 と 5-2D の OpenMC のモデル（feat/p3-p5-openmc）で見つけた。ライブラリの段階は別の作業で整えて
  いるので、このブランチでは足さない
