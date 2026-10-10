---
status: proposed
date: 2026-10-11
supersedes: []
superseded-by:
issues: ["#0058"]
---

# ADR-0120: 案内管に入れる棒（Pyrex の棒、制御棒）の中の He も真空にする

## 背景

[ADR-0054](0054-model-vera-materials.md) は、ペレットと被覆管のギャップの He を真空にすると決めた。VERA 2E〜2H では、
案内管に入れる棒の中にも He がある。Pyrex の棒では内管の中と Pyrex の内外の隙間、制御棒では吸収体と被覆管の間である。
仕様書（CASL-U-2012-0131-004、Problem 2 の MIPLIB の入力の注）は、参照解（KENO-VI）ではペレットと被覆管のギャップを
公称の密度の He（表 P2-4 の Gap、He-4 2.68714e-5 /(barn cm)）とし、「真空か空気としてもよい」、「制御棒と吸収体の棒の
ほかのギャップも同じに扱う」としている。

- 参照解（OpenMC）と同じ組成で比べられること
- ライブラリに足す核種（核データの処理と保守）

## 検討した案

- 案 A: He-4 を仕様の個数密度で入れる。KENO-VI と同じ形になるが、ライブラリと OpenMC のデータに He-4 を足す
- 案 B（採用）: ペレットのギャップと同じく真空にする

## 決定と理由

案 B にする。He-4 の熱中性子の散乱断面積は約 0.76 b で、仕様の個数密度では巨視的断面積が約 2e-5 /cm になる。
最も長い弦（Pyrex の内管の中、直径 0.43 cm）でも衝突の確率は約 1e-5、制御棒の隙間（厚さ 0.004 cm と 0.013 cm）では
約 1e-7 で、k への効果はどの比較の統計誤差よりも桁違いに小さいと見込む（計算では確かめていない）。仕様書も真空を
許している。xlc と OpenMC の参照解は同じに真空とする（[ADR-0031](0031-verify-against-openmc-with-same-data.md)）。

## 確認

xlc の入力（examples/vera/p2/assembly-2E.toml〜assembly-2H.toml）と OpenMC のモデル（examples/vera/openmc_model.py）が、
棒の中の He の領域を真空にする。確かめる数値はない約束事。

## 結果

- KENO-VI とは He の扱いが違う。上の見積もりでは、OpenMC と KENO-VI の差（2G で +129 pcm。examples/vera/verification.md、
  c0662d26、2026-10-11）を説明する大きさにはならない
- He-3 は天然の存在比（約 1.3e-6）が小さく、吸収（熱中性子で約 5300 b）を入れても巨視的断面積は約 2e-7 /cm
