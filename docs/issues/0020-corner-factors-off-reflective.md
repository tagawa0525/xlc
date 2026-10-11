---
status: open
labels: [lattice, core]
opened: 2026-10-08
closed:
milestone: M3
blocked-by: []
resolution:
---

# 均質な解が平らでない格子の角の不連続因子（CDF）

## 問題

CDF（角の不均質な中性子束 / 角の均質な中性子束）は、外周がすべて鏡面反射の単一集合体では、均質な解が平らなので
ノード平均の中性子束で割ればよい。反射体の格子（燃料と反射体を並べたもの）やカラーセットでは均質な解が平らでなく、
角の均質な中性子束を、ノードの中の 2D の均質な解から求める必要がある。面の ADF に使う NEM の1ノードの解は方向ごとの
1D の形しか持たず、角の値を与えない。ピン出力の再構成は、燃料と反射体の角でも角の中性子束の連続を使う。

関係する ADR: [ADR-0094](../adr/0094-define-adf-by-one-node-nem.md)、[ADR-0125](../adr/0125-derive-corner-factors-from-lattice-corner-balance.md)

## 完了の条件

非反射の格子で角の均質な中性子束を求める方法（例: ピン出力の再構成と同じ解析的な 2D の展開を、参照の面の正味の中性子流と
角の条件で解く。反射体の CDF を燃料の CDF で割って持つかを含む）を、文献（PARCS、SIMULATE-3、MPACT の角の扱い）と
比べて決め、ADR にする。REF1-2D と 5-2D の外周のピン出力を、同じ多群定数の xlc の MOC と比べて確かめる。

## 経緯

- 2026-10-08: ステップ 5 の ADR（ADR-0094〜0097）のレビューで、非反射の格子の CDF の定義がないことが分かった
- 2026-10-11: PARCS v3.0 Theory Manual（5.1 節）、MPACT Theory Manual 4.3（13.1.3 節）、RENUS と Scarabee の角の扱いを調べ、
  格子の中で core と同じ角の点のバランスを解いて CDF を作る ADR-0125 を下書きした（proposed）
- 2026-10-11: ADR-0125 を実装した。I 型の反射体の格子で、core の再構成が MOC のピン出力を RMS 0.080% で再現する
  （examples/reflector/verification.md、23e3a4a8）。REF1-2D と 5-2D の比較は残る（#0043、#0029、#0065）
