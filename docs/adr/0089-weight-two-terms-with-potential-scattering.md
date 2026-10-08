---
status: accepted
date: 2026-10-08
supersedes: [ADR-0023]
superseded-by:
issues: []
---

# ADR-0089: 等価原理の2項有理近似の重みに、吸収体自身のポテンシャル散乱 σp を入れる

## 背景

等価原理の2項有理近似（[ADR-0074](0074-shield-outside-range-by-equivalence.md)）で2つの σ0 の表を合わせる重みは、M2 では NR 近似の中性子束 σ0,i/(σt + σ0,i)
（[ADR-0023](0023-shield-by-bondarenko-and-equivalence.md)）とし、σp を落としていて、2つの項の重みがずれていた。

## 検討した案

- σp を落とした重み（M2 の形）: 2項を合わせた捕獲が最大 1.3% ずれる
- ライブラリに σ0 ごとの群の中性子束を持たせる: 厳密になるが、形式を変える必要がある
- NR 近似の群の中性子束 β_i (σp + σ0,i)/(σt(σ0,i) + σ0,i)（採用）

## 決定と理由

重みは、NR 近似の群の中性子束 β_i (σp + σ0,i)/(σt(σ0,i) + σ0,i) とする（σp は吸収体自身のポテンシャル散乱で、中性子束の
線源に入る。σp は、分離共鳴の範囲の s 波の散乱半径から求めた値（[ADR-0073](0073-store-pointwise-data-in-library.md)）をすべての群に使う）。ライブラリの分離共鳴の範囲の群は flux calculator の中性子束で重みづけしてあり（[ADR-0072](0072-weight-shielded-data-by-flux-calculator.md)）、この式は近似に
なる。xs の flux calculator で求めた群の中性子束と比べると、σp を入れた式は U-238（σ0 = 30〜100 b）で −12〜+5% 違う
ものの、VERA 1B の燃料の U-238 で2項を合わせた捕獲は 0.4% 以内で一致した（σp を落とすと最大 1.3% ずれる。
docs/plans/m3-core-hzp.md、99deb84、2026-10-08）。どちらも燃料では UFG の範囲の群で、実際には UFG の値を使う。この差では、
ライブラリの形式を変える理由にならない。

## 確認

`two_terms_are_weighted_by_the_narrow_resonance_flux`。k は 1B〜1D、2B（枠なし）で 1 pcm 以内しか変わらない
（examples/vera/verification.md）。

## 結果

燃料でこの重みが効くのは UFG の範囲の外の群だけである。
