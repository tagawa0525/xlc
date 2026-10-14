---
status: open
labels: [core]
opened: 2026-10-14
closed:
milestone: M3
blocked-by: []
resolution:
---

# 真空の面に接する薄いノードで、反復の途中に中性子束が負になって止まる

## 問題

Problem 5 の 3D の最小の入力（examples/vera/p5/hzp/small.toml、small-half.toml。3.10% の 3 × 3 と反射体の輪、径方向の
外周は真空）を core で解くと、NEM の 4 回目の外側の反復で止まる（8ff11f3d の bin/core、md5 a8d61ca7、2026-10-14）:

```text
error: state half: node 150 next to a vacuum face has a flux -0.000009789832978246111; it must be positive
```

- node 150 は、5 × 5 の位置の左下の角の反射体（角型の格子の反射体のノードの定数）で、z の層 70.761〜70.817 cm（仕様書の高さ
  15.761〜15.817 cm。Pyrex の棒の下の端栓の端と、KENO の軸方向の層の境界の間の、厚さ 0.056 cm の層）。径方向の幅は 21.5 cm
- 止まるのは、真空の面の横方向の漏れ（面の漏れ = ノードの漏れ × 面の中性子束 / ノード平均。ADR-0099）で、ノード平均が正で
  ないとき（crates/core/src/nem.rs の `transverse_shape`）。反復は 1〜3 回目で k が 0.666、0.671、0.735 と動いている途中で、
  収束した解でも負になるかは確かめていない
- 同じ入力で 70.761 cm の境界だけをメッシュから除くと（この 3.10% の集合体には、その高さの層の境界はない）、収束して
  k = 0.996097（state half）
- Problem 3B（同じメッシュ、同じ 0.056 cm の層）は径方向が鏡面反射で、真空の面に接するノードがないので止まらない

再現: worktree の直下で `cd examples/vera/p5/hzp && ../../../../bin/core small-half.toml`（格子と平板の段階の出力が要る。
`python3 trials/p5-hzp/run_stages.py`、または dvc の段階）。

Problem 5 は、Pyrex の棒のある集合体の層に合わせて全炉心のメッシュに 15.761 cm を入れるので、全炉心（core-hzp）も同じ形の
ノード（角と辺の反射体の列）を持つ。

関係する ADR: ADR-0099

## 完了の条件

真空の面に接する薄い（径方向に比べて z が非常に短い）ノードの横方向の漏れの扱いを決め（反復の途中で正でないときの扱い、
または厚さの比に応じた形）、ADR に書いて core に入れる。small-half.toml が 70.761 cm の境界を含むメッシュで収束し、境界を
除いた解との差を記録に書く。#0092 の Problem 5 の 3D の炉心（core-hzp、small）を、元の層（examples/vera/p5/p5_inputs.py の
`CORE_MERGED_BOUNDS` を空にし、15.761 cm の境界を戻したもの）で流し直し、examples/vera/verification.md の Problem 5 の 3D の
節を書き直す。戻すと、毒物の下端の支持格子の層（`*-poison-igrid`）の 8 つの格子（単一集合体 6、反射体 2）の段階が
生成のツールに戻るので、それも流す。

## 経緯

- 2026-10-14: #0092 の最小の入力の確認で見つけた。#0092 は、炉心のメッシュと層から 15.761 cm の境界を除いて進めた
  （ADR-0154。Problem 3B で k +0.03 pcm、下端の層の軸方向の出力 +0.17%）
