---
status: open
labels: [lattice, shielding, performance]
opened: 2026-11-03
closed:
milestone:
blocked-by: []
resolution:
---

# 減速計算の作業の配列の数え方の残り（#0212 のレビュー）

## 問題

issue #0212 の直し（ADR-0213 の「## 変更」、2026-11-03）は、1 回の減速計算の量から、`solve_cell` が手で挙げた格子点に比例する配列の確保を
除いて持つ区間の上限にする。fix/0212-lattice-depletion-memory の 1 周目の /code-review（2026-11-03、HEAD ac7fa356e）が挙げた次の点は、
落ちるテストか再現する手順が無いので直していない（間引かない格子の 2A の 1 ステップは MaxRSS 69.2 GiB で流れた。Slurm 652）。

- 数えていない配列: `Thermal` の今のスイープの累積の積率（`prefix`。リングごとに格子点（start から）あたり 48 B）、`finish` で前の
  スイープに移した積率と `flux.to_vec()` の中性子束の写し。自由ガスの範囲の広い入力では、数えた列と同じ桁になりうる（測っていない）
- 倍に伸びる確保: `Prefix.values` は空から伸び、`finish` のたびに伸ばし直す。自由ガスの σs の列は `collect::<Result<Vec<_>, _>>()` で
  伸びる。`lethargy_grid` の返す配列も足した点で伸びうる
- 数えすぎ: 自由ガスの枝の `energies` は `Thermal::new` の後に捨てるが数える。核種の格子の和、間引いた格子、`lethargy_grid` の
  `u` と `energies` は格子点を作った後は使わないが、関数の終わりまで残して数える（先に捨てれば、メモリも差し引きも減る）
- 形: 配列を足すたびに `solve_cell` の `work +=` に足す必要があり、漏れやすい。配列を持つ型（`Sweep`、`Thermal`、`Prefix`）が
  自分の確保を返す形にすると、数える場所と確保する場所が揃う。数え漏れを止める試験（小さな入力で held と work の和、または
  計数するアロケータ）が無い
- ログ: `UFG slowing down` の `work` は、`Scattering::Static` と `sweeps = 1`（差し引かない経路）でも出るので、コメントの「持つ区間の
  上限から除いた」が成り立たない
- 数えない段: `release` の後の移行の計数（`segment_groups`、`collapse`）と飛程の表は、どの量にも入らない。間引かない格子で測っていない

関係する ADR: ADR-0213

## 完了の条件

- 上の配列を数えるか捨てるかを決めて直し、間引かない格子の 2A の 1 ステップ（trials/0212/2A-fixed2.toml、コミット ac7fa356e）で、
  ログの held と work の和が量を超えないことと MaxRSS を確かめる

## 経緯

- 2026-11-03: #0212 のブランチの 1 周目の /code-review の指摘を、直す指摘の基準（再現の手順）に当たらないものとしてまとめた
