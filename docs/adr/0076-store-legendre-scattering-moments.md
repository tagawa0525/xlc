---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: ["#0009", "#0014"]
requires: []
---

# ADR-0076: ライブラリの散乱行列を、Legendre のモーメント P0〜P_L で持つ

## 背景

MOC を非等方散乱で解く（[ADR-0075](0075-solve-moc-with-anisotropic-scattering.md)）ので、ライブラリに散乱行列の高次のモーメントが要る。自己遮蔽した核種で散乱行列の行を
どう合わせるかを、すべての次数について決める必要があった。

## 検討した案

- 自己遮蔽した核種の散乱行列は、無限希釈の形のまま、行の和を自己遮蔽した弾性散乱断面積に合わせる。計画で一時採った
  （2026-10-06、7ebd447。M1 の計画）: 行の和には弾性散乱のほかに非弾性散乱や多重度も入るので、格子計算で遮蔽を
  使い始めたとき（d072037、2026-10-07）に下の形で実装した
- MF6 の角度のデータも使う（M2 の計画）: [ADR-0082](0082-leave-out-mf6-angular-data.md) の「検討した案」

ほかの持ち方を比べた記録はない。

## 決定と理由

ライブラリの散乱行列を、Legendre のモーメント P0〜P_L で持つ（L は核データの処理の入力。既定 3）。

- 2体の反応（弾性散乱と離散準位の非弾性散乱）は、出射のエネルギーで実験室系の散乱角が運動学で決まるので、エネルギーの
  移行の積分に P_l(μ_lab) を掛ける。重心系の角度分布（弾性散乱と離散準位の MF4）も同じ形で入る。まず重心系で等方
  として運動学の分（水素の効果の大半）を入れ、MF4 の角度分布はその後に入れる（[ADR-0081](0081-use-mf4-angular-distributions.md)）
- 熱中性子散乱則と自由ガスは、μ が α（またはエネルギー）の式で表せるので、核の μ の積分に P_l(μ) を掛ける
- 連続的な放出は、角度分布を持つ形式（MF5 の反応の MF4。LCT で実験室系か重心系かを見る）はその角度分布を使う。角度
  分布が等方のときだけ P1 以上が 0（Zircaloy-4 の Sn の MT16, 17, 22, 28, 91 などは MF4 が実験室系で非等方）
- 自己遮蔽した核種では、P0 の行の和（非弾性散乱や多重度を含む）に弾性散乱の変化を足した値と元の和の比で、すべての
  次数の行を変える。散乱行列は σ0 ごとに持たず、角度分布と群の間の移行の形は無限希釈のまま
  とする（散乱行列の自己遮蔽は効果が小さいので後回しにする）。高次のモーメントにも P0 と同じ比を
  使う点は MPACT と同じ（MPACT の比は遮蔽した散乱断面積と無限希釈の比。MPACT Theory Manual v4.3, 式 9.13〜9.14）

M2 の計画はこの項目に、MF6 の角度のデータ（LAW=2 の重心系の分布、連続的な放出の Kalbach-Mann など）を使うことも
書いていたが、実装の前に M3 で使わないことにした（[ADR-0082](0082-leave-out-mf6-angular-data.md)）。

## 確認

`hydrogen_moments_follow_the_kinematics`、`hydrogen_moments_keep_the_mean_cosine`、`heavy_target_moments_follow_the_kinematics`、
`free_gas_moments_match_the_s_alpha_beta_table`、`kernel_moments_match_the_free_gas_double_differential`、
`shielding_scales_every_scattering_moment`、`isotropic_moments_reproduce_p0`。NJOY GROUPR との比較（examples/vera/verification.md）。

## 帰結

- L は核データの処理の入力にしておらず、xs のコードに 3 で固定してある（`crates/xs/src/main.rs` の `SCATTERING_ORDER`。
  xs.toml に項目がない）。格子計算は、ライブラリにある次数までを入力で選べる。#0014
- MF4 の角度分布は [ADR-0081](0081-use-mf4-angular-distributions.md)、MF6 の角度のデータは [ADR-0082](0082-leave-out-mf6-angular-data.md)
- U-238 の最も高い群の弾性散乱が NJOY と 2e-3 違う（#0009）
