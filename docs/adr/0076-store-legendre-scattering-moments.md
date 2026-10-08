---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: ["#0009"]
---

# ADR-0076: ライブラリの散乱行列を、Legendre のモーメント P0〜P_L で持つ

## 背景

MOC を非等方散乱で解く（[ADR-0075](0075-solve-moc-with-anisotropic-scattering.md)）ので、ライブラリに散乱行列の高次のモーメントが要る。

## 検討した案

計画には、ほかの持ち方を比べた記録はない。

## 決定と理由

ライブラリの散乱行列を、Legendre のモーメント P0〜P_L で持つ（L は核データの処理の入力。既定 3）。

- 2体の反応（弾性散乱と離散準位の非弾性散乱）は、出射のエネルギーで実験室系の散乱角が運動学で決まるので、エネルギーの
  移行の積分に P_l(μ_lab) を掛ける。重心系の角度分布（弾性散乱と離散準位の MF4、MF6 の LAW=2）も同じ形で入る。まず
  重心系で等方として運動学の分（水素の効果の大半）を入れ、MF4 と LAW=2 の角度分布はその後に入れる（[ADR-0081](0081-use-mf4-angular-distributions.md)）
- 熱中性子散乱則と自由ガスは、μ が α（またはエネルギー）の式で表せるので、核の μ の積分に P_l(μ) を掛ける
- 連続的な放出は、角度分布を持つ形式（MF6 の Kalbach-Mann など、MF5 の反応の MF4。LCT で実験室系か重心系かを見る）は
  その角度分布を使う。角度分布が等方のときだけ P1 以上が 0（Zircaloy-4 の Sn の MT16, 17, 22, 28, 91 などは MF4 が
  実験室系で非等方）
- 自己遮蔽した核種では、P0 の行の和（非弾性散乱や多重度を含む）に弾性散乱の変化を足した値と元の和の比で、すべての次数の
  行を変える（角度分布と群の間の移行の形は無限希釈のまま。散乱行列の自己遮蔽を後回しにすることは [ADR-0050](0050-tabulate-sigma0-per-nuclide-and-temperature.md)）。高次の
  モーメントにも P0 と同じ比を使う点は MPACT と同じ（MPACT の比は遮蔽した散乱断面積と無限希釈の比。MPACT Theory Manual
  v4.3, 式 9.13〜9.14）

## 確認

`hydrogen_moments_follow_the_kinematics`、`hydrogen_moments_keep_the_mean_cosine`、`heavy_target_moments_follow_the_kinematics`、
`free_gas_moments_match_the_s_alpha_beta_table`、`kernel_moments_match_the_free_gas_double_differential`、
`shielding_scales_every_scattering_moment`、`isotropic_moments_reproduce_p0`。NJOY GROUPR との比較（examples/vera/README.md）。

## 結果

- MF6 の角度のデータの扱いは [ADR-0082](0082-leave-out-mf6-angular-data.md)
- U-238 の最も高い群の弾性散乱が NJOY と 2e-3 違う（#0009）
