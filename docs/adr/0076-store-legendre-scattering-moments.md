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
  移行の積分に P_l(μ_lab) を掛ける。重心系の角度分布も同じ形で入る。重心系で等方としても、運動学の分（水素の効果の
  大半）が入る。どの角度分布を使うかは [ADR-0081](0081-use-mf4-angular-distributions.md) と [ADR-0082](0082-leave-out-mf6-angular-data.md)
- 熱中性子散乱則と自由ガスは、μ が α（またはエネルギー）の式で表せるので、核の μ の積分に P_l(μ) を掛ける
- 連続的な放出は、出射エネルギーと角度が独立な実験室系の角度分布があれば、エネルギーの移行に P_l の係数を掛ける
  （MF5 の反応の MF4。[ADR-0081](0081-use-mf4-angular-distributions.md)）。角度分布を使わない放出（MF6 の LAW=1。[ADR-0082](0082-leave-out-mf6-angular-data.md)）と、等方の放出の P1 以上は 0
- 自己遮蔽した核種では、P0 の行の和（非弾性散乱や多重度を含む）に弾性散乱の変化を足した値と元の和の比で、すべての次数の
  行を変える（角度分布と群の間の移行の形は無限希釈のまま。散乱行列の自己遮蔽を後回しにすることは [ADR-0050](0050-tabulate-sigma0-per-nuclide-and-temperature.md)）。高次の
  モーメントにも P0 と同じ比を使う点は MPACT と同じ（MPACT の比は遮蔽した散乱断面積と無限希釈の比。MPACT Theory Manual
  v4.3, 式 9.13〜9.14）

## 確認

`hydrogen_moments_follow_the_kinematics`、`hydrogen_moments_keep_the_mean_cosine`、`heavy_target_moments_follow_the_kinematics`、
`free_gas_moments_match_the_s_alpha_beta_table`、`kernel_moments_match_the_free_gas_double_differential`、
`shielding_scales_every_scattering_moment`、`isotropic_moments_reproduce_p0`。NJOY GROUPR との比較（examples/vera/README.md）。

## 結果

- 角度分布の出所は [ADR-0081](0081-use-mf4-angular-distributions.md)（MF4）と [ADR-0082](0082-leave-out-mf6-angular-data.md)（MF6）
- U-238 の最も高い群の弾性散乱が NJOY と 2e-3 違う（#0009）
