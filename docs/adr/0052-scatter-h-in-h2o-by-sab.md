---
status: accepted
date: 2026-10-06
supersedes: [ADR-0038]
superseded-by:
issues: []
---

# ADR-0052: H in H2O の熱中性子散乱は S(α,β) で扱い、酸素は自由ガスのままとする

## 背景

M0 では熱中性子の散乱を自由ガスで扱っていた（[ADR-0038](0038-use-zero-kelvin-cross-sections-in-m0.md)）。水素は水の分子に束縛されている。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

H in H2O は ENDF の S(α,β)（MF7 の非干渉性非弾性散乱）から散乱行列を作り、4 eV 未満の水素の自由ガス模型を置き換える。
酸素は ENDF/B-VIII.0 に O in H2O がないので、自由ガスのままとする。

## 確認

`reads_incoherent_inelastic_scattering`、`sab_matrix_satisfies_detailed_balance`、`sab_ignores_mf3_elastic_below_cutoff`。
NJOY THERMR との比較（`crates/xs/examples/compare_njoy_thermr.rs`、examples/m1/README.md）。

## 結果

- S(α,β) の温度は 283.6〜800 K しかない（[ADR-0047](0047-choose-temperatures-per-nuclide.md)）。565 K は [ADR-0080](0080-interpolate-sab-in-temperature.md)
- 4 eV 未満では静止した自由原子の散乱核が成り立たないので、減速計算の下限を 4 eV にする（[ADR-0073](0073-store-pointwise-data-in-library.md)、#0001）
