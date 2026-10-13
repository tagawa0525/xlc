---
status: accepted
date: 2026-10-14
requires: [ADR-0164]
supersedes: []
superseded-by:
issues: []
---

# ADR-0167: 分離共鳴の MLBW のエネルギー依存の散乱半径（NRO=1）を ENDF-6 の式どおり位相シフトに使う

## 背景

燃焼鎖（[ADR-0164](0164-build-full-chain-from-endf.md)）の核種の断面積を xs で作るには（docs/plans/m4-depletion.md の「進め方」3）、xs が
ENDF/B-VIII.0 の中性子のデータを読めなければならない。xs は、分離共鳴の範囲で散乱半径がエネルギーの表で与えられたもの
（MF2 の範囲の見出しの NRO ≠ 0）を読めず、エラーにする（crates/xs/src/endf/mf2.rs）。

中性子のサブライブラリの 557 の評価を `crates/chain/examples/survey_neutrons.rs` で調べると（表は docs/plans/m4-depletion.md の
「核データ」）、分離共鳴の範囲で NRO = 1 のものは 52 評価で、すべて MLBW（LRF = 2）、NAPS = 0 だった。Reich-Moore と
R-matrix limited に NRO ≠ 0 はない。52 評価は、短寿命の核分裂生成物（I-133、I-134、Te-131、Nd-149、Pm-150 など）と、
Ag-108〜117、Gd-159、Pu-245 などで、範囲はどれも 1 keV 以下にある。非分離共鳴の NRO = 1（Au-197、U-238）は、すでに読める。

ENDF-6 の形式（ENDF-102、2.2.1 節）では、NRO = 1 のとき、範囲の見出しの次に散乱半径 AP(E) の TAB1 があり、NAPS = 0 では、
透過因子とシフト因子に a = 0.123 AWRI^(1/3) + 0.08（10^-12 cm）を、位相シフトに AP(E) を使う。NAPS = 1 は両方に AP(E)、
NAPS = 2 は透過因子とシフト因子に見出しの AP、位相シフトに AP(E) を使う。NJOY2016 の RECONR（`csmlbw`）は、NAPS = 0 では
この式どおりに計算する。Reich-Moore と R-matrix limited の NRO ≠ 0 は、SAMMY の方法の読み込み（`rdsammy`）がエラーにする。

52 評価で、範囲の中の AP(E) の変化は、範囲の両端の差で 0.01〜0.7%（OpenMC 0.15.2 の `openmc.data` で読み、範囲を
400 点で数えた。2026-10-14。最大は Sm-146 の 0.68%、Cu-64 の 0.37%）。

判断の軸:

- 評価済みデータを、形式の定めた式どおりに使うこと
- 必要な形だけを実装し、評価にない組み合わせは推測で補わずエラーにする

## 検討した案

- 案 A: MLBW の NRO = 1、NAPS = 0 だけを、ENDF-6 の式どおりに実装する（透過因子とシフト因子は a、位相シフトは AP(E)）。
  ほかの組み合わせはエラーのまま（採用）
- 案 B: 表を読み飛ばし、見出しの一定の AP を使う: ENDF-6 では、NRO = 1 の見出しの AP は NAPS = 2 の透過因子にだけ使う値で、
  位相シフトに使う根拠がない。差は 0.7% 以下で小さいが、評価済みデータを使わない理由にならない
- 案 C: NAPS = 1、2 と、Reich-Moore の NRO = 1 も実装する: VIII.0 に使う評価がない。NJOY の `rdf2bw` は、NAPS = 1 では共鳴
  エネルギーでの透過因子を見出しの一定の AP で、断面積を求めるエネルギーでは AP(E) で求めるので、Γn(E) = Γn(Er) P(E)/P(Er) の
  2 つの P が別の半径になる（NAPS = 2 も同じく食い違う）。比べる参照がなく、どちらに合わせるかを決める材料もないので、評価が
  使うときに決める
- 案 D: 52 評価を鎖の輸送から外す: Nd-149（Pm-149 → Sm-149 の経路の先行核）や I-133、Te-131 などの核分裂生成物を外すことになる。
  外すかどうかは縮約（[ADR-0165](0165-reduce-chain-by-contribution-and-folding.md)）が吸収への寄与で決めるもので、読めないことを理由にしない

## 決定と理由

案 A にする。

- MF2 の読み込み: 分離共鳴の範囲の NRO = 1 は、LRF = 2（MLBW）で NAPS = 0 のときだけ読み、AP(E) の TAB1 を持つ。表が範囲
  [EL, EH] を覆わないか、不連続（x の重複）を持てばエラー（範囲の外へ外挿しない。非分離共鳴の LRF = 1 の核分裂幅の表と
  同じ。不連続では、その点の半径と範囲の上端での内側の値が決まらない。VIII.0 の 52 評価の表は不連続を持たず、範囲の中に点も持たない。OpenMC 0.15.2 の `openmc.data` で読んで数えた。2026-10-14）。NRO = 1 で
  ほかの LRF か NAPS、NRO ≥ 2 はエラー
- 断面積: 透過因子とシフト因子（共鳴エネルギーでの値と、求めるエネルギーでの値の両方）は a（AWRI から）で、位相シフトは
  求めるエネルギー E での AP(E)（表の補間則で補間する。表の外はエラー）で求める。NJOY RECONR の `csmlbw` の NAPS = 0 と同じ。
  表の点で位相シフトが折れるので、範囲の中の表の点を、0 K の再構成の格子に必ず入れる点（共鳴のピークと半値幅の位置と同じ扱い）
  に加える
- 分離共鳴の上端とポテンシャル散乱（ライブラリの `resolved.potential`。範囲より上の減速計算の散乱源と、等価理論の重みに
  使う一定の σp）: 4π a² の a に、範囲の上端 EH での AP(EH) を使う。範囲より上の散乱源に使う値なので、上端の値で
  つなぐ。範囲の中の変化は 0.7% 以下なので、等価理論の重みへの差は小さい

## 確認

- MF2 の読み込みの単体テスト（crates/xs/src/endf/mf2/tests.rs）: `reads_energy_dependent_radius_of_mlbw_range`、
  `rejects_mlbw_radius_table_short_of_the_range`、`rejects_mlbw_radius_table_with_a_step`、
  `rejects_energy_dependent_mlbw_radius_unless_naps_is_zero`、`energy_dependent_radius_in_reich_moore_is_error`
- 断面積の単体テスト（crates/xs/src/resonance/mlbw/tests.rs）: `energy_dependent_radius_sets_the_s_wave_phase_shift_at_each_energy`
  （共鳴から遠いエネルギーで、弾性散乱が s 波の硬球散乱 (4π/k²) sin²(k AP(E)) に 1e-9 で一致する）、
  `energy_dependent_radius_keeps_the_awri_channel_radius_near_a_p_wave_resonance`（l = 1 の共鳴の近くで、各エネルギーの AP(E) を
  一定の AP にした範囲と一致する）、`energy_outside_the_radius_table_is_error`
- 再構成の単体テスト（crates/xs/src/pointwise/tests.rs）: `mlbw_tables_have_points_at_the_breaks_of_the_radius_table`
- ポテンシャル散乱の単体テスト（crates/xs/src/process/tests.rs）: `potential_uses_the_radius_at_the_top_of_an_energy_dependent_range`
- NJOY との比較（xlc 46ba7a41、2026-10-14）: 0 K の再構成を、NJOY2016 RECONR（許容誤差 1e-3）の PENDF と `compare_njoy_pendf` で
  比べた。Nd-149、I-133 と、範囲の中の AP(E) の変化が最大の Sm-146 で、全断面積、弾性散乱、捕獲の NJOY の全格子点での
  相対差の最大は 9.1e-4〜1.0e-3 で、許容誤差の 2 倍を超える点はない（NJOY の入力と出力は references/njoy-runs/nro/<核種>/。git で追跡しない。RECONR の入力は
  references/njoy-runs/U238/input と同じ形で MAT だけが違う）
- 棚卸しの道具（`crates/chain/examples/survey_neutrons.rs`、6548fca7、2026-10-14）で、NRO = 1 の 52 評価がすべて読める
  （読めない理由に NRO が現れない）

## 帰結

- 52 評価の 0 K の再構成と多群の処理ができるようになる。燃焼鎖の核種のうち xs が読めないものは、残りの形式（計画の
  「核データ」の表）だけになり、その扱いは燃焼鎖の核種のライブラリの ADR で決める
- NRO = 1 の NAPS = 1、2 と、Reich-Moore、R-matrix limited の NRO ≠ 0 は、評価が使うようになったら（ENDF/B-VIII.1 など）、
  NJOY の食い違い（案 C）を含めて決め直す
