---
status: accepted
date: 2026-11-07
requires: [ADR-0033, ADR-0039, ADR-0164, ADR-0172, ADR-0173, ADR-0183, ADR-0186, ADR-0187, ADR-0192, ADR-0200, ADR-0216]
supersedes: []
superseded-by:
issues: ["#0153", "#0264"]
---

# ADR-0238: 燃焼度の軸を持つテーブルに、ライブラリの族ごとの遅発中性子の量と群速度の逆数から作る動特性のデータと、輸送しない核種の放射化断面積を含む全核種の 1 群の微視的断面積を書く

## 背景

[ADR-0033](0033-store-kinetics-and-cooling-data-in-m4-tables.md) は、動特性と冷却計算が使うデータ（即発と遅発の ν、族ごとの遅発中性子の割合と崩壊定数、
遅発中性子のスペクトル、群速度の逆数、冷却計算の燃焼をやり直すための全核種の 1 群の微視的断面積）を、本番の分岐のテーブルより前に入れると
決めた（後から足すとテーブルをすべて作り直す）。[ADR-0187](0187-tabulate-burnup-branches-from-history.md) は置き場を決めた: 遅発中性子のデータと群速度の逆数は
状態ごと（少数群定数と同じく状態の格子で内挿する）、全核種の 1 群の微視的断面積は燃焼度の点ごとに公称の状態（履歴の輸送）の値だけ。
[ADR-0216](0216-tabulate-micro-depletion-cross-sections-per-state.md) はミクロ燃焼の側を入れ、ADR-0033 のデータを、ライブラリに族の量を足す決定と輸送しない核種の
断面積の出所の決定に分けて残した（#0153）。

前提:

- ライブラリ（crates/format/src/library.rs の `Library`、`NuclideData`、`TemperatureData`）は、遅発中性子の族ごとの量、崩壊定数、群速度の
  逆数を持たない。χ は即発と遅発の中性子を ν_p と ν_d で重みを付けた 1 本（[ADR-0186](0186-include-delayed-neutrons-in-chi.md)）で、ν は MF1 MT452 の全数
- xs は ENDF を自前で読み（NJOY を使わない）、χ のために MF1 MT455 の ν_d(E)（崩壊定数の LIST は読み飛ばす）と MF5 MT455 を読んでいる
  （crates/xs/src/process.rs の `fission_parts`、crates/xs/src/endf/mf1.rs の `parse_delayed_nubar`）
- ENDF/B-VIII.0 の中性子の評価（references/endf/ENDF-B-VIII.0_neutrons の 561 ファイル）のうち MF1 MT455 を持つのは 84 で、
  32 は 6 族の崩壊定数と MF5 MT455（6 つの小区分。族の割合 p_k(E) と族ごとの分布）を持ち、24 は 6 族の崩壊定数を持つが MF5 MT455 を
  持たず（Np-239 など）、28 は 1 族の形（Cm-244 は λ = 1 /s の 1 族）で MF5 MT455 を持たない（2026-11-07 に MF1 MT455 の 2 行目の NPL と
  MF5 MT455 の見出しの NK を数えた）。族の割合は MF5 MT455 にしかない。崩壊定数は核種ごとに違う（U-235 の第 1 族 0.013336 /s、Pu-239 の
  0.013271 /s）
- 燃焼チェーン（examples/depletion/chain.msgpack。[ADR-0164](0164-build-full-chain-from-endf.md)）の 3,820 核種のうち、中性子の反応か核分裂を持つのは ENDF/B-VIII.0 の
  中性子のデータのある 555 核種で、ほかの 3,265 核種は崩壊だけ。2A の多群ライブラリ（examples/depletion/p2/xs.toml）は 90 核種で、どれも
  チェーンの反応を持つ。反応を持つ 555 核種のうち 465 はライブラリにない（2026-11-07 にチェーンとライブラリの入力を数えた）
- 格子計算の燃焼（[ADR-0173](0173-deplete-lattice-pins-by-ce-cm.md)）は、ライブラリにない核種を崩壊だけで追う

判断の軸:

- 評価済みデータの量（ν_d、族の割合、崩壊定数、スペクトル）をそのまま使い、テーブルの値が格子で均質化した反応率と同じ重みで作られること
- core がテーブルを内挿した値が、格子点で格子計算の値そのものになること（ADR-0187 の内挿の性質）
- 重い計算（ライブラリと履歴の作り直し）を、項目の足し方で増やさないこと

## 検討した案

遅発中性子の族のデータの出所:

- 案 D1: xs が MF1 MT455 の崩壊定数 λ_k と ν_d(E)、MF5 MT455 の族の割合 p_k(E) と族ごとの分布を読み、核種ごとに λ_k と、温度ごとに族ごとの
  群の遅発中性子数 ν_d,k,g = ∫ν_d p_k σf w dE / ∫σf w dE と族ごとの遅発スペクトル χ_d,k,g（族 k の中性子の ν_d p_k σf w で重みを付けた出射の
  分布）を持つ（採用）。χ を作る読み込みと群への振り分け（ADR-0186）をそのまま使う
- 案 D2: NJOY の GROUPR（MF3 MT455 の族ごとの遅発の ν σf、MF5 MT455 の族ごとのスペクトル）でライブラリを作る: xs は ENDF を自前で処理し、
  NJOY は参照解（OpenMC のライブラリ）にだけ使う。遅発の量だけを NJOY から取ると、同じ核種の群定数の処理が 2 つの経路に分かれる
- 案 D3: 核種によらない族の定数（U-235 の 6 族など）を持つ: Pu の核分裂の割合が燃焼で増えて遅発中性子の割合が変わる（ADR-0033 が燃焼点ごとに
  持つとした理由）のを表せない
- 案 D4: 群ごとの割合 β_k,g = ν_d,k,g / ν_g を持ち、自己遮蔽した νσf に掛ける: 案 D1 と同じく群の中の形を無限希釈の重みで代える近似で、
  遅発の生成が ν_p によって変わる形になる。ν_d と σf はどちらも評価の量で、自己遮蔽した σf はライブラリの `Shielded` にある

MF5 MT455 のない評価（2A のライブラリの Np-239、Cm-244）:

- 案 M1: 族の量を持たず、遅発中性子は ADR-0186 の χ と同じく即発の中性子に数える（採用）。族の割合は MF5 MT455 にしかなく、ADR-0186 は
  χ でこの評価の遅発中性子を MF5 MT18 で数えて warn を出す
- 案 M2: 別の核種の族の割合を借りる: 評価にない量を作る。Cm-244 の MF1 MT455 は 1 族の形で、6 族の割合に分けられない

族の崩壊定数の均質化:

- 案 L1: テーブルに族ごとの遅発中性子の生成 ν_dΣf,k,G と、それを崩壊定数で割った ν_dΣf,k,G / λ（核種ごとの λ で割ってから足す）を少数群の
  断面積として持ち、core がノードの中性子束で λ_k = Σ_G ν_dΣf,k,G φ_G / Σ_G (ν_dΣf/λ)_k,G φ_G を作る（採用）。平衡の先行核は
  C_i,k = ν_d,i,k Σf,i φ / λ_i,k で、この λ_k は平衡の先行核の和と崩壊率の和（= 生成）の両方を保つ。どちらの量も組成と中性子束に線形なので、
  ほかの断面積と同じく均質化し、内挿できる
- 案 L2: 格子の中性子束で均質化した λ_k を状態ごとに持ち、内挿する: 比を内挿するので、生成と先行核の和の内挿と一致しない。ノードの中性子束で
  重みを付け直せない
- 案 L3: 核分裂する核種ごとの族（核種の数 × 6）を持つ: core が核種ごとの先行核を解くことになり、炉心の動特性の解法（6 族）に要らない量になる

群速度の逆数:

- 案 V1: xs がライブラリに群ごとの (1/v)_g = ∫_g w/v dE / ∫_g w dE（w は [ADR-0039](0039-weight-with-njoy-iwt4.md) の重み関数。NJOY GROUPR の MT259 と同じ定義）を持ち、
  lattice が MOC の中性子束で少数群にする（採用）
- 案 V2: lattice が群の境界から 1/E の重みで作る: ライブラリの群定数と違う重みになる（熱の群の重みは Maxwell 分布）
- 案 V3: 超微細群の中性子束で重みを付ける: 超微細群は分離共鳴の範囲だけで、1/v はその範囲で滑らか

冷却計算の全核種の 1 群の微視的断面積の出所:

- 案 A1: xs の入力で放射化の核種（`activation = true`）とした核種は、チェーンの反応（[ADR-0172](0172-store-chain-reactions-in-library.md) の MT と経路の並び）と核分裂の無限
  希釈の群の断面積だけを温度ごとに作り、ライブラリの `activation` に置く（輸送しない）。テーブルの 1 群の微視的断面積は、輸送する核種は
  ADR-0216 の微視的断面積の定義（自己遮蔽した断面積）を 1 群にしたもの、放射化の核種は無限希釈の断面積を同じ定義で 1 群にしたもの（採用）。
  微量の核種は自己遮蔽しない。全核種は、チェーンが反応を持つ ENDF/B-VIII.0 の中性子のデータのある核種（前提の 555）
- 案 A2: 輸送する核種だけ: 格子計算の燃焼（ADR-0173）と同じ扱いになるが、2A ではチェーンの反応を持つ 555 核種のうち 465 の反応が抜け、
  ADR-0033 の「全核種」にならない
- 案 A3: 輸送しない核種もライブラリの輸送する核種として処理する: 自己遮蔽の表、散乱行列、超微細群の断面積を作り（2A のライブラリの 90 核種で
  xs の CPU 時間 4,596 秒。examples/depletion/p2/stage-cpu-seconds.json）、格子計算の輸送する核種も変える。輸送する核種の選び方は ADR-0173 の決定で、
  冷却計算に要るのは反応率だけ
- 案 A4: 中性子のデータのない核種に、別の放射化ライブラリ（TENDL、JEFF の放射化ファイルなど）の断面積を足す: チェーンはその核種に反応を
  持たない（ADR-0164 の作り方）ので、チェーンの作り方から変え、評価済みデータが混ざる

放射化の断面積の置き場:

- 案 P1: 多群ライブラリの別の列 `activation`（採用）: 格子計算が読むライブラリは 1 つのまま。核種ごとの xs のステージ（ADR-0138、ADR-0193）を
  そのまま使える
- 案 P2: 別のファイル: lattice の入力に 2 つ目のライブラリを足し、群構造と温度が同じことを別に確かめる

## 決定と理由

案 D1、M1、L1、V1、A1、P1 にする。

### ライブラリ

- `Library` に群ごとの群速度の逆数 `inverse_velocity`（s/cm）を持つ。(1/v)_g = ∫_g w(E)/v(E) dE / ∫_g w(E) dE、v = √(2E/m_n)（ADR-0039 の
  重み関数。NJOY GROUPR の MT259 と同じ定義）。核種によらない
- 遅発中性子の族のデータを持つのは、MF1 MT455 と MF5 MT455 の両方を持つ評価だけ。`NuclideData` に族ごとの崩壊定数 λ_k（MF1 MT455 の LIST。
  1/s）、`TemperatureData` に族ごとの群の遅発中性子数 ν_d,k,g = ∫ν_d p_k σf w dE / ∫σf w dE と、族ごとの遅発スペクトル χ_d,k,g（和は 1）を
  持つ。p_k(E) は MF5 MT455 の族 k の小区分の割合。族の数は MF1 MT455 の LIST と MF5 MT455 の小区分の数で同じであること（違えばエラー）。
  ほかの評価は持たず（nil）、その遅発中性子は ADR-0186 の χ と同じく即発に数える
- ν_d,k,g は無限希釈の重み（σf w）で作り、格子計算は自己遮蔽した σf に掛ける（遅発の生成 = ν_d,k,g σf,g）。ν_d(E) と p_k(E) は数 MeV までほぼ
  一定で、群の中の共鳴の形は σf にあるので、群の中の重みの違いは σf の自己遮蔽が受け持つ
- 放射化の核種（`activation`）は、名前、ZA、AWR と、温度ごとのチェーンの反応の群の断面積（`TemperatureData` の `reactions` と同じ形）と核分裂の
  断面積 σf（無限希釈、重み w）だけを持つ。輸送する核種と同じ名前を持てない。格子計算は輸送しない（ADR-0173 の輸送する核種は `nuclides` だけ）

### テーブル

ADR-0187 の置き場のとおり、`[history]` の実行が書く燃焼度の軸を持つテーブルに 2 つの項目を足す。燃焼度の軸を持たないテーブルは持たない。

- `kinetics`: 状態ごと（`points` の順）に、族ごとの遅発中性子の生成 `delayed_nu_fission`（ν_dΣf,k,G。1/cm）、族ごとの `delayed_precursors`
  （Σ_i ν_d,i,k Σf,i / λ_i,k。s/cm）、族ごとの遅発スペクトル `delayed_chi`（族 k の生成で重みを付けた χ_d,k,G。和は 1）、群速度の逆数
  `inverse_velocity`（s/cm）。少数群定数（`fewgroup::homogenize`）と同じく、MOC の中性子束と体積で重みを付けて SPH 因子を掛けた量を足し、
  少数群の中性子束で割る（格子全体の均質化）。遅発の量は燃焼の領域の輸送する核種の数密度と自己遮蔽した σf から作る。燃焼の領域の核種の
  νσf の和が少数群定数の νΣf と一致すること（燃焼させない物質が核分裂しないこと）を確かめる。即発の ν と即発の χ は、少数群定数の νΣf と χ
  からこの遅発の量を引いて求まる
- `cooling`: 燃焼度の点ごとに、公称の状態（ADR-0200 の制御棒のない組の x₀。点が 1 つの軸はその点なので、入力の物質がどの軸でも
  格子点にあるときに履歴の輸送と同じ問題）の、核種ごとの 1 群の微視的断面積（barn）。
  核種は履歴の `transportable` と、ライブラリの放射化の核種のうち履歴の核種にあるもの（チェーンの順）。核種ごとに、チェーンの反応の経路ごとの
  断面積 `reactions`（MT と経路の並び）、核分裂の断面積 `fission`、核分裂の収率の境界より下の割合 `thermal_fission`（[ADR-0183](0183-weight-fission-yields-by-thermal-fast-cutoff.md)。
  核分裂しない核種は nil）。定義は ADR-0216 の微視的断面積の定義を 1 群にしたもの（燃焼の領域の数密度で重みを付け、均質化した数密度が 0 の
  核種は燃焼の領域で 1 として重みを付ける）。放射化の核種の群の断面積は燃焼の領域の温度の値（ない温度はエラー。ADR-0187 と同じく近い温度で
  代用しない）

族の数は、族のデータを持つすべての核種で同じであること（違えばエラー）。核分裂があるのに族のデータを持つ核種がなければエラーにする。

### core

- 動特性のデータは、少数群定数と同じ重み（ADR-0200 の状態の組み立てと、[ADR-0192](0192-interpolate-burnup-by-cubic-hermite.md) の燃焼度の内挿）で内挿する。族の崩壊定数は、
  内挿した `delayed_nu_fission` と `delayed_precursors` から、与えた少数群の中性子束で L1 の式で作る
- 冷却計算の 1 群の微視的断面積は、燃焼度だけで ADR-0192 の内挿をする（`thermal_fission` は ADR-0216 と同じく、境界より下の σf を足して比に戻す）
- 解法（動特性、冷却計算）は、それぞれを作るときに入れる（ADR-0033）

## 確認

- ライブラリ（crates/format/src/library/tests.rs）: `delayed_data_need_every_temperature_and_the_families_of_the_decay_constants`、
  `inverse_velocities_are_positive_for_every_group`、`activation_cross_sections_are_looked_up_at_their_exact_temperature`、
  `activation_nuclides_are_checked_like_the_chain_reactions_of_the_library`、`round_trips_bit_exactly`、`merge_joins_the_nuclides_in_order`
- xs（crates/xs）: `reads_the_decay_constants_and_the_delayed_nubar`（MF1 MT455）、`delayed_families_split_the_delayed_nubar_by_their_fractions`、
  `delayed_nubar_varying_in_a_group_is_weighted_by_the_fission_rate`、`inverse_velocity_averages_one_over_v_with_the_weight_function`（1/E の
  範囲の解析値と相対 1e-9）、`delayed_families_take_the_decay_constants_and_the_family_spectra`、
  `delayed_families_need_one_spectrum_per_decay_constant`、`evaluations_without_the_delayed_spectrum_have_no_families`、
  `activation_nuclides_have_the_chain_reactions_of_a_processed_nuclide`、`reads_activation_nuclides`、`an_activation_result_must_match_its_input`、
  `the_pointwise_range_leaves_out_activation_nuclides`、`libraries_carry_the_inverse_velocities_of_their_groups`
- テーブル（crates/format/src/table/tests.rs）: `kinetics_and_cooling_data_round_trip_bit_exactly_and_are_written_only_when_held`、
  `kinetics_data_need_a_burnup_axis_one_sub_node_and_matching_sizes`、`cooling_data_have_one_point_per_burnup_and_the_same_reactions_at_every_point`
- lattice（crates/lattice/src/history/tests.rs）: `kinetics_split_the_delayed_production_by_family_and_divide_it_by_the_decay_constants`
  （遅発の生成が ν_d,k × Σf、崩壊定数で割った量、族のスペクトル、群速度の逆数が相対 1e-12）、`a_fissile_lattice_needs_delayed_neutron_families`、
  `cooling_cross_sections_are_the_one_group_micro_cross_sections_of_transported_and_activation_nuclides`、
  `activation_cross_sections_need_the_temperature_of_the_depletion_regions`、`history_branches_need_the_nominal_state_for_the_cooling_cross_sections`、
  入力（crates/lattice/src/input/tests.rs）の `history_tables_are_not_split_without_micro_nuclides_either`。1 つの入力の制御棒の分岐（crates/lattice/tests/history_rod_branch.rs の
  `history_rod_branches_match_separate_lattices`）で、状態ごとの動特性のデータが値ごとの別の格子計算とビットで一致し、冷却計算の断面積が
  棒のない値の入力と同じ
- core（crates/core/src/state/tests.rs）: `kinetics_data_at_the_grid_points_are_the_lattice_values`（格子点でテーブルの値そのもの）、
  `kinetics_data_take_the_weights_of_the_few_group_constants`、`decay_constants_keep_the_delayed_production_and_the_precursors_under_the_given_flux`、
  `cooling_cross_sections_are_interpolated_in_burnup_only`（燃焼度の点ではその点の値）、`kinetics_and_cooling_data_need_the_data_and_the_burnup_range`
- 2A（examples/depletion/verification.md の「動特性と冷却計算のデータ（2A）」、636d3b1d3、2026-11-07、Slurm 726）: 放射化の 455 核種を足した
  ライブラリと作り直した履歴から作ったテーブルで、燃焼度の点ごとの公称の状態の遅発中性子の割合は燃焼度 0 の 712.62 pcm から 60 MWd/kg の
  464.33 pcm に減る。冷却計算の断面積を持つ核種は 544（輸送する 89 と放射化の 455）

## 帰結

- ライブラリの形式に必須の項目（`inverse_velocity`、`activation`、核種の `decay_constants`、温度の `delayed`）を足したので、今あるライブラリは
  どれも読めない。depletion/p2 のほかの xs のステージを持つ例（examples の m0、m1、m3-nuclides、mox、reflector、core、vera、vera/endfb70、
  depletion）のライブラリは、それを読む作業が次に来たときに作り直す。核種ごとの CPU 時間の記録（`stage-cpu-seconds.json`）の和は、vera 6,031 秒、
  vera/endfb70 2,286 秒、depletion 11,616 秒。depletion/p2 のライブラリは作り直した（「確認」の 2A）
- p2 の dvc.lock の格子計算の燃焼の 7 ステージ（lattice@2A-density、2A-rods、2F、2M、2N、2O、2P）は、deps の p2 のライブラリ
  （examples/depletion/p2/depletion.msgpack。p2 の 8 本の燃焼が共有し、作り直した）が前の版のままで、changed になる。deps の bin/lattice も前の版
  （48387d65）だが、main の dvc.lock の bin/lattice が 2927836cb（2026-11-04）で別の値に移ったときから古く、ADR-0238 の帰結ではない。前の実行のログの `total_seconds`（Slurm 666〜672、
  2026-11-03〜04。2F と 2N の他は割り込まれうるジョブ）は 2A-density 2,270 秒、2A-rods 3,510 秒、2F 3,345 秒、2M 11,139 秒、2N 6,299 秒、
  2O 4,393 秒、2P 5,004 秒で、和は 35,960 秒（約 10 時間）。それらを読む作業が次に来たときに再実行する
- テーブルの `kinetics` と `cooling` は、燃焼の履歴から作ったテーブルだけのキーで、持たないテーブル（今ある燃焼度の軸のテーブルを含む）は
  そのまま読める。core の動特性と冷却計算の解法はまだなく（ADR-0033）、内挿の関数だけがある
- `[history]` の入力は、`micro` がなくても部分ノードに分けないテーブルと、公称の状態（ADR-0200 の x₀）を求め、そうでなければ
  輸送の前に止まる。燃焼の履歴からの分岐は、族のデータを持つ核種のないライブラリ（今あるライブラリ）ではエラーになる。チェックポイントの状態にも
  2 つの項目を足した
- 2A のライブラリの放射化の核種は、チェーンの反応を持つ輸送しない 465 核種のうち、xs が処理できない 10 核種（荷電粒子の対を持つ分離共鳴、
  生成物の状態の経路の和が反応を超える評価）を除いた 455 核種（examples/depletion/p2/xs.toml。#0264）
- 放射化の核種の断面積は無限希釈で、自己遮蔽しない。格子計算の燃焼（ADR-0173）は放射化の断面積を使わず、輸送しない核種を
  崩壊だけで追うままにした。チェーンと照らした核種の確かめ、中性子のデータのない長寿命の核種（案 A4）、格子計算の燃焼で放射化の
  反応を数えるかは #0264
- 群速度の逆数の群の中の重みは重み関数で、熱の群の中の実際の中性子束の形との差は測っていない
