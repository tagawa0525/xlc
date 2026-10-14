---
status: accepted
date: 2026-10-14
requires: [ADR-0011, ADR-0073, ADR-0164]
supersedes: []
superseded-by:
issues: []
---

# ADR-0172: ライブラリに燃焼鎖の反応の群の断面積を持たせる（捕獲は σ0 ごとに遮蔽しポイントワイズも持ち、しきい値の反応と状態ごとの経路は無限希釈）

## 背景

格子計算の燃焼（docs/plans/m4-depletion.md の「進め方」5）は、燃焼鎖（[ADR-0164](0164-build-full-chain-from-endf.md)）の反応ごとの反応率を要る: (n,γ)、(n,2n)、
(n,3n)、(n,4n)、(n,p)、(n,α) と核分裂、MF8 が状態を挙げる反応では状態ごとの経路（割合は MF9 の多重度か MF10 の状態ごとの生成
断面積。鎖は割合を持たず、xs が状態ごとの群の断面積を作ると決めた）。ライブラリ（crates/format/src/library.rs）は、σ0 ごとの
全、弾性散乱、核分裂、ν 核分裂と、無限希釈の散乱行列、UFG のポイントワイズの全、弾性散乱、核分裂、ν 核分裂しか持たない。
xs は σ0 ごとの捕獲（MT102）を群定数を作る途中で求めているが、書き出していない。

判断の軸:

- 燃焼の反応率の大部分を占める捕獲は、核分裂と同じく自己遮蔽する（Pu-240 の 1.06 eV、Gd、Sm-149 など）
- 核種ごとの反応と経路の集合を、鎖と同じ規則で決める（鎖と食い違うと、経路の反応率が欠けるか余る）
- 形式を変えるので、後方互換は作らない（利用者はまだいない）

## 検討した案

捕獲:

- 案 C1: 無限希釈だけを持つ: 共鳴の自己遮蔽が抜け、Pu-240 や U-238 の捕獲を大きく誤る
- 案 C2: σ0 ごとの群の捕獲（`shielded` に足す）と、UFG の範囲のポイントワイズの捕獲を持つ（採用）。核分裂と同じ扱いで、lattice は
  UFG と等価原理で捕獲を遮蔽できる

しきい値の反応（(n,2n)、(n,3n)、(n,4n)、(n,p)、(n,α)）と状態ごとの経路:

- 案 T1: σ0 ごとに持つ: しきい値の反応は高いエネルギーにあり、自己遮蔽は小さい。記憶が σ0 の数倍に増える
- 案 T2: 無限希釈の群の断面積を、鎖の反応ごと、経路ごとに持つ（採用）。捕獲の状態ごとの経路は、無限希釈の経路の断面積の比で、
  遮蔽した捕獲を分ける（MF9 の多重度はエネルギーに緩やかなので、群の中の遮蔽で比が大きく変わらない見込み。確かめていない）

反応の集合の決め方:

- 案 R1: xs が MF3 にあるすべての MT を書き出し、lattice が鎖の反応にまとめる: まとめる規則（(n,2n) の 16 と 875〜891 など）が
  鎖と lattice の 2 か所になる
- 案 R2: 鎖の反応の MT の範囲（ADR-0164）を xlc-format（`CHAIN_REACTIONS`）に、MF8 の状態の読み込みを xlc-xs に置き、鎖
  （xlc-chain）と xs の両方がそれを使う（採用）。規則が 1 か所になる。xs は xlc-depletion（faer）に依存しない（ADR-0166）ので、
  範囲は鎖の型の置き場でなく、両方が依存する xlc-format に置く

## 決定と理由

案 C2、T2、R2 にする。

- 形式: σ0 ごとの群定数（`shielded`）に `capture`（MT102）を、UFG のポイントワイズの断面積に `capture` を足す。温度ごとに `reactions`
  を足す: 鎖の反応（代表の MT。102、16、17、37、103、107）ごとに、経路（`final_state`。MF8 が状態を挙げなければ nil、挙げれば
  MF9 か MF10 と LFS。基底状態を挙げない反応の残りは LFS が nil）ごとの無限希釈の群の断面積 (barn)。経路の並びと `final_state` は
  鎖の `Channel` と同じ
- 反応の断面積: 再構成した断面積の代表の MT（部分の MT（875〜891、600〜649、800〜849）があれば、再構成が部分の和で作り直した
  もの。NJOY RECONR と同じ）。経路の断面積は、MF9 なら多重度 × 反応の断面積、MF10 なら状態の生成断面積、残りは反応の断面積から
  挙げた経路を引いたもの。残りが丸めの程度（相対 1e-9）より負ならエラー（挙げた状態が反応を超えている）。群の平均は、重み
  関数だけで取る（無限希釈）。MF10 は 0 K の値をそのまま使う（NJOY BROADR も MF10 を広げない）
- 反応を持つかは、鎖と同じく、代表か部分の MT のどれかが MF3 にあること

## 確認

- 形式の単体テスト（crates/format/src/library/tests.rs）: `keeps_the_chain_reactions_and_the_capture`、
  `reaction_channels_need_one_value_per_group`、`reactions_must_not_repeat_and_need_channels`、
  `interpolation_mixes_the_capture_and_the_reactions`、`interpolates_shielded_cross_sections_in_ln_sigma0`（捕獲）
- xs の単体テスト（crates/xs/src/process/tests.rs）: `capture_is_kept_per_sigma0_and_pointwise`、
  `chain_reactions_sum_the_partial_mts_and_split_the_mf9_states`、`capture_states_without_the_ground_state_leave_the_remainder_to_it`、
  `listed_states_beyond_the_reaction_are_an_error`
- lattice の単体テスト: `capture_is_shielded_with_the_same_flux_as_fission`（UFG）、`absorber_rods_are_shielded_by_slowing_down`
  （棒のリングの平均の捕獲）
- 鎖の出力（examples/depletion/chain.msgpack）の md5 が変わらない（120d7435…）
- examples/vera のライブラリを作り直すと、既存の項目は main のライブラリ（キャッシュ）と最大で相対 1.9e-15 違い（ADR-0170 の J, K の
  表の原点の移動による丸め）、6,026 の配列のうち 4,785 が同じ。VERA 1B の k は 1.182383 のまま（xlc 7c80cfc8、
  分岐点 c9c6e901、2026-10-14。ADR-0144 の判定）

## 帰結

- 既存のライブラリの形式が変わるので、xs の段階をすべて作り直す。bin/xs と bin/lattice が変わるので、下流の段階も走る
- 捕獲の状態ごとの経路の比を群の中で一定とする近似が残る（Am-241 の (n,γ) の Am-242m など）。燃焼の比較（「進め方」5）で
  OpenMC の状態ごとの反応率と比べる
- 核分裂の反応率は、既存の `fission` を使う
