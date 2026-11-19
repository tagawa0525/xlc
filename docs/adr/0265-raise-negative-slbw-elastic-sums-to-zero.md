---
type: ADR
status: accepted
date: 2026-11-20
requires: [ADR-0167]
supersedes: []
superseded-by:
issues: ["#0305"]
---

# ADR-0265: 分離共鳴の SLBW の弾性散乱が、MF3 のバックグラウンドとの和で負になるところを 0 にする

## 背景

xs は分離共鳴の SLBW（MF2 の LRU = 1、LRF = 1）を読めず、ENDF/B-VII.0 と VII.1 の燃焼のライブラリを作れない（#0305）。ユーザーは、
本番の核データの版を決める比較（#0303、ADR-0264）のために、4 つの版の全部からライブラリを作れるようにしておくことを求めた
（チェーンが VII.0 と VII.1 で止まるという報告を受けた場面。ユーザー 2026-11-20）。

SLBW のパラメータは MLBW と同じ並びで、断面積は ENDF-102（ENDF-6 Formats Manual、BNL-224854-2023-INRE、2023-09-28）の付録 D.1.1 の
式 (D.1)〜(D.8) で求める。準位の間の干渉を持たないので、捕獲と核分裂は負にならないが、弾性散乱は、共鳴の下の負の干渉の裾が重なると
負になる。ENDF-102 の 2.4.9.1 は、これは自然には無く、評価者は負のエネルギーの準位か、MF3 のバックグラウンドで打ち消すべきとし、
多群定数を作るなら「弾性散乱が負になったら 0 にする」程度の手当てでよく、そのときは全断面積もそれに合わせて増やすとする。

NJOY2016 の RECONR は、範囲ごとの共鳴の部分（`csslbw` の結果）の各反応を、負なら 0 にしてから MF3 のバックグラウンドを足し（`sigma`）、
足したあとの弾性散乱が 1e-8 b 以下なら 1e-8 b にする（`emerge`）。全断面積は部分反応の和で作り直す（`emerge` の redundant reactions）。

判断の軸:

- 評価済みデータ（MF2 と MF3 の和）を、形式の定めた式どおりに使うこと
- 物理として無い負の断面積を、後の処理（ドップラー拡がり、群定数）に渡さないこと
- 参照コードの扱いが評価の意図と食い違うなら、真似ずに正しい方法を実装し、差を warn で知らせること（リポジトリの CLAUDE.md の「方針」）

## 検討した案

- 案 A: 共鳴の部分と MF3 のバックグラウンドの和が負になるところを 0 にする（採用）
- 案 B: NJOY と同じく、共鳴の部分を負なら 0 にしてからバックグラウンドを足す: 負の裾を MF3 で打ち消した評価（ENDF-102 の 2.4.9.1 が
  挙げる手当ての 1 つ）では、打ち消しの分だけ弾性散乱が大きくなる。評価の和と違う値になる根拠が無い
- 案 C: 負のまま残す: 物理として無い値で、全断面積も負の分だけ小さくなる。ENDF-102 は手当てを求める
- 案 D: SLBW のパラメータを MLBW の式で計算する（ENDF-102 の 2.4.9.1 は、点ごとの断面積に SLBW を使う理由は無いとする）: 評価者は
  SLBW の式と MF3 の和で断面積を合わせたので、式を替えると評価の値から外れる

## 決定と理由

案 A にする。

- MF2 の読み込み: LRF = 1 を MLBW と同じ並びで読む（同じ構造体 `BreitWigner`）。散乱半径のエネルギーの表（NRO ≠ 0）と競合反応の幅
  （LRX ≠ 0）は MLBW と同じくエラー（ADR-0167）。標的のスピンと l から作れない J（VII.0 の Am-241、Am-243 の全部の共鳴と、VII.0 と
  VII.1 の Pm-148m の 1 つの共鳴は J = I）は、MLBW ではエラーのままにし、SLBW では読んで、その数を warn で出す。SLBW は J を統計因子
  g = (2J + 1) / (2 (2I + 1)) にだけ使い、準位を J の系列に分けないので、J = I は g = 1/2 として式が決まる（ENDF-102 の 2.4.7 は
  J の分からない共鳴に J = I を置いた古い評価を挙げ、すべての共鳴が J = I なら非負になるとする。NJOY の `csslbw` も同じ g を使う）
- 断面積: ENDF-102 の (D.1)〜(D.7)。透過因子、シフト因子、位相シフト、チャネル半径は MLBW と同じ（NJOY の `csslbw` と同じ）。全幅は
  MLBW と同じく Γn(E) + Γγ + Γf（LRX = 0 では競合反応の幅が無く、GT との差はパラメータの丸め）
- 負の弾性散乱: 0 K の再構成で、範囲の中の弾性散乱を max(共鳴の部分 + MF3, 0) にする。和を 0 で切ると MF3 の表の点で折れるので、範囲の
  中の MF3 の弾性散乱の表の点を、再構成の格子に必ず入れる点に加える。MF3 が段差を持つ点では、和は段差の左右の両方の値を持つので、
  範囲の下端では右、上端では左の値、範囲の中では小さい方に対して切り、どちらの側も負にしない。全断面積は部分反応の和で作り直すので、0 にした分だけ増える
  （ENDF-102 の手当てのとおり）。共鳴の部分が負の点があれば、和を 0 にした点と、MF3 が打ち消した点（NJOY と値が違う点）の数を
  warn で出す
- NJOY の 1e-8 b の下限は真似ない: 0 にした弾性散乱を 1e-8 b にする理由が評価に無い

## 帰結

- VII.0 と VII.1 の SLBW の核種の 0 K の再構成と多群の処理ができるようになる
- MF3 が負の裾を打ち消す評価では、弾性散乱が NJOY より小さい。VII.0 と VII.1 の燃焼のライブラリの SLBW の核種で、弾性散乱が負になる
  Pu-238、Pu-242、Cm-242 は、範囲の中の MF3 の弾性散乱が 0 なので、NJOY と同じ値になる（examples/depletion/verification.md の
  「版ごとのライブラリの分離共鳴の SLBW」）
- 0 にした弾性散乱を NJOY は 1e-8 b にするので、その点の相対差は 1 になる。NJOY と比べるときは、0 に近い点の差を絶対値で見る

## 確認

- 単体テスト（crates/xs/src/resonance/slbw/tests.rs）: `single_s_wave_level_matches_exact_breit_wigner`、
  `several_levels_sum_the_endf_single_level_terms`、`cooperating_levels_give_the_negative_elastic_of_the_formula`
- 再構成の単体テスト（crates/xs/src/pointwise/tests.rs）: `slbw_range_adds_the_single_level_formula_to_the_background`、
  `negative_slbw_elastic_with_its_background_is_raised_to_zero`、`negative_slbw_elastic_is_raised_against_both_sides_of_a_background_step`
- MF2 の読み込み（crates/xs/src/endf/mf2/tests.rs）: `reads_slbw_range_in_the_form_of_mlbw`、`energy_dependent_radius_in_slbw_is_error`
- ポテンシャル散乱（crates/xs/src/process/tests.rs）: `resolved_range_has_its_top_and_potential_scattering`
- MF2 の読み込みの J（crates/xs/src/endf/mf2/tests.rs）: `slbw_reads_resonances_whose_spin_equals_the_target_spin`、
  `mlbw_rejects_resonances_whose_spin_equals_the_target_spin`
- NJOY との比較（xlc 9a5525c5、2026-11-20）: VII.1 の Xe-135、Pm-148m、VII.0 の Pu-238、Pu-242、Am-241 の 0 K の再構成を、NJOY2016 の
  RECONR（許容誤差 1e-3）の PENDF と `compare_njoy_pendf` で比べた。分離共鳴の範囲の全断面積、核分裂、捕獲の相対差の最大は 9.7e-4〜1.0e-3 で、
  許容誤差の 2 倍を超える点は無い。弾性散乱の 2 倍を超える点は、Pu-238 と Pu-242 の NJOY の値が 3.9e-4 b 以下の点だけで、差の絶対値は
  1.7e-6 b 以下（NJOY の丸めと 1e-8 b の下限。値と入力は examples/depletion/verification.md の「版ごとのライブラリの分離共鳴の SLBW」）
- VIII.0 の既存の核種の出力が変わらないこと: Zr-90、Xe-135（MLBW）、U-238（Reich-Moore）の `xs nuclide` の出力のバイトが、この変更の前の
  main（3e134ac6）と同じ
