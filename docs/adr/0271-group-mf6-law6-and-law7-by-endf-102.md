---
type: ADR
status: accepted
date: 2026-11-20
requires: [ADR-0057, ADR-0082]
supersedes: []
superseded-by:
issues: ["#0311", "#0321", "#0322", "#0323", "#0324"]
---

# ADR-0271: MF6 の中性子の LAW=6 は入射エネルギーごとに ENDF-102 の式で、LAW=7 は余弦で積分したエネルギー分布で群に振り分け、NJOY の表の作り方を真似ない

## 背景

xs は MF6 の中性子の LAW=6（N 体の位相空間分布）と LAW=7（実験室系の角度とエネルギーの表）を読めず、H-2（4 つの版）、
ENDF/B-VIII.1 の Li-6、Be-9（4 つの版）で止まっていた（#0311）。4 つの版のライブラリで全部の核種が流れることが求められている
（ユーザー 2026-11-20。別のセッションが伝えたもの）。読み方と、参照の NJOY2016 GROUPR（groupr.f90 の `getmf6`、`f6psp`、
`ll2lab`）の扱いのどちらに合わせるかを決める。

- 正確さ: ENDF-102 の定義のとおりに群に振り分けること
- NJOY GROUPR との比較で、差の原因を示せること
- 角度のデータの扱いは、LAW=1 と同じく ADR-0082 に揃える

## 検討した案

- NJOY GROUPR のとおりに表を作る: LAW=6 は、多重度の最初の点 / 1.1999 から 1.2 倍ずつの入射エネルギーで、E' = 0 から始まる
  実験室系の表にし、その間を単位基底で補間する。LAW=7 は、余弦で積分した点の並びを `ll2lab` の判定で間引く。`ll2lab` の判定
  （前の点が要るか）は、前の点での直線の値を f1·y(前の前) + (1 − f1)·y(今)（f1 = (x(前) − x(前の前))/(x(今) − x(前の前))）
  としていて、正しい (1 − f1)·y(前の前) + f1·y(今) と重みが入れ替わっている。ヒストグラムの E' の表では段の直前の点が消え、段が
  傾きに替わる。Be-9 の ENDF/B-VIII.0 の 2 MeV の入射で、群ごとの割合が ENDF のヒストグラムから差の絶対値の和で 31% ずれ、和が
  14% 大きくなる（重みを正しくすると 0.3%。`ll2lab` を真似た examples/depletion/issue-0311/ll2lab_thinning.py で測った。examples/depletion/verification.md）。
  採らない
- LAW=7 の余弦が開いた範囲を覆うかを、運動学から前方の円錐の μ_min を求めて余弦の範囲と比べて確かめる: 開いた範囲は出射の
  エネルギーごとに違い、式 6.2 の規格化で足りる。採らない（#0323）
- LAW=7 の余弦の範囲を確かめずに通す: 範囲の一部を欠いた分布が中性子を黙って失う。採らない（#0323）
- ENDF-102 の定義のとおりに振り分ける（採用）

## 決定と理由

ENDF-102 の定義のとおりに群に振り分け、NJOY の表の作り方は真似ない。読み込みで `tracing::warn!` を出す。式の出典は ENDF-6 Formats
Manual（ENDF-102）の 6.2.8（LAW=6。式 6.21〜6.24、6.27、6.28）と 6.2.9（LAW=7）、f の定義は 6.1 の式 6.1、規格化は式 6.2。

- LAW=6 は、重心系で等方な μ と E' の密度 C_n √E' (E_max − E')^(3n/2 − 4)（式 6.21。E' だけで積分すると 1/2 で、μ の −1〜1 と合わせて 1。C_3 = 4/(π E_max²)、C_4 = 105/(32 E_max^(7/2))、
  C_5 = 256/(14π E_max⁵) は式 6.22〜6.24）を、入射エネルギーの積分点ごとに実験室系の群に振り分ける。E_max = (M − m_i)/M · E_a（式 6.27）、
  E_a = m_T/(m_p + m_T) E + Q（式 6.28）で、M は APSX（この分布で扱う n 個の粒子の質量の和、中性子の質量を単位）、中性子では m_i = m_p = 1、
  m_T = A なので E_max = (APSX − 1)/APSX (A/(A + 1) E + Q)。Q は MF3 の QI（式 6.28 の定義）。分布は重心系で与えられる（式 6.21）ので、LCT に
  よらず重心系とする。n は ENDF-102 が係数を与える 3〜5 に限る。式 6.28 は 1 段の反応の E_a で、2 段の反応では 1 段目の反跳のエネルギー
  になる（6.2.8）。APSX が生成粒子の質量の和（およそ A + 1）より中性子の質量の半分を超えて小さい節は 2 段の反応の 2 段目として、エラーにする
  （1 段の H-2 と Li-6 の APSX は A + 1 との差が Q/(m_n c²) ほど、0.002〜0.004）。重心系の E' = E_max sin²θ について積分すると密度の両端の √ の特異性が消えるので、実験室系の範囲が群の境界を
  横切るところで θ を分けて Gauss-Legendre で積分する。入射エネルギーの積分は、しきい値と、実験室系の範囲の端が群の境界を横切る
  入射エネルギーで区切る。入射エネルギーの格子を置かないので、NJOY の 1.2 倍の格子とその補間の誤差が無い
- LAW=7 は、同じ E' での値を余弦の補間則（lin-lin かヒストグラム）で補間して、表の余弦の範囲で積分したエネルギー分布 g(E') にする
  （余弦が lin-lin なら NJOY の `ll2lab` の P0 と同じ積分。`ll2lab` は余弦の補間則によらず lin-lin で積分する。4 つの版の Be-9 は lin-lin）。
  f(μ, E, E') は式 6.2 の規格化（μ と E' で積分して 1）で、余弦と出射エネルギーはどちらも実験室系、
  余弦は TAB2 のふつうの補間則で与えられ、両方の変数が出射の粒子に開いた範囲の全部を覆う（6.2.9）。E' の表は lin-lin かヒストグラムに限る。そのとき g は全部の表の点の間で線形なので、
  ヒストグラムの段を同じ E' の 2 点として残した lin-lin の表で厳密に表せる。LAW=7 は実験室系の分布なので、重心系の LCT と組んだら
  エラーにする。開いた範囲は重心が中性子より速いと前方の円錐 μ_min〜1 になるので、余弦の範囲は −1〜1 を求めず、表の範囲で積分する。
  余弦が 1 つだけ（積分の幅が 0）と、−1〜1 の外の余弦はエラーにする。範囲の全部を覆うかは、運動学から開いた範囲を求めずに、式 6.2 の規格化で
  確かめる: 入射エネルギーごとの ∫∫ f dμ dE'（xs が積分に使う補間則で）が 1 から 1e-5 を超えてずれたら、入射エネルギーと積分の値を示して
  止める（黙って規格化し直さない）。しきい値（最初の入射エネルギー）の 0 の分布は、断面積が 0 なので通す。正しい前方の円錐は 1 になり、範囲の
  一部を欠いた分布は欠いた割合だけ足りない。1e-5 は見込みからの値で、ENDF の 11 桁の欄の実数（指数が 1 桁なら 7 桁、2 桁なら 6 桁）の値と
  格子の丸めの相対誤差 5e-7〜5e-6 から和のずれを 1e-6 の桁と見込み、余裕をとった（格子の近い点の幅の丸めは見込みに入れていない）。4 つの版の
  Be-9 の MT16 のずれは最大 7.0e-7（VII.0、VII.1 の 2.7 MeV。VIII.0、VIII.1 はしきい値の 1.74883 MeV の 5.2e-7。xs と同じ補間則で積分する
  examples/depletion/issue-0311/law7_norm.py で、72b8297f、2026-11-20）。xs の検査も 4 つの版の Be-9 を通す（読み込みの走査）。LAW=1 と
  MF5 の表の規格化は確かめていない（#0324）
- 入射エネルギーの間の補間は、LAW=1 と同じくファイルのとおり（ADR-0057。NJOY は INT=2 を単位基底に替える）
- P1 以上は、LAW=1 と同じく 0 とする（ADR-0082）。H-2、Li-6、Be-9 は軽く、重心の運動による実験室系の P1 は小さくない
  （NJOY の P1 の行の和は P0 の 4〜9 割。examples/depletion/verification.md）が、どれも 1.7 MeV 以上のしきい反応である
  （PWR の計算への影響は測っていない）

### 帰結

- H-2、Li-6、Be-9 の 4 つの版の評価を、輸送する核種として処理できる
- NJOY との散乱行列の差は、LAW=6 では NJOY の入射エネルギーの格子と実験室系の単位基底の補間（Li-6 で確かめた。H-2 の残りの原因は
  確かめていない）、LAW=7 では単位基底（ADR-0057）と `ll2lab` の間引きの重み、P1 以上は ADR-0082 による（差の大きさは
  examples/depletion/verification.md）
- n が 3〜5 の外の LAW=6、2 段の反応の LAW=6 と、lin-lin とヒストグラムの外の E' の表と余弦の補間則、2 つ未満の余弦、−1〜1 の外の余弦の
  LAW=7 と、∫∫ f dμ dE' が 1 から 1e-5 を超えてずれる LAW=7 は、エラーで止まる（#0323）
- H-2 の (n,2n) のしきい値の群の差の残りの原因は #0321、軽い核の P1 以上を 0 にしていること（ADR-0082）は #0322

### 確認

- テスト: `law6_neutrons_are_phase_space_in_the_centre_of_mass`、`phase_space_max_energy_follows_endf_102`、
  `phase_space_density_integrates_to_one`、`phase_space_spreads_by_the_centre_of_mass_motion`（実験室系の群ごとの割合を、
  重心系のエネルギーの細かい中点則と 2e-7 で比べる）、`law7_neutrons_integrate_the_cosines_into_a_laboratory_spectrum`、
  `law7_histogram_cosines_short_of_one_are_error`、`law7_single_cosine_is_error`、`law7_zero_distribution_at_threshold_passes`、`law7_histogram_cosines_of_a_normalized_distribution_take_the_left_table`、`law7_forward_cone_normalized_to_one_passes`、`law7_distribution_short_of_one_is_error`、`two_step_phase_space_is_error`、
  `phase_space_rows_carry_the_neutrons_above_threshold`、`phase_space_rows_follow_the_kinks_above_threshold`（しきい値の群の行を、
  入射エネルギーの細かい中点則と 1e-5 で比べる）、
  `a_phase_space_emission_takes_its_q_value_from_mf3`
- NJOY2016 GROUPR との比較と差の原因は、examples/depletion/verification.md の「MF6 の LAW=6、7 の散乱行列（NJOY との比較）」
