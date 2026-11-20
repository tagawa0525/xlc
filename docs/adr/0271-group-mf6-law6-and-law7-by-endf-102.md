---
type: ADR
status: accepted
date: 2026-11-20
requires: [ADR-0057, ADR-0082]
supersedes: []
superseded-by:
issues: ["#0311"]
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
- ENDF-102 の定義のとおりに振り分ける（採用）

## 決定と理由

ENDF-102 の定義のとおりに群に振り分け、NJOY の表の作り方は真似ない。読み込みで `tracing::warn!` を出す。

- LAW=6 は、重心系で等方な密度 C_n √E' (E_max − E')^(3n/2 − 4)、E_max = (A_p − 1)/A_p (A/(A + 1) E + Q) を、入射エネルギーの
  積分点ごとに実験室系の群に振り分ける。Q は MF3 の QI。LCT によらず重心系とする（ENDF-102 の 6.2.8）。n は ENDF-102 が係数を
  与える 3〜5 に限る。重心系の E' = E_max sin²θ について積分すると密度の両端の √ の特異性が消えるので、実験室系の範囲が群の境界を
  横切るところで θ を分けて Gauss-Legendre で積分する。入射エネルギーの格子を置かないので、NJOY の 1.2 倍の格子とその補間の誤差が無い
- LAW=7 は、同じ E' での値を余弦の補間則（lin-lin かヒストグラム）で補間して −1〜1 で積分したエネルギー分布 g(E') にする
  （NJOY の `ll2lab` の P0 と同じ積分）。E' の表は lin-lin かヒストグラムに限る。そのとき g は全部の表の点の間で線形なので、
  ヒストグラムの段を同じ E' の 2 点として残した lin-lin の表で厳密に表せる。LAW=7 は実験室系の分布なので、重心系の LCT と組んだら
  エラーにする
- 入射エネルギーの間の補間は、LAW=1 と同じくファイルのとおり（ADR-0057。NJOY は INT=2 を単位基底に替える）
- P1 以上は、LAW=1 と同じく 0 とする（ADR-0082）。H-2、Li-6、Be-9 は軽く、重心の運動による実験室系の P1 は小さくない
  （NJOY の P1 の行の和は P0 の 4〜9 割。examples/depletion/verification.md）が、どれも 1.7 MeV 以上のしきい反応である
  （PWR の計算への影響は測っていない）

### 帰結

- H-2、Li-6、Be-9 の 4 つの版の評価を、輸送する核種として処理できる
- NJOY との散乱行列の差は、LAW=6 では NJOY の入射エネルギーの格子と実験室系の単位基底の補間（Li-6 で確かめた。H-2 の残りの原因は
  確かめていない）、LAW=7 では単位基底（ADR-0057）と `ll2lab` の間引きの重み、P1 以上は ADR-0082 による（差の大きさは
  examples/depletion/verification.md）
- n が 3〜5 の外の LAW=6 と、lin-lin とヒストグラムの外の E' の表と余弦の補間則の LAW=7 は、エラーで止まる

### 確認

- テスト: `law6_neutrons_are_phase_space_in_the_centre_of_mass`、`phase_space_max_energy_follows_endf_102`、
  `phase_space_density_integrates_to_one`、`phase_space_spreads_by_the_centre_of_mass_motion`（実験室系の群ごとの割合を、
  重心系のエネルギーの細かい中点則と 2e-7 で比べる）、`law7_neutrons_integrate_the_cosines_into_a_laboratory_spectrum`、
  `law7_histogram_cosines_weight_the_left_tables`、`phase_space_rows_carry_the_neutrons_above_threshold`、
  `a_phase_space_emission_takes_its_q_value_from_mf3`
- NJOY2016 GROUPR との比較と差の原因は、examples/depletion/verification.md の「MF6 の LAW=6、7 の散乱行列（NJOY との比較）」
