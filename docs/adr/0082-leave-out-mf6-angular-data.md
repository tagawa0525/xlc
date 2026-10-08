---
status: accepted
date: 2026-10-08
supersedes: [ADR-0056]
superseded-by:
issues: []
---

# ADR-0082: MF6 の角度のデータは散乱行列に入れない

## 背景

MF4 の角度分布は散乱行列に入れる（[ADR-0081](0081-use-mf4-angular-distributions.md)）。MF6 の角度のデータ（LAW=1 の Kalbach-Mann と Legendre、LAW=2 の
Legendre）も入れるかを、対象の核種で調べて決めた（2026-10-08）。M2 の計画では、MF6 の角度のデータを使うことに一時して
いたが、実装していなかった（下の「検討した案」）。

- PWR への影響
- 入射エネルギーの間の単位基底の補間に角度を通す手間

## 検討した案

- MF6 の角度のデータも入れる。計画で一時採った（2026-10-07、df9171a）: 連続的な放出のうち角度分布を持つ形式（MF6 の
  Kalbach-Mann など）はその角度分布を使い、MF6 の LAW=2 の離散準位は重心系の角度分布を使うとしていた。実装の前に、
  下の調査で採らないことにした。入射エネルギーの間の単位基底の補間に角度を通す手間がかかる
- 入れない（採用）

## 決定と理由

MF6 の角度のデータは入れない。

- MF6 の連続的な放出（LAW=1）の P1 以上は 0 とする
- MF6 の LAW=2, 3 で与えられる離散準位の非弾性散乱は、MF3 だけの離散準位と同じく、Q 値による重心系で等方な2体反応と
  する（[ADR-0056](0056-treat-mf6-law2-law3-as-two-body.md) の扱いを引き継ぐ）。LAW=3 は重心系で等方なので厳密で、LAW=2 は角度分布を無視する近似になる。P1 以上には
  重心系で等方な2体反応の運動学の分だけが入る

角度のデータを持つ MF6 の中性子は、重心系の Kalbach-Mann（O-16、B-11、U-234、U-236 の (n,2n) や連続の非弾性散乱、
構造材の MT5）と、実験室系の Legendre（Si-28 の MT91）だけで、WIMS 69 群の上端 10 MeV より下に入るのは O-16 の
(n,nα)（8 MeV 以上）、Si-28 の MT91（9.5 MeV 以上）と、微量の U-234, U-236 の反応だけである。重心系で等方な放出
（U-235、U-238、Zr など）の実験室系の P1 も、重い核では重心の速さが小さくわずか。PWR への影響は無視できる大きさと
見込み、手間に見合わない。

## 確認

上の調査（2026-10-08）。影響の大きさは測っていない。2体反応の扱いは `two_body_mf6_level_uses_q_value`、
`two_body_neutron_law_is_marked`、`isotropic_two_body_neutron_law_is_marked`、`lab_frame_two_body_angles_are_error`。

## 結果

MF6 の連続的な放出の P1 以上は 0 である（`crates/xs/src/multigroup/scatter.rs` のコメント。LAW=2 の角度のデータを
読み飛ばすことは `crates/xs/src/endf/mf6.rs`）。
