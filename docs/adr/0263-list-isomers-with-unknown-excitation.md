---
type: ADR
status: accepted
date: 2026-11-20
requires: [ADR-0164]
supersedes: []
superseded-by:
issues: ["#0303"]
---

# ADR-0263: 燃焼チェーンは、負の ELIS を持つ準安定状態を、入力が挙げたときだけ励起エネルギーの分からない準位（0）として受け入れる

## 背景

燃焼チェーン（ADR-0164）は、核種の励起エネルギーを decay の MF1/MT451 の ELIS から取り、MF8 の生成物の状態を、励起エネルギーの
相対の差の小さい組から準安定状態に対応させる。チェーンの検査（`Chain::check`）は、負の励起エネルギーと、0 でない基底状態の
励起エネルギーを止める。準安定状態の励起エネルギーを評価が与えないときは 0 と書く（同じ検査の注釈）。

本番の版を決め直す検討（#0303）で ENDF/B-VII.0 からチェーンを作ると、decay の Nb-88m の ELIS が -1000 eV で、この検査で止まる。
ENSDF で励起エネルギーの分からない準位の値で、同じ Nb-88m を VII.1、VIII.0、VIII.1 は 1e6 eV と書く。4 つの版の decay で、
負の ELIS はこの 1 件だけである（VII.0 の 3,838 の物質と、ほかの版の全ファイルの MF1/MT451 の 2 行目を見た。2026-11-20）。

判断の軸:

- 黙ったフォールバックを足さない（怪しい入力はエラーで止める）
- 4 つの版のどれからもチェーンを作れる

## 検討した案

- 案 A: 負の ELIS を止めたままにする。VII.0 のチェーンを作れない
- 案 B: 準安定状態の負の ELIS を、すべて 0（与えられない）にして警告に出す。最初の実装の形で、レビューで、ファイルの符号の誤り
  （Am-242m の -48630 など）も 0 にして通し、MF8 の状態の対応で 0 にした準安定状態が負けて別の状態に対応させうる、と指摘された
- 案 C: 後の版と同じ 1e6 eV に置き換える。後の版の値が分からない準位の約束の値なのかを確かめていない。版をまたいだ値を持ち込む
- 案 D（採用）: 入力 `unknown-excitations` が挙げた準安定状態の負の ELIS だけを 0（与えられない）にし、警告に出す。挙げない負の
  ELIS は止める

## 決定と理由

案 D にする。評価の書き方の癖を、核種を名前で挙げて受け入れ、規則で黙って値を変えない（天然の元素の評価と同じ形。ADR-0262）。

- 入力 `unknown-excitations`（核種の名前の並び）を足す。省くと空にする（VIII.0 の入力 `examples/depletion/chain.toml` を変えないため）
- 挙げた核種の ELIS を 0 にし、MF8 の状態の対応の前に置き換える。置き換えたことは核種ごとの警告に出す
- 次はエラーにする: 挙げない準安定状態の負の ELIS、decay にない核種、基底状態、ELIS が負でない核種、同じ核種の 2 度の記載。
  基底状態の 0 でない励起エネルギーは、今までどおり `Chain::check` が止める
- VII.0 の入力 `examples/depletion/endfb70/chain.toml` は Nb-88m を挙げる

### 帰結

- 0 にした準安定状態は、MF8 の状態の対応で、どの状態とも相対の差が 1 になる。VII.0 の中性子の評価に、Nb-88 を生成物の状態に
  挙げる MF8 はないので、VII.0 のチェーンでは対応に効かない（`bin/chain` の MF8 の対応の警告に Nb-88 は出ない）
- 入力の型は、天然の元素の評価を当てる核種（ADR-0262）と合わせて、評価にないものを入力で補う 1 つの値（`Supplements`）にする

### 確認

- 規則: crates/chain/src/build/tests.rs（`records_a_listed_negative_isomer_excitation_as_not_given`、
  `rejects_an_unlisted_negative_isomer_excitation`、`rejects_a_listed_isomer_whose_excitation_is_not_negative`、
  `rejects_a_ground_state_with_an_excitation`）
- 実データ: VII.0 の入力で `bin/chain` が終了コード 0 でチェーンを作る（#0303 の経緯）
