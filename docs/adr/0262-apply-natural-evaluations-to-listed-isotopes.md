---
type: ADR
status: accepted
date: 2026-11-20
requires: [ADR-0079, ADR-0164, ADR-0175]
supersedes: []
superseded-by:
issues: ["#0303"]
---

# ADR-0262: 燃焼チェーンは、天然の元素の中性子の評価（A = 0）を、入力が挙げた同じ元素の核種の反応として当てる

## 背景

燃焼チェーン（ADR-0164）は、中性子の評価を MF1/MT451 の (ZA, LISO) で decay の核種に対応させる。本番の版を決め直す検討
（#0303）で、ENDF/B-VII.0 と VII.1 からもチェーンを作る。この 2 つの版の中性子のファイルには、天然の元素の評価（A = 0）が
ある。VII.0 は C、V、Zn（`n-006_C_000.endf`、`n-023_V_000.endf`、`n-030_Zn_000.endf`）、VII.1 は C だけで、どちらもその元素の
核種ごとの評価を持たない。VIII.0 と VIII.1 には天然の元素の評価がない。decay には A = 0 の核種がないので、`bin/chain` は
「neutron data for (Z = 6, A = 0, state 0), which is not in the decay data」で止まる。

xs の側では、VIII.0 に天然の C がないので、VERA の天然の C の個数密度を C-12 と C-13 に分ける（ADR-0079）。VII.0 の
検証用ライブラリは、C-12 と C-13 の両方に天然の C のファイルを当てる（ADR-0175。入力 `examples/vera/endfb70/xs.toml` の
`[nuclide.C-12]` と `[nuclide.C-13]` が `n-006_C_000.endf` を指す）。どの核種に天然の元素の評価を当てるかは、入力が決めている。

判断の軸:

- 黙ったフォールバックを足さない（怪しい入力はエラーで止める）
- xs の側の天然の C の扱い（ADR-0079、ADR-0175）と食い違わない。xs が断面積を作る核種に、チェーンが反応を持たないことがない
- 新しいデータの表（存在比など）を持たない

## 検討した案

- 案 A: 天然の元素の評価をチェーンのどの核種にも当てず、外したことを警告に出す。C-12 と C-13 は反応を持たないので、
  xs のライブラリが C-12 と C-13 に天然の C の断面積を持っても、燃焼で (n,γ) と (n,α) による減りと生成が起きない。xs の側と
  食い違う
- 案 B: decay で安定（NST = 1）な同じ元素の核種のうち、自分の評価のないものに自動で当てる。VII.0 の decay では、天然に
  ある V-50 と Zn-70 が不安定（NST = 0）なので外れる（同じ版の K-40 も不安定）。安定と天然にあるとは違い、規則が黙って核種を落とす
- 案 C: 天然の存在比の表を持ち、存在比が 0 でない核種に当てる。表を作り、保つ費用が要る。xs の側は存在比を持たず、入力で
  核種とファイルを結ぶので、2 つの側で当てる核種が食い違いうる
- 案 D（採用）: チェーンの入力（`natural-isotopes`）が当てる核種を名前で挙げ、天然の元素の評価を、挙げた同じ元素の核種
  それぞれの評価として使う。xs の入力（ADR-0175）と同じく、入力が核種とファイルを結ぶ

## 決定と理由

案 D にする。xs の側と同じく入力で当てる核種を決め、規則で黙って核種を足したり落としたりしない。

- 入力 `natural-isotopes`（核種の名前の並び。ADR-0015 の形）を足す。省くと空にする（天然の元素の評価のない VIII.0 の
  入力 `examples/depletion/chain.toml` を変えないため）
- 中性子の評価のうち Z > 0、A = 0 のものを天然の元素の評価とし、`natural-isotopes` のうち同じ Z の核種それぞれに、その MF3 の
  MT、MF8 の生成物の状態、核分裂の Q 値をそのまま当てる。反応の生成物は、当てた核種の (Z, A) から ADR-0164 の規則で決める
  （C-12 の (n,γ) は C-13、C-13 の (n,γ) は C-14）
- 当てたことは、元素ごとに 1 つの警告に出す（データにない対応を規則で補ったので。ADR-0164 の警告と同じ扱い）
- 次はエラーにする: どの核種も挙げない天然の元素の評価、天然の元素の評価のない元素の核種、decay にない核種、準安定状態、
  同じ核種の 2 度の記載、同じ元素の 2 つの天然の元素の評価。自分の評価のある核種を挙げたときは、同じ核種の評価が 2 つになり、
  ADR-0164 のとおり止まる

版ごとの入力で挙げる核種は、天然に存在する同位体（VII.0: C-12、C-13、V-50、V-51、Zn-64、Zn-66、Zn-67、Zn-68、Zn-70。
VII.1: C-12、C-13）にする。xs の側で天然の元素の評価を当てる核種を足すときは、両方の入力を揃える。VII.0 の入力は、decay の
Nb-88m の負の励起エネルギーの扱いを決めてから置く（#0303）。

### 帰結

- 天然の元素の評価の MF8 が生成物の状態を挙げると、その ZA は天然の元素の評価のもので、当てた核種の生成物と食い違い、
  ADR-0164 の検査で止まる。VII.0 と VII.1 の天然の元素の評価では止まらない（下の「確認」。VII.0 は Nb-88m の励起エネルギーを
  ADR-0263 で扱った入力で確かめた）
- 当てた核種の反応の断面積は、xs の側で同じ評価から作る。チェーンと xs の入力の挙げる核種が食い違っても、どちらのツールも
  気づかない（チェーンは xs の入力を読まない）。版ごとのライブラリの入力を作るとき（#0303 の次の作業）に揃える

### 確認

- 規則: crates/chain/src/build/tests.rs（`applies_a_natural_evaluation_to_each_listed_isotope`、
  `rejects_a_natural_evaluation_that_no_listed_isotope_takes`、`rejects_a_listed_isotope_without_a_natural_evaluation_of_its_element`、
  `rejects_a_listed_isotope_that_has_its_own_evaluation`、
  `rejects_a_listed_isotope_outside_the_decay_data_or_in_an_excited_state`、`rejects_two_natural_evaluations_of_one_element`）
- 実データ: VII.1 の入力（`examples/depletion/endfb71/chain.toml`）と、Nb-88m の励起エネルギーを ADR-0263 で扱った VII.0 の
  入力（`examples/depletion/endfb70/chain.toml`）で `bin/chain` が終了コード 0 でチェーンを作る（#0303 の経緯）
