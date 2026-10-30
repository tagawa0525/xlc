---
status: open
labels: [core, lattice, format]
opened: 2026-10-30
closed:
milestone: M4
blocked-by: []
resolution:
---

# 履歴の反応率の補正（ADR-0215 の案 S2）のレビューで残したもの

## 問題

ブランチ feat/0190-history-micro-sigma（#0190）のローカルのレビュー（1 周目、HEAD 23863f7f。2 周目、3 周目、記録の周の分は後の項）で、再現できる誤りでない指摘を残した。
どれも出力の数値は変えない。

関係する ADR: ADR-0215

- 型: crates/core/src/case.rs の `node_micro` は、線形の補正の係数を読む `correct` と同じ型 `MicroAt` を返し、`fission`（Σf の係数でもある）と
  収率の境界より下の割合だけを補正する。補正した値が線形の補正に渡ると、ADR-0215 の案 S2 の「係数は補正しない」が崩れ、履歴の分岐の点に
  (cᵢ^j − cᵢ(j)) ΔN̄ᵢ が残る。今は `node_constants` が `interpolate_micro` を、反応率が `node_micro` を別に呼ぶので起きないが、型が止めない
  （反応率の断面積を別の型にする候補）
- テスト: crates/core/src/depletion/tests.rs の `micro_depletion_rates_take_the_weighted_rate_residuals_of_the_history` の w = 0 の半分は、
  Δσ を 0 にした `histories` を持つテーブルと比べ、`histories` を消したテーブルとは比べない（どちらも w = 0 で同じ経路を通る）。注釈は
  「補正のないテーブル」と書く。case の `nominal_histories_leave_the_rate_cross_sections_bit_identical` は `histories` を消して比べている
- エラーの文: 核分裂する追う核種の σf の群の和が内挿で 0 以下になると、`interpolate_micro` は割合を `None` にし、w が 0 でないノードでは
  `add_rate_residual` が「核分裂の有無が違う」と報告する（本当の原因は σf の和）
- 重複: crates/lattice/src/residual.rs の `residuals` で、Δσ と R が同じ燃焼度の点の走査と条件の点の引き直しを 2 回する（1 回の走査で組を返せる）。
  crates/core/src/state.rs の `scaled_rate`、`add_rate`、`add_rate_residual` の反応の経路の足し方が `scaled_micro`、`accumulate_micro` と同じ形で、
  (MT、経路の数) の形の比べ方が lattice の `rate_difference`、形式の `validate_rate_residuals`、`validate_micro` にある
- 効率: crates/core/src/depletion.rs の `Solve` が解ごとにノードの履歴の変数のコピーを持つ（平衡 Xe の反復でも解ごと）。反応率を作るときに
  読むだけなので、`n` と同じく引数で渡せる。炉心で効くかは測っていない

- 2 周目（HEAD 375cf059）: `add_rate_residual` は負の σ と割合の範囲を変数ごとに確かめるので、変数が 2 つ以上（制御棒、Pyrex）になると、
  途中の和で負になり後の変数で戻る σ を止める（すべての変数を足してから確かめる候補）。`node_micro` の `table.micro` の検査は、先に呼ぶ
  `interpolate_micro` が同じ条件で止めるので通らない。`add_rate_residual` の長さのエラーが `names` の長さを出さない。crates/core/src/rods.rs の
  `mix_micro` の注釈は σf の和が正である理由を `interpolate_micro` とだけ書き、履歴の補正のノードでは `add_rate_residual` の割合の検査が支える
- 3 周目（HEAD 15093f32）: examples/depletion/verification.md の途中で条件を変えた燃焼の表の見出しと本文の「w × 18.50」は、丸める前の
  18.50496（core-history.json の `applied_max_pcm`）で掛けた値で、2/3 × 18.50 = 12.333 と 12.34 が合わない（見出しを丸める前の値にするか、
  丸めた値で掛ける。表は main の版から同じ形）
- 記録の周（HEAD 4be63bb3。歯止めの 3 周の後で、指摘は直さない。ADR-0215 は proposed で、凍結する部分に当たりうる 3 件は親の判断に回した）:
  - ADR-0215 の検討した案（ミクロ燃焼の微視的断面積の段落）の補正前の数密度の差の出典「#0190 の経緯」は、Pu-241 −0.828% を載せていない
    （main の c421d15a の examples/depletion/verification.md にある）。差の燃焼度の点（60 MWd/kg）も書いていない
  - ADR-0215 の「理由」の S1 との差の文に Σᵢ と f(S2) − f(S1) の向きがない（検討した案と帰結には書いた）
  - examples/depletion/verification.md の途中で条件を変えた燃焼: 注記の「公称の状態にある同じ w のノードの目安」は Δσ の状態への依存（測って
    いない）を含まないことを書かず、ADR-0215 の確認の文と範囲が食い違う。「この燃焼の誤差は…比例の分」は、履歴の補正のない 1 の水準
    （−3.808〜+7.714 pcm）を除いた、履歴の補正の分に限って書くべき。表の列「点の誤差の上限」と本文の「誤差の上限」の名前が、この燃焼に
    当たらない目安になった注記と合わない（列名は compare_core_history.py と verification_tables.py が出す）

## 完了の条件

- 上の項目ごとに、直すか、直さない理由を書く

## 経緯

- 2026-10-30: feat/0190-history-micro-sigma のレビューで起こした
