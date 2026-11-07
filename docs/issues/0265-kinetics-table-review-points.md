---
status: closed
labels: [lattice, core, xs]
opened: 2026-11-07
closed: 2026-11-08
milestone: M4
blocked-by: []
resolution: wontfix: レビューで直さなかった指摘を集めた issue は作らない決定（nixfiles の ADR-0066）に合わせ、再現のある項目が無いので閉じた
---

# 動特性と冷却計算のデータの実装に残したレビューの指摘

## 問題

ADR-0238 の実装（ブランチ feat/0153-kinetics-table-data）のレビューで出た、直さなかった指摘（再現を書けない文書、効率、重複）:

- crates/lattice/src/history.rs: `nominal_state` の前に足した関数が、`derivative_nuclides` の説明の行（「数密度の微分を持つ核種の名前…」）と
  関数の間に入り、その行が `nominal_state` の説明になっている。`derivative_nuclides` に説明がない
- crates/core/src/state.rs の `interpolate_kinetics`: 説明は遅発スペクトルの和が 1 から外れるとエラーとするが、コードは負の重みのときの負の値
  だけを確かめ、和を確かめない
- crates/xs/src/main.rs の `run`: ポイントワイズの範囲のために、輸送する核種の入力とテープを複製する。足した `transported` と
  `transported_tapes` を使っていない
- crates/xs/src/main.rs の `check_shared`: 群速度の逆数を `Weight::lwr()` で核種ごとに求め直す（`settings_of` の重みと 2 か所になる）
- crates/lattice/src/main.rs: 冷却計算の断面積の燃焼度の点を `solved_points[i].first()`（燃焼度が最初の軸）で取り、`HistoryBranches::point_of`
  を使わない。history.rs の `kinetics`、`cooling`、`homogenize_nuclides` も `index / self.branch.len()` を繰り返す
- crates/lattice/src/history.rs の `cooling`: 輸送する核種の 1 群の断面積に `homogenize_nuclides` を使い、使わない散乱行列も作る
- crates/lattice/src/history.rs: `micro`、`kinetics`、`cooling` が同じ状態の SPH 因子と少数群の中性子束をそれぞれ求める
- crates/lattice/src/history.rs の `history_branches`: 公称の状態の検査は入力だけで決まるが、ライブラリと履歴を読み、値ごとの形状の
  トラックを作った後に行う（`LatticeInput::parse` で止められる）
- crates/lattice/src/history.rs: 放射化の核種の燃焼の領域の温度のデータ（`cooling`）と、遅発中性子の族のデータの有無と族の数（`kinetics`）は
  ライブラリと分岐の温度で決まるが、最初の公称の状態を輸送した後に確かめる
- crates/lattice/src/input.rs の `[history]` の split の検査: 直前の `matches!` で `output` が `Some` と分かった後に `if let Some` で
  取り直し、`None` の枝が起きない形になっている
- crates/lattice/src/history.rs の公称の状態のないときのエラーの文: 点が 1 つの軸は入力と違っても公称になるのに、どの軸も入力の物質を残すことを
  求める文で、どの軸が原因かを示さない（branches.rs の部分格子のエラーは軸の名前を出す）。ADR-0238 の帰結も、止まる条件（点が 2 つ以上の軸に
  入力の物質の格子点がない）を書かない
- 点が 1 つの軸が入力の物質と違う値のとき、公称の状態 x₀ は履歴の輸送と違う状態で、冷却計算の断面積はその状態で黙って作られる（ADR-0238 の
  決定は x₀ と書く）。止めるか、許す理由を書くかを決めていない
- crates/lattice/src/history.rs: `nominal_state` の説明の先頭に `derivative_nuclides` の説明の行が残る（上の 1 項目と同じ）。`HistoryBranches.nominal` は
  検査の後には必ず値があるのに `Option` のままで、説明の「なければ None」は起きない
- crates/lattice/src/history/tests.rs: `history_branches_need_the_nominal_state_for_the_cooling_cross_sections` は文に "nominal" があることだけを確かめ、
  ほかの検査の文でも通る。`off_nominal_states_are_solved_at_their_burnups_only` の近くの説明は、間引いたときだけ公称の状態が要ると読める
- crates/lattice/src/history.rs の解く燃焼点のコメント: 「間引いたかによらず」は前の挙動と比べた書き方で、どこで何を求めるかもあいまい

- crates/lattice/src/main.rs のチェックポイントの適合（`fits`）: 保存した状態の動特性のデータの有無が違うときのエラーの文に、その有無が出ない
- crates/lattice/src/history.rs の `cooling`: 放射化の核種の燃焼の領域の温度のデータは、数密度が 0 で重みが 0 の領域では引かないので、温度の
  検査が数密度に依存する

- examples/depletion/kinetics_table.py: `--report` が位置引数のテーブルと `--nuclide` を黙って無視する。`cooling` か `kinetics` のないテーブルで
  KeyError になり、何が足りないかを言わない

- examples/depletion/kinetics_table.py の `nominal_points`: 公称の状態の点の番号を探す `index` の ValueError が、何が無いかを言わない
- crates/lattice/src/history.rs の `kinetics`: 族のデータを持つ核種がないときと、族の生成が 0 のときを同じエラーの文で言う

- examples/depletion/kinetics_table.py: 生成の検査は 4 つの条件（k、νΣf、k の有限、k_inf）を 1 つの文で言い、k_inf が 0 のときも「no production」と
  出る。k_inf が 0 の場合と k が有限でない場合のテストがない
- ADR-0238 の帰結は、作り直した ../p2/depletion.msgpack を読む入力（examples/depletion/branches の branch-2A.toml、branch-2F-all.toml、
  branch-2O-all.toml、shielding-2A-every.toml と、examples/depletion/history/residual-rods-branch.toml など）を数えない。前のライブラリで作った
  履歴と組むものもある
- ADR-0238 の帰結は、ほかの例のライブラリの作り直しに xs の時間だけを数え、その下流の格子計算のステージ（examples/depletion の lattice-1C、
  lattice-2A、lattice-mox1、examples/vera の格子計算など）の再実行の時間を数えない
- examples/depletion/verification.md の「動特性と冷却計算のデータ（2A）」と #0264 の H-1 の根拠（`thermal_scattering` と `resolved`）は、`lookup` の
  内挿が拒む残りの条件（表の形、チェーンの反応の経路、遅発中性子の族、S(α,β) の格子、ポイントワイズ）に触れない。H-1 の 565 K の内挿を
  支える出典には、p2 の格子計算の燃焼が 565 K の減速材の H-1 を内挿して解いていること（例えば lattice@2A のログ）を使える
- ADR-0238 の帰結の p2 の 7 本の燃焼の項: 「main の dvc.lock の bin/lattice が 2927836cb（2026-11-04）で別の値に移った」は事実と違う。2927836cb は
  別のブランチのコミット（bin/lattice を a277b17b にした）で、main の lock が 48387d65 から 7e959698 に移ったのはマージ e78c01e63（2026-11-05）。
  同じ項は、再実行が下流の compare@ を実行して記録の timed の時間を上書きしうること（#0270）に触れない
- 上の H-1 の項が出典に挙げた「lattice@2A のログ」には、減速材の H-1 の内挿を示す行が無い（燃焼の領域の H-1 の decay only の WARN だけ）。
  根拠になるのは、565 K の減速材に H-1 を持つ lattice@2A が lookup のエラーなしに最後まで走ったこと
- examples/depletion/kinetics_table.py: `summary` の生成 νΣf · φ は、`fundamental` の k を φ の和で割った値と同じ内積を計算し直し、生成の検査の
  k > 0 と生成 > 0 は同じ符号を 2 度見る

関係する ADR: なし

## 完了の条件

- 上の点を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-11-07: ADR-0238 の実装のレビューの 1〜3 周目で出た指摘のうち、直さなかったものをまとめた
- 2026-11-08: レビューで直さなかった指摘を集めた issue は作らない決定（nixfiles の ADR-0066 の「## 変更」の 2026-11-08）に合わせて閉じた。再現のある項目は無い
