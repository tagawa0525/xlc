---
status: open
labels: [lattice, core, xs]
opened: 2026-11-07
closed:
milestone: M4
blocked-by: []
resolution:
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

関係する ADR: なし

## 完了の条件

- 上の点を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-11-07: ADR-0238 の実装のレビューの 1 周目と 2 周目で出た指摘のうち、直さなかったものをまとめた
