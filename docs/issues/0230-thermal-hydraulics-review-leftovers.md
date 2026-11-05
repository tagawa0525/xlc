---
status: open
labels: [core]
opened: 2026-11-06
closed:
milestone: M4
blocked-by: []
resolution:
---

# 熱水力のフィードバックのレビューの残り（効率、重複、将来のエラーの形）

## 問題

ブランチ feat/0228-thermal-hydraulics（ADR-0229）の 1 周目のローカルのレビュー（/code-review high）の指摘のうち、今の入力で誤動作を
再現できないもの（/git-merge の「直す指摘の基準」の non-blocking）:

- crates/core/src/case.rs の `check_densities` は、すべてのテーブルの減速材の軸の格子点を IF97 の領域 1 で確かめるので、熱水力の圧力の
  飽和温度（15.5132 MPa で 618.010 K。IF97 の式 (31)）より高い格子点を持つテーブルは、ノードがその点に届かなくても止まる。今のライブラリの温度の上限は
  600 K で届かない。#0228 の「残り」でライブラリに高い温度を足すとき、格子点を飽和温度の下に置くか、飽和温度より上の点の扱いを決める
  - 2026-11-06: ADR-0230 で格子点を飽和温度（IF97 で 618.010 K）の下に置き（Problem 6 の減速材の上端 610 K）、`check_densities` は
    変えないと決めた（案 M4 を採らない）。飽和温度より上の点を持つテーブルの扱いは、要るテーブルが出るまで残す
- crates/core/src/input.rs の `CaseState::at` は、入口の温度の状態を減速材と燃料が入口の温度の一様な状態にする。今は燃焼のケースと
  熱水力の組み合わせを入力で止めているので、燃焼の側（crates/core/src/depletion.rs の `at` の呼び出し）には届かない。#0228 の燃焼の
  熱水力を入れるとき、`at` を入口の温度の状態でエラーにし、熱水力の解の基準の状態を別に作る
- crates/core/src/thermal/rod.rs の `pellet_points` は `gauss_legendre(16)` のエラーを捨て、点がないという別のメッセージにする
  （16 点の Newton 法は決まった入力で失敗しない）
- 効率: `film_coefficient` が IF97 の領域 1 を 2 度計算し（`conductivity` の中で再び）、`water::temperature` が飽和温度と飽和の
  エンタルピーをノードごとに計算し直す。`solve_thermal` は探索のホウ素濃度ごとに `Channels::new` と `check_densities` をやり直す
- 重複: `ThermalReport` は `ChannelSolution` の項目を写し、`EffectiveEntry` は `Effective` の写し。核分裂のあるテーブルの判定
  （`thermal::fissile`）は crates/core/src/state.rs の重みの判定と同じ式で、テーブルか state.rs に置く方が依存の向きに合う

2 周目（同じ日）の non-blocking の指摘:

- `Channels::new` は反射体の位置を核分裂の判定の前に外すので、反射体の位置が核分裂のあるテーブルを使う入力は、2 回目の繰り返しで
  「流路のない位置の出力」のエラーになる（止まるが、入力の誤りを指さない）。作るときに止める
- `solve_thermal` の「流路に核分裂のあるノードがない」のエラーは、`Channels::new` が核分裂のある流路だけを作るので届かない。
  `Channels::present` は呼ばれていない
- 物質の番号から部分挿入のノードの区間のテーブルを引く書き方が、`density_error`、`Channels::new`、燃焼の `node_uses` の 3 か所にある
  （`Layout` のメソッドにまとめる）
- `report_thermal` の層ごとの燃料の温度の平均は、出力が正のノードで数え、`fuel_average` は核分裂のあるテーブルで数える
- 探索のある熱水力の状態の標準出力の「searched in N solves」はホウ素濃度の数で、各濃度の熱水力の繰り返しの固有値計算の回数を含まない。
  `thermal-hydraulics iterations` は最後の濃度の回数だけ
- 1 周目の直しのコミット（2963bf427）は、テストの入力の `[state.hfp]` の追加を実装と同じコミットに入れた（TDD のコミットの分け方）

関係する ADR: ADR-0230

## 完了の条件

上の項目ごとに直すか、直さない理由を書く。前の 2 つは #0228 の続き（ライブラリの温度の点、燃焼の熱水力）の中で扱う

## 経緯

- 2026-11-06: feat/0228-thermal-hydraulics の 1 周目のレビューの non-blocking の指摘から起こした
