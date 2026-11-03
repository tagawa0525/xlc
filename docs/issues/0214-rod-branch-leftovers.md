---
status: open
labels: [lattice, examples]
opened: 2026-11-03
closed:
milestone: M3
blocked-by: []
resolution:
---

# 制御棒の分岐の残り（examples/core の再実行、テストの範囲、始める解の軸）

## 問題

この issue は、#0093 の制御棒の分岐（ADR-0224）のレビュー（/code-review high、2026-11-03）で出た、直さずに残した指摘。

- examples/core の `dvc.lock` は、ステージを替えた（値ごとの 3 つの格子計算と rod_table.py の結合 → 1 つの格子計算）後に
  再実行していないので、新しい `lattice` の記載が無く、古い `lattice@out`、`lattice@aic`、`lattice@b4c`、`table` の記載が残る。
  `scripts/dvc-checkout.py` は `table-rods.msgpack` を戻さない。examples/core のステージは main の時点（f86d7fda）ですでに
  bin の変更で 57 本すべてが changed で（`scripts/count-changed-stages.py examples/core/dvc.yaml`）、`rods.msgpack` は今の
  ライブラリの形式で読めない。記録（examples/core/verification.md）の数値は ADR-0209 の前のもの
- 結合テスト `rod_branches_match_separate_lattices`（crates/lattice/tests/rod_branch.rs）は自己遮蔽なし、直積、チェックポイント
  なしで、値ごとの自己遮蔽、部分格子（値の最初の状態が状態 0 でない）、値の途中からの再開を通らない。自己遮蔽のある一致は
  VERA の集合体の 1 回の実行（examples/core/verification.md の「制御棒の分岐と値ごとの格子計算の結合の一致」）だけが支える
- 値の最初の状態で始める解を捨てるのは、`solve` の中の特例（`present_states.first()` との比較）で、`Starts` は制御棒の軸を
  普通の軸として扱い、前の値の解（領域の数が違う）を返し、チェックポイントにも残す。引き継がない軸を `Starts` に持たせれば
  特例が要らない
- 値ごとのトラッキングはライブラリを読むのと重ねず、実行の終わりまでメモリに残る（自己遮蔽に渡す領域の物質の名前を値ごとの
  形状の燃焼の領域の名前にすることは、#0213 の ADR-0225 で済んだ）
- 標準出力の `state <番号>` はテーブルの軸の直積の順だが、エラーの文脈とログ（`lattice` の `state=`）は解く順（値が最初の軸）の
  番号で、同じ状態に 2 つの番号が付く。標準出力の番号のテストもない
- 制御棒の軸の有無を、main.rs の 3 か所（`solve_counts` と `solved_points` の `rods > 1`、ログとエラーの `values` の
  `rod_values.len() > 1`）、`Branches::table_axes`、input.rs のノードの出力の検査（`!b.rods().is_empty()`）で決め、点の並べ替えは番号の列
  （`lexicographic`）と `Option` の取り出しで行う（点と状態の組を並べ替える形のほうが短い）
- ノードごとの出力と制御棒の分岐を組む入力は、どのノードにも置き換えるセルがある格子（すべての集合体の案内管に棒を入れるカラーセットなど）でもエラーに
  なる（ADR-0224）。平板（#0030）と径方向の反射体の格子（#0213）でノードごとに軸を付け分ける仕組みを作るとき、この場合も通す
- ノードごとに制御棒の軸を付け分ける仕組みを作るときに決めること: カラーセットのノードと、径方向の反射体の格子の燃料のノード
  （棒が入る）は ADR-0127 と ADR-0128 の外。平板の燃料の媒質のノードは ADR-0224 で「決めていない」としたが、ADR-0128 の
  「制御棒の値の定数は、その層全体に棒があり、燃料には棒がない平板から作る」からは棒の入らないノードとも読める（input.rs の
  検査のコメントは棒の入らないノードに数える）。#0030 で決める
- 媒質を出力する格子で、棒のある値の状態も細かい群の媒質を作ってチェックポイントに書き、最後に捨てる（計算と書き出しの無駄）

- #0213 の履歴からの分岐と制御棒の分岐の組み合わせ（ADR-0225）のレビュー（/code-review high、2026-11-03）で、直さずに残した指摘:
  `Starts::start` は解いた状態をすべて走査して添字の列を作り直す（呼ぶたびに解いた数に比例。番号の順に後ろから見て最初の
  前の軸が同じ状態で止めれば同じ点になる）。値ごとの均質化した数密度は、割る格子の体積だけが違うのに分子の和を値ごとに求め直す。
  値の燃焼点（`rod-burnup`）だけが燃焼点を間引いて公称の状態がないときのエラーの文が `off-nominal-burnup` を名指しする。
  数密度の微分のログ（`derivative`）の `state=` は値の中の番号で、`lattice` のログ（解く順の番号）と標準出力（テーブルの順の番号）と
  合わない。main.rs の値ごとの自己遮蔽と摂動の自己遮蔽を 2 つの `Option` で持ち、値の切り替えを 2 か所で書く（1 つにまとめれば
  食い違わない）。`HistoryBranches::rods` はテストだけが使う。2 周目: `rod-burnup` を書かない値は `off-nominal-burnup` を使うので、棒のない値の
  公称でない状態を間引くと制御棒の値も黙って間引く（ADR-0200 は制御棒を入れた組を測るまですべての燃焼点で解く。ADR-0225 の既定）。
  1 つの入力と別の入力のビットの一致は、`rod-burnup` が `off-nominal-burnup` より多くの点を持つ場合と、2 軸以上を変える部分格子では
  比べていない。`HistoryBranches::micro` は値の格子を引数でもう一度受け、値の均質化した数密度の体積と同じ格子かを確かめない。
  値ごとの形状の `ValueGeometry`（history.rs）は main.rs の `Geometry` のコピーで、`burnups` の値 0 の分岐（`off_nominal` を選ぶ特例）は、`HistoryInput::rod_burnup` が書かない値に `off_nominal` を返すのと同じ規則を
  2 か所に書いている。
  3 周目: history.rs のモジュール、`ValueShape`、`history_branches` の説明は、燃焼の領域の照合を（物質、組の代表のセル、体積）と
  書くが、コードは ADR-0226 の `region_mismatch`（組のセルの位置と MOC の領域の数も比べる）を使う。
  `rod_burnups_are_history_burnups_of_rod_values` は `is_err()` だけを見て、どの検査で止まったかを確かめない。

関係する ADR: ADR-0224、ADR-0225

## 完了の条件

examples/core を今の bin とライブラリで再実行して `dvc.lock` と記録を書き直す（読む作業が来たとき）。結合テストに自己遮蔽の
ある小さい格子か、部分格子とチェックポイントの再開の場合を足す。始める解と媒質の扱いは、直すかを決める。

## 経緯

- 2026-11-03: #0093 のレビューで出た指摘のうち、落ちるテストか誤った数値を示せないものをまとめた
