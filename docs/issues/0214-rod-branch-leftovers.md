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
  `scripts/dvc-checkout.py` は `table-rods.msgpack` を戻さない。examples/core のステージは main の時点（b29e4a53）ですでに
  bin の変更で 57 本すべてが changed で（`scripts/count-changed-stages.py examples/core/dvc.yaml`）、`rods.msgpack` は今の
  ライブラリの形式で読めない。記録（examples/core/verification.md）の数値は ADR-0209 の前のもの
- 結合テスト `rod_branches_match_separate_lattices`（crates/lattice/tests/rod_branch.rs）は自己遮蔽なし、直積、チェックポイント
  なしで、値ごとの自己遮蔽、部分格子（値の最初の状態が状態 0 でない）、値の途中からの再開を通らない。自己遮蔽のある一致は
  VERA の集合体の 1 回の実行（examples/core/verification.md の「制御棒の分岐と値ごとの格子計算の結合の一致」）だけが支える
- 値の最初の状態で始める解を捨てるのは、`solve` の中の特例（`present_states.first()` との比較）で、`Starts` は制御棒の軸を
  普通の軸として扱い、前の値の解（領域の数が違う）を返し、チェックポイントにも残す。引き継がない軸を `Starts` に持たせれば
  特例が要らない
- 自己遮蔽に渡す領域の物質の名前は、値 0 だけ履歴の燃焼の領域の名前（`<物質>@<セル>`）で、ほかの値は形状の名前（main.rs の `solve`）。
  `[history]` と `[branches.rod]` が排他なので今は正しいが、#0213 で組むときは値ごとの形状に燃焼の領域の名前を持たせる。値ごとの
  トラッキングはライブラリを読むのと重ねず、実行の終わりまでメモリに残る
- 標準出力の `state <番号>` はテーブルの軸の直積の順だが、エラーの文脈とログ（`lattice` の `state=`）は解く順（値が最初の軸）の
  番号で、同じ状態に 2 つの番号が付く。標準出力の番号のテストもない
- 制御棒の軸の有無を、main.rs の 3 か所（`solve_counts` と `solved_points` の `rods > 1`、ログとエラーの `values` の
  `rod_values.len() > 1`）、`Branches::table_axes`、input.rs のノードの出力の検査（`!b.rods().is_empty()`）で決め、点の並べ替えは番号の列
  （`lexicographic`）と `Option` の取り出しで行う（点と状態の組を並べ替える形のほうが短い）
- ノードごとの出力と制御棒の分岐を組む入力は、棒の入らないノードを持たない格子（燃料の集合体だけのカラーセット）でもエラーに
  なる（ADR-0224）。平板（#0030）と径方向の反射体の格子（#0213）でノードごとに軸を付け分ける仕組みを作るとき、この場合も通す
- 媒質を出力する格子で、棒のある値の状態も細かい群の媒質を作ってチェックポイントに書き、最後に捨てる（計算と書き出しの無駄）

関係する ADR: ADR-0224

## 完了の条件

examples/core を今の bin とライブラリで再実行して `dvc.lock` と記録を書き直す（読む作業が来たとき）。結合テストに自己遮蔽の
ある小さい格子か、部分格子とチェックポイントの再開の場合を足す。始める解と媒質の扱いは、直すかを決める。

## 経緯

- 2026-11-03: #0093 のレビューで出た指摘のうち、落ちるテストか誤った数値を示せないものをまとめた
