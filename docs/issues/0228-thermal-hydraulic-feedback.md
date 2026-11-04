---
status: open
labels: [core, lattice]
opened: 2026-11-05
closed:
milestone: M4
blocked-by: []
resolution:
---

# 熱水力のフィードバック（Problem 6、7、9）の ADR と実装

## 問題

VERA の Problem 6〜9 は HFP で、仕様書（CASL-U-2012-0131-004、Problem 9 の節、85 ページ）は「Thermal-hydraulic feedback is required
for at-power conditions」とし、燃料の温度を、温度、照射、燃焼度とともに変わる燃料棒のモデルで求めることを求める。Problem 9 の入力は
定格の流量 131.7 Mlbs/hr（バイパスの 9% を除いた値で、サイクルを通じて一定）、圧力 2250 psia（表 P9-1）と、ステップごとの入口の温度
（表 P9-2、554.9〜558.8 °F）で、熱水力の物性は与えない（docs/plans/m4-depletion.md の「熱水力のフィードバック」）。

今の xlc には熱水力のモデルが無く、決める ADR も無い（[ADR-0029](../adr/0029-defer-thermal-hydraulic-feedback.md) は初期のマイルストーンで
扱わないと決めた）。core はケースの状態ごとに減速材の温度と燃料の温度を 1 つずつ持ち、炉心のすべてのノードで同じ値を使う
（crates/core/src/input.rs の `CaseState`、crates/core/src/state.rs の `State`。減速材の密度はテーブルの格子点の組で決まる）。
分岐のテーブルの温度の格子点は HZP の 550〜580 K（ADR-0097）で、HFP の燃料の温度と減速材の密度の範囲を持たない。

これを待つもの（docs/plans/m4-depletion.md の「進め方」）:

- 6: 本番の分岐のテーブルの軸の格子点（HFP の範囲）。ADR-0033 のデータ（#0153）はその前に入れる。HFP の範囲の 3 次の項と M×F の誤差の
  測定（#0176）
- 7: core の燃焼の入口の温度（ADR-0218 は、熱水力の ADR の後まで状態の減速材の温度で書く形にとどめた）
- 9: Problem 9 の燃焼と臨界ホウ素濃度の比較

関係する ADR: なし

## 完了の条件

- 計画の「熱水力のフィードバック」の「ADR で決める」の項目（M4 で入れるか、流路の粒度、燃料の温度のモデル（熱伝導度の燃焼度依存、
  ギャップのコンダクタンス））を ADR で決める
- core がノードごとの減速材の温度と密度、燃料の温度を熱水力の解から求め、テーブルをノードごとの状態で内挿して、出力との繰り返しを
  解く。Problem 6 か 7 の入力で解ける

## 経緯

- 2026-11-05: 計画の「進め方」8 が未着手で、追う issue が無かったので起こした（ブランチ docs/m4-plan-progress）
