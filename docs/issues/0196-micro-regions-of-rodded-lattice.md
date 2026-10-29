---
status: closed
labels: [lattice, depletion]
opened: 2026-10-29
closed: 2026-10-29
milestone: M4
blocked-by: []
resolution: コミット 139ddd8a5
---

# 制御棒を入れた格子の履歴からの分岐で、微視的断面積が履歴の MOC の領域の番号を引く

## 問題

履歴からの分岐（crates/lattice/src/history.rs の `HistoryBranches`）は、燃焼の領域を入力の格子から作り直し（`depletion_regions`）、
履歴の燃焼の領域と（物質、組の代表のセル、体積）が同じことを確かめる。コメントは「MOC の領域の番号は入力のものを使う（案内管の中身を
変えた格子でも、燃料のセルが同じなら同じ領域になる）」と書くが、`HistoryBranches` に入れるのは `regions: history.regions.clone()`
（履歴の格子の MOC の領域の番号）で、作り直した `regions` は照合にだけ使う。組成（`region_material`）は物質とセルの名前で引くので
番号によらないが、微視的断面積（`HistoryBranches::micro`。ADR-0216）は `d.regions` の番号で自己遮蔽と MOC の中性子束を引くので、案内管の
セルの領域の数が履歴の格子と違う格子（制御棒を入れた格子）では、別の領域を引く。

再現（2026-10-29、main 596c2437 の `bin/lattice`）: examples/depletion/branches/branch-2A-micro.toml の案内管を AIC の制御棒にした入力
（`branch-2A-micro-rods.toml`。`cr` のセルは 12 のリング、`gt` は 7）を、燃焼度 20 MWd/kg、1 状態に絞って解くと、
``micro cross sections: region 27140 of material `fuel`: 4 shielded nuclides for 89 transported`` で止まる（領域 27140 は、制御棒の
格子では燃料でない領域。4 核種は水の候補）。同じ絞り方の棒のない入力は通る。領域の数が違っても、ずれた先が燃料の領域なら、
エラーにならずに別の領域の中性子束で重みを付けうる（確かめていない）。

この不具合は #0192 の確認（部分挿入のノードのミクロ燃焼。制御棒を入れた組の微視的断面積つきのテーブルが要る）を止める。

制御棒を入れた組は、制御棒を入れた格子の入力で同じ履歴を読んで作る（ADR-0200）。微視的断面積の定義は ADR-0216。

関係する ADR: ADR-0217

## 完了の条件

- 微視的断面積が、入力の格子の燃焼の領域の MOC の領域の番号を引く（修正の候補は、`HistoryBranches` の `regions` に作り直した
  `regions` を入れること。確かめていない）
- 案内管のセルの領域の数が履歴と違う格子で、微視的断面積が同じ燃料の領域を引くことを試験する

## 経緯

- 2026-10-29: #0192 の確認（feat/0192-partial-rod-micro-depletion）で、2A の制御棒を入れた組のテーブルを作る前に 1 状態を手元で測って見つけた
- 2026-10-29: `HistoryBranches` に、入力の格子から作り直した燃焼の領域（MOC の領域の番号は入力の格子のもの）を入れた（139ddd8a5）。
  試験（`micro_cross_sections_take_the_regions_of_the_branch_lattice_when_the_guide_tubes_differ`）は、2 × 2 の格子で案内管の
  セルの領域を 2 から 5 にし、重みを付ける領域がどれもその燃焼の領域の物質（`fuel@<セル>`）であることと、Σᵢ σᵢ N̄ᵢ が分岐の格子を
  直接解いた燃料の反応率と合うことを確かめる。ずれた先が燃料の領域になる場合も、前者の照合が止める
