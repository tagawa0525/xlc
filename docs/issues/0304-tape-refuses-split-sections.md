---
type: Issue
status: open
labels: [xs, depletion]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# ENDF のテープの読み込み（Tape）が、離れた場所に再び現れるセクションを止める

## 問題

- xlc-xs の `endf::Tape::section` は、同じ (MAT, MF, MT) の行のうち最初の連続した行だけを返す。同じ (MAT, MF, MT) が別の行を
  挟んで再び現れるテープ（同じ MAT の物質が 2 つあるなど）では、後の行を黙って読まない。`Tape::mts` は連続した重なりだけを
  除くので、同じ MT を 2 度返す。xs の `read_evaluation`、`final_states` はこのテープをエラーなく読む
- 燃焼チェーン（crates/chain/src/endf.rs の `materials`）は、#0303 で複数の物質のテープを読むときに、行の制御欄を自分で
  読み直して、このテープを止める。制御欄の読み方（非 ASCII の文字の置き換え、桁、空欄を 0 とすること）を `Tape::parse` から
  コピーしているので、`Tape::parse` を変えると食い違いうる
- 検査を `Tape::parse` か `Tape` の索引の側に置けば、xs を含むすべての読み手が止まり、チェーンのコピーを消せる。#0303 のチェーンの作業では
  xs を変えなかった。`bin/xs` を呼ぶステージは changed になるが、作り直しの残高で、それを理由に待たない（ADR-0258 の「changed の扱い」）

関係する ADR: なし

## 完了の条件

- [ ] `Tape` が、離れた場所に再び現れる (MAT, MF, MT) をエラーにする（テストは小さな合成のテープ）
- [ ] crates/chain/src/endf.rs の `control` と `materials` の重なりの検査を消し、`Tape` の検査に任せる

## 経緯

- 2026-11-20: #0303 の燃焼チェーンのブランチのレビューで、根本の原因が xs の `Tape::section` にあると指摘されて起こした
