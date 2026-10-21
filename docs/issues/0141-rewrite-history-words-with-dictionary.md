---
status: closed
labels: [docs]
opened: 2026-10-21
closed: 2026-10-21
milestone:
blocked-by: []
resolution: コミット 0f118e06
---

# 履歴の語を辞書で書き換える（nixfiles の history-rewrite）

## 問題

pre-commit の辞書（nixfiles の ADR-0014、prh.yml）に 2026-10-21 に規則（実行する、燃焼させる、ツール、チェーン、
中性子収支、バランス、支持格子、ステージ、確率テーブル、シード、クロージャ、再利用、監視する、非アクティブ/アクティブバッチ
など）が入った。hook はステージしたファイルの全体を見るので、辞書に当たる語が残る既存の文書とコードは、語を全部直さないと
どの行も直せない（2026-10-21、main edda0385 で docs/plans/m4-depletion.md に 52 か所、CLAUDE.md に 16、
crates/lattice/src/depletion.rs に 20）。履歴の中の語も同じ辞書で直し、記録の hash を付け直す。

ツールは nixfiles の history-rewrite（ADR-0045。辞書で語を直し、番号は `commit-refs rewrite-history`（ADR-0013）で
付け直す。手順は nixfiles の .claude/skills/history-rewrite/SKILL.md）。前回の書き換え（2026-10-11）のツール
references/history-rewrite/（git で追跡しない）は、これで置き換える。

関係する ADR: なし（ツールの決定は nixfiles 側）

## 完了の条件

- すべてのブランチとタグの語が辞書に合い、記録の hash が付け直され、すべての worktree のブランチが新しい履歴を指す
- 未コミットの変更が壊れていない。行う直前に、すべての worktree の `git status` で列挙する
- references/history-rewrite/ を消し、#0077 の参照を history-rewrite に変える

## 手順

1. nixfiles の PR（ブランチ feat/history-rewrite-tool）をマージして rebuild する
2. 誰も作業していないとき（worktree のブランチの `Claude-Session:`、`squeue`、AWS の監視）に 1 度だけ行う。
   先にミラーを保存する
3. 除くパス: .textlintignore に挙げたファイル（#0142 の調査の表は、直す対象の語をわざと含む）
4. 書き換えで変わった箇所を見て、辞書をチューニングする（例: 同じ語を文脈で別の語に分けたい規則）。辞書で直せない語は #0142
5. 古いオブジェクトの掃除は #0077

## 経緯

- 2026-10-21: 辞書に規則が入り、既存の文書が直せなくなった。ツールは nixfiles の ADR-0045 で作成中
- 2026-10-21: 書き換えた（3478 コミットのうち 2780 が変わり、28006 か所の番号を付け直した）。除くパスは
  docs/issues/0142-fix-context-dependent-terms.md だけ。書き換えた後の main（0f118e06）で、追跡している 1194 ファイルが辞書を
  通る。feat/core-cmfd の未コミットの変更は、退避したパッチと同じ中身に戻った。コメントと文書を直した deps のある DVC の
  ステージ（build、OpenMC のライブラリなど）は changed になるが、実行せずにそのままにする（ADR-0130）。辞書の調整は
  要らなかった（文脈で変わる語は #0142）。references/history-rewrite/ を消した
