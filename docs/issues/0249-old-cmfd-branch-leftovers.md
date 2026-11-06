---
status: open
labels: [docs]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# 旧ブランチ feat/core-cmfd の片付けの残り

## 問題

ブランチ docs/0076-close-old-cmfd-branch（#0076 を閉じた）のレビュー（medium、1 周目）の non-blocking の指摘。旧ブランチ feat/core-cmfd と
worktree xlc-feat-core-cmfd は 2026-11-06 に消し、先端 68f762520 にタグ archive/feat-core-cmfd を付けた（docs/branches.md、#0076）。

- 未コミットだった `crates/numeric/src/cmfd.rs` と `lib.rs` の変更は、タグに入らない。残るコピーは履歴の書き換え（#0141）の控え
  `~/.local/share/history-rewrite/xlc-2026-10-21/` の `cmfd.rs.before`（10,077 バイト）と `worktree-2.after.patch`（lib.rs の
  `pub mod cmfd;` の 1 行）で、2026-10-21 の時点のもの（消した worktree の最後の中身と同じかは確かめていない）。#0077 の完了の条件は、
  この作業ディレクトリを消す前に控えのパッチの変更がコミットされたかを確かめることを求めるが、変更はそのままはコミットしておらず
  （ADR-0232 の案 S1 に書き、`xlc-cmfd` に取り込んだ）、コピーを消してよいかはユーザーが決めていない。#0076 と docs/branches.md も
  このコピーの場所を書いていない
- ADR-0232 の検討した案（凍結する部分）は、ブランチ feat/core-cmfd と旧 worktree のパスを指したままで、タグを指す文が ADR に無い
  （`git show 177d63bc5:…` はタグがあるので引ける）
- #0156 の 2026-10-29 の経緯は「別のブランチ feat/core-cmfd がある。中身は確かめていない」のままで、タグと ADR-0232 を指す行が無い
- タグ archive/feat-core-cmfd は注釈なしのタグで、付けた理由を持たない（注釈つきのタグなら理由をタグに残せる）。「消さない」は
  docs/branches.md の文だけが支える
- 閉じた #0076 の問題の節は、worktree に未コミットの実装があるという当時の状態のまま（経緯の最後の行で片付けを書いた）

関係する ADR: なし

## 完了の条件

控えのコピーを残すか消すかをユーザーが決め、#0077 の完了の条件をそれに合わせる。ほかの項目は、直すか直さないかを決めて閉じる。

## 経緯

- 2026-11-06: docs/0076-close-old-cmfd-branch のレビューの指摘から起こした
