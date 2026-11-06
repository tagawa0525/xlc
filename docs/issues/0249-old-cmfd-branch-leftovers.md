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
worktree xlc-feat-core-cmfd は 2026-11-06 に消し、先端 4d82b4ce8 にタグ archive/feat-core-cmfd を付けた（docs/branches.md、#0076）。

- 未コミットだった `crates/numeric/src/cmfd.rs` と `lib.rs` の変更は、タグに入らない。残るコピーは履歴の書き換え（#0141）の控え
  `~/.local/share/history-rewrite/xlc-2026-10-21/` の `cmfd.rs.before`（10,077 バイト）、`worktree-2.patch` と `worktree-2.after.patch`（lib.rs の
  `pub mod cmfd;` の 1 行。書き換えの前と後）で、旧 worktree の HEAD（13fddc36e）は同じ場所の `worktree-state` にある。13fddc36e は書き換えの前の hash で、同じ場所の `transaction` の
  `update refs/heads/feat/core-cmfd 4d82b4ce8… 13fddc36e…`（新しい値、古い値）のとおり、書き換えの後は 4d82b4ce8（タグ archive/feat-core-cmfd の先）。どれも 2026-10-21 の時点のもの（消した worktree の最後の中身と同じかは確かめていない）。#0077 の完了の条件は、
  この作業ディレクトリを消す前に控えのパッチの変更がコミットされたかを確かめることを求めるが、変更はそのままはコミットしておらず
  （ADR-0232 の案 S1 に書き、`xlc-cmfd` に取り込んだ）、コピーを消してよいかはユーザーが決めていない。#0076 と docs/branches.md も
  このコピーの場所を書いていない
- ADR-0232 の検討した案（凍結する部分）は、ブランチ feat/core-cmfd と旧 worktree のパスを指したままで、タグを指す文が ADR に無い
  （`git show de5b93436:…` はタグがあるので引ける）
- #0156 の 2026-10-29 の経緯は「別のブランチ feat/core-cmfd がある。中身は確かめていない」のままで、タグと ADR-0232 を指す行が無い
- タグ archive/feat-core-cmfd は注釈なしのタグで、付けた理由を持たない（注釈つきのタグなら理由をタグに残せる）。「消さない」は
  docs/branches.md の文だけが支える
- 2 周目以降の non-blocking: #0077 と閉じた #0076 の経緯から #0249 を指す行が無い。#0076 の閉じる行の「上の 2026-11-06 の項」は
  2026-11-06 の項が 3 つあり、どれか分かりにくい。#0076 と docs/branches.md は下書きを de5b93436 で引くが、ブランチの先端 4d82b4ce8 も
  0119 番の ADR を直している（書式の更新）。docs/branches.md はブランチの一覧の文書で、タグの「消さない」は表の外の段落にある

関係する ADR: なし

## 完了の条件

控えのコピーを残すか消すかをユーザーが決め、#0077 の完了の条件をそれに合わせる。ほかの項目は、直すか直さないかを決めて閉じる。

## 経緯

- 2026-11-06: docs/0076-close-old-cmfd-branch のレビューの指摘から起こした
