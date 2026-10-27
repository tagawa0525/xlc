# マージしないで残すブランチ

ユーザーの判断で、マージも削除もしないブランチ。理由は参照先の ADR か issue にある。ブランチの整理（`git branch -d` など）で消さない。

| ブランチ                 | 理由の参照先                                                   |
| ------------------------ | -------------------------------------------------------------- |
| `feat/essm`              | ADR-0068（採らなかった ESSM の実装）                           |
| `feat/free-gas-ufg`      | ADR-0070（断面積一定の自由ガスの散乱核の実装）                 |
| `refactor/tally-slots`   | #0078（MOC のタリーの入れ物を固定の 16 個にする案）            |
| `fix/0044-graded-chords` | #0044（弦の求積の改善）                                        |
| `feat/dancoff-potential` | #0007（自動レビューの指摘 10 件が未対応のまま残っている）      |
| `feat/core-cmfd`         | #0076（ADR の下書き。worktree に未コミットの実装の途中がある） |

`feat/core-cmfd` の worktree（`xlc-feat-core-cmfd`）にある未コミットの `crates/numeric/src/cmfd.rs` と `lib.rs` の変更は、消さない。
