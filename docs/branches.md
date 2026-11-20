---
type: Reference
description: ユーザーの判断でマージも削除もしないブランチと、その理由の参照先の一覧。
---

# マージしないで残すブランチ

ユーザーの判断で、マージも削除もしないブランチ。理由は参照先の ADR か issue にある。ブランチの整理（`git branch -d` など）で消さない。

| ブランチ                 | 理由の参照先                                              |
| ------------------------ | --------------------------------------------------------- |
| `feat/essm`              | ADR-0068（採らなかった ESSM の実装）                      |
| `feat/free-gas-ufg`      | ADR-0070（断面積一定の自由ガスの散乱核の実装）            |
| `refactor/tally-slots`   | #0078（MOC のタリーの入れ物を固定の 16 個にする案）       |
| `fix/0044-graded-chords` | #0044（弦の求積の改善）                                   |
| `feat/dancoff-potential` | #0007（自動レビューの指摘 10 件が未対応のまま残っている） |

片付けたブランチ `feat/core-cmfd` の先端（68f762520）は、タグ `archive/feat-core-cmfd` で残す。ADR-0232 の「検討した案」が引く下書き（177d63bc5 の 0119 番の ADR）を
履歴に残すため（#0076）。タグを消さない。

記録と ADR が引く、ブランチの上で rebase する前のコミットは、タグ `record/<ADR か issue>-<内容>` で残す（文書が引くコミットは
HEAD、ブランチ、タグのどれかから辿れることを、crates/docs の検査が確かめる）。タグを消さない。

| タグ                               | 引く文書                                   |
| ---------------------------------- | ------------------------------------------ |
| `record/adr-0157-fast-placement`   | ADR-0157、#0105、examples/build-alignment/ |
| `record/adr-0157-slow-placement`   | ADR-0157、#0105、examples/build-alignment/ |
| `record/adr-0185-pointwise-trials` | ADR-0185 の検討した案の案 L                |
| `record/adr-0225-rod-axis-check`   | ADR-0225、#0213、examples/depletion/       |
| `record/adr-0251-coverage-848-854` | ADR-0251 の検討した案の案 A、#0286         |
| `record/adr-0268-dilute-absorbers` | ADR-0268 の検討した案の案 E、#0292         |
