---
type: ADR
status: accepted
date: 2026-11-19
requires: []
supersedes: []
superseded-by:
issues: ["#0296", "#0297"]
---

# ADR-0256: docs/ を OKF v0.3 のバンドルにし、design-records check で type と索引を検査する

## 背景

OKF（Open Knowledge Format）は、マークダウンと YAML のフロントマターだけで知識を表す形式で、概念の文書に `type` を求め、
各ディレクトリの `index.md` で内容を一覧にする（仕様の原典は https://github.com/open-knowledge-format/open-knowledge-format の
`SPEC.md`、v0.3、コミット 527e81af、2026-11-19 取得。節の番号は v0.3 のもの）。ユーザーは issue の整理について、xlc の docs/ が
OKF の形である前提で尋ねたが、xlc は OKF を適用していなかった。ユーザーは、仕様の原典を web から取ること、claude-harness にも
適切な処置をすることを求め、OKF を v0.3 に上げると決めた（ユーザー 2026-11-19）。

claude-harness の ADR-0092 は、`docs/index.md` のフロントマターに `okf_version` を持つリポジトリだけを OKF のバンドルとして
`design-records check` で検査すると決めた（lock 3d1a4f85、2026-11-19 の配備。`okf_version` は `"0.2"` と `"0.3"` を通す）。
検査の対象は docs/ の全体で、バンドルのルートは docs/ に固定されている。xlc の docs/ には `type` も `index.md` も無かった
（2026-11-19 の main 9f1a69e72 で、ADR 249 件、issue 291 件、計画 7 件、ほかの .md 4 件）。

- 既存の文書の本文と、ADR と issue の `status` の値を変えない
- `type` は文書の種類を表す値にし、検査を通すためだけの一律の値にしない
- 新しい文書を足すときの手作業を増やさない

## 検討した案

- 案 1: docs/ の全体を OKF v0.3 のバンドルにする。docs/ の .md のすべてに `type` を書き、文書を持つディレクトリごとに
  `index.md` を置いて、`design-records check` に検査させる
- 案 2: 案 1 を v0.2 で宣言する。v0.2 はルートの `index.md` に `okf_version` のほかのキーを許さず、リンクの解き方が違う
- 案 3: docs/ の .md のすべてに、種類を分けない 1 つの `type`（`Document` など）を書く
- 案 4: 索引に載る説明のために、ADR と issue の全件に `description` を今書く

## 決定と理由

案 1 を採る。v0.3 はユーザーの決定（2026-11-19）で、案 2 は採らない。

- バンドルのルートの docs/index.md のフロントマターは `okf_version: "0.3"`、`title`、`description`。`title` と `description` は、
  ほかのバンドルと並べて読む consumer がバンドルを選べるように置く（§8.1）。`license` は置かない。xlc にはライセンスの定めが
  無く、決めるのはユーザーなので、推測で書かない
- `type` の値（OKF は値を登録しない、§4.1）:
  - ADR は `ADR`、docs/issues/ の issue は `Issue`（`design-records check` がこの値を求める）
  - docs/plans/ の計画は `Plan`（ステージの目的、範囲、作業の順序を書く文書で、決定の ADR と分けて引ける）
  - docs/environment.md と、記録の書き方と運用を書く docs/adr/README.md、docs/issues/README.md は `Guide`（手順と決まりを説く
    文書。claude-harness の docs/windows.md と同じ値）。OKF は README を `type` の無い概念でないファイルとして許すが（§3.2）、
    `design-records check` は `type` の無い .md を違反にする（claude-harness の ADR-0092）
  - docs/branches.md は `Reference`（引くための一覧。§4.1 の例の値）
- 案 3 は、`type` で文書を振り分ける consumer（§4.1）に何も伝えないので採らない
- 計画とほかの 4 件には、索引の行が引く `description` を 1 文ずつ書く。ADR と issue には今は書かない（案 4 を採らない）。
  `design-records check` は説明の無い索引の行を違反にも注意にもせず（2026-11-19 に確かめた）、ADR の `description` は
  `design-records` が知らないキーの違反にする。数百件を書く前に要否を決めるため、#0296 に残す
- 索引は `design-records index` で作り、足す。行の題は文書の題（ADR は 1 段目の見出し）。#0012 は題の見出しにリンクを含み、
  `design-records` の索引の行の読み取りが題の中の `](` で切れて「載っていない」の違反になるので、docs/issues/index.md の
  その行だけリンクを外した題を手で書く（題の違いの注意が 2 件残る）
- docs/issues/README.md と docs/adr/README.md の「一覧のファイルは作らない」は、手で書く一覧が各ファイルとずれることを理由に
  していた。索引は `design-records check` が文書と揃っていることを検査するので、索引（`index.md`）は置くと改める
- crates/docs の検査（#0013）は、ADR と issue の `type` を任意のキーとして受け付ける。有無と値は `design-records check` が
  検査するので、重ねて検査しない
- ADR と issue の `status` のキーと値は今のままにする（claude-harness の ADR-0092 は、ADR の `status` を OKF の `status` と別の
  ものとして扱い、キーの改名を別のトピックにした。claude-harness の #261 は 2026-11-19 に open）

### 帰結

- 新しい ADR は `design-records new` が `type: ADR` を書き、索引の行も足す。issue は `design-records new --kind issue` の後に
  フロントマターを手で埋めるとき `type: Issue` を書く（docs/issues/README.md のひな形）。ほかの新しい文書は、`type` を書き、
  `design-records index` で索引に足す。どちらも欠ければ pre-commit の `design-records check --changed` が止める
- OKF の consumer は、`status` の知らない値（ADR の `accepted`、issue の `open` と `closed`）を `stable` でないものとして、今の
  知識に絞った結果から外す（§5.4）。ADR は claude-harness の #261 の改名を待つ。issue は #0297
- 索引の ADR と issue の行は題だけ（#0296）

### 確認

- `design-records check` が xlc の docs/ で違反 0 で通る（注意は、main の 385 件と同じものと、#0012 の題の違いの 2 件）
- crates/docs の `accepts_the_okf_type_key_in_adrs_and_issues` と `treats_okf_index_and_log_files_as_not_records`
- `type` の無い文書と、索引に載っていない文書を足したコミットを、pre-commit の hook が止めることを手で確かめた（2026-11-19。
  「a concept needs a YAML front matter with `type`」と「is not listed in the index」の違反）
