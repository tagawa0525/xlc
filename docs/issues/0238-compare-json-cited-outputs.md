---
status: open
labels: [dvc, verification, review]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# 比較の JSON が md5 で引く生成物の棚卸しと ADR-0222 の対象の残り

## 問題

ADR-0222 の対象に、比較のツールが書いて git で追跡する JSON が md5 で引く生成物を足した（ADR-0222 の「## 変更」の
2026-11-06）。ブランチ docs/0222-compare-json-generated のローカルのレビュー（/code-review high）で、再現を書けず直さなかった
指摘（/git-merge の「直す指摘の基準」の non-blocking）を残す。

関係する ADR: なし（ADR-0222 のフロントマターの `issues` には足していない。下の 7）

1. 対象を広げる前に作った比較の JSON が md5 で引く生成物は棚卸ししておらず、`.dvc` が無い。git で追跡する `.json` の
   `path` と `md5` の組のうち、git で追跡せず、`.dvc` が無く、`dvc.lock` の outs に無く、`references/` の下でもないパスを
   粗く数えると 47 ファイルで 233 件（2026-11-06、main d8eb8a668。JSON の相対パスの基準を解いていないので、ステージの出力を
   含む上限で、確かめていない）。例: `examples/depletion/branches/compare-core-small-guess.json` が引く
   `core-small-search-flat-1.msgpack`（md5 058973f1）ほかの小さな炉心の履歴は、main の作業ツリーに `.msgpack` が無い
   （`.toml` だけ）。ADR-0222 の「帰結」の「失われたもの」にも載っていない
2. 比較の JSON は `.log` と `.out` も md5 で引く（例: `compare-core-small-guess.json` の `core-small-search-flat-1.log`
   e1794956、`compare-core-2A-boron.json` の `core-2A-boron-uniform.out`）。ADR-0222 の対象の「コマンドの出力先としてだけ名前が
   出るログと k の出力は追わない」との境が文から読みにくい（JSON が引くものは「出力先としてだけ」ではない読み方になる）。
   `trial-0231/core-p9.log` は git で追跡している
3. 「比較のツールが書いて git で追跡する JSON」の範囲（どのツールの JSON か。例: `examples/vera/p6/compare-6.json` は
   追跡しない `.toml` と `.out` を md5 で引く）を定めていない
4. 記録の md5 と JSON の `inputs` を、`.dvc`、git、`dvc.lock` の outs と突き合わせる検査が無い（#0209 と同じ。JSON は機械で
   読めるので、検査の対象にしやすい）
5. スキル verification-record の「記録に書く試行」の 4 は試行の出力だけを指し、試行でない比較の JSON（`branches/` など）と
   記録が直接引く生成物の `dvc add` の手順はスキルに無い（#0209）
6. ADR-0222 の表題と「背景」は記録と ADR が引く生成物だけを書く（凍結する部分。対象を広げたことは「## 変更」にある）
7. ADR-0222 のフロントマターの `issues` に #0234（対象を広げたきっかけ）が無く、#0234 の「関係する ADR:」にも ADR-0222 が無い

## 完了の条件

1 の棚卸しをして、残っている生成物を `dvc add` し、失われたものを ADR-0222 の「帰結」に足す。2 と 3 の境を決める（要るなら
ADR-0222 の「## 変更」）。4 は #0209 で扱う。

## 経緯

- 2026-11-06: ブランチ docs/0222-compare-json-generated のレビューの 1 周目の指摘から起こした
