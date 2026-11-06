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

関係する ADR: なし

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
3. 「比較のツールが書いて git で追跡する JSON」の範囲（どのツールの JSON か。例: `examples/vera/p6/compare-6.json` も
   `path` と `md5` の組を持つが、引くのは git で追跡する `.toml` と、`examples/vera/p6/dvc.yaml` のステージの出力の `.out` だけ）を
   定めていない
4. 記録の md5 と JSON の `inputs` を、`.dvc`、git、`dvc.lock` の outs と突き合わせる検査が無い（#0209 と同じ。JSON は機械で
   読めるので、検査の対象にしやすい）
5. スキル verification-record の「記録に書く試行」の 4 は試行の出力だけを指し、試行でない比較の JSON（`branches/` など）と
   記録が直接引く生成物の `dvc add` の手順はスキルに無い（#0209）
6. ADR-0222 の表題と「背景」は記録と ADR が引く生成物だけを書く（凍結する部分。対象を広げたことは「## 変更」にある）
7. ADR-0222 のフロントマターの `issues` に #0234（対象を広げたきっかけ）と #0238 が無く、#0234 の「関係する ADR:」にも
   ADR-0222 が無い（この issue の「関係する ADR:」も、そのため「なし」）
8. 時機が「記録に書く前」と「JSON だけが引くものは JSON をコミットする前」の 2 つに分かれ、JSON と記録の両方が引く生成物は、
   JSON を先にコミットして記録を書くまでの間 `.dvc` が無くてよい読み方になる（「引く記録、ADR、JSON の最初のコミットの前」の
   1 つにまとめる案がある）
9. `trial-0231/.gitignore` の注は、記録が同じ置き場の `timing-p9.json`（と、それを出す `core-p9.log`）を読むことを書かなくなった
10. ADR-0222 の「## 変更」は、ユーザーの判断を要約で書く（発言は「両方推奨にしたがオマス」（ユーザー 2026-11-06）で、推奨の文は
    発言に無い）。どこまでの範囲（ログ、既にある JSON の分）に同意したかは、発言からは読めない

11. 1 の小さな炉心の履歴は、ADR-0231 と同じブランチ（perf/0156-warm-start-nem）で 2026-11-06 に作った
    `branches/compare-core-small-guess.json` が引く 4 件（`core-small-search-flat-1` 058973f1、`core-small-search-1` adc9e50d、
    `core-small-equilibrium-flat-1` d6762f89、`core-small-equilibrium-1` 0c3828ae）である。2026-11-07（main afd692dda）に、
    adc9e50d は main の `core-small-search-none-1.msgpack.dvc`（50a53f237 で追加）、0c3828ae は `core-small-equilibrium-nested-1.msgpack.dvc`
    （ce4befbd9 で追加）が指し、DVC のキャッシュにある。058973f1 と d6762f89 は main の作業ツリー（`.toml` だけ）にも DVC のキャッシュにも
    無い。main の `core-small-search-1.msgpack.dvc` と `core-small-equilibrium-1.msgpack.dvc` は、JSON が引くものと別の版（4675c7f9、
    cab0f27e）を指す。ADR-0222 の「帰結」の「失われたもの」に載せるかを、棚卸しで決める
12. スキル verification-record の「記録に書く試行」で、記録に書かず 2 で消す試行の出力を比較の JSON が引くときの扱いが無い
13. スキル verification-record の「記録に書く試行」の 4 は `dvc add` を 3 の移す時に置くが、1 で `trials/` のパスのまま記録が
    試行の出力の md5 を書くと、CLAUDE.md の「計算」の `dvc add` の項の時機（記録に書く前）と食い違う（`trials/` の中では
    `dvc add` できないので、両方には従えない）
14. このブランチの分岐点（d8eb8a668）の後に、main は比較の JSON が引く履歴の `.dvc` を足した（`trial-0076`、`trial-0156`、
    `trial-0234`、`trial-0244` の `core-p9.msgpack.dvc`、`branches/core-small-equilibrium-nested-1.msgpack.dvc` など）。ADR-0222 の
    「帰結」の棚卸しの行と、1 の数（main d8eb8a668 で数えたもの）は、棚卸しのときに今の main で数え直す
15. main の `examples/depletion/p9-timing/trial-0156/.gitignore` の注は、隣に `core-p9.msgpack.dvc` があるのに「追跡しない」と書く
    （このブランチは `trial-0231/.gitignore` の注だけを直した）
16. 1 の例の「main の作業ツリーに `.msgpack` が無い（`.toml` だけ）」は main d8eb8a668 の時点のもので、2026-11-07（main afd692dda）
    には `core-small-search-1.msgpack` と `core-small-equilibrium-1.msgpack`（別の版。11）がある。`.toml` だけなのは `-flat-1` の 2 件
17. 11 は `compare-core-small-guess.json` が引く `.msgpack` だけを挙げ、同じ JSON が md5 で引くログ 4 件（`core-small-search-flat-1.log`
    e1794956 など）を挙げていない（2 の境が決まるまで、棚卸しに入れるかは決めていない）
18. 2026-11-07 の経緯の行は 13〜15 を足したことだけを書き、3 と 11 を直したこと（135a705c6）を書かない
19. 15 の `trial-0156/.gitignore` の注は、`trial-0076` と `trial-0244` の `compare-acceleration.json` が `trial-0156/core-p9.msgpack` を
    md5 で引くことも書かない（9 の `trial-0231` と同じ抜け）

## 完了の条件

1 の棚卸しをして、残っている生成物を `dvc add` し、失われたものを ADR-0222 の「帰結」に足す。2 と 3 の境を決める（要るなら
ADR-0222 の「## 変更」）。4 は #0209 で扱う。

## 経緯

- 2026-11-06: ブランチ docs/0222-compare-json-generated のレビューの 1 周目の指摘から起こし、2 周目の指摘（8〜10）と 3 周目の指摘（11、12）を足した
- 2026-11-07: 記録の周の後の例外の周の指摘で 3 と 11 を直し、直さない指摘（13〜15）を足した。例外の 2 周目の指摘（16〜19）を足した
