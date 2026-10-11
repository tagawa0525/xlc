---
name: dvc-relock
description: DVC の段階（examples/*/dvc.yaml）を含むブランチのマージの前後の手順。rebase のあとに dvc.lock を整え、マージのあとに main の出力を置き換える。/git-merge の rebase の段の後とマージの後に使う。
allowed-tools:
  - Bash(nix develop -c dvc checkout examples/*)
  - Bash(nix develop -c dvc checkout)
  - Bash(nix develop -c dvc status)
  - Bash(nix develop -c dvc status examples/*)
  - Bash(nix develop -c dvc status --json*)
  - Bash(nix develop -c dvc repro -s *)
  - Bash(realpath*)
  - Bash(git reflog*)
  - Bash(git range-diff*)
  - Bash(git status*)
  - Bash(git diff*)
  - Bash(git log*)
  - Bash(git show*)
  - Bash(git merge-base*)
  - Bash(commit-refs remap*)
  - Bash(git add*)
  - Bash(git commit*)
---

# dvc.lock の整理（rebase のあとと、マージのあと）

`dvc.yaml` のある例（examples/m0、examples/vera など）を含むブランチの、マージの前後の手順。決定は
ADR-0130 の「帰結」、規則は CLAUDE.md の「計算」にある。

- 記録（verification.md、ADR）の数値には、走らせたコミットと分岐点（main のコミット）を書く
- rebase で deps が変わっても、重い段階は走らせ直さない。分岐点にそこまでのブランチのコミットを当て直して復元する
- 数秒で終わる段階だけは走らせ直してよい

手順は前半（rebase のあと。ブランチの worktree）と後半（マージのあと。main の worktree）に分かれる。呼ばれたときの
worktree と、マージ済みかどうかで、どちらかを選ぶ。

/git-merge のレビューの指摘を直して再び rebase したときは、前半をもう一度行う（分岐点が変わるので、下の 5 も直す）。

`dvc` は `default` の devShell にある。リポジトリ（worktree）の直下で `nix develop -c dvc ...` と書く（direnv で
入っていれば `nix develop -c` は要らない）。`references` の devShell の中では走らせない（ADR-0130）。
`examples/<問題>/dvc.yaml` の指定は、リポジトリの直下からのパスで書く。

## 前半: rebase のあと、チェックの前（ブランチの worktree で）

1. main 側で変わった上流の出力を戻す:

   ```bash
   nix develop -c dvc checkout
   ```

   `dvc status` と `dvc repro` の前に必ず行う。しないと、古い上流の出力から走らせ直して、誤った md5 を `dvc.lock` に
   書く。失敗したら（キャッシュにない出力、`references` の symlink がない、など）そこで止まり、原因を調べる。
   `references` が main の `references/` を指すかは `realpath references` と main の `references/` の `realpath`
   を比べる（symlink なので文字列の比較は合わない）。出力がキャッシュにないときは、手で `dvc repro` で埋めない。
   上流の段階を所有するブランチか main の記録を確かめてから、ユーザーに聞く。

2. deps の変わった段階を挙げる:

   ```bash
   nix develop -c dvc status
   ```

   対象は `changed deps` の段階だけ。`changed outs` や `.dvc` ファイルの変更が混ざるので、`--json` のキーを
   そのまま段階の一覧にしない（`deps` の項目を持つ段階を読む）。`changed outs` だけの段階は、出力が手元で変わっている
   ので走らせず、原因を調べる。ブランチが触れた例だけを見るなら、`dvc status examples/<問題>/dvc.yaml`。
   `dvc status` が失敗したら（`dvc.lock` が読めない、キャッシュの参照が切れた、など）止める。

3. 段階ごとに、走らせ直すかを決める。目安（所要時間は ADR-0130 の「確認」の値）:

   - 走らせ直してよい: 数秒で終わる lattice の問題（lattice-1B）と、その出力だけを入力にする小さな段階（`compare-*`）
   - 走らせない: xs（約 9 分）、OpenMC（40〜60 分）、ライブラリ（核種ごと、逐次で約 30 分）、数分かかる lattice の問題
     （Slurm で流す重さ。CLAUDE.md の「計算」）
   - `references-shell` は `always_changed` なので `dvc status` に常に出る。対象から外す

   数秒で終わると分かっている段階だけ走らせる。迷えば走らせない。走らせるときは Slurm を通さず直接、段階を指定する:

   ```bash
   nix develop -c dvc repro -s <段階>
   ```

   CLAUDE.md の「計算」は手で段階を選ばないと定めるが、rebase のあとの relock は ADR-0130 の「数秒で終わる段階だけ走らせ直してよい」の手段として、重い段階を
   巻き込まないために `-s` を使う。`-s` は上流を走らせない。上流が `dvc status` で changed の段階には使わない
   （古い出力から作ることになる）。lattice-1B の上流の xs が changed のときは、xs を走らせないので lattice-1B も
   走らせない（xs を含む約 9 分は、数秒の段階にならない）。走らせた段階の出力を deps にする下流の段階（lattice-1B なら `compare-1B`）も
   changed になるので、同じように数秒なら続けて走らせ、重いなら走らせない。

4. 走らせたら、`dvc.lock` の差分だけをコミットする（出力の追跡ファイルが変わっていればそれも含める）:

   ```bash
   git add examples/<問題>/dvc.lock   # git status で、dvc が書き換えたほかの追跡ファイル（.gitignore など）も足す。rebase の hook が付け替えた記録は足さない（5 で別のコミットにする）
   git commit -F <メッセージのファイル>
   ```

   件名は `chore: relock <段階> after the rebase`。本文に、なぜ走らせ直したか（deps の何が変わったか）、分岐点
   （`git merge-base HEAD main`）を書き、`Refs: ADR-0130` と、/git-commit の規則の末尾（Co-Authored-By など）を
   添える。走らせたコミットはこのコミット自身なので、記録には、このコミットか、その結果を出したコード（道具を含む）を
   含む最初のコミットを書く。

5. 記録（verification.md、ADR）の hash と分岐点を直す。段階で分ける:

   - 走らせ直した段階: 走らせたコミットを 4 のコミット（またはその結果を出したコードを含む最初のコミット）に、分岐点を今の `git merge-base HEAD main` にする
   - 走らせなかった段階: 分岐点は、結果を出したときの main のコミットのまま変えない（今の merge-base に変えると、
     復元が合わなくなる）。走らせたコミットが rebase で付け替わっていたら、hash だけを付け替える。記録に走らせた
     コミットと分岐点が書かれていなければ、書く（ADR-0130 の「確認」の形）

   hash の付け替えは、手で書き換えず、`commit-refs remap` で行う。rebase と amend では、`remap.commitRefs` が true なら
   post-rewrite hook が自動で呼ぶ（付け替えはコミットされないので、自分でコミットする）。hook が失敗したときは、
   `$(git rev-parse --absolute-git-dir)/commit-refs-remap.input` の対応表で `commit-refs remap <そのファイル>` を呼び直す
   （詳細は verification-record の「rebase のあと」の節）。手で対応表を作るのは、この対応表が無いときだけ。rebase の前の tip
   （`git reflog` の `ORIG_HEAD` など）と今の HEAD を `git range-diff` で対応させ、中身を確かめてから渡す（全ファイルを
   書き換えるので、推測した対応を渡さない）。測り直しはしない（記録の hash は再現のため）。付け替えのコミットは `docs:` で分ける。

この後は /git-merge の「チェック」に進む。

## 後半: マージのあと（main の worktree で）

ブランチの worktree ではなく、main の worktree の直下で行う（出力を置き換える先が main）。

main の出力は、DVC を入れる前のもの（`vera.msgpack` など）か、前のマージの時点のものが残っている。dvc.lock の出力に
置き換える:

```bash
nix develop -c dvc checkout examples/<問題>/dvc.yaml
nix develop -c dvc status examples/<問題>/dvc.yaml
```

`-P` ではなく、ブランチが `dvc.yaml` を足した、または変えた例ごとに行う。`dvc status` の changed が、前半で走らせ
なかった段階（deps が `dvc.lock` と合わない段階）と、その下流の段階だけであることを確かめる。`references-shell`
（`always_changed`）は常に出るので数えない。ほかの段階が changed なら、
出力が足りない、またはキャッシュにない。`dvc checkout` が失敗したら止めて、原因を調べる（推測で `dvc repro` を走らせない）。

走らせなかった段階の changed は、そのままでよい。`dvc.lock` は、どの deps から出力を作ったかの記録であり、結果が
要るときに HEAD で `dvc repro` すれば走り直す（ADR-0130）。
