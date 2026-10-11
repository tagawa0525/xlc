---
name: dvc-relock
description: DVC の段階（examples/*/dvc.yaml）を含むブランチのマージの前後の手順。rebase のあとに dvc.lock を整え、マージのあとに main の出力を置き換える。/git-merge の rebase の段の後とマージの後に使う。
argument-hint: [before|after]
allowed-tools:
  - Bash(nix develop -c dvc checkout examples/*)
  - Bash(nix develop -c dvc checkout)
  - Bash(nix develop -c dvc status examples/*)
  - Bash(nix develop -c dvc status --json*)
  - Bash(nix develop -c dvc repro -s *)
  - Bash(realpath*)
  - Bash(git status*)
  - Bash(git diff*)
  - Bash(git log*)
  - Bash(git show*)
  - Bash(git merge-base*)
  - Bash(git add*)
  - Bash(git commit*)
---

# dvc.lock の整理（rebase のあとと、マージのあと）

`dvc.yaml` のある例（examples/m0、examples/vera など）を含むブランチの、マージの前後の手順。決定は
ADR-0130 の「帰結」、規則は CLAUDE.md の「計算」にある。

- 記録（verification.md、ADR）の数値には、走らせたコミットと分岐点（main のコミット）を書く
- rebase で deps が変わっても、重い段階は走らせ直さない。分岐点にそこまでのブランチのコミットを当て直して復元する
- 数秒で終わる段階だけは走らせ直してよい

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

| 段階                         | 所要時間       | 判断             |
| ---------------------------- | -------------- | ---------------- |
| lattice の 1 問題            | 数秒〜数分     | 走らせ直してよい |
| xs                           | 約 9 分        | 走らせない       |
| OpenMC                       | 40〜60 分      | 走らせない       |
| ライブラリ（NJOY、核種ごと） | 逐次で約 30 分 | 走らせない       |

   数秒で終わると分かっている段階だけ走らせる。迷えば走らせない。走らせるときは Slurm を通さず直接、段階を指定する:

   ```bash
   nix develop -c dvc repro -s <段階>
   ```

   `-s` は指定した段階だけを走らせ、上流の段階を巻き込まない（上流が古いときは、古い出力から作ることになるので、
   上流が `dvc status` で changed の段階には使わない）。

4. 走らせたら、`dvc.lock` の差分だけをコミットする（出力の追跡ファイルが変わっていればそれも含める）:

   ```bash
   git add examples/<問題>/dvc.lock
   git commit   # 件名: chore: relock <段階> after the rebase
   ```

   本文に、なぜ走らせ直したか（deps の何が変わったか）、走らせたコミット、分岐点（`git merge-base HEAD main`）を書く。
   `Refs: ADR-0130` を添える。

5. 走らせなかった段階は、記録に、走らせたコミットと分岐点が書かれていることを確かめる。分岐点は
   `git merge-base HEAD main` で、rebase の直後の main のコミットになる。記録が rebase 前の hash を指していれば、
   指している先が今のブランチに残っているかを `git log` で見て、書き直す（記録の hash は再現のためで、checkout すれば
   当時の結果が出ること）。

   - verification.md の結果ごとの「測ったときの xlc のコミットと日付」
   - ADR の「確認」の、走らせたコミットと分岐点

   測り直さずに、道具や文書の変更で付け替える。付け替えたコミットは `chore:` や `docs:` で分ける。

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
なかった段階（deps が `dvc.lock` と合わない段階）だけであることを確かめる。ほかの段階が changed なら、
出力が足りない、またはキャッシュにない。`dvc checkout` が失敗したら止めて、原因を調べる（推測で `dvc repro` を走らせない）。

走らせなかった段階の changed は、そのままでよい。`dvc.lock` は、どの deps から出力を作ったかの記録であり、結果が
要るときに HEAD で `dvc repro` すれば走り直す（ADR-0130）。
