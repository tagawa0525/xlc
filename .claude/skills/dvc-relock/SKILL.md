---
name: dvc-relock
description: DVC のステージ（examples/*/dvc.yaml）を含むブランチのマージの前後の手順。rebase のあとに dvc.lock を整え、マージのあとに main の出力を置き換える。/git-merge の rebase の段の後とマージの後に使う。
allowed-tools:
  - Bash(nix develop -c dvc checkout examples/*)
  - Bash(nix develop -c dvc checkout)
  - Bash(nix develop -c dvc checkout dvc.yaml)
  - Bash(nix develop -c dvc status)
  - Bash(nix develop -c dvc status --json*)
  - Bash(nix develop -c dvc status --with-deps examples/*)
  - Bash(nix develop -c scripts/check-dvc-lock.sh*)
  - Bash(nix develop -c scripts/compare-dvc-lock.py*)
  - Bash(nix develop -c scripts/count-changed-stages.py*)
  - Bash(scripts/check-dvc-lock.sh*)
  - Bash(scripts/check-references-shell.sh*)
  - Bash(head -1 examples/*/dvc.lock)
  - Bash(head -1 dvc.lock)
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

- 記録（verification.md、ADR）の数値には、実行したコミットと分岐点（main のコミット）を書く
- rebase で deps が変わっても、重いステージは再実行しない。分岐点にそこまでのブランチのコミットを当て直して復元する
- 数秒で終わるステージだけは再実行してよい

手順は前半（rebase のあと。ブランチの worktree）と後半（マージのあと。main の worktree）に分かれる。呼ばれたときの
worktree と、マージ済みかどうかで、どちらかを選ぶ。

/git-merge のレビューの指摘を直して再び rebase したときは、前半をもう一度行う（分岐点が変わるので、下の 5 も直す）。

`dvc` は `default` の devShell にある。リポジトリ（worktree）の直下で `nix develop -c dvc ...` と書く（direnv で
入っていれば `nix develop -c` は要らない）。`references` の devShell の中では実行しない（ADR-0130）。
`examples/<問題>/dvc.yaml` の指定は、リポジトリの直下からのパスで書く。

## 前半: rebase のあと、チェックの前（ブランチの worktree で）

1. main 側で変わった上流の出力を戻す:

   ```bash
   nix develop -c dvc checkout
   ```

   `dvc status` と `dvc repro` の前に必ず行う。しないと、古い上流の出力から再実行して、誤った md5 を `dvc.lock` に
   書く。

   引数なしの `dvc checkout` は、main の上でも終了コード 255 で終わる。`dvc.lock` に hash がまだ無いステージの出力（まだ
   流していない examples/depletion/p2 の本番と examples/vera/endfb70/openmc の OpenMC。`No file hash info found`）が
   戻らないためで、戻らなかった出力は `Checkout failed for following targets:` の後に並ぶ（2026-10-23、43a0af25d で
   p2 の 21 件と endfb70 の 5 件）。並んだ出力ごとに、そのステージの `dvc.lock`（`examples/depletion/p2/dvc.lock`、
   `examples/vera/endfb70/openmc/dvc.lock` など）を `git diff main -- <その dvc.lock>` で見る。差が無ければ main でも戻らない
   既知の出力なので、続ける。ブランチで変えた `dvc.lock` の出力が並んだら、そこで止まり、原因を調べる。`references` の symlink が無いときも止まる。
   `references` が main の `references/` を指すかは `realpath references` と main の `references/` の `realpath`
   を比べる（symlink なので文字列の比較は合わない）。出力がキャッシュにないときは、手で `dvc repro` で埋めない。
   上流のステージを所有するブランチか main の記録を確かめてから、ユーザーに聞く。

2. 実行ファイルと Python の正規形を作り直してから、deps の変わったステージを挙げる（rebase でリポジトリの直下の `dvc.lock`
   か `python-ast/` が衝突したら、main の版を取って rebase を終え、ここで `build` と `python-ast` を実行して書き直す。手で行を
   選んで混ぜない）:

   ```bash
   nix develop -c dvc repro -s build
   nix develop -c dvc repro -s dvc.yaml:python-ast
   nix develop -c dvc status
   ```

   例のステージの deps の .py は、`python-ast`（リポジトリの直下の `dvc.yaml`。数秒）が書く正規形で（ADR-0210）、`build` と
   同じく、先に実行しないとコードの変わった .py の下流が changed に出ない。docstring とコメントだけの変更なら正規形は変わらず、
   下流は changed にならない。`python-ast` を実行したら、リポジトリの直下の `dvc.lock` と `python-ast/` も 4 でコミットする

   計算のステージの deps は `bin/` の実行ファイルで（ADR-0159）、1 の `dvc checkout` は `dvc.lock` が指す古い実行ファイルを
   戻す。先に `build`（リポジトリの直下の `dvc.yaml`。ビルドは数秒〜十数秒で、手元でよい）を実行しないと、ソースが変わった
   ステージも changed に出ず、3 の `-s` が古い実行ファイルで走る。`build` が「didn't change」なら、ソースは変わっていない。
   作り直した実行ファイルが `dvc.lock` と同じバイトなら、計算のステージは changed にならない。`build` を実行したら、
   リポジトリの直下の `dvc.lock` も 4 でコミットする

   対象は `changed deps` のステージだけ。`changed outs` や `.dvc` ファイルの変更が混ざるので、`--json` のキーを
   そのままステージの一覧にしない（`deps` の項目を持つステージを読む）。`changed outs` だけのステージは、出力が手元で変わっている
   ので実行せず、原因を調べる。ブランチが触れた例だけを見るなら、`dvc status --with-deps examples/<問題>/dvc.yaml`（`--with-deps` がないと、上流の `build` を示さない）。
   `dvc status` が失敗したら（`dvc.lock` が読めない、キャッシュの参照が切れた、など）止める。

   changed のステージの数は、`dvc status` の文字列を grep で数えず（タブで始まる `changed deps:` の行を拾って誤る。2026-10-28 に
   468 を 942 と数えた）、`nix develop -c scripts/count-changed-stages.py [examples/<問題>/dvc.yaml...]` で数える（yaml ごとの
   changed のステージ数と、deps、outs、command、always changed の内訳を出す）。

3. ステージごとに、再実行するかを決める。目安（所要時間は ADR-0130 の「確認」の値）:

   - 再実行してよい: 数秒で終わる lattice の問題（lattice-1B）と、その出力だけを入力にする小さなステージ（`compare-*`）
   - 実行しない: xs（約 9 分）、OpenMC（40〜60 分）、ライブラリ（核種ごと、逐次で約 30 分）、数分かかる lattice の問題
     （Slurm で流す重さ。CLAUDE.md の「計算」）
     ライブラリの核種ごとのステージを作り直すときは、relock の中でなく、Slurm で `examples/stage_runner.py` を流す（ADR-0193）
   - `references-shell` は `always_changed` なので `dvc status` に常に出る。対象から外す。ただし出力の `references-shell.txt` は
     `dvc repro` まで書き直されず、devShell が変わっても `dvc status` に出ないので、中身を `scripts/check-references-shell.sh` で
     確かめる（コミットのたびに pre-commit の hook が同じものを実行するが、rebase で main の flake の変更と混ざった版は hook を
     通っていない）。食い違えば止め、出力の文言の修正方法に従う（ADR-0212 の (a)）
   - `build` と `python-ast` は 2 で走らせ終えている

   数秒で終わると分かっているステージだけ実行する。迷えば実行しない。実行するときは Slurm を通さず直接、ステージを指定する:

   ```bash
   nix develop -c dvc repro -s <ステージ>
   ```

   CLAUDE.md の「計算」は `-s` を /dvc-relock だけに限る（ADR-0184）。rebase のあとの relock は ADR-0130 の「数秒で終わるステージだけ再実行してよい」の手段として、重いステージを
   巻き込まないために `-s` を使う。`-s` は上流を実行しない。上流が `dvc status` で changed のステージには使わない
   （古い出力から作ることになる）。lattice-1B の上流の xs が changed のときは、xs を実行しないので lattice-1B も
   実行しない（xs を含む約 9 分は、数秒のステージにならない）。実行したステージの出力を deps にする下流のステージ（lattice-1B なら `compare-1B`）も
   changed になるので、同じように数秒なら続けて実行し、重いなら実行しない。

4. 実行したら、`dvc.lock` の差分だけをコミットする（出力の追跡ファイルが変わっていればそれも含める）:

   ```bash
   git add examples/<問題>/dvc.lock dvc.lock python-ast   # dvc.lock と python-ast/ はリポジトリの直下の build と python-ast の分。git status で、dvc が書き換えたほかの追跡ファイル（.gitignore など）も足す。rebase の hook が付け替えた記録は足さない（5 で別のコミットにする）
   git commit -F <メッセージのファイル>
   ```

   件名は `chore: relock <ステージ> after the rebase`。本文に、なぜ再実行したか（deps の何が変わったか）、分岐点
   （`git merge-base HEAD main`）を書き、`Refs: ADR-0130` と、/git-commit の規則の末尾（Co-Authored-By など）を
   添える。実行したコミットはこのコミット自身なので、記録には、このコミットか、その結果を出したコード（ツールを含む）を
   含む最初のコミットを書く。

5. 記録（verification.md、ADR）の hash と分岐点を直す。ステージで分ける:

   - 再実行したステージ: 実行したコミットを 4 のコミット（またはその結果を出したコードを含む最初のコミット）に、分岐点を今の `git merge-base HEAD main` にする
   - 走らせなかったステージ: 分岐点は、結果を出したときの main のコミットのまま変えない（今の merge-base に変えると、
     復元が合わなくなる）。実行したコミットが rebase で付け替わっていたら、hash だけを付け替える。記録に実行した
     コミットと分岐点が書かれていなければ、書く（ADR-0130 の「確認」の形）

   hash の付け替えは、手で書き換えず、`commit-refs remap` で行う。rebase のあとは、`remap.commitRefs` が true なら
   post-rewrite hook が自動で呼ぶ（rebase の途中の amend では呼ばず、終わりに一度）（付け替えはコミットされないので、自分でコミットする）。hook が失敗したときは、
   `$(git rev-parse --absolute-git-dir)/commit-refs-remap.input` の対応表で `commit-refs remap <そのファイル>` を呼び直す
   （詳細は verification-record の「rebase のあと」の節）。手で対応表を作るのは、この対応表が無いときだけ。rebase の前の tip
   （`git reflog` の `ORIG_HEAD` など）と今の HEAD を `git range-diff` で対応させ、中身を確かめてから渡す（全ファイルを
   書き換えるので、推測した対応を渡さない）。測り直しはしない（記録の hash は再現のため）。再実行したステージの新しいコミットと分岐点の値は、rebase の対応表に載らないので、手で書く（書いた hash は `git rev-parse` で確かめる）。付け替えのコミットは `docs:` で分ける。

6. 作り直した `dvc.lock` を確かめる（`build` か `python-ast` を実行したら、リポジトリの直下の `dvc.lock` も同じように）。`head -1 examples/<問題>/dvc.lock` が `schema: '2.0'` を出すこと（dvc は、この行の無い lock を
   読めず、書くこともない。欠けたら、手で書き換えた箇所を探し、dvc が書いた版に戻す。戻したら、4 のコミットを `git commit --fixup` で直して `git rebase -i --autosquash` でまとめ、5 の付け替えのコミットは rebase のあとに積み直す）。続けて
   `nix develop -c scripts/check-dvc-lock.sh examples/<問題>/dvc.lock`（direnv で入っていれば `nix develop -c` は要らない）が終了コード 0 で終わること（dvc が lock を読めるかの検査は、このツールが
   唯一の定義。pre-commit の hook `dvc-lock` も同じツールを呼ぶ。issue #0090）。

出力を作り直さず `dvc commit` で lock を付け直してよいのは、ADR-0212 の (b) の条件（重いステージで、変更の親で changed でなく、
変わった deps が値の同じ TOML だけで、親の版のバイトが出力に入っていない）を満たすステージと、ADR が決めた移行だけである。実行ファイルと
`references-shell.txt` の差では付け直さない。
付け直したとき（deps の置き場や正規形を替えた ADR-0210 の移行など）は、deps の hash
だけが変わったことを `nix develop -c scripts/compare-dvc-lock.py <付け直す前のコミット> HEAD` で確かめる（2 つのコミットの間で変わった
すべての `dvc.lock` について、ステージの集合と、各ステージの `cmd`、`outs`、`params` を比べ、違いはファイルとステージと欄で出して
終了コード 1 で終わる）。ステージの集合の差も違いとして出る（python-ast のステージの追加のように増えて当然のときは、出力を読んで確かめる）。
ステージを再実行した lock は `outs` が変わるので、この確認は付け直しだけのコミットの範囲に使う。

この後は /git-merge の「チェック」に進む。

### マージの前に: main の作業ツリーに残る古い出力

ブランチが DVC の outs を `.gitignore` から外したとき（`cache: false` に変えて git で追跡するようにしたときなど）、main の作業ツリーには
以前の `dvc checkout` / `dvc repro` が置いた古い出力が、まだ無視されたファイルとして残っている。マージでブランチ側の追跡されたファイルが
同じ場所に入ると、main では未追跡のファイルが上書きされる形になり、pre-merge の hook（未追跡の衝突）が止まる（2026-10-15 に
p5 と p3 の `compare-*.out` で 2 回起きた）。マージの前に、ブランチで `.gitignore` から外れた出力（`git diff main...HEAD -- '*.gitignore'`）が
main の作業ツリーに残っていないかを `git -C <main の worktree> status --short <パス>` と `ls` で見て、あれば main 側で消す。消すのは
DVC のキャッシュから戻せる古い出力だけで、マージのあとに後半の `dvc checkout` が置き直す。追跡されている版に差があるファイルは消さず、原因を調べる。

## 後半: マージのあと（main の worktree で）

ブランチの worktree ではなく、main の worktree の直下で行う（出力を置き換える先が main）。

main の出力は、DVC を入れる前のもの（`vera.msgpack` など）か、前のマージの時点のものが残っている。dvc.lock の出力に
置き換える。`bin/` の実行ファイル（リポジトリの直下の `dvc.yaml` の build。ADR-0159）も戻す:

```bash
nix develop -c dvc checkout dvc.yaml
nix develop -c dvc checkout examples/<問題>/dvc.yaml
nix develop -c dvc status --with-deps examples/<問題>/dvc.yaml
```

changed のステージの数は、前半の 2 と同じく `nix develop -c scripts/count-changed-stages.py examples/<問題>/dvc.yaml` で数える。

`-P` ではなく、ブランチが `dvc.yaml` を足した、または変えた例ごとに行う。`dvc status` の changed が、前半で走らせ
なかったステージ（deps が `dvc.lock` と合わないステージ）と、その下流のステージだけであることを確かめる。`references-shell`
（`always_changed`）は常に出るので数えない。ブランチが `build` を実行せずにコードを変えたなら、`build` とその下流も出る
（ADR-0159 の帰結。`dvc repro` で作り直すまで、そのままでよい）。ほかのステージが changed なら、
出力が足りない、またはキャッシュにない。`dvc checkout` が失敗したら止めて、原因を調べる（推測で `dvc repro` を実行しない）。

走らせなかったステージの changed は、そのままでよい。`dvc.lock` は、どの deps から出力を作ったかの記録であり、結果が
要るときに HEAD で `dvc repro` すれば再実行する（ADR-0130）。
