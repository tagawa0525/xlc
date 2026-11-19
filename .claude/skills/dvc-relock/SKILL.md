---
name: dvc-relock
description: DVC のステージ（examples/*/dvc.yaml）を含むブランチのマージの前後の手順。rebase のあとに dvc.lock を整え、マージが .gitignore から外れた古い出力で止まったら main の出力を消してよいかを決め、マージのあとに main の出力を置き換える。/git-merge の rebase の段の後、マージの hook がその一覧で止まったとき、マージの後に使う。
allowed-tools:
  - Bash(nix develop -c scripts/dvc-relock.py*)
  - Bash(nix develop -c scripts/dvc-checkout.py*)
  - Bash(nix develop -c dvc status)
  - Bash(nix develop -c dvc status --json*)
  - Bash(nix develop -c dvc status --with-deps examples/*)
  - Bash(nix develop -c dvc status examples/*)
  - Bash(nix develop -c scripts/check-dvc-lock.sh*)
  - Bash(nix develop -c scripts/compare-dvc-lock.py*)
  - Bash(nix develop -c scripts/count-changed-stages.py*)
  - Bash(scripts/check-dvc-lock.sh*)
  - Bash(scripts/check-references-shell.sh*)
  - Bash(head -1 examples/*/dvc.lock)
  - Bash(head -1 dvc.lock)
  - Bash(nix develop -c scripts/dvc-repro.py -s *)
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
ADR-0130 の「帰結」と ADR-0258 の「rebase の後の付け直し」、規則はスキル compute にある。

- 記録（verification.md、ADR）の数値には、実行したコミットと分岐点（main のコミット）を書く
- rebase で deps が変わったことを理由には再実行しない。分岐点にそこまでのブランチのコミットを当て直して復元する
- 再実行するかは、その結果が要るかと時間の見積もりで、流す側が決める

手順は前半（rebase のあと。ブランチの worktree）と後半（マージのあと。main の worktree）に分かれる。呼ばれたときの
worktree と、マージ済みかどうかで、どちらかを選ぶ。マージが pre-merge-commit の hook の「.gitignore から外した生成物」の
一覧で止まったときは、その間の「マージが止まったら」（main の worktree）を行う。

/git-merge のレビューの指摘を直して再び rebase したときは、前半をもう一度行う（分岐点が変わるので、下の 4 も直す）。

`dvc` は `default` の devShell にある。リポジトリ（worktree）の直下で `nix develop -c dvc ...` と書く（direnv で
入っていれば `nix develop -c` は要らない）。`references` の devShell の中では実行しない（ADR-0130）。
`examples/<問題>/dvc.yaml` の指定は、リポジトリの直下からのパスで書く。

## 前半: rebase のあと、チェックの前（ブランチの worktree で）

rebase でリポジトリの直下の `dvc.lock`、`python-ast/`、`binary-forms/`、`toml-values/` が衝突したら、main の版を取って
rebase を終える（手で行を選んで混ぜない。1 のコマンドが書き直す）。

1. 付け直しのコマンドを実行する（何をどの順に行い、どこで止まるかは、スクリプトの先頭にある）:

   ```bash
   nix develop -c scripts/dvc-relock.py [--run <yaml>:<ステージ>]... [--allow-reference <yaml>:<ステージ>]... [examples/<問題>/dvc.yaml...]
   ```

   - `--run` には、結果が要るステージだけを挙げる。挙げなければ、計算のステージは実行しない。挙げたステージの出力を deps に
     する下流（lattice-1B なら `compare-1B`）も changed になるので、要るなら上流の後に挙げる
   - 参照解のステージ（cmd が references の devShell で包まれたもの）は、`--allow-reference` にも挙げる。重いステージを挙げる
     なら、コマンドごと Slurm に投げる（スキル compute）。ライブラリの核種ごとのステージは、ここでなく Slurm で
     `examples/stage_runner.py` を流す（ADR-0193）
   - `examples/<問題>/dvc.yaml` を渡すと、`dvc status` をその例と上流に限る（ブランチが触れた例だけを見るとき）
   - `CHANGED:` は残した changed のステージで、そのままでよい（ADR-0130）

2. 終了コード 1 なら、`STOP:` の行の場所で原因を調べ、直したら 1 をやり直す。手で `dvc repro` して埋めない:

   - 戻らない出力（`MISSING:`）: 上流のステージを所有するブランチか main の記録を確かめてから、ユーザーに聞く。`references` が
     main の `references/` を指すかは `realpath references` と main の `references/` の `realpath` を比べる（symlink なので文字列の
     比較は合わない）。dvc が「Can't remove the following unsaved files」で止まったときは、並んだファイルが `.dvc` で追う生成物
     （ADR-0222）で、このブランチで作り直したものなら、`dvc add <ファイル>` で `.dvc` を付け直してコミットしてからやり直す。
     `--force` は作り直した版を失うので付けない
   - `changed outs` だけのステージ: 出力が手元で変わっている。実行せず、何が変えたかを調べる
   - `references-shell.txt` の食い違い: 出力の文言の修正方法に従う（ADR-0212 の (a)）
   - `dvc status` の失敗: dvc の文（`dvc.lock` が読めない、キャッシュの参照が切れた、など）から原因を調べる
   - 読めない lock（始めの確認。作業ツリーは書き換えていない）: 手で書き換えた箇所を探し、dvc が書いた版に戻す。戻したら、その
     lock を入れたコミットを `git commit --fixup` で直して `GIT_SEQUENCE_EDITOR=true git rebase -i --autosquash main` でまとめ、4 の
     付け替えのコミットは rebase のあとに積み直す。終わりの確認で止まったら、dvc が書いた lock なので、コミットせずに原因を調べる
   - `--run` のステージ: `scripts/dvc-repro.py -s` の `BLOCKED:` と `STOP:` の行に従う。前に挙げたステージは `RAN:` のとおり実行済み

3. `COMMIT:` のファイルだけをコミットする（rebase の hook が付け替えた記録は足さない。4 で別のコミットにする）。
   件名は `chore: relock <ステージ> after the rebase`。本文に、なぜ再実行したか（deps の何が変わったか）、分岐点
   （`git merge-base HEAD main`）を書き、`Refs: ADR-0130` と、/git-commit の規則の末尾（Co-Authored-By など）を
   添える。実行したコミットはこのコミット自身なので、記録には、このコミットか、その結果を出したコード（ツールを含む）を
   含む最初のコミットを書く。

4. 記録（verification.md、ADR）の hash と分岐点を直す。ステージで分ける:

   - 再実行したステージ: 実行したコミットを 3 のコミット（またはその結果を出したコードを含む最初のコミット）に、分岐点を今の `git merge-base HEAD main` にする。
     上流が changed のまま実行したら、その `dvc status` に changed で出た上流のすべての出どころ（それぞれを実行したコミットと
     分岐点）を別に書く（ADR-0184）
   - 走らせなかったステージ: 分岐点は、結果を出したときの main のコミットのまま変えない（今の merge-base に変えると、
     復元が合わなくなる）。実行したコミットが rebase で付け替わっていたら、hash だけを付け替える。記録に実行した
     コミットと分岐点が書かれていなければ、書く（ADR-0130 の「確認」の形）

   hash の付け替えは、手で書き換えず、`commit-refs remap` で行う。rebase のあとは、`remap.commitRefs` が true なら
   post-rewrite hook が自動で呼ぶ（rebase の途中の amend では呼ばず、終わりに一度）（付け替えはコミットされないので、自分でコミットする）。hook が失敗したときは、
   `$(git rev-parse --absolute-git-dir)/commit-refs-remap.input` の対応表で `commit-refs remap <そのファイル>` を呼び直す
   （詳細は verification-record の「rebase のあと」の節）。手で対応表を作るのは、この対応表が無いときだけ。rebase の前の tip
   （`git reflog` の `ORIG_HEAD` など）と今の HEAD を `git range-diff` で対応させ、中身を確かめてから渡す（全ファイルを
   書き換えるので、推測した対応を渡さない）。再測定はしない（記録の hash は再現のため）。再実行したステージの新しいコミットと分岐点の値は、rebase の対応表に載らないので、手で書く（書いた hash は `git rev-parse` で確かめる）。付け替えのコミットは `docs:` で分ける。

出力を作り直さず `dvc commit` で lock を付け直してよいのは、ADR が決めた移行だけである（値の同じ TOML の差は、ADR-0254 の
正規形で deps に表れない）。実行ファイルと
`references-shell.txt` の差では付け直さない。
付け直したとき（deps の置き場や正規形を替えた ADR-0210 の移行など）は、deps の hash
だけが変わったことを `nix develop -c scripts/compare-dvc-lock.py <付け直す前のコミット> HEAD` で確かめる（2 つのコミットの間で変わった
すべての `dvc.lock` について、ステージの集合と、各ステージの `cmd`、`outs`、`params` を比べ、違いはファイルとステージと欄で出して
終了コード 1 で終わる）。ステージの集合の差も違いとして出る（python-ast のステージの追加のように増えて当然のときは、出力を読んで確かめる）。
ステージを再実行した lock は `outs` が変わるので、この確認は付け直しだけのコミットの範囲に使う。

この後は /git-merge の「チェック」に進む。

## マージが止まったら: ブランチが .gitignore から外した古い出力（main の worktree で）

ブランチがステージを消すか outs の名前を変えると、DVC が書いた `.gitignore` の行も消え、main の作業ツリーに残る古い出力
（以前の `dvc checkout` / `dvc repro` が置いたもの）がマージの結果では未追跡のファイルになる。pre-merge-commit の hook は
それを「次のファイルは、マージするブランチが .gitignore から外した生成物」の見出しの別の一覧に出して止まり、記録された版
（DVC のキャッシュなど）と中身が同じで戻せるなら消す、違うか記録が無ければ原因を調べる、と示す（2026-11-03 に #0093 が
examples/core の `lattice` のステージ（`foreach: [out, aic, b4c]`）を消し、その出力の `fg-{out,aic,b4c}.msgpack` と
`lattice-{out,aic,b4c}.out` の 6 本で止まった）。この節は、その一覧の DVC の出力について、消してよいかの決め方を書く。

/git-merge が「ブランチで直す」と書くのは品質チェックで拒まれたときだけで、作業ツリーにマージ結果以外の変更があるという
拒否はこの節が扱う。/topic-start は main の作業ツリーで dvc を実行しないとするが、この一覧のファイルは main の作業ツリーにだけ
ある追跡しないファイルで、ブランチでは直せない。下の手順（main の worktree での
`git merge --abort`、`dvc status`、消すこと）はその例外で、ブランチは変わらないので rebase、レビュー、チェックはやり直さない。
一覧のほかの行（見出しの前の変更と未追跡のファイル）は hook の修正方法と /git-merge のとおりに扱う。

コマンドは main の worktree の直下で実行する（`cd <main の worktree> && …`。ブランチの worktree の `dvc.lock` には消えたステージが
無く、出力でないと判定される）。

1. main の worktree で `git merge --abort` する。この後は main の HEAD の `dvc.lock` と `.dvc` で判断する（マージの途中の
   作業ツリーにはブランチの `dvc.lock` があり、消えたステージが載っていない）。hook が別の一覧を出さず、1 つの一覧だけを出した
   とき（「判定に失敗した」と出たとき、rebuild の前の古い hook）は、その一覧のパスのうち、`git check-ignore -- <パス>...` が
   出すもの（main の `.gitignore` で無視されている。ディレクトリのパターンの中のファイルも出る）を同じく扱う
2. 一覧のパスごとに、main の記録と比べる（パスはリポジトリの直下から）:

   ```bash
   nix develop -c dvc status <パス>
   ```

   パスを持つステージ（`.dvc` で追う生成物なら、その `.dvc`）の状態だけが出る。見方（dvc 3.67.1 の `Output.status`）:

   - 終了コードが 0 でない: `does not exist as an output or a stage name` なら、どのステージの出力でも `.dvc` でもなく、記録が
     無い。ほかの失敗（dvc.yaml が読めない、`nix develop` の失敗など）も含め、消さずに原因を調べる
   - 終了コード 0 で、`changed outs` にそのパスが無い（ディレクトリの出力の中のファイルは、そのディレクトリの行で見る。
     `changed deps` は見ない）: 中身が main の `dvc.lock` の md5 と同じで、キャッシュにある。消してよい。ただし `dvc.yaml` の
     outs に `cache: false` のある出力は、dvc がキャッシュを確かめず、キャッシュに無くても何も出ないので、消さずに原因を調べる
   - `modified:`: 中身が `dvc.lock` の md5 と違う（DVC を入れる前の出力、手で流した出力、後半を行わなかった前のマージの版など）。
     main の記録からは戻らないので、消さずに何の出力でなぜ違うかを調べる
   - `not in cache:`: `dvc.lock` の md5 の版がキャッシュに無い。中身が同じでも消すと戻らないので、消さずに原因を調べる

   1 本のファイルの md5 で確かめる方法（`md5sum` とキャッシュの置き場）は、`cache: false` の出力と前の版の区別ができないので
   使わない
3. 消してよいと決めたパスだけを main の worktree で消し、同じ `git merge` をやり直す。消さないパスが残れば、マージせずに
   原因と一緒にユーザーに聞く。マージをやめたら、消したパス（ディレクトリの出力なら、そのディレクトリ）だけを指定して
   `nix develop -c scripts/dvc-checkout.py <パス>...` で戻す（例の `dvc.yaml` を指定すると、ほかの出力も置き換える）

## 後半: マージのあと（main の worktree で）

ブランチの worktree ではなく、main の worktree の直下で行う（出力を置き換える先が main）。

main の出力は、DVC を入れる前のもの（`vera.msgpack` など）か、前のマージの時点のものが残っている。dvc.lock の出力に
置き換える。`bin/` の実行ファイル（リポジトリの直下の `dvc.yaml` の build。ADR-0159）も戻す:

```bash
nix develop -c scripts/dvc-checkout.py dvc.yaml
nix develop -c scripts/dvc-checkout.py examples/<問題>/dvc.yaml
nix develop -c dvc status --with-deps examples/<問題>/dvc.yaml
```

changed のステージの数は `nix develop -c scripts/count-changed-stages.py examples/<問題>/dvc.yaml` で数える（`dvc status` の文字列を grep で数えると誤る）。

`-P` ではなく、ブランチが `dvc.yaml` を足した、または変えた例ごとに行う。`dvc status` の changed が、前半で走らせ
なかったステージ（deps が `dvc.lock` と合わないステージ）と、その下流のステージだけであることを確かめる。`references-shell`
（`always_changed`）は常に出るので数えない。ブランチが `build` を実行せずにコードを変えたなら、`build` とその下流も出る
（ADR-0159 の帰結。`dvc repro` で作り直すまで、そのままでよい）。ほかのステージが changed なら、
出力が足りない、またはキャッシュにない。`scripts/dvc-checkout.py` が失敗したら止めて、原因を調べる（推測で `dvc repro` を実行しない）。

ブランチが `.dvc` ファイル（記録が引く生成物。ADR-0222）を足すか変えたなら、main の置き場もキャッシュに合わせる（`HEAD` は
マージのコミット。名前を変えた `.dvc` は新しい名前で並ぶ）。並べたものが空なら実行しない（引数なしの checkout は、すべての
dvc.yaml の出力を対象にする）:

```bash
git diff --name-only --diff-filter=AMR HEAD^1 HEAD -- '*.dvc'
nix develop -c scripts/dvc-checkout.py <上で並んだ .dvc>
```

main の作業ツリーに同じ名前の、キャッシュに無い古い版が残っていると、dvc は「Can't remove the following unsaved
files」で止まる。`--force` を付けると古い版を失うので付けず、記録がその版を md5 で引くかを確かめてユーザーに聞く。

走らせなかったステージの changed は、そのままでよい。`dvc.lock` は、どの deps から出力を作ったかの記録であり、結果が
要るときに HEAD で `dvc repro` すれば再実行する（ADR-0130）。
