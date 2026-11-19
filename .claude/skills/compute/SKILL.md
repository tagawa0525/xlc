---
name: compute
description: xlc で計算を流す前に読む。Slurm の投げ方と区分、スレッドの数、OpenMC と KENO の使い分け、DVC のステージの実行（dvc repro、stage_runner）、走るステージの数え方と見積もり、worktree での出力の戻し方、deps の正規形、dvc add、記録の数値とコミットの対応。
---

# 計算を流す

本番の作り直し（ライブラリ、燃焼の履歴、テーブル、参照解）の前に、`scripts/check-exclusions.py --list` で印の一覧を出し、その計算の入力が拠る印
（例の入力の excluded と provisional）を docs/plans/upstream-error-budget.md の項目に対応させる（#0294）。

- 重い計算（OpenMC、NJOY、ライブラリの作り直し、xlc の集合体以上、時間を記録する測定）は、Slurm で 1 本ずつ流す:
  `~/.claude/scripts/slurm-run.sh --class <区分> --ref '<参照>' --purpose '<目的>' [--top] [--timed | -c <CPU 数>] -t <見積もり> <名前> '<コマンド>'`（排他で、終わるまで待つ）。区分は timed（時間を記録する測定。`--timed` を付けなくてもパーティション timed）、record（記録に載せる本番）、decision（決定に要る確認）、trial（検討）で、この順に優先する。`--ref` は `#NNNN`、`ADR-NNNN`（複数はカンマ）、または `- <理由>` で、`#` がシェルの注釈になるので引用符で囲む。`--purpose` は結果が何を決めるかの 1 文。待ち行列は `~/.claude/scripts/slurm-queue.sh` で見る。時間を記録する測定は `--timed` を付ける（割り込まれないパーティション timed で、走っている main のジョブを止めて先に走る）。付けないジョブは main で、timed のジョブに止められうるので、その時間は記録に書かない（`TIME:` 行に `preemptible:` の印が付く）。`--top` は同じ区分の中で前に出す。走っているジョブの終わりを待てない数分の計算は、区分を変えずに `--timed` で投げる。本数の上限は置かないので、投げる前に 1 本の見積もりと合計の時間を報告に書き、再実行するときは理由を書く。直接の `sbatch` は hook が止める。直接実行せず、`pgrep` で他の計算を監視しない。時間を記録しない長い計算は `-c <CPU 数>`（ハードウェアスレッドで数える）を付けて排他なしに流せるが、その時間は記録に書かない。
- OpenMC のスレッド数は docs/environment.md の値を使う
- lattice のスレッドの数は、`RAYON_NUM_THREADS` が無ければプロセスに許された論理 CPU（SMT を含む）の数で、r995 の排他のジョブでは
  32（1 コア 2 スレッドの 16 コア。`lscpu`）。始めにログに出す（`thread pool threads=`）。xs の `--threads` の既定は
  機械の物理コアの数をプロセスに許された CPU の数で抑えたもの、runner の予算は許された CPU の中の物理コアの数で、r995 の排他のジョブでは
  どちらも 16、`-c 16` では xs 16、runner 8（ADR-0214、ADR-0193）。xs の数はログに出ないので、記録と
  依頼に書くときは、この規則と、ジョブで付けた `--threads`、`-c` の数から書く（#0186）
- 数時間の OpenMC の参照解も、xlc の計算と同じく r995 の Slurm で流す。AWS は枠が無いので今は使わない（ユーザー 2026-11-19）。
  短い OpenMC はローカル。投げ方は上の Slurm の項とスキル openmc-reference に従う。AWS の手順は docs/environment.md の「外部の計算機（AWS）」
- changed は作り直しの残高で、即時の義務ではない。修正は、結果を変えるものも変えないものも待たずに入れる。
  本番の作り直しは、docs/plans/upstream-error-budget.md の「進め方の決まり」のとおり「要」の項目を見積もり、結果を変える修正を
  まとめてから 1 回で行う（ADR-0258）
- 重い計算を足す前に、既存の結果で足りるか、結論が変わりうるかを見る。形式（履歴、テーブル、ライブラリ）を変える修正は、
  下流が古い出力を読めずに止まるので、読めなくなる出力と作り直す時間を数えて依頼と報告に書き、入れる時期を決める。
  読めなくなった出力の再実行は、それを読む作業が来るまで後回し（ユーザー 2026-10-11、2026-10-24）
- 高速化が進んでいる間に本番の計算を流すかは、場面ごとに決める。決まりではなく、高速化を終わらせてから流した方が時間を
  節約できるだろうという程度の提案（ユーザー 2026-10-24、2026-11-20）で、遅い bin で流した計算が再実行になるかを見て決める（#0300 の 5）
- 見込みの外れたジョブは、理由を書いて止める
- 仕様書に KENO-VI の参照解がある問題は ENDF/B-VII.0 の検証用ライブラリで KENO と直接比べ、OpenMC の参照解は仕様書に参照のない
  問題（燃焼、ブランチ計算、HFP）と診断だけに作る（ADR-0175）。3D の全炉心と Problem 4 の OpenMC は作らない（ADR-0124、docs/plans/m3-core-hzp.md。ユーザー 2026-10-15）
- 結果を変えない deps の変更で、重いステージを再実行せず `dvc commit` で付け直してよいのは、ADR が決めた移行
  （ADR-0159、ADR-0210。ADR-0254 の移行は lock の deps の書き換え）だけ（値の同じ TOML の差は、ADR-0254 の正規形で deps に表れない）。ファイルやディレクトリを動かしただけの変更は、`dvc commit` でなく
  `scripts/rewrite-lock-paths.py [--base <動かす前のコミット>]`（既定は HEAD。移動をコミットする前に実行する）で lock のパスを書き換え、ADR を書かずに
  コミットしてよい。`UNVERIFIED:` に出たステージは書き換わらず changed に残り（終了コード 1）、`COMPARISON:` の比較のステージは名前で再実行する（ADR-0260）。参照のステージ（NJOY、OpenMC）が使わない Python のパッケージは、flake.nix の `pythonPackages` にだけ足す（`references` の
  devShell の store のパスが参照解の deps。食い違いは pre-commit の hook `references-shell` が止める。ADR-0212）
- `examples/<問題>/dvc.yaml` のある例は、`dvc repro examples/<問題>/dvc.yaml`（すべては `dvc repro -P`）でステージを
  実行する。重いときは全体を上の形で Slurm に流し、ジョブのコマンドは `scripts/dvc-repro.py --max-stages <見込む本数>
  [--allow-reference <yaml>:<ステージか組>]... <目標>` にする（走りうるステージを数え、上限を超えるか許していない参照解が入れば
  何も実行せずに止まる。生の `dvc repro` を `slurm-run.sh` に渡すと hook が止める。ADR-0257）。重いステージを 1 本ずつ流す、または参照解の待ちで全体が失敗するときは、
  目標を指定する `dvc repro <yaml>:<ステージ>`（上流の changed は実行される。ADR-0184）。`-s` は上流を飛ばすので、rebase のあとに
  名前で挙げたステージを再実行する /dvc-relock（`scripts/dvc-relock.py --run`）と、上流の changed が deps だけで outs が `dvc.lock` のとおりのときに、目標だけを
  `dvc.lock` の上流の出力で実行するとき（重さによらない）だけに使う。例のステージの `-s` は、手元でも Slurm でも
  `scripts/dvc-repro.py -s [--allow-reference <yaml>:<ステージ>] <yaml>:<ステージ>` で実行する（上流が `changed outs`、
  `changed command`、`changed checksum` を持てば何も実行せずに止まり、changed の正規形は先に実行する。生の `dvc repro -s` は
  上流を確かめない。ADR-0184、ADR-0258）。記録には上流の出どころを別に書く。実行したら `dvc.lock` をコミットする（上流の `build`、`python-ast`、`binary-forms`、`toml-values` が書き換えたリポジトリの直下の `dvc.lock`、`python-ast/`、`binary-forms/`、`toml-values/` も。ADR-0159、ADR-0210、ADR-0255、ADR-0254）
- 核種ごとの xs のステージ（`xs-nuclide`）は、`dvc repro` でなく `python3 examples/stage_runner.py <yaml> xs-nuclide` で流す（changed の
  ものを CPU 時間の長い順に同時に実行して `dvc commit` する。上流の changed も実行する。下流の `xs` は続けて `scripts/dvc-repro.py --max-stages 1 <yaml>:xs`。ADR-0193）。
  runner が書き直す CPU 時間の記録 `<yaml のディレクトリ>/stage-cpu-seconds.json` も `dvc.lock` と一緒にコミットする
- `dvc repro` を投げる前に、出力を戻した、ジョブを実行する worktree の直下で `scripts/dvc-repro.py` に同じ引数と `--count-only` を
  付けて数え、`--max-stages` と見積もりを決める（ジョブも直下から投げる。ほかの場所ではスクリプトが止まる。手元で数えるのが長ければ（#0268）、見込みの本数を上限にして投げ、ジョブの中で数えさせる）。
  時間は、過去のログで実行したステージの数の一覧（`grep -o "Running stage '[^'@]*" <ログ> | sort | uniq -c`）の全体が同じジョブから取る。
  2026-10-22 に、65 核種の作り直しを多群の作成のないジョブから 45 分と見積もって 2 時間 5 分かかった（Slurm 452。ラダーの和をベクトル化した ADR-0188 の後は 38 分、Slurm 473）
- worktree の `references` は main の `references/` への symlink（`.worktreelink`。worktree-add.sh が作る）で、DVC の
  キャッシュと核データを共有する（ADR-0126）
- `examples/vera/vera.msgpack` と VERA の参照解（OpenMC のライブラリ、statepoint、k）は examples/vera/dvc.yaml のステージの
  出力で、git で追跡しない。worktree では `scripts/worktree-restore-outputs.sh`（中身は `scripts/dvc-checkout.py`）でキャッシュから戻す（ADR-0130）。
  一部だけ戻すときも `scripts/dvc-checkout.py [<dvc.yaml か .dvc の指定>...]` を使い、生の `dvc checkout` は `.claude/settings.json`
  の hook が止める（未実行の参照解でいつも 255 で終わり、キャッシュに無い出力が埋もれる。ADR-0223）。ステージが呼ぶ
  `bin/` の実行ファイルは、リポジトリの直下の `dvc.yaml` の `build` の出力で、同じスクリプトで戻る（ADR-0159）。ステージの
  deps は `bin/` でなく、`binary-forms` が書く正規形（行と列の番号を 0 にしたコピーの md5。`binary-forms/<名前>.md5`、git で追跡）に
  する（cmd は `bin/` を呼ぶ。`scripts/test_normalize_binary.py` が食い違いを止める。ADR-0255）。
  ステージの deps の .py は、元のファイルでなく、直下の `dvc.yaml` の `python-ast` が書く正規形（docstring とコメントを除いた形。
  `python-ast/<パス>.ast`、git で追跡）にする（cmd は元の .py を呼ぶ。ADR-0210）。入力の TOML（ステージの出力を除く）も同じく、`toml-values` が書く
  値の正規形（`toml-values/<パス>.json`）を deps にする（ADR-0254）。foreach と正規形は pre-commit の hook が deps から書き直す（ADR-0258）。
  NJOY と OpenMC を使うステージは cmd が `nix develop ../..#references -c` で包んであるので、`dvc repro` は default の
  devShell で実行する
- 記録か ADR が md5 か名前で引く、または比較のツールが書いて git で追跡する JSON が md5 で引く、git で追跡しない生成物
  （DVC のステージの出力でないもの）は、記録に書く前（JSON だけが引くものは JSON をコミットする前）に置き場で `dvc add` し、
  `.dvc` ファイルをコミットする（同じ名前で作り直したら `dvc add` をやり直す）。worktree の後片付けで消えず、
  `scripts/worktree-restore-outputs.sh` で戻り、前の版はキャッシュに残る（`dvc gc` は `--all-commits` で行う。ADR-0222、ADR-0123）
- 記録の数値には、それを出したコミットと日付を添える（checkout すれば再現できるため）。ブランチの上で実行した計算は、実行した
  コミットと分岐点（main のコミット）を書く。rebase で deps が変わったことを理由には再実行せず、分岐点から復元する（ADR-0130）。
  rebase で hash が変わっても
  再測定せず、ツール（`commit-refs remap`。rebase と amend では、`remap.commitRefs` が true なら post-rewrite hook が自動で呼ぶ（失敗したときだけ手で呼ぶ）。付け替えはコミットされないので、自分でコミットする）で付け替える。文書やテストの変更で hash が最新でなくなっても再測定しない。
  再測定するのは計算（コード、入力、環境）を変えたときだけ（ADR-0122）
