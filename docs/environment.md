# 検証用の環境

xlc の計算は、NJOY（核データ処理）と OpenMC（モンテカルロ計算）の結果と比べて検証する。ツール（NJOY、OpenMC と
その Python API）と評価済み核データは、`flake.nix` の packages で版とハッシュごと固定する（ADR-0124）。
計算の置き場 `references/` は git で追跡しない。そのうち外から取る入力とツールは、消しても下の手順で作り直せる
（取る先の分からない仕様書の PDF を除く）。ツールで作る参照の結果は DVC のステージで作り直す（「VERA の参照解」。段階に
していない問題が残る。#0018）。
direnv を使う場合は、リポジトリの直下に `use flake` と書いた `.envrc` を置き（手元の設定として git では追跡しない）、
一度 `direnv allow` すると、`default` の devShell に自動で入る。Claude Code の LSP のサーバー（rust-analyzer）は Claude Code の
起動時の PATH を引き継ぐので、Claude Code は devShell の中で起動する（`.envrc` のない worktree では
`nix develop -c claude`）。devShell が変わったら起動し直す（ADR-0135）。

| devShell     | 用途                                                                                                         |
| ------------ | ------------------------------------------------------------------------------------------------------------ |
| `default`    | 開発（Rust と rust-analyzer。ADR-0135）、段階の駆動（DVC）と、検証用の Python スクリプト（numpy, h5py など） |
| `references` | 参照の計算: `njoy`、`openmc` と、OpenMC の Python API（openmc、endf）の入った `python3`                      |

## 段階の駆動（DVC）

入力から出力を作る段階は DVC（`default` の devShell に入っている）で実行する（ADR-0123）。段階は
`examples/<問題>/dvc.yaml` に書き、`dvc repro examples/<問題>/dvc.yaml` で、入力かコードが変わった段階だけを実行する。
各段階の入力・コード・出力の md5 は同じディレクトリの `dvc.lock` に書かれ、コミットする。出力は
`references/dvc-cache`（`.dvc/config`）にキャッシュし、`dvc checkout` で戻せる。worktree では `references` を main の
`references/` への symlink にして、キャッシュと核データを共有する。`~/.claude/scripts/worktree-add.sh` が、チェックアウトした
ブランチの `.worktreelink` を読んで作る（出力の `SHARED: references`）。`.worktreelink` のないブランチ（main に rebase する前の
もの）や、ほかの方法で作った worktree では出ないので、worktree の直下で手で張る:
`ln -s "$(git rev-parse --path-format=absolute --git-common-dir)/../references" references`。symlink のない linked
worktree は、その worktree の `cargo test`（crates/docs の `linked_worktrees_share_references`。run-checks.sh が実行する）が
止める（ADR-0126）。

## 評価済み核データ（ENDF/B-VIII.0）

NNDC の zip（neutrons、thermal_scatt、decay、nfy）を URL と sha256 で固定し（`.#endf-neutrons` などの packages）、
展開した木を `.#endf` にする。`references/endf` は、その木への symlink（GC の root を兼ねる）にする:

```sh
nix build .#endf -o references/endf
```

main の checkout の直下で実行する。worktree から実行すると、GC の root が worktree の中のパス（symlink の
`references` を通したパス）で登録され、worktree を消すと root が外れて、ENDF の木が GC で消えうる。
`references/endf` に手で取った木（ディレクトリ）があると `nix build` は symlink を作れないので、先にどける。worktree の
`references` は main の `references/` を指すので、ほかの worktree で ENDF を読む計算が走っていないときに行う。

入力の TOML と dvc.yaml は `references/endf/ENDF-B-VIII.0_<part>/...` を指す。

## 仕様書の PDF

VERA、C5G7、MOX の仕様書は再配布できないので、URL と sha256 だけを docs/specs.sha256 に書く。URL から
`references/` の同じパスに取り、リポジトリの直下で照らす:

```sh
sha256sum -c docs/specs.sha256
```

「URL 不明」の PDF（C5G7 と MOX の仕様）は取り直す先がないので、`references/` から消さない。

## NJOY2016（2016.79）

`.#njoy2016` は、NJOY2016 2016.79 に `docs/patches/njoy2016-acer-law0.patch` を当ててビルドしたもの。
`references` の devShell の `PATH` に `njoy` として入る。

`docs/patches/njoy2016-acer-law0.patch` は ACER のバグを直す。MF6 の生成物の分布が LAW=0（分布が
与えられていない）だと、ACER は荷電粒子の生成データにその項目を LNW, LAW, IDAT の3語だけで書くが、
書き出し（`change`）は分布の見出しが続くものとして読み進め、止まる（B-10 の MT700 の α 粒子）。
パッチは LAW=0 の項目でその3語のあと次の反応に進む。荷電粒子の生成データだけにかかわり、
中性子の断面積は変わらない（OpenMC は荷電粒子の生成データを読まない）。

NJOY の入力と出力は `references/njoy-runs/` に置き、`crates/xs/examples/compare_njoy_*` で比べる。

## OpenMC（0.15.2）

`.#openmc` は OpenMC 0.15.2 を OpenMP つきでビルドしたもの（同梱の fmt が GCC 16 でコンパイルできないので GCC 13）。
`.#openmc-python` は、同じソースから作った openmc の Python API と endf 0.1.12 を、検証のスクリプトの依存に足した
Python。どちらも `references` の devShell に入る。スクリプトは `openmc` と `njoy` を `PATH` から呼ぶ。

コマンドはリポジトリ（worktree）の直下で実行し、作業のディレクトリを引数で渡す。worktree の `references` は main の
`references/` への symlink なので、`references/openmc` に `cd` してから `../..` で flake やスクリプトを指すと、
worktree でなく main の flake とスクリプトが使われる。たとえば C5G7 の参照解は:

```sh
nix develop .#references -c sh -c \
  'OMP_NUM_THREADS=16 python3 examples/c5g7/openmc_model.py references/openmc/work 1000000 400 100 2'
```

## 外部の計算機（AWS）

r995 で待ち時間が長くなる OpenMC の参照解は、AWS EC2 の c7a のスポットで流す（ADR-0151。検討段階）。c7a は vCPU が物理コアで、
x86_64 の nix のクロージャをそのまま送れる。起動、実行、回収、削除のツールは未実装（feat/0073-openmc-library のブランチで作業中）。

ユーザーが用意するもの:

- AWS アカウントと支払い方法
- 予算アラート（AWS Budgets）
- スポットの vCPU クォータ（承認に日数がかかるので最初に申請する。4 本並列で 64 vCPU。"All Standard Spot Instance Requests"。リージョンごと。既定値は確認していないので、足りなければ引き上げを申請する）
- リージョン（us-east-1 か us-east-2 を既定にする）
- 次の方式のどちらか
  - 方式 A（SSH の接続先を渡す）: Ubuntu の c7a を起動し（1 本に c7a.4xlarge を 1 台）、SSH 鍵と、自分の IP からの SSH を許可するセキュリティグループを用意して、ホスト名とユーザー名を渡す。VM の削除と、スポット中断後の再起動はユーザーが行う
  - 方式 B（認証情報を渡す）: EC2 の起動、停止、終了、describe に絞った IAM ユーザー（またはロール）のアクセスキーを渡す。起動には `iam:PassRole` を使う場合のほか、キーペア、セキュリティグループ、スポットの要求の権限も要るので、最小の権限はツールを作るときに決める。予算上限を設定し、使い終わったら失効させる。起動、再起動、削除までこちらで行える

こちらで行うこと:

- VM に nix を入れる
- `nix copy --to ssh-ng://<host> .#openmc .#openmc-python` でクロージャを送る
- 核データ（591 MB）と入力、スクリプトを送る
- 1 本ずつ別の VM で、並列に実行する（スレッド数は VM の物理コア数に合わせる）
- 結果（数十 MB）を回収する
- VM を削除する（方式 B）

## VERA の参照解

VERA Problem 1B〜1D, 2A〜2K, 2Q の参照解は、ENDF/B-VIII.0 から手元の NJOY で作った連続エネルギーの
データで OpenMC を動かして作る（温度は examples/vera/openmc-library.toml。565 K は燃料、冷却材、被覆管の核種を 565 K でも
作り、H in H2O の S(α,β) は 550 K と 573.6 K の内挿。ADR-0080）。段階は examples/vera/dvc.yaml にあり（ADR-0130）、`default` の devShell で、
リポジトリ（worktree）の直下から実行する:

```sh
dvc repro examples/vera/dvc.yaml
```

NJOY と OpenMC を使う段階は、cmd が `nix develop ../..#references -c` で `references` の devShell に入るので、
`dvc repro` を `references` の devShell の中で実行しない（cargo と dvc がない）。1 分を超えるので、`dvc repro` 全体を
Slurm で流す（CLAUDE.md の「計算」）。

| 段階                                               | 入力                  | 出力（examples/vera/ の下）                                              |
| -------------------------------------------------- | --------------------- | ------------------------------------------------------------------------ |
| `references-shell`                                 | flake                 | `references-shell.txt`（`references` の devShell の store のパス）       |
| `openmc-neutron@<核種>`、`openmc-thermal@<散乱則>` | `openmc-library.toml` | `openmc-library/neutron/<核種>.h5`、`openmc-library/thermal/<散乱則>.h5` |
| `openmc-index`                                     | `openmc-library.toml` | `openmc-library/cross_sections.xml`                                      |
| `openmc-1B`                                        | `p1/openmc-1B.toml`   | `p1/openmc-1B/`（statepoint、summary、model.xml）                        |
| `openmc-1B-extract`                                | `p1/openmc-1B.toml`   | `p1/openmc-1B-results/k.json`                                            |
| `xs-range`、`xs-nuclide@<核種>`、`xs`              | `xs.toml`             | `xs-range.toml`、`xs-nuclides/<核種>.msgpack`、`vera.msgpack`            |
| `lattice-1B`                                       | `p1/pincell-1B.toml`  | `p1/lattice-1B.out`                                                      |
| `compare-1B`                                       |                       | `p1/compare-1B.out`（xlc と OpenMC の k の差）                           |

計算の条件（核種と温度、粒子数、世代、非アクティブバッチ、乱数のシード、スレッド数）は入力の TOML に書く。スクリプトに既定値は
なく、足りないキーと知らないキーで止まる。出力は `references/dvc-cache` にキャッシュされ、`dvc checkout
examples/vera/dvc.yaml` で戻せる。`dvc.lock` は、記録した実行の deps（スクリプト、入力、ツール）から出力を作ったことの記録で、
あとのコミットで deps が変わると `dvc status` は changed を示す。記録の結果を出し直すときは、記録に書いたコミットを
checkout して `dvc repro` する（ADR-0130）。ほかの問題（1C〜2K、2E-empty、2Q）は、`p1/openmc-1B.toml` と同じ形の入力と、`openmc-1B`、
`openmc-1B-extract` と同じ形の段階を足して実行する（#0018）。段階にするまでは、同じ入力を作って
`python3 examples/vera/openmc_model.py run <入力> <作業ディレクトリ>` と `extract <入力> <作業ディレクトリ> <結果の
ディレクトリ>` を `references` の devShell で実行する。

同じ多群定数で OpenMC の多群モードと比べる（examples/vera/verification.md の「同じ多群定数での OpenMC との比較」）には、
負にならない散乱核のライブラリと、そこから多群定数を作ってから解く。格子の入力は `library` をそのライブラリに、
`self-shielding = false` と `scattering-order` を同じ次数にしたもの:

```sh
nix develop .#references -c python3 examples/vera/positive_library.py examples/vera/vera.msgpack \
  references/openmc/positive-P3.msgpack 3
cargo run --release -p xlc-lattice --example mg_materials -- <格子の入力> references/openmc/mg-2B.msgpack
nix develop .#references -c sh -c \
  'OMP_NUM_THREADS=16 python3 examples/vera/openmc_mg.py references/openmc/mg-2B.msgpack 2B \
     references/openmc/work-mg-2B 1000000 300 50 1'
```
