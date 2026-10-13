# 検証用の環境

xlc の計算は、NJOY（核データ処理）と OpenMC（モンテカルロ計算）の結果と比べて検証する。道具（NJOY、OpenMC と
その Python API）と評価済み核データは、`flake.nix` の packages で版とハッシュごと固定する（ADR-0124）。
計算の置き場 `references/` は git で追跡しない。そのうち外から取る入力と道具は、消しても下の手順で作り直せる
（取る先の分からない仕様書の PDF を除く）。道具で作る参照の結果は DVC の段階で作り直す（「VERA の参照解」。段階に
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

入力から出力を作る段階は DVC（`default` の devShell に入っている）で走らせる（ADR-0123）。段階は
`examples/<問題>/dvc.yaml` に書き、`dvc repro examples/<問題>/dvc.yaml` で、入力かコードが変わった段階だけを走らせる。
各段階の入力・コード・出力の md5 は同じディレクトリの `dvc.lock` に書かれ、コミットする。出力は
`references/dvc-cache`（`.dvc/config`）にキャッシュし、`dvc checkout` で戻せる。worktree では `references` を main の
`references/` への symlink にして、キャッシュと核データを共有する。`~/.claude/scripts/worktree-add.sh` が、チェックアウトした
ブランチの `.worktreelink` を読んで作る（出力の `SHARED: references`）。`.worktreelink` のないブランチ（main に rebase する前の
もの）や、ほかの方法で作った worktree では出ないので、worktree の直下で手で張る:
`ln -s "$(git rev-parse --path-format=absolute --git-common-dir)/../references" references`。symlink のない linked
worktree は、その worktree の `cargo test`（crates/docs の `linked_worktrees_share_references`。run-checks.sh が走らせる）が
止める（ADR-0126）。

## 評価済み核データ（ENDF/B-VIII.0）

NNDC の zip（neutrons、thermal_scatt、decay、nfy）を URL と sha256 で固定し（`.#endf-neutrons` などの packages）、
展開した木を `.#endf` にする。`references/endf` は、その木への symlink（GC の root を兼ねる）にする:

```sh
nix build .#endf -o references/endf
```

main の checkout の直下で走らせる。worktree から走らせると、GC の root が worktree の中のパス（symlink の
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

コマンドはリポジトリ（worktree）の直下で走らせ、作業のディレクトリを引数で渡す。worktree の `references` は main の
`references/` への symlink なので、`references/openmc` に `cd` してから `../..` で flake やスクリプトを指すと、
worktree でなく main の flake とスクリプトが使われる。たとえば C5G7 の参照解は:

```sh
nix develop .#references -c sh -c \
  'OMP_NUM_THREADS=16 python3 examples/c5g7/openmc_model.py references/openmc/work 1000000 400 100 2'
```

## 外部の計算機（AWS）

r995 で待ち時間が長くなる OpenMC の参照解は、AWS EC2 の c7a のスポットで流す（ADR-0151。検討段階）。c7a は vCPU が物理コアで、
x86_64 の nix の閉包をそのまま送れる。起動、終了、送る、実行、回収の道具は次の節（方式 A の変種で、起動と終了も道具で行う）。

ユーザーが用意するもの:

- AWS アカウントと支払い方法
- 予算アラート（AWS Budgets）
- スポットの vCPU クォータ（承認に日数がかかるので最初に申請する。4 本並列で 64 vCPU。"All Standard Spot Instance Requests"。リージョンごと。既定値は確認していないので、足りなければ引き上げを申請する）
- リージョン（us-east-1 か us-east-2 を既定にする）
- 次の方式のどちらか
  - 方式 A（SSH の接続先を渡す）: Ubuntu の c7a を起動し（1 本に c7a.4xlarge を 1 台）、SSH 鍵と、自分の IP からの SSH を許可するセキュリティグループを用意して、ホスト名とユーザー名を渡す。VM の削除と、スポット中断後の再起動はユーザーが行う
  - 方式 B（認証情報を渡す）: EC2 の起動、停止、終了、describe に絞った IAM ユーザー（またはロール）のアクセスキーを渡す。起動には `iam:PassRole` を使う場合のほか、キーペア、セキュリティグループ、スポットの要求の権限も要るので、最小の権限は道具を作るときに決める。予算上限を設定し、使い終わったら失効させる。起動、再起動、削除までこちらで行える

こちらで行うこと:

- VM に nix を入れる
- `nix copy --to ssh-ng://<host> .#openmc .#openmc-python` で閉包を送る
- 核データ（591 MB）と入力、スクリプトを送る
- 1 本ずつ別の VM で、並列に実行する（スレッド数は VM の物理コア数に合わせる）
- 結果（数十 MB）を回収する
- VM を削除する（方式 B）

## 外部の計算機での参照解（Problem 3、5-2D）

Problem 3 と 5-2D の本番（数時間の走行を 5 本）は、r995 でなく外部の VM（AWS、GCP など）で流す。VM には Nix だけを入れ、
OpenMC と Python API は手元でビルドした閉包を送る（VM で再ビルドしない）。Slurm は使わない。VM は中断されうる（スポット）ので、
入力の `statepoint-interval` で途中の statepoint を書き、中断のあとは `--restart` で再開する。

### 道具を使う手順

道具は `examples/vera/` にある。リージョンは us-east-1（`AWS_DEFAULT_REGION` で変えられる）、認証は `aws sts get-caller-identity` が通る状態。
キーペア `xlc-aws`（`~/.ssh/xlc-aws.pem`）とセキュリティグループ `xlc-ssh`（22/tcp を手元の IP の /32 だけ）は作成済みで、残してある。

```sh
examples/vera/aws_launch.sh xlc3 --count 2          # 起動。c7a.4xlarge のスポット（one-time、gp3 40 GB）を xlc3-1、xlc3-2 の名前で
examples/vera/remote_send.sh --allow-ip xlc3-1 /home/ubuntu/xlc   # 送る（各 VM に。下の「VM に Nix を入れる」を先に）
examples/vera/aws_status.sh                         # 走っている VM、起動からの時間、スポット価格
examples/vera/remote_fetch.sh --allow-ip xlc3-1 /home/ubuntu/xlc references/openmc/remote/xlc3-1   # 戻す
examples/vera/aws_terminate.sh xlc3                 # 終了（terminated を待ち、~/.ssh/config の Host を消す）
```

- `aws_launch.sh` は、AMI（Ubuntu 24.04 の最新。SCP で SSM の公開パラメータが引けないので `describe-images --owners 099720109477`）を引き、
  スポットで起動して public IP を待ち、`~/.ssh/config` に `# BEGIN xlc-aws <名前>` から `# END` までの Host（`User ubuntu`、`IdentityFile`、
  `StrictHostKeyChecking accept-new`）を書く（同名は置き換え）。起動後に `lscpu` と `free -g` を出すので、物理コア数（Thread(s) per core が 1 であること）と
  メモリを `--cores` と記録に使う。同名のインスタンスが残っていれば止まる
- IP の許可: 手元のグローバル IP が変わると SSH が通らなくなる。`aws_allow_ip.sh` が `xlc-ssh` の 22/tcp を今の IP の /32 だけにする（古い規則は消す。一致なら何もしない）。
  `aws_launch.sh` は起動の前に必ず呼ぶ（新しい VM に入れないと意味がないため）。`remote_send.sh` と `remote_fetch.sh` は `--allow-ip` を付けたときだけ呼ぶ。
  この 2 つは AWS 以外の VM（GCP など）にも使うこと、セキュリティグループの変更は外から見える副作用なので、黙って呼ばず明示にした。
  ssh が通らなくなったら、`aws_allow_ip.sh` を単体で走らせるか、`--allow-ip` を付ける
- 中断（スポットの回収）: `aws_status.sh` で消えたことを見て、同じ名前で `aws_launch.sh` を走らせ直す（IP が変わるので Host も更新される）。
  Nix を入れ直し、`remote_send.sh` で送り直してから、`run_remote.sh --restart` で最新の statepoint から再開する。statepoint は VM の EBS にあり、
  one-time のスポットが回収されると消える。再開できるのは、中断前に `remote_fetch.sh` で戻した statepoint を送り返した場合（作業ディレクトリの `results/<名前>/run/`）だけで、
  戻していなければその 1 本をやり直す
- 終わったら必ず `aws_terminate.sh`。`aws_status.sh` で空であることを確認する（走らせたままだと課金が続く）

費用の見積もり: 台数 × 起動から終了までの時間（時間）× スポット価格（`aws_status.sh` が出す $/時間。c7a.4xlarge は 2026-10 に確認していない。
ADR-0151 の c7a.16xlarge の $0.015〜0.023/物理コア時間 × 16 コア ≒ $0.24〜0.37/時間を目安にする）。EBS は gp3 40 GB × 起動時間（月 $0.08/GB なので 1 日あたり約 $0.1）で、終了で消える。
本番は 2 台 × 2 回 × 約 6 時間なら約 $6〜9。試走（c7a.large、`aws_launch.sh` から `aws_terminate.sh` まで約 1 分）は 1 セント未満。

#### 道具の確認（2026-10-13、xlc のコミット 66077c0c のスクリプト。以降のレビューの修正のあと、同じ往復をもう一度通した）

c7a.large を 1 台、`aws_launch.sh`（約 30 秒。AMI の検索から ssh が通って `lscpu` を出すまで）→ `aws_allow_ip.sh --ip 203.0.113.7`（わざと別の IP にすると ssh は
タイムアウト）→ `aws_allow_ip.sh`（戻すと ssh が通る）→ `aws_terminate.sh`（約 37 秒、terminated を確認、Host が消えて `~/.ssh/config` が元に戻る）の往復を通した。
c7a.4xlarge は起動していない。

### VM に Nix を入れる、送る、流す、戻す

1. VM に Nix を入れる。試走（2026-10-13、Ubuntu 24.04、c7a.large）で通した手順は Determinate Systems の installer（multi-user）。
   約 10 秒。`nix copy --to ssh://` は未署名のパスを送るので、VM の ssh のユーザーを信頼するユーザーにする。非対話の ssh でも
   `nix-store` が PATH にある（`/nix/var/nix/profiles/default/bin`）ので、PATH の追加は要らない。`rsync` は Ubuntu の AMI に入っている:

   ```sh
   ssh <接続先> 'curl -fsSL https://install.determinate.systems/nix | sh -s -- install --no-confirm'
   ssh <接続先> 'echo "trusted-users = root $USER" | sudo tee /etc/nix/nix.custom.conf; sudo systemctl restart nix-daemon'
   ```

   接続先は `~/.ssh/config` の Host にしておくと、`nix copy` と rsync にも同じ設定が効く（`aws_launch.sh` が書く）
2. 手元（リポジトリの直下）から送る。`--dry-run` は、閉包のパス（62 個、約 800 MB）と rsync の一覧だけ出して何も送らない:

   ```sh
   dvc checkout examples/vera/dvc.yaml          # examples/vera/openmc-library（約 630 MB）
   examples/vera/remote_send.sh --dry-run <ssh の接続先> <VM の作業ディレクトリ>
   examples/vera/remote_send.sh <ssh の接続先> <VM の作業ディレクトリ>
   ```

   閉包（flake の `.#openmc` と `.#openmc-python`）は `nix copy --to ssh://<接続先>`、ライブラリ、入力（p3、p5）、スクリプトは
   rsync で `<VM の作業ディレクトリ>` の下に同じ配置で送る。VM 上の道具の場所は `tools.env` に書く。VM の CPU は x86_64 の Linux
   （閉包が x86_64-linux のため）。試走の時間: 閉包 800 MB と核データ 630 MB の送信は合わせて 2 分 36 秒（手元から us-east-1 へ）
3. VM で、作業ディレクトリの直下から流す。1 本あたりのスレッド数は `--cores`（物理コア数）を入力の数で割った値。省略すると `lscpu` から数える。
   c7a は SMT がなく vCPU が物理コア（c7a.large は 2、`lscpu` の Thread(s) per core が 1。試走で確認）。SMT のある機種では物理コアは
   vCPU 数の半分になる。切り離して走るので、ssh を切ってよい:

   ```sh
   examples/vera/run_remote.sh --cores <物理コア数> examples/vera/p3/openmc-3A.toml examples/vera/p3/openmc-3B.toml \
     examples/vera/p5/openmc-5B-2D.toml examples/vera/p5/openmc-5C-2D.toml
   tail -n 3 results/*/log; ls results/*/done     # 進み具合と完了
   ```

   run が終わると続けて extract も走り、`results/<名前>/` に k.json、pins.msgpack、layers.msgpack などができる
4. 中断されたら、VM を起動し直し（上の「中断」）、同じコマンドに `--restart` を足す。最新の statepoint（`results/<名前>/run/statepoint.<世代>.h5`。最新の
   2 個を残す）から `openmc -r` で再開する。入力は同じものを使う（`--threads` だけ変えてよい）
5. 結果を手元に戻す。k.json、pins.msgpack などと最後の世代の statepoint と summary.h5（一次の成果物。手元で extract し直す。途中の statepoint は戻さない）:

   ```sh
   examples/vera/remote_fetch.sh --dry-run <ssh の接続先> <VM の作業ディレクトリ> references/openmc/remote
   examples/vera/remote_fetch.sh <ssh の接続先> <VM の作業ディレクトリ> references/openmc/remote
   ```

### 試走（2026-10-13）

c7a.large のスポット（2 vCPU = 2 物理コア）で、上の 1〜5 を通した（xlc のコミット 6b39ad13 のスクリプト）。入力は本番の写しで粒子数 1 万・20 世代・捨てる 5・
`statepoint-interval = 5`（`cross-sections` はリポジトリの `examples/vera/openmc-library` を指す）。

- 起動（今は `aws_launch.sh` が行う。試走のときの手順）: キーペア `xlc-aws`（秘密鍵 `~/.ssh/xlc-aws.pem`）、セキュリティグループ `xlc-ssh`（手元のグローバル IP の 22/tcp だけ）を作り、
  AMI は `aws ec2 describe-images --owners 099720109477 --filters "Name=name,Values=ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*"`
  に `--query 'sort_by(Images,&CreationDate)[-1].ImageId'` を付けた最新（このアカウントは SCP で SSM の公開パラメータが引けない）を `aws ec2 run-instances --instance-type c7a.large
  --instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time}'` で起動。EBS は gp3 40 GB。起動から ssh が通るまで約 40 秒。終わったら `aws ec2 terminate-instances` で終了する（EBS は DeleteOnTermination で消える。
  キーペアとセキュリティグループは残す）
- 3A と 5A-2D を `run_remote.sh`（`--cores` を省略。物理コア 2 を数え、2 本で 1 スレッドずつ）で同時に流し、約 1 分で k.json などまで
  できた。statepoint は 5 世代ごとに書かれ、最新の 2 個（15、20）が残る
- 3A を `--threads 2` で流して世代 10 の statepoint の後に殺し、`--restart` で再開すると、世代 10 から続けて完了し、最初の走行と k が
  一致した（1.172117 ± 0.002783）。殺すと `done` は書かれない
- `remote_fetch.sh` で戻した statepoint から手元で extract すると、k.json と pins.msgpack がバイト一致した（`openmc-python` の
  python で走らせる。`nix develop .#references` でもよい）

記録には、VM の機種、vCPU 数と物理コア数、`threads` を書き、時間は書かない（共有の環境で、スポットの中断もありうる）。

## VERA の参照解

VERA Problem 1B〜1D, 2A〜2K, 2Q の参照解は、ENDF/B-VIII.0 から手元の NJOY で作った連続エネルギーの
データで OpenMC を動かして作る（温度は examples/vera/openmc-library.toml。565 K は燃料、冷却材、被覆管の核種を 565 K でも
作り、H in H2O の S(α,β) は 550 K と 573.6 K の内挿。ADR-0080）。段階は examples/vera/dvc.yaml にあり（ADR-0130）、`default` の devShell で、
リポジトリ（worktree）の直下から走らせる:

```sh
dvc repro examples/vera/dvc.yaml
```

NJOY と OpenMC を使う段階は、cmd が `nix develop ../..#references -c` で `references` の devShell に入るので、
`dvc repro` を `references` の devShell の中で走らせない（cargo と dvc がない）。1 分を超えるので、`dvc repro` 全体を
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

計算の条件（核種と温度、粒子数、世代、捨てる世代、乱数の種、スレッド数）は入力の TOML に書く。スクリプトに既定値は
なく、足りないキーと知らないキーで止まる。出力は `references/dvc-cache` にキャッシュされ、`dvc checkout
examples/vera/dvc.yaml` で戻せる。`dvc.lock` は、記録した実行の deps（スクリプト、入力、道具）から出力を作ったことの記録で、
あとのコミットで deps が変わると `dvc status` は changed を示す。記録の結果を出し直すときは、記録に書いたコミットを
checkout して `dvc repro` する（ADR-0130）。ほかの問題（1C〜2K、2E-empty、2Q）は、`p1/openmc-1B.toml` と同じ形の入力と、`openmc-1B`、
`openmc-1B-extract` と同じ形の段階を足して走らせる（#0018）。段階にするまでは、同じ入力を作って
`python3 examples/vera/openmc_model.py run <入力> <作業ディレクトリ>` と `extract <入力> <作業ディレクトリ> <結果の
ディレクトリ>` を `references` の devShell で走らせる。

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
