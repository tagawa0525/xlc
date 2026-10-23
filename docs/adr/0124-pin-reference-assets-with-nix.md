---
status: accepted
date: 2026-10-11
supersedes: [ADR-0031]
superseded-by:
issues: ["#0018"]
requires: []
---

# ADR-0124: 参照解のツールと核データは flake で版とハッシュごと固定し、references/ は作り直せる置き場にする

## 背景

格子計算と炉心計算の参照解は、同じ ENDF/B-VIII.0 から手元の NJOY で作った連続エネルギーのデータで OpenMC を解いて
求めている（ADR-0031）。ADR-0031 はツールを「ソースから `references/` にビルドする」と決め、手順は docs/environment.md
にあった: ENDF の zip は `curl` で取ってハッシュを確かめず、NJOY2016 と OpenMC は版を指定した `git clone` と手での
ビルド、OpenMC の Python API は venv に `pip install ./src` と PyPI の `endf`（版の指定なし）を入れていた。仕様書の
PDF は取る手順の記録がなかった。

ADR-0122 は、記録の hash を checkout すれば、固定した環境のもとで当時の結果を出し直せることを求め、その環境を #0018 で
整えるとした。ADR-0123 は計算のステージを DVC で駆動すると決め、examples/m0/dvc.yaml は `../../references/endf/...` の
ENDF のファイルと `flake.lock` を deps にしている。

判断の軸:

- 外から取るものが差し替わったら気づけること（URL とハッシュ）
- ツールの版、パッチ、コンパイラ、依存のライブラリを、コミットした記述から作り直せること
- `references/` を消しても、手順なしに作り直せること
- 環境の固定（Nix）と計算のステージの駆動（DVC、ADR-0123）の役割を分け、DVC の deps と食い違わないこと

## 検討した案

固定の手段:

- 案 A（採用）: flake.nix の packages にする。ENDF の zip は `fetchurl`（URL と sha256）、NJOY と OpenMC はソースを
  `fetchFromGitHub`（コミットと NAR ハッシュ）で取り、ビルドの手順を derivation に書く。コンパイラと依存の版は
  `flake.lock` の nixpkgs で決まる。nixpkgs には NJOY2016 も OpenMC もない（2026-10-11 に確かめた）ので、自前で書く
- 案 B: 取得とビルドのスクリプトに、URL と sha256 の一覧を添える。ハッシュは照らせるが、コンパイラと依存のライブラリは
  スクリプトを実行したときの環境に任され、作ったものがどの記述から来たかが残らない
- 案 C: DVC のステージで取り、ビルドする（`dvc import-url` など）。DVC は内容の md5 で出力を追うが、ツールのビルドの環境
  （コンパイラ、ライブラリ）を固定しない。環境は Nix、ステージの駆動は DVC と分ける（ADR-0123 の案 E と同じ切り分け）

ENDF の zip の取り方:

- 案 D（採用）: `fetchurl` で zip を取り、展開は別の derivation（`unzip`）にする。固定するハッシュは zip そのものの
  sha256 なので、配布元のファイルに `sha256sum` で照らせる
- 案 E: `fetchzip` で展開した木を固定する。ハッシュは展開した木の NAR ハッシュになり、配布物と直接照らせない

OpenMC の Python API（openmc、endf）の固定:

- 案 F（採用）: `buildPythonPackage` で作る。openmc は C++ と同じソース（同じ `fetchFromGitHub`）から、endf は PyPI の
  sdist（`fetchPypi`、0.1.12、sha256）からビルドする。numpy、h5py などの依存は nixpkgs のもので、venv が
  `--system-site-packages` で使っていたものと同じである。`python3.withPackages` で 1 つの Python にまとめる
- 案 G: venv のまま、`pip install --require-hashes` で endf の版とハッシュを手順に書く。ハッシュは固定できるが、venv は
  nix store の Python のパスを焼き込んだ `references/` の中の手作りの置き場で、Python が store から消えると黙って
  壊れ、どの flake から作ったかも残らない
- 案 H: PyPI の wheel（manylinux）を `fetchurl` で取って入れる。今の venv とバイナリが同じになるが、nix store の外の
  ライブラリを前提にした wheel を patchelf で直すことになり、ソースからのビルドより脆い
- 案 I: uv2nix、poetry2nix などのロックのツールを flake の入力に足す。依存 2 つのために入力が増え、`flake.lock` が
  変わる（examples/m0/dvc.yaml は `flake.lock` を deps にしているので、M0 のステージが再実行する）

固定した資産を `references/` に見せる方法:

- 案 J（採用）: ENDF は `nix build .#endf -o references/endf` で、store の木への symlink（out-link）にする。入力の TOML と
  dvc.yaml が指すパス（`references/endf/ENDF-B-VIII.0_neutrons/...`）が変わらず、中身が同じなら DVC の md5 も
  変わらない。out-link は GC の root になるので、`nix-collect-garbage` で消えない。ツール（NJOY、OpenMC、Python API）は
  `references/` に置かず、flake の `references` の devShell で `PATH` に入れる（スクリプトは `njoy`、`openmc` を
  `PATH` から呼ぶ）
- 案 K: 入力の TOML と dvc.yaml に store のパスを直接書く。nixpkgs や derivation が変わるたびにパスが変わり、すべての
  入力を書き換えることになる
- 案 L: store の木を `references/` にコピーする（`cp -r`）。容量が倍になり、コピーがどの derivation から来たかが残らない
- 案 M: ツールも `nix build -o references/njoy` などの out-link にする。GC の root になるが、使うのは `PATH` からなので、
  devShell で足りる

仕様書の PDF:

- 案 N（採用）: 再配布できないので、URL と sha256 を docs/specs.sha256（`sha256sum -c` の形）に書く。取る手順の記録が
  ないものは「URL 不明」と書く
- 案 O: nixpkgs の `requireFile` で固定出力にする。ハッシュは照らせるが、手で `nix-store --add-fixed` することになり、
  store と `references/` に同じ PDF が 2 つ置かれる。ビルドの入力ではないので、照らせれば足りる

## 決定と理由

案 A、D、F、J、N にする。

参照解の作り方（ADR-0031 から引き継ぐ）:

- 参照解は、ENDF/B-VIII.0 から手元の NJOY で作った連続エネルギーのデータで、同じ組成・同じ温度で OpenMC を解いて
  求める。核データの違いを除き、手法の誤差だけを取り出すためである。これはコード同士の検証（verification）であり、
  実測との比較による妥当性確認（validation）ではない。VERA の公開値（KENO-VI、ENDF/B-VII.0）とも比べるが、核データの
  差を含むので参考とする
- 参照解を作るツール（NJOY、OpenMC）は検証にだけ使い、本体の依存に含めない。flake の packages と `references` の
  devShell にだけ置き、`default` の devShell と Cargo のワークスペースには入れない
- 3D の全炉心の OpenMC は計算時間が大きいので作らない。Problem 5 は仕様書の実測値（妥当性確認）と KENO の値で比べる
- 参照解の組成は計算と同じにする（[ADR-0054](0054-model-vera-materials.md)）

flake.nix の packages:

| package                                                         | 中身                                                                                                                               |
| --------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `endf-neutrons`、`endf-thermal-scatt`、`endf-decay`、`endf-nfy` | NNDC の ENDF/B-VIII.0 の zip（`fetchurl`、URL と sha256）                                                                          |
| `endf`                                                          | 4 つの zip を展開した木（`ENDF-B-VIII.0_<part>/`）                                                                                 |
| `njoy2016`                                                      | NJOY2016 2016.79（タグのコミット）に docs/patches/njoy2016-acer-law0.patch を当て、GCC 13 と gfortran 13 でビルドしたもの          |
| `openmc`                                                        | OpenMC 0.15.2（タグのコミットとサブモジュール）を OpenMP つき、GCC 13（同梱の fmt が GCC 16 でコンパイルできない）でビルドしたもの |
| `openmc-python`                                                 | 検証のスクリプトの依存（numpy、h5py など）に、openmc と endf 0.1.12 の Python API を足した Python                                  |

- OpenMC の CMake は版を `git describe` か `.git_archival.txt` から読む。取ったソースには `.git` がなく、
  `.git_archival.txt` は `git archive` が埋めるテンプレートのままなので、タグのコミットと版を書いた形に置き換える
- openmc の Python のパッケージは `ipython` と `matplotlib` に依存すると書くが、使うのは対話と作図だけで、検証の
  スクリプトは使わないので外す（venv にも入っていなかった）
- `default` の devShell には dvc の実行ファイルだけを入れる。`pkgs.dvc` をそのまま入れると、Python の setup hook が
  dvc の依存を `PYTHONPATH` に足し、検証用の Python のパッケージより先に読まれる。`references` の devShell は
  `PYTHONPATH` を外し、外の shell の Python のパッケージを読まない
- ツールを足しても flake の入力は増やさない。`flake.lock` は M0 のステージの deps なので（ADR-0123）、変えると M0 が再実行する

`references/` の扱い:

- `references/` は、消しても作り直せる置き場にする。外から取るものとツールは flake から作る: ENDF は
  main の checkout で `nix build .#endf -o references/endf`（worktree から実行すると、GC の root が worktree の中の
  パスで登録され、worktree を消すと外れる）、ツールは `nix develop .#references`（`PATH` に `njoy`、`openmc`、OpenMC の
  Python API の入った `python3`）。仕様書の PDF は docs/specs.sha256 の URL から手で取り、`sha256sum -c` で照らす。
  ただし取る先の分からない PDF（「URL 不明」）は作り直せないので、`references/` から消さない
- ツールで作る参照の結果（NJOY の出力、OpenMC の HDF5 と statepoint）は、#0018 の次のブランチで DVC のステージにする。
  それまでは `references/` に手で置いたものを使う

Nix と DVC の役割:

- Nix は環境（外から取る入力、ツール、そのコンパイラとライブラリ）を固定する。何から作ったかは flake.nix と
  `flake.lock` の内容で決まり、store のパスがそれを表す
- DVC は計算のステージを駆動し、ステージの deps と outs を内容の md5 で追う（ADR-0123）。ENDF のファイルは、out-link を通した
  パスで deps にする。DVC は symlink をたどって中身の md5 を取るので、store のパスが変わっても中身が同じなら走らない

## 確認

flake の packages を作り、手で取ったもの・手でビルドしたもの（main の `references/`）と比べた（d52538b5 の
trials/pin-reference-assets/compare_tools.sh、Slurm の上で 16 スレッド、2026-10-11）。flake の側の NJOY と openmc.data は、
`.#endf`（nix store の読み出し専用のファイル）を読んだ。

- `nix build` は 4 つの zip を取り、NJOY、OpenMC、Python API をビルドして約 1 分 40 秒で終わった（r995）。`nix flake check`
  が通る
- `.#endf` の木は `references/endf` の展開した木と、4 つのディレクトリ（561、91、3825、35 ファイル）のすべてで
  `diff -r` が一致した。入力の TOML と dvc.yaml が指すファイルの中身は変わらない
- docs/specs.sha256 は `sha256sum -c` が 4 つとも通る。VERA と MOX の pin_cell は、書いた URL から取ったファイルの
  sha256 が一致した
- NJOY: references/njoy-runs の B10-acer（moder、reconr、broadr、パッチを通る acer）と O17-groupr（reconr、groupr）で、
  ACE（B10 の tape24、tape25）と GENDF（O17 のバイナリの tape23、テキストの tape24）はバイト単位で一致した。
  バイナリの ENDF と PENDF（B10 の tape21〜23、O17 の tape21、22）は、最初のレコードの 97〜155 バイト目の 6〜9 バイト
  だけが違う。その値はメモリのアドレスに見え（`0x7ffff7eac4bc` など）、初期化されない領域がそのまま書かれているのが候補である
  （NJOY のソースでは確かめていない）。ASLR の
  効いた手元では、同じ手でビルドした `njoy` を 2 回実行しても同じ位置が違った。ほかのバイトは一致した
- openmc.data（Python API）が NJOY を通して作る H-1 の 600 K の HDF5 は、手で取った ENDF から作ったものと、18 の
  データセットのすべてが一致した
- OpenMC: `openmc --version` は、どちらも 0.15.2、コミット e23760b0 を示した。C5G7（多群、2 万粒子、120 世代、
  種 1）の k は両方とも 1.186062 ± 0.000477、VERA 1B（連続エネルギー、同じ粒子数）は 1.181836 ± 0.000647 で、
  一致した。世代ごとの k と tally の差は相対で 1e-14 以下、組み合わせた k の推定の差は 3e-11 以下で、同じツールを 2 回
  実行した差と同じ大きさだった（スレッドの和の順序による丸めが候補。確かめていない）
- `default` の devShell の `PYTHONPATH` は空で、`sys.path` に dvc の依存が入らず、`dvc --version` は 3.67.1 を示す。
  `references` の devShell は、`default` の devShell の中から入っても、`sys.path` が自分の Python の環境だけになる

## 帰結

- main の `references/` には、手でビルドしたツール（`NJOY2016/`、`openmc/src`、`openmc/build`、`openmc/install`、
  `openmc/venv`）と、手で取った `endf/`（zip と展開した木）が残っている。`references/endf` を out-link に置き換えるのは、
  このブランチをマージしたあとに手で行う（ディレクトリをどけてから `nix build .#endf -o references/endf`）。中身が同じ
  ことは「確認」で確かめたので、M0 の `dvc.lock` は変わらない
- ツールは devShell の中だけにあり、GC の root がない。`nix-collect-garbage` のあとは、`nix develop .#references` が作り
  直す（同じ記述から作る。2026-10-11 の r995 で、4 つの zip の取得を含めすべての packages の `nix build` が約 1 分 40 秒）
- OpenMC の Python API の版（`openmc.__version__`）は、venv では 0.0.0 だった（原因は確かめていない）が、0.15.2 になる
- ADR-0031 はツールを `references/` に手でビルドすると決めていたので、この ADR が置き換える
- 次のブランチ（#0018）で参照の結果を DVC のステージにするとき、ツールを表す deps（flake.nix、`flake.lock`、
  docs/patches/）の並べ方を決める。flake.nix はこの ADR の packages のほか devShell も含むので、そのまま deps にすると
  関係のない変更でも走る
