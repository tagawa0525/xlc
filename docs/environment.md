# 検証用の環境

xlc の計算は、NJOY（核データ処理）と OpenMC（モンテカルロ計算）の結果と比べて検証する。道具（NJOY、OpenMC と
その Python API）と評価済み核データは、`flake.nix` の packages で版とハッシュごと固定する（ADR-0124）。
計算の置き場 `references/` は git で追跡しない。そのうち外から取る入力と道具は、消しても下の手順で作り直せる
（取る先の分からない仕様書の PDF を除く。道具で作る参照の結果を作り直す段階は #0018 で整える）。
direnv を使う場合は、リポジトリの直下に `use flake` と書いた `.envrc` を置き（手元の設定として git では追跡しない）、
一度 `direnv allow` すると、`default` の devShell に自動で入る。

| devShell     | 用途                                                                                    |
| ------------ | --------------------------------------------------------------------------------------- |
| `default`    | 開発（Rust）、段階の駆動（DVC）と、検証用の Python スクリプト（numpy, h5py など）       |
| `references` | 参照の計算: `njoy`、`openmc` と、OpenMC の Python API（openmc、endf）の入った `python3` |

## 段階の駆動（DVC）

入力から出力を作る段階は DVC（`default` の devShell に入っている）で走らせる（ADR-0123）。段階は
`examples/<問題>/dvc.yaml` に書き、`dvc repro examples/<問題>/dvc.yaml` で、入力かコードが変わった段階だけを走らせる。
各段階の入力・コード・出力の md5 は同じディレクトリの `dvc.lock` に書かれ、コミットする。出力は
`references/dvc-cache`（`.dvc/config`）にキャッシュし、`dvc checkout` で戻せる。worktree では、
worktree の直下で `ln -s "$(git rev-parse --path-format=absolute --git-common-dir)/../references" references` とし、main の
`references/` を指して、キャッシュと核データを共有する。

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

## VERA の参照解

VERA Problem 1B〜1D, 2B〜2J の参照解は、ENDF/B-VIII.0 から手元の NJOY で作った連続エネルギーの
データで OpenMC を動かして作る（`examples/vera/`）。データは `references/openmc/endfb80/` に置く:

```sh
nix develop .#references -c python3 examples/vera/openmc_library.py references/endf references/openmc/endfb80 8
nix develop .#references -c sh -c \
  'OMP_NUM_THREADS=16 python3 examples/vera/openmc_model.py references/openmc/endfb80/cross_sections.xml 2B \
     references/openmc/work-2B 1000000 400 100 1'
```

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
