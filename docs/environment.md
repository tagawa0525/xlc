# 検証用の環境

xlc の計算は、NJOY（核データ処理）と OpenMC（モンテカルロ計算）の結果と比べて検証する。これらのツールと
評価済み核データは git で追跡しない `references/` に置く。環境は `flake.nix` の devShell で固定する。
direnv を使う場合は、リポジトリの直下に `use flake` と書いた `.envrc` を置き（手元の設定として git では追跡しない）、
一度 `direnv allow` すると、`default` の devShell に自動で入る。

| devShell     | 用途                                                                            |
| ------------ | ------------------------------------------------------------------------------- |
| `default`    | 開発（Rust）と、検証用の Python スクリプト（numpy, h5py など）                  |
| `references` | NJOY と OpenMC のビルド（GCC 13。OpenMC 0.15.2 は GCC 16 でコンパイルできない） |

## 評価済み核データ（ENDF/B-VIII.0）

```sh
mkdir -p references/endf && cd references/endf
for s in neutrons thermal_scatt decay nfy; do
  curl -LO https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_$s.zip
  unzip -q ENDF-B-VIII.0_$s.zip
done
```

## NJOY2016（2016.79）

```sh
git clone --branch 2016.79 https://github.com/njoy/NJOY2016.git references/NJOY2016
git -C references/NJOY2016 apply ../../docs/patches/njoy2016-acer-law0.patch
nix develop .#references -c sh -c \
  'cmake -S references/NJOY2016 -B references/NJOY2016/build -DCMAKE_BUILD_TYPE=Release &&
   cmake --build references/NJOY2016/build -j'
```

NJOY の入力と出力は `references/njoy-runs/` に置き、`crates/xs/examples/compare_njoy_*` で比べる。

`docs/patches/njoy2016-acer-law0.patch` は ACER のバグを直す。MF6 の生成物の分布が LAW=0（分布が
与えられていない）だと、ACER は荷電粒子の生成データにその項目を LNW, LAW, IDAT の3語だけで書くが、
書き出し（`change`）は分布の見出しが続くものとして読み進め、止まる（B-10 の MT700 の α 粒子）。
パッチは LAW=0 の項目でその3語のあと次の反応に進む。荷電粒子の生成データだけにかかわり、
中性子の断面積は変わらない（OpenMC は荷電粒子の生成データを読まない）。

## OpenMC（0.15.2）

```sh
mkdir -p references/openmc && cd references/openmc
git clone --depth 1 --branch v0.15.2 --recurse-submodules --shallow-submodules \
  https://github.com/openmc-dev/openmc.git src
nix develop ../..#references -c sh -c \
  'cmake -S src -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PWD/install -DOPENMC_USE_OPENMP=on &&
   cmake --build build -j && cmake --install build &&
   python3 -m venv --system-site-packages venv &&
   venv/bin/pip install --no-deps --no-build-isolation ./src && venv/bin/pip install endf'
```

実行するときも、ビルドと同じ `references` の devShell に入る（OpenMC は GCC 13 の OpenMP（libgomp）に
リンクしていて、外では読み込めずに止まる）。`install/bin` を `PATH` に入れ、`venv/bin/python` でスクリプトを
動かす。たとえば C5G7 の参照解は:

```sh
cd references/openmc
nix develop ../..#references -c sh -c \
  'PATH=$PWD/install/bin:$PATH OMP_NUM_THREADS=16 venv/bin/python ../../examples/c5g7/openmc_model.py work 1000000 400 100 2'
```

## VERA の参照解

VERA Problem 1B〜1D, 2B〜2H の参照解は、ENDF/B-VIII.0 から手元の NJOY で作った連続エネルギーの
データで OpenMC を動かして作る（`examples/vera/`）。データは `references/openmc/endfb80/` に置く:

```sh
cd references/openmc
nix develop ../..#references -c sh -c \
  'PATH=$PWD/../NJOY2016/build:$PATH venv/bin/python ../../examples/vera/openmc_library.py ../endf endfb80 8'
nix develop ../..#references -c sh -c \
  'PATH=$PWD/install/bin:$PATH OMP_NUM_THREADS=16 venv/bin/python ../../examples/vera/openmc_model.py \
     endfb80/cross_sections.xml 2B work-2B 1000000 400 100 1'
```

同じ多群定数で OpenMC の多群モードと比べる（examples/vera/verification.md の「同じ多群定数での OpenMC との比較」）には、
負にならない散乱核のライブラリと、そこから多群定数を作ってから解く。格子の入力は `library` をそのライブラリに、
`self-shielding = false` と `scattering-order` を同じ次数にしたもの:

```sh
cd references/openmc
nix develop ../..#references -c sh -c \
  'venv/bin/python ../../examples/vera/positive_library.py ../../examples/vera/vera.msgpack positive-P3.msgpack 3'
cargo run --release -p xlc-lattice --example mg_materials -- <格子の入力> mg-2B.msgpack
nix develop ../..#references -c sh -c \
  'PATH=$PWD/install/bin:$PATH OMP_NUM_THREADS=16 venv/bin/python \
     ../../examples/vera/openmc_mg.py mg-2B.msgpack 2B work-mg-2B 1000000 300 50 1'
```
