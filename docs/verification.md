# 検証用の環境

xlc の計算は、NJOY（核データ処理）と OpenMC（モンテカルロ計算）の結果と比べて検証する。これらの道具と
評価済み核データは git で追跡しない `references/` に置く。環境は `flake.nix` の devShell で固定する。

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
nix develop .#references -c sh -c \
  'cmake -S references/NJOY2016 -B references/NJOY2016/build -DCMAKE_BUILD_TYPE=Release &&
   cmake --build references/NJOY2016/build -j'
```

NJOY の入力と出力は `references/njoy-runs/` に置き、`crates/xs/examples/compare_njoy_*` で比べる。

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

実行するときは `install/bin` を `PATH` に入れ、`venv/bin/python` でスクリプトを動かす。たとえば C5G7 の参照解は:

```sh
cd references/openmc
PATH=$PWD/install/bin:$PATH OMP_NUM_THREADS=16 \
  venv/bin/python ../../examples/c5g7/openmc_model.py work 200000 400 100
```
