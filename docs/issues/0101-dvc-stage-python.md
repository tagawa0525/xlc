---
status: closed
labels: [verification]
opened: 2026-10-13
closed: 2026-10-13
milestone:
blocked-by: []
resolution: コミット 2733b747
---

# DVC のステージの python3 が検証用の Python でなく dvc のラッパーの素の python3 になる

## 問題

default の devShell の `dvc`（flake.nix。pkgs.dvc の `bin/dvc` へのリンク）は、nixpkgs の makeWrapper のラッパーで、
`PATH` の先頭に dvc の依存の Python のパッケージの bin、nettle、git-minimal と素の python3（`python3-3.14.7/bin`）を足してから
本体を呼ぶ。段階の cmd は
この `PATH` で走るので、`python3 ../compare_p5.py ...` のような段階は、devShell の検証用の Python（`withPackages`）でなく
素の python3 で走り、`ModuleNotFoundError: No module named 'msgpack'` で止まる。

feat/0086-periphery-constants で examples/vera/p5/dvc.yaml の比較の段階を実行して見つけた（2026-10-13。手元の
`dvc repro`、`nix develop -c dvc repro`、Slurm のジョブ 256 のどれでも同じ）。同じ形の段階は 2026-10-12 のジョブ 247
（同じブランチ、`dvc repro examples/vera/p5/dvc.yaml`）では通った。その違いの理由は確かめていない（候補は、そのときの
シェルの環境に検証用の Python のパッケージを指す PYTHONPATH が残っていたこと）。

`nix develop .#references -c` で包んだ段階（NJOY、OpenMC）は、包んだ shell の python3 を使うので、この影響を受けない。

## 修正方法

devShell の `dvc` を、ラッパーの中の `.dvc-wrapped` を直接呼ぶ小さなスクリプトにする。`.dvc-wrapped` は自分の依存の site の
ディレクトリを自分で足すので、`PATH` を変えなくても動く。ラッパーが付けていた `PYTHONNOUSERSITE=true` は残す。
`.dvc-wrapped` は nixpkgs の Python のラッパーの中の名前なので、無ければ devShell のビルドで止める。git-minimal は
`PATH` から外れるが、空の環境で git を `PATH` に置かずに `dvc status` と `dvc repro --dry` が動くことを確かめた（2026-10-13、2733b747 の devShell）。
flake の check `dvc-stage-python` が、devShell の dvc で msgpack と numpy を読む段階を実行する（nixpkgs のラッパーの dvc では
ModuleNotFoundError で止まることを確かめた）。

ADR-0124 の「`default` の devShell には dvc の実行ファイルだけを入れる」は変えない（入れるのは今も実行ファイル 1 つで、
dvc の setup hook は shell に入らない）。default の devShell が外の `PYTHONPATH` を外さないことは #0103 に分けた。

採らなかった案:

- 段階の cmd を `nix develop ..#default -c` で包む: すべての段階の cmd と deps が変わり、段階ごとに shell の評価の時間が
  かかる。原因（dvc の `PATH`）が残る
- pkgs.dvc を `overridePythonAttrs` の `makeWrapperArgs` で `--suffix PATH` にする: dvc を作り直すことになり、nixpkgs の
  キャッシュが使えない

関係する ADR: なし

## 完了の条件

default の devShell の `dvc repro` で、`python3` を呼ぶ段階が検証用の Python で走る。

## 経緯

- 2026-10-13: 見つけて、devShell の `dvc` を `.dvc-wrapped` を呼ぶスクリプトにした（2733b747）。`.dvc-wrapped` が無ければ
  ビルドで止め（cb97b8f0）、flake の check を足した（847af602）
