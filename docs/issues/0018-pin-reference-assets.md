---
status: open
labels: [verification]
opened: 2026-10-08
closed:
milestone:
blocked-by: []
resolution:
---

# 検証用の資産（references/）を版とハッシュで固定する

## 問題

`references/`（git の管理外、約 10 GB）には、性格の違う3種類が混ざっている（比べるのに使う小さな結果は
examples/ にコミットしていて、対比のために表の最後に並べる）。手順は docs/environment.md にあるが、再現性は
手作業に頼っている。

| 中身                                                              | 今の扱い                                                                                                         |
| ----------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| 外から取る入力（ENDF/B-VIII.0、VERA・C5G7・MOX の仕様書の PDF）   | ENDF は flake の `fetchurl` で URL と sha256 ごと固定した。PDF は docs/specs.sha256 に URL と sha256（ADR-0124） |
| 道具（NJOY2016 2016.79 とパッチ、OpenMC 0.15.2）                  | flake の packages と `references` の devShell（ADR-0124）                                                        |
| 道具で作る参照の結果（NJOY の出力、OpenMC の HDF5 と statepoint） | 一部はスクリプトで作れるが、`njoy-runs/` の多くは作り方がリポジトリから辿れない                                  |
| （対比）比べるのに使う小さな結果                                  | examples/ にコミット                                                                                             |

参照の結果は、作り方の一部がリポジトリから辿れず、何がどこにあるかの一部はリポジトリの外のメモに頼っている。

一般には、取得のスクリプトと URL・sha256 の一覧、Nix の固定出力の取得（`fetchurl`）、Git LFS、DVC、submodule、
アーカイブ（DOI）などで扱う。仕様書の PDF は再配布できないので、URL とハッシュだけを記録する。

外からの入力と道具の固定、`references/` への見せ方は [ADR-0124](../adr/0124-pin-reference-assets-with-nix.md) で
決めた（ADR-0031 を置き換えた）。

関係する ADR: [ADR-0122](../adr/0122-run-trials-from-committed-inputs.md)、[ADR-0123](../adr/0123-drive-stages-with-dvc.md)、[ADR-0124](../adr/0124-pin-reference-assets-with-nix.md)、[ADR-0126](../adr/0126-share-references-across-worktrees-by-worktreelink.md)、[ADR-0127](../adr/0127-drive-reference-results-with-dvc.md)

## 完了の条件

- 済（ADR-0124）: 外からの入力を固定する。ENDF の zip は flake の `fetchurl` で URL と sha256 ごと、再配布できない
  仕様書の PDF は docs/specs.sha256 に URL と sha256
- 済（ADR-0124）: NJOY2016（2016.79 とパッチ）と OpenMC 0.15.2 とその Python API を、版とパッチごとビルドする flake の
  パッケージにする（nixpkgs にはどちらもない）
- 済（ADR-0127）: 決めた方式を ADR にする（外からの入力と道具は ADR-0124、参照の結果の作り方は ADR-0127）
- 一部済（ADR-0127）: 参照の結果（NJOY → OpenMC のライブラリ → OpenMC、xlc との比較）を DVC の段階にする。
  examples/vera/dvc.yaml で、OpenMC のライブラリ（核種ごと）と VERA 1B（OpenMC → 抽出、xs → lattice → 比較）を
  端から端まで通した。残り（同じ形で足す）:
  - VERA 1C〜2J の OpenMC と抽出（2B〜2J はピン出力の `pins.msgpack`）、xlc の lattice、比較（`compare_pins.py` を
    含む）の段階。ADR-0121 の `p2/openmc-<問題>.msgpack` の置き場を、抽出の結果のディレクトリに合わせるかを決める
  - MOX（examples/mox/openmc_model.py）、C5G7（examples/c5g7/openmc_model.py）、#0005 の均質な媒質
    （`openmc_homogeneous.py`）、`openmc_balance.py`、`openmc_mg.py` を、計算の条件を入力の TOML で受ける形に直し、
    段階にする（どれもまだ位置の引数と既定値で受ける）
- 参照の結果は、作るスクリプト（入力、コマンド、粒子数、乱数の種、道具の版）をリポジトリに置き、`references/` は
  いつ消しても作り直せる作業の置き場にする。`njoy-runs/` も作るスクリプトにまとめ直す。すべてを段階にしたら、
  `references/openmc/` の手で作った作業ディレクトリ（`endfb80`、`vera-*`、`work-*` など）と手でビルドした道具を消す

## 経緯

- 2026-10-08: ユーザーとの議論で一般の扱いを調べ、Nix の flake で固定する方針を推奨した。ADR と検証の記録の整理の
  後に回す（ユーザーの判断）
- 2026-10-11: ブランチ feat/pin-reference-assets で、ENDF の zip、NJOY2016、OpenMC とその Python API を flake の
  packages にし、`references/endf` を `nix build .#endf` の out-link にすると決めた（ADR-0124）。flake の NJOY と
  OpenMC が手でビルドしたものと同じ結果を出すことを、ADR-0124 の「確認」で確かめた。main の `references/endf` を
  out-link に置き換えるのは、マージのあとに手で行う。参照の結果を DVC の段階にするのは次のブランチ
- 2026-10-11: `references/parcs/theory.pdf` は PARCS の理論マニュアルでなく、NRC のサーバーが返した「Access Denied」の
  HTML（402 バイト）だった。どの記録も引いていないので docs/specs.sha256 に入れていない。`references/papers/` の
  論文も記録から引かれていないので入れていない
- 2026-10-11: 仕様書から抜いたテキスト（`references/{vera,c5g7,mox}/spec.txt`）は、記録が引き（examples/vera/compare_keno.py、
  #0058、examples/m3-nuclides/verification.md）、作る手順（pdftotext の選択肢）がどこにもない。作るスクリプトにするか、
  ハッシュを記録するかを、参照の結果と一緒に決める。C5G7 と MOX の仕様の PDF は取る先が分からないので、消すと
  作り直せない
- 2026-10-11: worktree の `references` は main の `references/` への symlink なので、`references/openmc` に `cd` して
  `../..` で flake とスクリプトを指すと、worktree でなく main のものが使われる（`nix flake metadata ../..` が main を
  示した）。docs/environment.md と examples/c5g7/verification.md の手順はリポジトリの直下から走らせる形にした。
  examples/mox/verification.md と examples/vera/verification.md の手順（「references/openmc で」）はまだこの形で、
  参照の結果を DVC の段階にするときに直す
