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

| 中身                                                              | 今の扱い                                                                                                             |
| ----------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| 外から取る入力（ENDF/B-VIII.0、VERA・C5G7・MOX の仕様書の PDF）   | ENDF は手順の `curl` で、ハッシュを確かめない。PDF は取る手順の記録がない（URL があるのは MOX の pin_cell.pdf だけ） |
| 道具（NJOY2016 2016.79 とパッチ、OpenMC 0.15.2）                  | 版を指定した `git clone` と手でのビルド（ビルド環境は flake の `references` の devShell）                            |
| 道具で作る参照の結果（NJOY の出力、OpenMC の HDF5 と statepoint） | 一部はスクリプトで作れるが、`njoy-runs/` の多くは作り方がリポジトリから辿れない                                      |
| （対比）比べるのに使う小さな結果                                  | examples/ にコミット                                                                                                 |

配布物が差し替わっても気づけず、何がどこにあるかの一部はリポジトリの外のメモに頼っている。

一般には、取得のスクリプトと URL・sha256 の一覧、Nix の固定出力の取得（`fetchurl`）、Git LFS、DVC、submodule、
アーカイブ（DOI）などで扱う。仕様書の PDF は再配布できないので、URL とハッシュだけを記録する。

道具の作り方は [ADR-0031](../adr/0031-verify-against-openmc-with-same-data.md) が「ソースから `references/` に
ビルドする」と決めている。

関係する ADR: [ADR-0122](../adr/0122-run-trials-from-committed-inputs.md)

## 完了の条件

- 外からの入力（ENDF の zip、仕様書の PDF）を flake の `fetchurl` で URL と sha256 ごと固定する
- NJOY2016（2016.79 とパッチ）と OpenMC 0.15.2 を、版とパッチごとビルドする flake のパッケージにする（nixpkgs に
  あれば使えるか確かめる）
- 参照の結果は、作るスクリプト（入力、コマンド、粒子数、乱数の種、道具の版）をリポジトリに置き、`references/` は
  いつ消しても作り直せる作業の置き場にする。`njoy-runs/` も作るスクリプトにまとめ直す
- 決めた方式を ADR にする。道具の作り方を変えるなら、その ADR は ADR-0031 を置き換える（supersedes）

## 経緯

- 2026-10-08: ユーザーとの議論で一般の扱いを調べ、Nix の flake で固定する方針を推奨した。ADR と検証の記録の整理の
  後に回す（ユーザーの判断）
