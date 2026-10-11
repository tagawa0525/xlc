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

| 中身                                                                | 今の扱い                                                                                                         |
| ------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| 外から取る入力（ENDF/B-VIII.0、VERA・C5G7・MOX の仕様書の PDF）     | ENDF は flake の `fetchurl` で URL と sha256 ごと固定した。PDF は docs/specs.sha256 に URL と sha256（ADR-0124） |
| ツール（NJOY2016 2016.79 とパッチ、OpenMC 0.15.2）                  | flake の packages と `references` の devShell（ADR-0124）                                                        |
| ツールで作る参照の結果（NJOY の出力、OpenMC の HDF5 と statepoint） | 一部はスクリプトで作れるが、`njoy-runs/` の多くは作り方がリポジトリから辿れない                                  |
| （対比）比べるのに使う小さな結果                                    | examples/ にコミット                                                                                             |

参照の結果は、作り方の一部がリポジトリから辿れず、何がどこにあるかの一部はリポジトリの外のメモに頼っている。

一般には、取得のスクリプトと URL・sha256 の一覧、Nix の固定出力の取得（`fetchurl`）、Git LFS、DVC、submodule、
アーカイブ（DOI）などで扱う。仕様書の PDF は再配布できないので、URL とハッシュだけを記録する。

外からの入力とツールの固定、`references/` への見せ方は [ADR-0124](../adr/0124-pin-reference-assets-with-nix.md) で
決めた（ADR-0031 を置き換えた）。

関係する ADR: [ADR-0122](../adr/0122-run-trials-from-committed-inputs.md)、[ADR-0123](../adr/0123-drive-stages-with-dvc.md)、[ADR-0124](../adr/0124-pin-reference-assets-with-nix.md)、[ADR-0126](../adr/0126-share-references-across-worktrees-by-worktreelink.md)、[ADR-0130](../adr/0130-drive-reference-results-with-dvc.md)

## 完了の条件

- 済（ADR-0124）: 外からの入力を固定する。ENDF の zip は flake の `fetchurl` で URL と sha256 ごと、再配布できない
  仕様書の PDF は docs/specs.sha256 に URL と sha256
- 済（ADR-0124）: NJOY2016（2016.79 とパッチ）と OpenMC 0.15.2 とその Python API を、版とパッチごとビルドする flake の
  パッケージにする（nixpkgs にはどちらもない）
- 済（ADR-0130）: 決めた方式を ADR にする（外からの入力とツールは ADR-0124、参照の結果の作り方は ADR-0130）
- 一部済（ADR-0130）: 参照の結果（NJOY → OpenMC のライブラリ → OpenMC、xlc との比較）を DVC のステージにする。
  examples/vera/dvc.yaml で、OpenMC のライブラリ（核種ごと）と VERA 1B（OpenMC → 抽出、xs → lattice → 比較）を
  端から端まで通した。残り（同じ形で足す）:
  - VERA 1C〜2K の OpenMC と抽出（2B〜2K はピン出力の `pins.msgpack`）、xlc の lattice、比較（`compare_pins.py` を
    含む）の段階。ADR-0121 の `p2/openmc-<問題>.msgpack` の置き場を、抽出の結果のディレクトリに合わせるかを決める
  - 燃料の実効断面積の tally（`effective-xs-groups`）を使う段階。xlc のライブラリを deps にすると xs のコードを変えるたびに
    OpenMC が再実行するので、群構造だけを入力に書くかを決める（ADR-0130 の「帰結」）
  - OpenMC の段階を「入力とスクリプト → model.xml」と「model.xml とライブラリ → statepoint」に分け、スクリプトを作り直しても
    model.xml が同じなら OpenMC が再実行しないようにする（ライブラリの段階も同じく分ける）。分けたあと 1B を一度走らせ直せば、
    記録の「実行したコミットと分岐点からの復元」（ADR-0130）は要らなくなる
  - `openmc-library.toml` の H-1 の温度は、冷却材と被覆管の温度と S(α,β) の温度を合わせたものを手でコピーしている（main の
    `openmc_library.py` は和集合を計算していた）。食い違いを止める検査を足す
  - `compare_k.py` は lattice の出力が k.json と同じ問題のものかを確かめない。1C 以降の比較の段階を足すときに、問題を
    照らす
  - MOX（examples/mox/openmc_model.py）、C5G7（examples/c5g7/openmc_model.py）、#0005 の均質な媒質
    （`openmc_homogeneous.py`）、`openmc_balance.py`、`openmc_mg.py` を、計算の条件を入力の TOML で受ける形に直し、
    段階にする（どれもまだ位置の引数と既定値で受ける）
- 参照の結果は、作るスクリプト（入力、コマンド、粒子数、乱数のシード、ツールの版）をリポジトリに置き、`references/` は
  いつ消しても作り直せる作業の置き場にする。`njoy-runs/` も作るスクリプトにまとめ直す。すべてを段階にしたら、
  `references/openmc/` の手で作った作業ディレクトリ（`endfb80`、`vera-*`、`work-*` など）と手でビルドしたツールを消す

## 経緯

- 2026-10-08: ユーザーとの議論で一般の扱いを調べ、Nix の flake で固定する方針を推奨した。ADR と検証の記録の整理の
  後に回す（ユーザーの判断）
- 2026-10-11: ブランチ feat/pin-reference-assets で、ENDF の zip、NJOY2016、OpenMC とその Python API を flake の
  packages にし、`references/endf` を `nix build .#endf` の out-link にすると決めた（ADR-0124）。flake の NJOY と
  OpenMC が手でビルドしたものと同じ結果を出すことを、ADR-0124 の「確認」で確かめた。main の `references/endf` を
  out-link に置き換えるのは、マージのあとに手で行う。参照の結果を DVC のステージにするのは次のブランチ
- 2026-10-11: `references/parcs/theory.pdf` は PARCS の理論マニュアルでなく、NRC のサーバーが返した「Access Denied」の
  HTML（402 バイト）だった。どの記録も引いていないので docs/specs.sha256 に入れていない。`references/papers/` の
  論文も記録から引かれていないので入れていない
- 2026-10-11: 仕様書から抜いたテキスト（`references/{vera,c5g7,mox}/spec.txt`）は、記録が引き（examples/vera/compare_keno.py、
  #0058、examples/m3-nuclides/verification.md）、作る手順（pdftotext の選択肢）がどこにもない。作るスクリプトにするか、
  ハッシュを記録するかを、参照の結果と一緒に決める。C5G7 と MOX の仕様の PDF は取る先が分からないので、消すと
  作り直せない
- 2026-10-11: worktree の `references` は main の `references/` への symlink なので、`references/openmc` に `cd` して
  `../..` で flake とスクリプトを指すと、worktree でなく main のものが使われる（`nix flake metadata ../..` が main を
  示した）。docs/environment.md と examples/c5g7/verification.md の手順はリポジトリの直下から実行する形にした。
  examples/mox/verification.md と examples/vera/verification.md の手順（「references/openmc で」）はまだこの形で、
  参照の結果を DVC のステージにするときに直す
- 2026-10-11: ブランチ feat/reference-stages で、参照の結果を DVC のステージにする方式を決め（ADR-0130）、examples/vera/dvc.yaml
  で OpenMC のライブラリ（核種ごと）と VERA 1B を端から端まで通した。OpenMC 1B の k は手で実行した値と一致し、ライブラリは
  `references/openmc/endfb80-565`（同じ温度で手で作ったもの）とバイト単位で一致した。examples/vera/verification.md と examples/mox/verification.md の手順は
  リポジトリの直下から実行する形に直した
- 2026-10-12: 手でビルドした `references/NJOY2016` と `references/openmc/{src,build,install,venv}` は、参照の結果をすべて段階に
  するまで残す。消すのは派生物（build、install、venv）だけにし、git clone（src）は原本なので残すかをユーザーに確かめる。
  `references/openmc/endfb80-565`（565 K を含む endfb80 の上位集合）で `endfb80` を置き換えるかは決めていない
