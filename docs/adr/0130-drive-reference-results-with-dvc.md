---
status: accepted
date: 2026-10-11
supersedes: []
superseded-by:
issues: ["#0018"]
requires: [ADR-0122, ADR-0123, ADR-0124, ADR-0126]
---

# ADR-0130: 参照の結果は DVC の段階で作り、OpenMC の段階は references の devShell で包む

## 背景

参照の結果（NJOY で作る OpenMC の連続エネルギーのライブラリ、OpenMC の問題ごとの結果、xlc との比較）は、
docs/environment.md の手順で人が走らせ、`references/openmc/` の手で名付けた作業ディレクトリ（`endfb80`、`vera-1B`、
`work-2F` など）に置いてきた。どの結果がどのスクリプトと入力から作られたかは、記録のコミットの hash と手順からしか
辿れない（#0018）。

前提になる決定:

- 計算の段階は DVC で駆動し、`examples/<問題>/dvc.yaml` に素のコマンドで書き、`dvc.lock` をコミットする。OpenMC の
  statepoint は outs にしてキャッシュし、比べる小さな結果は別の段階で抽出する（ADR-0123）
- 道具（`njoy`、`openmc`、OpenMC の Python API の入った `python3`）は flake の `references` の devShell に、核データは
  `references/endf`（`nix build .#endf` の out-link）に固定した。`default` の devShell には Rust の道具と dvc がある
  （ADR-0124）
- 記録に書く計算の条件は入力のファイルに書き、本番のコードに与える（ADR-0122）

OpenMC のスクリプト（`openmc_library.py`、`openmc_model.py`、`openmc_homogeneous.py`）は、粒子数、世代、捨てる世代、
乱数の種を位置の引数で受け、省くと既定値を使っていた。`openmc_library.py` は核種と温度をスクリプトの中の定数に持ち、
出力のディレクトリにある HDF5 を使い回していた。

判断の軸:

- `references/` を消しても `dvc repro` で作り直せること（ADR-0124 の「作り直せる置き場」を参照の結果にも及ぼす）
- 道具の版が変わったときに走り直し、関係のない変更（flake.nix の別の devShell、文書）では走らないこと
- 重い段階（ライブラリ、OpenMC）を、関係のない入力の変更で走らせ直さないこと
- 計算の条件が入力のファイルにあり、既定値で黙って走らないこと

## 検討した案

段階の道具の環境（Rust の段階は `default` の devShell、OpenMC の段階は `references` の devShell が要る）:

- 案 A（採用）: `dvc repro` は `default` の devShell で走らせ、OpenMC と NJOY を使う段階の cmd を
  `nix develop ../..#references -c` で包む。段階ごとに環境が cmd に書かれ、`dvc.lock` に残る。`nix develop` の評価が
  段階ごとに 1 秒ほどかかる
- 案 B: `dvc repro` 全体を `references` の devShell の中で走らせる。`references` の devShell には cargo も dvc もないので、
  Rust の段階が走らない。両方を足した devShell を作ると、`default` と `references` の 2 つの `python3` が `PATH` で
  ぶつかる（ADR-0124 が `PYTHONPATH` を外した理由と同じ）
- 案 C: 段階を問題ごとに 2 つの dvc.yaml に分け、片方を `references` の devShell で走らせる。xlc と OpenMC を比べる段階が
  両方の出力にまたがり、どちらの devShell で走らせるかを人が覚えることになる

道具を表す deps（ADR-0124 の「結果」で次のブランチに回したもの）:

- 案 D（採用）: `references` の devShell の nix store のパスを書き出す段階（`always_changed`）を置き、その出力を道具の
  deps にする。store のパスは、道具の版、パッチ、コンパイラ、依存のライブラリ、devShell の中身のどれが変わっても変わり、
  flake.nix のほかの部分（`default` の devShell など）や worktree のパスでは変わらない（main と worktree で同じパスを
  示した。2026-10-11）。評価だけでビルドしない
- 案 E: flake.nix、`flake.lock`、docs/patches/ をそのまま deps にする。`default` の devShell に道具を足しても、参照の
  ライブラリと OpenMC（1 問題 16〜37 分）が走り直す
- 案 F: 道具を deps に書かない。道具の版を上げても走らず、古い結果が最新とみなされる

ライブラリの段階の粒度:

- 案 G（採用）: 核種ごと（中性子の核種と熱中性子散乱則のそれぞれ）を 1 段階にし（`foreach`）、`cross_sections.xml` を
  書く段階を別に置く。核種と温度はライブラリの入力（`openmc-library.toml`）に書き、各段階はその核種の項目だけを
  `params` で deps にする。核種を足すと足した段階だけが走る（これまでスクリプトが HDF5 を使い回していたことを、DVC の
  判定に置き換える）。DVC は段階を逐次に走らせるので、`ProcessPoolExecutor` で 8 並列にしていたより壁時計は長い
- 案 H: ライブラリ全体を 1 段階にする。並列に作れるが、核種を 1 つ足しても全部を作り直すか、スクリプトが DVC の外で
  出力を使い回すことになる。DVC は段階の前に outs を消すので、使い回しは効かない

入力の TOML を deps にする単位: ADR-0123 は入力の TOML のファイル全体を deps にすると決めた。ライブラリの核種ごとの
段階は、入力のうちその核種の項目だけを読むので、項目を `params` にする（ほかの核種の変更で走らせない）。
`cross_sections.xml` を書く段階と、OpenMC と比べる段階は、ADR-0123 のとおりファイル全体を deps にする。

OpenMC の段階が読むライブラリの deps:

- 案 I（採用）: ライブラリのディレクトリ全体。どの核種が変わっても走るが、漏れがない
- 案 J: 問題が使う核種の HDF5 だけを並べる。核種を足しても走らないが、並びは組成の写しになり、漏れると古い結果が最新と
  みなされる（ADR-0123 の「結果」の deps の写しと同じ危うさ）。`cross_sections.xml` は核種を足すと変わるので、
  これも deps から外すことになる

既存の出力（`references/openmc/endfb80`、`vera-1B` など）の扱い:

- 案 K（採用）: `dvc repro` で走らせ直す。`dvc.lock` が、記録した cmd と deps から実際に作った出力を指す。走らせ直した
  結果を既存の結果と比べることが、道具と段階の確認になる
- 案 L: 既存の出力を outs の場所に写し、`dvc commit` で取り込む。計算は要らないが、既存の出力は手でビルドした道具と
  別の引数で作ったもので、`dvc.lock` が、記録した cmd と deps から作っていない出力を指す

statepoint と抽出の分け方は ADR-0123 の案 M のとおりにする。

## 決定と理由

案 A、D、G、I、K にする。

- 段階は `examples/vera/dvc.yaml` に書く。ライブラリは VERA と MOX（examples/mox）が使うが、作るスクリプトと入力が
  examples/vera にあるので、そこに置く
- `dvc repro` は `default` の devShell で走らせる。NJOY と OpenMC を使う段階（ライブラリ、OpenMC、statepoint からの
  抽出）は cmd を `nix develop ../..#references -c` で包み、deps に `references-shell.txt`（`references` の devShell の
  store のパス。`cache: false` でコミットする）を入れる。Rust の段階の deps は ADR-0123 のまま
- 計算の条件は入力の TOML に書き、スクリプトは入力のパスと、出力のパスだけを引数で受ける。入力にないキーや、知らない
  キーがあればエラーで止まる（既定値を持たない）。入力のパスは入力のファイルからの相対で読む
  - ライブラリ（`openmc-library.toml`）: 中性子の核種ごとに ENDF のファイルと温度、熱中性子散乱則ごとに中性子と熱の
    ENDF のファイルと温度
  - OpenMC の問題（`p1/openmc-1B.toml` など）: 問題、`cross_sections.xml`、粒子数、世代、捨てる世代、乱数の種、
    スレッド数。燃料の核種ごとの実効断面積の tally は、群構造を読む xlc のライブラリを書いたときだけ足す
- OpenMC の段階は作業ディレクトリ（statepoint、summary、model.xml）を outs にする。抽出の段階は statepoint から k
  （`k.json`）と、問題に応じてピン出力と実効断面積を、結果のディレクトリに書く
- xlc の xs の段階は、ENDF の中性子と熱中性子散乱則のディレクトリを deps にする。ENDF は sha256 で固定した
  （ADR-0124）ので変わらず、入力の TOML の 55 核種のファイルを dvc.yaml に写さずに済む

## 確認

examples/vera/dvc.yaml の全段階を、`default` の devShell の `dvc repro examples/vera/dvc.yaml` で Slurm の上で走らせた
（2026-10-11、全体で 57 分）。走らせたコミットは 53d92a76（rebase の remap で付け替わった hash。以後も付け替わる）、そのときの
分岐点は main の 1aab726a。`dvc.lock` の openmc-* と xs の段階の deps は、その状態の md5 である。復元は、分岐点にそこまでの
ブランチのコミットを当て直し、`dvc status` で md5 が合うことを確かめてから行う（衝突を手で解いたコミットは差分が変わり
うるので、md5 で確かめる）。lattice-1B（数秒）は rebase のあとに 6c0d133c（分岐点 635563b2）で
走らせ直し、その lock をコミットした。ライブラリの温度は、main が 2A と 2Q のために足した 565 K と S(α,β) の温度を
含む（`openmc-library.toml`）。

- OpenMC 1B の k は 1.181783 ± 0.000063 で、examples/vera/verification.md の「参照解」の 1B（869bd88f に手で走らせた
  `references/openmc/vera-1B`、同じ粒子数・世代・種）と一致した。k の推定の差は相対 4e-14、世代ごとの k も最初の世代から
  同じで、ADR-0124 の「確認」の、同じ道具を 2 回走らせた差と同じ大きさである。1B は最も近い温度のモードで 600 K の表だけを
  使うので、ライブラリに 565 K などの温度が増えても k は変わらない。model.xml は、手で走らせたものと `cross_sections` の
  パスだけが違った
- 核種ごとの段階で作ったライブラリの 62 の HDF5（中性子 61、熱中性子散乱則 1）は、手で作った
  `references/openmc/endfb80-565`（`openmc_library.py` に温度を引数で渡して作った、同じ温度のもの）とバイト単位で一致した（`cmp`）。
  rebase の前に 600 K だけの温度で走らせたとき（2 回）は、`references/openmc/endfb80` とバイト単位で一致した
- xlc の 1B の k∞ は 1.184735 で、記録の 64e069de の値と同じ。`compare-1B` は +295 ± 6 pcm を書いた。rebase のあとに
  6c0d133c で走らせ直した lattice-1B も k∞ 1.184735 で、`compare-1B` の入力の md5 が変わらないので走らなかった
- 走らせたあとにもう一度 `dvc repro` すると、`references-shell`（評価だけ）のほかは走らず、`references-shell.txt` も
  変わらなかった（11 s）。`references` の devShell の store のパスは main と worktree で同じだった
- 入力のキーが足りない、または知らないキーがあると、`openmc_model.py` が両方を挙げて止まることを、手元の小さな入力で
  確かめた（trials/ で、記録しない）
- 時間（Slurm、32 コア）: ライブラリの 63 段階が約 30 分（逐次）、OpenMC 1B が 17.8 分（statepoint の total）、xs と
  lattice-1B が約 9 分

## 帰結

- ブランチの上で走らせた計算は、走らせたコミットと分岐点（main のコミット。rebase で書き換わらない）を記録に書く。
  rebase で deps が変わっても走らせ直さず、分岐点から復元する（分岐点にそこまでのブランチのコミットを当て直し、
  `dvc status` の md5 で確かめる）。数秒で終わる段階だけは走らせ直してよい。根本の直しは、OpenMC とライブラリの段階を
  「入力とスクリプト → model.xml」と「model.xml とライブラリ → statepoint」に分け、スクリプトの作り直しで計算が
  走り直さないようにすること（#0018）
- rebase やあとのコミットで deps（スクリプト、入力、コード）が変わると、`dvc.lock` は HEAD と合わなくなり、`dvc status` は
  changed を示す。`dvc.lock` は、どの deps から出力を作ったかの記録なので、そのままでよい。HEAD で `dvc repro` すれば
  変わった段階が走り直すので、結果が要るときだけ流す
- ライブラリと OpenMC の段階はディレクトリを deps にする（案 I）ので、どの段階の outs でもないファイル（入力から消した
  核種の HDF5、試しに置いたファイル）がディレクトリに残ると、そのハッシュが変わって走り直す。消した核種の HDF5 は手で消す
- 段階にしたのは、ライブラリと VERA 1B だけである。1C〜2K、MOX、C5G7、#0005 の均質な媒質（`openmc_homogeneous.py`）、
  `openmc_balance.py`、`openmc_mg.py` は、同じ形（入力の TOML と段階）で足す（#0018）。`openmc_mg.py` は、
  `openmc_model.py` から移した位置の引数の読み方（既定値あり）を、それまで持つ
- 2B〜2K の抽出の段階は、ピン出力を結果のディレクトリの `pins.msgpack` に書く。比べる道具と ADR-0121 の置き場
  （`p2/openmc-<問題>.msgpack`）との合わせ方は、2B〜2K を段階にするブランチで決める
- 入力の `effective-xs-groups`（燃料の実効断面積の tally の群構造を読む xlc のライブラリ）は、まだどの段階も使わない。
  使う段階は、そのファイルを deps に足す（dvc.yaml のコメント）。xlc のライブラリは xs のコードを変えるたびに作り直されるので、
  そのたびに OpenMC が走り直す。群構造だけを入力に書くかは、使う段階を足すときに決める（#0018）
- ライブラリのディレクトリ全体を OpenMC の段階の deps にするので（案 I）、核種を 1 つ足すとすべての OpenMC の段階が
  走り直す（2B〜2K は 1 問題 30 分ほど）。核種を足すことが増えたら見直す
- main の examples/vera/ には DVC を入れる前の `vera.msgpack` が残る。マージのあとに `dvc checkout
  examples/vera/dvc.yaml` で dvc.lock の出力に置き換える（ADR-0123 の M0 と同じ）
- `references/openmc/` の手で作った作業ディレクトリ（`endfb80`、`endfb80-565`、`vera-*`、`work-*`）は消していない。すべてを段階に
  するまでは、まだ段階にしていない記録の出し直しに使う
- `nix develop` と `nix eval` の評価は段階ごとに 1 秒ほどで、ライブラリの 63 段階では 1 分ほどになる
- `x86_64-linux` を cmd に書いた（flake.nix の `system` と同じ）。ほかの system で走らせるときは直す
