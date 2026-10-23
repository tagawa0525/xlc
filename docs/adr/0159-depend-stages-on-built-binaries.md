---
status: accepted
date: 2026-10-14
requires: [ADR-0122, ADR-0126, ADR-0130]
supersedes: [ADR-0123]
superseded-by:
amended-by: [ADR-0184]
issues: ["#0080", "#0110"]
---

# ADR-0159: DVC のステージは、ソースでなくビルドした実行ファイルに依存させる

## 背景

ADR-0123 は、計算のステージを DVC で駆動し、Rust のステージの deps に、そのステージのバイナリを作るコード（ステージのクレートと、パスで
依存するクレートの `src` と `Cargo.toml`、ワークスペースの `Cargo.toml` と `Cargo.lock`、`rust-toolchain.toml`、
`.cargo/config.toml`、`flake.lock`）を並べた（ADR-0123 の案 H）。ビルドしたバイナリを deps にする案（ADR-0123 の案 I）は、
バイナリに worktree の絶対パスが埋まり、同じコードでも worktree ごとに md5 が変わるとして退けた。

ソースを deps にすると、計算の出力を変えない変更でもステージが古くなる:

- テストの移動（#0089、ADR-0143）は `crates/*/src` の md5 を変え、xs、lattice、core のステージをほぼすべて古くする。#0089 は
  移動の前後で strip した release の実行ファイルがバイト一致することを確かめたが、ADR-0143 の文面では、190 ステージの
  `dvc repro -P`（5〜8 時間）を再実行することになる
- `flake.nix` の default の devShell の変更（環境変数、ツールチェーンの override）は、ビルドを変えても deps に表れない
  （#0080）。逆に `flake.nix` を deps にすると、開発のツールを足すだけでステージが古くなる（ADR-0135）
- リポジトリの外の設定（ホームの `~/.cargo/config.toml` のリンカ、`RUSTC_WRAPPER`）は deps にできない

判断の軸:

- ステージが走るかどうかを、ステージが実際に呼ぶプログラムで決めること（コードの変更が実行ファイルを変えなければ走らない）
- `dvc.lock` が、出力を作ったプログラムを正しく指すこと（ADR-0122 の「hash から出し直せる」）
- 依存のクレートの並びのような、手でコピーする一覧を持たないこと（ADR-0123 の「帰結」が挙げたコピーの危うさ）
- worktree、Slurm、`nix develop` の今の運用で、同じソースから同じ実行ファイルができること

## 検討した案

ステージの deps に何を書くか:

- 案 A（採用）: ビルドのステージ `build` をリポジトリの直下の `dvc.yaml` に置き、xs、lattice、core の release の実行ファイルを
  `bin/` に作る（outs。キャッシュする）。計算のステージは deps にその実行ファイルを書き、cmd は `cargo run` でなくその
  実行ファイルを呼ぶ。`build` の deps は実行ファイルを作りうるもの全体（`crates`、`Cargo.toml`、`Cargo.lock`、
  `rust-toolchain.toml`、`.cargo/config.toml`、`flake.nix`、`flake.lock`、ビルドのスクリプト）にする。どれかが変わると
  `build` だけが走り、実行ファイルのバイトが変わらなければ、計算のステージは走らない
- 案 B: ソースを deps にしたまま、出力を変えない変更のあとで `dvc commit` して deps だけを取り込む。ADR-0143 が退けた形で、
  出力を変えないかの判断が人に残り、`dvc.lock` が、記録した deps のコードで実行していない出力を指す（ADR-0130 の案 L と同じ）
- 案 C: ビルドを DVC の外に置き、計算のステージの deps を `target/release/xs` などにする。ソースが変わっても、人がビルドする
  まで DVC は古いと分からない。`target/release` は手で試すビルド（`cargo build --config` でコンパイルの設定を変えたものなど）とも共有するので、
  ステージがそのビルドを呼びうる。`target/` は `.dvcignore` で外してある
- 案 D: 実行ファイルの md5 をテキストに書くステージを置き、それを deps にする。判定は案 A と同じだが、実行ファイルを
  キャッシュしないので、`dvc checkout` で出力を作った実行ファイルを戻せない。間に 1 つファイルが増えるだけで得るものがない
- 案 E: 実行ファイルごとに `build` を分け（xs、lattice、core）、それぞれの deps に依存のクレートだけを並べる。lattice を
  変えても xs を作り直さずに済むが、クレートの並びのコピーが残る。案 A でも、lattice だけを変えたビルドで xs と core の
  実行ファイルはバイト一致する（「確認」の (c)）ので、xs と core のステージは走らない。分けて得るのはビルドの数秒だけである
- 案 F: `build` を `always_changed` にする。リポジトリの外の設定の変更も毎回の再ビルドで拾うが、ビルドするまで `dvc status`
  がソースの変更を示さない。リポジトリの外の設定を変えたときは、`dvc repro -f -s build` で作り直せば足りる

同じソースから同じ実行ファイルを作る方法（「確認」の (a) で、worktree ごとに違ったのは RUNPATH だけだった）:

- 案 P1（採用）: default の devShell に `NIX_NO_SELF_RPATH = "1"` を設定する。nixpkgs の stdenv の setup は、devShell の
  `$out`（`nix develop` を実行したディレクトリの `outputs/out`）の `lib` を `NIX_LDFLAGS` の `-rpath` に足す
  （`_addRpathPrefix`）。この変数はそれを止める。そこに置くライブラリはないので、実行ファイルの動きは変わらない
- 案 P2: 実行ファイルを strip する。strip は RUNPATH を消さないので、strip した版も worktree ごとに違った（「確認」の (a)）
- 案 P3: `--remap-path-prefix`。パニックの位置などのソースのパスは、ワークスペースのクレートはワークスペースからの相対、
  依存のクレートは `~/.cargo/registry` の下で、もともと worktree によらない。リンカの RUNPATH には効かない
- 案 P4: リンクのあとに patchelf で RUNPATH を書き換える。原因（stdenv が足す `-rpath`）を残したまま、ビルドの外のツールで
  結果を直すことになる

ビルドの置き場:

- 案 T1（採用）: 手で使う `target/release` と分けた `target/dvc` で、xlc-xs、xlc-lattice、xlc-core を毎回
  `cargo clean -p` してから `cargo build --release --locked` する。cargo はリンクの環境（`NIX_LDFLAGS`、ホームの設定の
  リンカ）を fingerprint に含めない（devShell の変数を変えても `cargo build` はリンクし直さなかった。2026-10-14）ので、
  消さなければ、前の環境でリンクした実行ファイルが残る。3 つのパッケージのコンパイルとリンクは毎回今の環境でやり直し、
  依存のクレートは再利用する（sccache なしで 5 秒、CPU 45 秒）。分けた置き場なので、手で試したビルド（`cargo build --config`
  でコンパイルの設定を変えたものなど）も拾わない。ワークスペース全体をビルドするのは、`-p` で選ぶと機能の統合が変わりうる
  ので、手で `cargo build --release` したものと同じ実行ファイルにするため
- 案 T2: `target/release` をそのまま使う。ビルドは 1 つで済むが、上の理由で、別の環境でリンクした実行ファイルや、手で
  試したビルドがステージに入りうる
- 案 T3: 毎回空のターゲットのディレクトリから全体をビルドする。依存のクレートも今の環境でコンパイルし直すが、ビルドが
  毎回 12 秒（CPU 100 秒。sccache なし）になる。依存のクレートの rlib はリンクの環境を使わない

## 決定と理由

案 A、P1、T1 にする。ADR-0123 の案 H（deps にコードのソースを並べる）を置き換え、ADR-0123 のほかの決定は次のとおり
引き継ぐ。

- 計算のステージは DVC で駆動する。ステージは `examples/<問題>/dvc.yaml` に書き（問題ごとの置き場所は ADR-0116）、cmd は
  作業ディレクトリをそのディレクトリにした素のコマンドにする。結果を標準出力に出すステージは、標準出力をファイルに受けて outs に
  する。実行するのは `dvc repro examples/<問題>/dvc.yaml`（すべては `dvc repro -P`）で、Slurm で流すときは `dvc repro`
  全体を流す。ステージは逐次に実行する（ADR-0123 の案 A。ツールの比較は ADR-0123 の「検討した案」）
- deps には、ステージの入力（入力の TOML はファイル全体と、そこから読むファイル）と、ステージが呼ぶ実行ファイル
  （`<リポジトリの直下>/bin/<xs|lattice|core>`）を並べる
- `dvc.yaml` と `dvc.lock` をコミットする。キャッシュは `references/dvc-cache` に置き、reflink でコピーする（copy に黙って戻らない。
  ADR-0123 の案 J）。worktree では `references` を main の `references/` への symlink にして共有する（ADR-0126）
- OpenMC の statepoint は outs にしてキャッシュし、比べる小さな結果は別のステージで抽出する（ADR-0123 の案 M）
- DVC の `core.analytics` と `core.check_update` は切る

`build` のステージ（リポジトリの直下の `dvc.yaml`）は `scripts/build-binaries.sh` を呼ぶ。スクリプトは、default の devShell の
中で走っていること（`NIX_NO_SELF_RPATH` が 1）を確かめ、`target/dvc` に案 T1 のとおり release でビルドし、worktree のパスが埋まって
いないことを確かめてから `bin/` にコピーする。どちらかが成り立たなければエラーで止まる（黙って md5 の違う実行ファイルを
作らない）。`bin/` の実行ファイルは DVC の outs で、git では追跡しない（DVC が `bin/.gitignore` を書く）。

案 A にするのは、ステージが呼ぶプログラムそのものの md5 で判定するので、出力を変えない変更（テストの移動、`src` の外の
変更、開発のツールの追加）では計算のステージが走らず、`dvc.lock` が出力を作った実行ファイルを指すからである。リポジトリの外の
設定やツールチェーンが変わったときも、実行ファイルが変わればステージが走る（#0080 の、ビルドの環境が deps に表れない問題が
なくなる）。依存のクレートの並びの写しもなくなる。案 A を取れるのは、案 P1 で同じソースから worktree によらず同じ
実行ファイルができる（「確認」の (a)、(b)）からで、ADR-0123 が案 I を退けた理由（worktree ごとに md5 が変わる）は、
その原因が stdenv の `-rpath` だったと分かって取り除けた。

### 移行

この変更で、すべての Rust のステージの deps と cmd が変わるので、`dvc status` はすべてを changed として示す。流し直さずに
`dvc.lock` を揃える。この 1 回だけ、`dvc commit` で deps を取り込む:

- 対象は、この変更の直前の main（0268e5aa）で `dvc status` が changed を示さなかった Rust のステージだけである。そのステージの
  `dvc.lock` は、今のコードの md5 を deps に記録している。つまり出力は今のコードから作られたもので、`build` が今のコードから
  作る実行ファイル（「確認」の (a)、(b) で、ソースとビルドの環境が同じなら、パスとビルドの仕方によらずバイト一致する）は、
  その出力を作ったプログラムと同じ計算をする。ADR-0130 の案 L と ADR-0143 が `dvc commit` を退けたのは、記録した deps の
  コードで実行していない出力を `dvc.lock` が指すからで、ここではコードが変わっていない
- この根拠は、ADR-0123 の deps が実行ファイルを決めていたという前提に乗る。古い deps には `flake.nix` とリポジトリの外の
  設定（ホームの `~/.cargo/config.toml` など）がなかったので（#0080）、出力を作ったときからそれらが変わっていないことは
  確かめていない。ADR-0123 の `dvc.lock` が持っていた限界を、そのまま引き継ぐ
- 変更の前から changed だったステージは、`dvc.lock` を変えない。古いソースの deps を記録したまま残るので、changed のままで、
  ふだんの規則（ADR-0130）で必要なときに実行する
- Rust でないステージで、この変更が deps を変えたもの（`dvc.yaml` を書くツールを deps にするステージ）は、数秒で終わるなら再実行する

## 確認

測ったのは main の 0268e5aa のソース（ツールチェーン 1.99.0、ホームの設定は sccache、clang、mold）で、2026-10-14。
`NIX_NO_SELF_RPATH` を足したあとの測定は、0268e5aa にこのブランチの 5d3c25ca の devShell の変更を当てたもの。

(a) 同じコミットを、パスの違う 2 つの worktree（`xlc-feat-0159-dvc-binary-deps` と、一時の `xlc-tmp-0159-repro`）で、
それぞれの直下から `nix develop -c cargo build --release` した。xs はバイト 25、lattice と core はバイト 577 から違った。
strip した版も同じ場所から違った。違いは RUNPATH の先頭の `<worktree>/outputs/out/lib` で、ほかに worktree のパスを含む文字列は
なかった。devShell に `NIX_NO_SELF_RPATH = "1"` を足してリンクし直すと、RUNPATH は nix store のパスだけになり、strip しない
実行ファイルが 2 つの worktree でバイト一致した（md5: xs `8376dbec…`、lattice `d9811966…`、core `2f204012…`）。

(b) 同じ md5 になったビルド:

- sccache を外し（`RUSTC_WRAPPER=`）、空のターゲットのディレクトリから全体（55 のクレート）をビルドしたもの（12 秒、CPU
  100 秒）
- 一時の worktree で e67c79cb の 1 つ前のソースにしてビルドしたあと、0268e5aa のソースに戻して incremental にビルドした
  もの（sccache あり、なしの両方。なしでは、ソースの変わった xlc-core と xlc-lattice だけをコンパイルし直した）
- `dvc repro -s build` と `dvc repro -f -s build`（`scripts/build-binaries.sh` の `target/dvc`。3 つのパッケージを消して
  コンパイルし直し、リンクし直したもの。sccache なしでも同じ）

(c) 実行ファイルが変わるかの感度（`dvc repro -s build` で作り直し、戻して元の md5 に戻ることも見た）:

- テストだけのファイル（crates/lattice/src/axial/tests.rs）の末尾にコメントの行を足す: 3 つとも変わらない
- crates/lattice/src/moc.rs の先頭にコメントの行を足す: lattice だけが変わり（`99b8a12c…`）、xs と core は変わらない。
  原因の候補は、`tracing` の呼び出しの位置とパニックの位置として、ソースのパスと行番号が実行ファイルに入ること
  （lattice の実行ファイルは文字列 `crates/lattice/src/moc.rs` を含む。行番号が入ることは確かめていない。#0110）
- default の devShell の `packages` に `pkgs.libxml2`（`NIX_LDFLAGS` に `-L` が増える）と `pkgs.hello` を足し、`dvc repro -f -s build`
  で作り直す: 3 つとも変わらない（5747d4fa の上で。下の (f)）

(d) DVC の振る舞い（DVC 3.67.1。一時のリポジトリで、ルートの `build` とサブディレクトリの計算のステージの組）:
`dvc repro <サブディレクトリ>/dvc.yaml` は、上流のルートの `build` も実行した。`build` のソースを、出力を変えないように
変えると、`build` は走り、計算のステージは「didn't change」で走らなかった。`dvc status <サブディレクトリ>/dvc.yaml` は
上流の `build` の changed を示さず、`dvc status`（引数なし）と `--with-deps` は示した。`bin/` の実行ファイルは
`dvc.lock` に `isexec: true` で記録され、消して `dvc checkout` すると実行できる形で戻った。

(e) 移行（「決定と理由」の「移行」）: 0268e5aa で `dvc status` が changed を示したのは 127 項目（`references-shell` の
`always changed` と、126 のステージ: examples/core 55、vera/p5 31、vera/p3 19、reflector 15、m0 3、vera 3）だった。切り替えの
あとに新たに changed になった 64 のステージのうち、63 の Rust のステージ（examples/vera の xs-range、xs-nuclide、xs など 62 と、
examples/reflector の xs）を `dvc commit` した。`dvc.lock` の差分は、その 63 のステージの cmd と、deps のコードの項目が
実行ファイルの 1 項目に替わったことだけで、outs と params とほかの deps は変わらなかった（両方の `dvc.lock` を YAML として
比べた）。残る 1 つの examples/core の `compare`（deps の rods_cases.py が変わった。数秒）は `dvc repro -s` で走らせ直し、
`compare.out` の md5 は変わらなかった（`a4a1c8e6…`）。そのあとの `dvc status` は、切り替えの前と同じ 127 項目を示した。

(f) マージの前に main の 5747d4fa（#0039 の自己遮蔽の高速化で crates/lattice が変わった）に rebase し、`build` を再実行した。
lattice は `bfe45232…` に変わり、xs と core は (a) と同じだった。lattice のステージはすべて切り替えの前から changed だったので、
changed のステージは増えなかった（127 項目のまま）。移行の `dvc commit` は 0268e5aa の上で行ったもので、そのステージは xs の
ステージと reflector の xs だけである。

## 帰結

- 出力を変えない変更のうち、実行ファイルのバイトを変えないもの（テストの変更、`src` の外の変更、開発のツールの追加、
  依存しないクレートの変更）では、`build` だけが走り（数秒〜十数秒）、計算のステージは走らない。#0089 のテストの移動は、
  実行ファイルがバイト一致すれば、`build` の `dvc.lock` を更新するだけで済む
- 実行ファイルのバイトを変える変更は、出力を変えなくても計算のステージを実行する。行の位置を変える整理（ADR-0144 の関数の
  分割など）、ループの整列（ブランチ perf/0105-code-alignment の `-align-loops`）、ツールチェーン、`flake.lock` の glibc、ホームの設定のリンカ（mold の版は
  `.comment` に入る）の変更がこれにあたる。行の位置で実行ファイルが変わる件は #0110 で扱う
- `build` の deps にないもの（ホームの `~/.cargo/config.toml`、環境変数）を変えたときは、`dvc repro -f -s build` で作り
  直す。実行ファイルが変わらなければ、計算のステージは走らない
- `dvc status examples/<問題>/dvc.yaml` は上流の `build` を示さない。ソースの変更で古くなったかは、`dvc status`（引数なし）か
  `--with-deps` で見る。どの計算のステージが古くなったかは、`build` を実行するまで分からない（実行ファイルが変わったかで
  決まる）。rebase のあとの手順（/dvc-relock）は、先に `build` を実行する
- リポジトリの直下の `dvc.lock` は、最後に `build` を実行したときの実行ファイルを記録する。計算のステージを実行したら、
  例の `dvc.lock` と一緒にコミットする（例の `dvc.lock` が記録する実行ファイルの md5 と揃えるため）。コードだけを変える
  ブランチは、ほかのステージと同じく（ADR-0130）`build` を古いまま残してよい。その main で `dvc checkout dvc.yaml` が戻すのは
  古いソースの実行ファイルで、`dvc repro` すれば `build` が作り直す
- `bin/` の実行ファイルは worktree ごとにあり、`dvc checkout` でキャッシュから戻る。ある commit の `dvc.lock` が指す実行ファイルが
  キャッシュにあれば、出力を作った実行ファイルそのものを戻せる。実行ファイルは動的にリンクしていて、RUNPATH と ELF の
  インタープリタが nix store の glibc などを指すので、実行するにはその store のパスが要る（nix の GC で消えたら、その
  commit の `flake.lock` で devShell に入れば戻る）。キャッシュは、作った実行ファイルの版ごとに約 13 MB
  （3 つの合計）ずつ増える。消し方は ADR-0123 の帰結のとおり `dvc gc --all-commits`
- 実行ファイルのバイトは、リポジトリの外の設定（ホームの `~/.cargo/config.toml` の clang と mold、sccache）とホストにも
  よる。手元と Slurm のジョブはどちらも r995 で同じホームの設定を使う。別のホストでビルドすると md5 が変わり、計算の
  ステージが走り直しうる
- Slurm で流すときは、`nix develop -c dvc repro …` のジョブの中で `build` が走る（ビルドの時間がジョブに入る）。`build` は
  default の devShell の中でしか走らない。`NIX_NO_SELF_RPATH` は devShell に入ったときに決まるので、この変更より前に入った
  shell（direnv を含む）は入り直す
- examples に dvc.yaml を足すブランチは、Rust のステージを同じ形（cmd で `bin/` の実行ファイルを呼び、deps にその実行ファイルを
  書く）にする。`dvc.yaml` を書くツール（rods_cases.py、p3_inputs.py、p5_inputs.py）は、この形を書く
- `build` は `crates` 全体を deps にするので、`crates/docs`（文書の検査）の変更でも走る。ビルドが数秒なので、並びを
  絞らない
