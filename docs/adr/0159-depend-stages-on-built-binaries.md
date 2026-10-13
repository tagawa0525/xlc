---
status: accepted
date: 2026-10-14
requires: [ADR-0122, ADR-0126, ADR-0130]
supersedes: [ADR-0123]
superseded-by:
issues: ["#0080", "#0110"]
---

# ADR-0159: DVC の段階は、ソースでなくビルドした実行ファイルに依存させる

## 背景

ADR-0123 は、計算の段階を DVC で駆動し、Rust の段階の deps に、その段階のバイナリを作るコード（段階のクレートと、パスで
依存するクレートの `src` と `Cargo.toml`、ワークスペースの `Cargo.toml` と `Cargo.lock`、`rust-toolchain.toml`、
`.cargo/config.toml`、`flake.lock`）を並べた（ADR-0123 の案 H）。ビルドしたバイナリを deps にする案（ADR-0123 の案 I）は、
バイナリに worktree の絶対パスが埋まり、同じコードでも worktree ごとに md5 が変わるとして退けた。

ソースを deps にすると、計算の出力を変えない変更でも段階が古くなる:

- テストの移動（#0089、ADR-0143）は `crates/*/src` の md5 を変え、xs、lattice、core の段階をほぼすべて古くする。#0089 は
  移動の前後で strip した release の実行ファイルがバイト一致することを確かめたが、ADR-0143 の文面では、190 段階の
  `dvc repro -P`（5〜8 時間）を流し直すことになる
- `flake.nix` の default の devShell の変更（環境変数、ツールチェーンの override）は、ビルドを変えても deps に表れない
  （#0080）。逆に `flake.nix` を deps にすると、開発の道具を足すだけで段階が古くなる（ADR-0135）
- リポジトリの外の設定（ホームの `~/.cargo/config.toml` のリンカ、`RUSTC_WRAPPER`）は deps にできない

判断の軸:

- 段階が走るかどうかを、段階が実際に呼ぶプログラムで決めること（コードの変更が実行ファイルを変えなければ走らない）
- `dvc.lock` が、出力を作ったプログラムを正しく指すこと（ADR-0122 の「hash から出し直せる」）
- 依存のクレートの並びのような、手で写す一覧を持たないこと（ADR-0123 の「帰結」が挙げた写しの危うさ）
- worktree、Slurm、`nix develop` の今の運用で、同じソースから同じ実行ファイルができること

## 検討した案

段階の deps に何を書くか:

- 案 A（採用）: ビルドの段階 `build` をリポジトリの直下の `dvc.yaml` に置き、xs、lattice、core の release の実行ファイルを
  `bin/` に作る（outs。キャッシュする）。計算の段階は deps にその実行ファイルを書き、cmd は `cargo run` でなくその
  実行ファイルを呼ぶ。`build` の deps は実行ファイルを作りうるもの全体（`crates`、`Cargo.toml`、`Cargo.lock`、
  `rust-toolchain.toml`、`.cargo/config.toml`、`flake.nix`、`flake.lock`、ビルドのスクリプト）にする。どれかが変わると
  `build` だけが走り、実行ファイルのバイトが変わらなければ、計算の段階は走らない
- 案 B: ソースを deps にしたまま、出力を変えない変更のあとで `dvc commit` して deps だけを取り込む。ADR-0143 が退けた形で、
  出力を変えないかの判断が人に残り、`dvc.lock` が、記録した deps のコードで走らせていない出力を指す（ADR-0130 の案 L と同じ）
- 案 C: ビルドを DVC の外に置き、計算の段階の deps を `target/release/xs` などにする。ソースが変わっても、人がビルドする
  まで DVC は古いと分からない。`target/release` は手で試すビルド（`cargo build --config` でコンパイルの設定を変えたものなど）とも共有するので、
  段階がそのビルドを呼びうる。`target/` は `.dvcignore` で外してある
- 案 D: 実行ファイルの md5 をテキストに書く段階を置き、それを deps にする。判定は案 A と同じだが、実行ファイルを
  キャッシュしないので、`dvc checkout` で出力を作った実行ファイルを戻せない。間に 1 つファイルが増えるだけで得るものがない
- 案 E: 実行ファイルごとに `build` を分け（xs、lattice、core）、それぞれの deps に依存のクレートだけを並べる。lattice を
  変えても xs を作り直さずに済むが、クレートの並びの写しが残る。案 A でも、lattice だけを変えたビルドで xs と core の
  実行ファイルはバイト一致する（「確認」の (c)）ので、xs と core の段階は走らない。分けて得るのはビルドの数秒だけである
- 案 F: `build` を `always_changed` にする。リポジトリの外の設定の変更も毎回の再ビルドで拾うが、ビルドするまで `dvc status`
  がソースの変更を示さない。リポジトリの外の設定を変えたときは、`dvc repro -f -s build` で作り直せば足りる

同じソースから同じ実行ファイルを作る方法（「確認」の (a) で、worktree ごとに違ったのは RUNPATH だけだった）:

- 案 P1（採用）: default の devShell に `NIX_NO_SELF_RPATH = "1"` を設定する。nixpkgs の stdenv の setup は、devShell の
  `$out`（`nix develop` を走らせたディレクトリの `outputs/out`）の `lib` を `NIX_LDFLAGS` の `-rpath` に足す
  （`_addRpathPrefix`）。この変数はそれを止める。そこに置くライブラリはないので、実行ファイルの動きは変わらない
- 案 P2: 実行ファイルを strip する。strip は RUNPATH を消さないので、strip した版も worktree ごとに違った（「確認」の (a)）
- 案 P3: `--remap-path-prefix`。パニックの位置などのソースのパスは、ワークスペースのクレートはワークスペースからの相対、
  依存のクレートは `~/.cargo/registry` の下で、もともと worktree によらない。リンカの RUNPATH には効かない
- 案 P4: リンクのあとに patchelf で RUNPATH を書き換える。原因（stdenv が足す `-rpath`）を残したまま、ビルドの外の道具で
  結果を直すことになる

ビルドの置き場:

- 案 T1（採用）: `cargo build --release --locked --target-dir target/dvc`。cargo はリンクの環境（`NIX_LDFLAGS`、ホームの設定の
  リンカ）を fingerprint に含めない（devShell の変数を変えても `cargo build` はリンクし直さなかった。2026-10-14）。手で
  使う `target/release` と分けると、前の環境でリンクした実行ファイルや手で試したビルドを `build` が拾わない（この決定の
  あとに作る `target/dvc` は、今の環境でリンクされる）。ワークスペース全体をビルドするのは、`-p` で選ぶと機能の統合が
  変わりうるので、手で `cargo build --release` したものと同じ実行ファイルにするため
- 案 T2: `target/release` を使う。ビルドは 1 つで済むが、上の理由で、別の環境でリンクした実行ファイルや、手で試した
  ビルドが段階に入りうる

## 決定と理由

案 A、P1、T1 にする。ADR-0123 の案 H（deps にコードのソースを並べる）を置き換え、ADR-0123 のほかの決定は次のとおり
引き継ぐ。

- 計算の段階は DVC で駆動する。段階は `examples/<問題>/dvc.yaml` に書き（問題ごとの置き場所は ADR-0116）、cmd は
  作業ディレクトリをそのディレクトリにした素のコマンドにする。結果を標準出力に出す段階は、標準出力をファイルに受けて outs に
  する。走らせるのは `dvc repro examples/<問題>/dvc.yaml`（すべては `dvc repro -P`）で、Slurm で流すときは `dvc repro`
  全体を流す。段階は逐次に走らせる（ADR-0123 の案 A。道具の比較は ADR-0123 の「検討した案」）
- deps には、段階の入力（入力の TOML はファイル全体と、そこから読むファイル）と、段階が呼ぶ実行ファイル
  （`<リポジトリの直下>/bin/<xs|lattice|core>`）を並べる
- `dvc.yaml` と `dvc.lock` をコミットする。キャッシュは `references/dvc-cache` に置き、reflink で写す（copy に黙って戻らない。
  ADR-0123 の案 J）。worktree では `references` を main の `references/` への symlink にして共有する（ADR-0126）
- OpenMC の statepoint は outs にしてキャッシュし、比べる小さな結果は別の段階で抽出する（ADR-0123 の案 M）
- DVC の `core.analytics` と `core.check_update` は切る

`build` の段階（リポジトリの直下の `dvc.yaml`）は `scripts/build-binaries.sh` を呼ぶ。スクリプトは、default の devShell の
中で走っていること（`NIX_NO_SELF_RPATH` が 1）を確かめ、`target/dvc` に release でビルドし、worktree のパスが埋まって
いないことを確かめてから `bin/` に写す。どちらかが成り立たなければエラーで止まる（黙って md5 の違う実行ファイルを
作らない）。`bin/` の実行ファイルは DVC の outs で、git では追跡しない（DVC が `bin/.gitignore` を書く）。

案 A にするのは、段階が呼ぶプログラムそのものの md5 で判定するので、出力を変えない変更（テストの移動、`src` の外の
変更、開発の道具の追加）では計算の段階が走らず、`dvc.lock` が出力を作った実行ファイルを指すからである。リポジトリの外の
設定やツールチェーンが変わったときも、実行ファイルが変われば段階が走る（#0080 の、ビルドの環境が deps に表れない問題が
なくなる）。依存のクレートの並びの写しもなくなる。案 A を取れるのは、案 P1 で同じソースから worktree によらず同じ
実行ファイルができる（「確認」の (a)、(b)）からで、ADR-0123 が案 I を退けた理由（worktree ごとに md5 が変わる）は、
その原因が stdenv の `-rpath` だったと分かって取り除けた。

### 移行

この変更で、すべての Rust の段階の deps と cmd が変わるので、`dvc status` はすべてを changed として示す。流し直さずに
`dvc.lock` を揃える。この 1 回だけ、`dvc commit` で deps を取り込む:

- 対象は、この変更の直前の main（387743fc）で `dvc status` が changed を示さなかった Rust の段階だけである。その段階の
  `dvc.lock` は、今のコードの md5 を deps に記録している。つまり出力は今のコードから作られたもので、`build` が今のコードから
  作る実行ファイル（「確認」の (a)、(b) で、ソースとビルドの環境が同じなら、パスとビルドの仕方によらずバイト一致する）は、
  その出力を作ったプログラムと同じ計算をする。ADR-0130 の案 L と ADR-0143 が `dvc commit` を退けたのは、記録した deps の
  コードで走らせていない出力を `dvc.lock` が指すからで、ここではコードが変わっていない
- 変更の前から changed だった段階は、`dvc.lock` を変えない。古いソースの deps を記録したまま残るので、changed のままで、
  ふだんの規則（ADR-0130）で必要なときに走らせる
- Rust でない段階で、この変更が deps を変えたもの（`dvc.yaml` を書く道具を deps にする段階）は、数秒で終わるなら走らせ直す

## 確認

測ったのは main の 387743fc のソース（ツールチェーン 1.99.0、ホームの設定は sccache、clang、mold）で、2026-10-14。
`NIX_NO_SELF_RPATH` を足したあとの測定は、387743fc にこのブランチの adf77759 の devShell の変更を当てたもの。

(a) 同じコミットを、パスの違う 2 つの worktree（`xlc-feat-0159-dvc-binary-deps` と、一時の `xlc-tmp-0159-repro`）で、
それぞれの直下から `nix develop -c cargo build --release` した。xs はバイト 25、lattice と core はバイト 577 から違った。
strip した版も同じ場所から違った。違いは RUNPATH の先頭の `<worktree>/outputs/out/lib` で、ほかに worktree のパスを含む文字列は
なかった。devShell に `NIX_NO_SELF_RPATH = "1"` を足してリンクし直すと、RUNPATH は nix store のパスだけになり、strip しない
実行ファイルが 2 つの worktree でバイト一致した（md5: xs `8376dbec…`、lattice `d9811966…`、core `2f204012…`）。

(b) 同じ md5 になったビルド:

- sccache を外し（`RUSTC_WRAPPER=`）、空のターゲットのディレクトリから全体（55 のクレート）をビルドしたもの（12 秒、CPU
  100 秒）
- 一時の worktree で 16781fc9 の 1 つ前のソースにしてビルドしたあと、387743fc のソースに戻して incremental にビルドした
  もの（sccache あり、なしの両方。なしでは、ソースの変わった xlc-core と xlc-lattice だけをコンパイルし直した）
- `dvc repro -s build`（`scripts/build-binaries.sh` の `target/dvc`）

(c) 実行ファイルが変わるかの感度（`dvc repro -s build` で作り直し、戻して元の md5 に戻ることも見た）:

- テストだけのファイル（crates/lattice/src/axial/tests.rs）の末尾にコメントの行を足す: 3 つとも変わらない
- crates/lattice/src/moc.rs の先頭にコメントの行を足す: lattice だけが変わり（`99b8a12c…`）、xs と core は変わらない。
  原因の候補は、`tracing` の呼び出しの位置とパニックの位置として、ソースのパスと行番号が実行ファイルに入ること
  （lattice の実行ファイルは文字列 `crates/lattice/src/moc.rs` を含む。行番号が入ることは確かめていない。#0110）

(d) DVC の振る舞い（DVC 3.67.1。一時のリポジトリで、ルートの `build` とサブディレクトリの計算の段階の組）:
`dvc repro <サブディレクトリ>/dvc.yaml` は、上流のルートの `build` も走らせた。`build` のソースを、出力を変えないように
変えると、`build` は走り、計算の段階は「didn't change」で走らなかった。`dvc status <サブディレクトリ>/dvc.yaml` は
上流の `build` の changed を示さず、`dvc status`（引数なし）と `--with-deps` は示した。`bin/` の実行ファイルは
`dvc.lock` に `isexec: true` で記録され、消して `dvc checkout` すると実行できる形で戻った。

(e) 移行（「決定と理由」の「移行」）: 387743fc で `dvc status` が changed を示したのは 127 項目（`references-shell` の
`always changed` と、126 の段階: examples/core 55、vera/p5 31、vera/p3 19、reflector 15、m0 3、vera 3）だった。切り替えの
あとに新たに changed になった 64 の段階のうち、63 の Rust の段階（examples/vera の xs-range、xs-nuclide、xs など 62 と、
examples/reflector の xs）を `dvc commit` した。`dvc.lock` の差分は、その 63 の段階の cmd と、deps のコードの項目が
実行ファイルの 1 項目に替わったことだけで、outs と params とほかの deps は変わらなかった（両方の `dvc.lock` を YAML として
比べた）。残る 1 つの examples/core の `compare`（deps の rods_cases.py が変わった。数秒）は `dvc repro -s` で走らせ直し、
`compare.out` の md5 は変わらなかった（`a4a1c8e6…`）。そのあとの `dvc status` は、切り替えの前と同じ 127 項目を示した。

## 帰結

- 出力を変えない変更のうち、実行ファイルのバイトを変えないもの（テストの変更、`src` の外の変更、開発の道具の追加、
  依存しないクレートの変更）では、`build` だけが走り（数秒〜十数秒）、計算の段階は走らない。#0089 のテストの移動は、
  実行ファイルがバイト一致すれば、`build` の `dvc.lock` を更新するだけで済む
- 実行ファイルのバイトを変える変更は、出力を変えなくても計算の段階を走らせる。行の位置を変える整理（ADR-0144 の関数の
  分割など）、ループの整列（ブランチ perf/0105-code-alignment の `-align-loops`）、ツールチェーン、`flake.lock` の glibc、ホームの設定のリンカ（mold の版は
  `.comment` に入る）の変更がこれにあたる。行の位置で実行ファイルが変わる件は #0110 で扱う
- `build` の deps にないもの（ホームの `~/.cargo/config.toml`、環境変数）を変えたときは、`dvc repro -f -s build` で作り
  直す。実行ファイルが変わらなければ、計算の段階は走らない
- `dvc status examples/<問題>/dvc.yaml` は上流の `build` を示さない。ソースの変更で古くなったかは、`dvc status`（引数なし）か
  `--with-deps` で見る。どの計算の段階が古くなったかは、`build` を走らせるまで分からない（実行ファイルが変わったかで
  決まる）。rebase のあとの手順（/dvc-relock）は、先に `build` を走らせる
- `bin/` の実行ファイルは worktree ごとにあり、`dvc checkout` でキャッシュから戻る。ある commit の `dvc.lock` が指す実行ファイルが
  キャッシュにあれば、出力を作った実行ファイルそのものを戻せる。キャッシュは、作った実行ファイルの版ごとに約 13 MB
  （3 つの合計）ずつ増える。消し方は ADR-0123 の帰結のとおり `dvc gc --all-commits`
- Slurm で流すときは、`nix develop -c dvc repro …` のジョブの中で `build` が走る（ビルドの時間がジョブに入る）。`build` は
  default の devShell の中でしか走らない。`NIX_NO_SELF_RPATH` は devShell に入ったときに決まるので、この変更より前に入った
  shell（direnv を含む）は入り直す
- examples に dvc.yaml を足すブランチは、Rust の段階を同じ形（cmd で `bin/` の実行ファイルを呼び、deps にその実行ファイルを
  書く）にする。`dvc.yaml` を書く道具（rods_cases.py、p3_inputs.py、p5_inputs.py）は、この形を書く
- `build` は `crates` 全体を deps にするので、`crates/docs`（文書の検査）の変更でも走る。ビルドが数秒なので、並びを
  絞らない
