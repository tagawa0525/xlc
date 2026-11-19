---
type: ADR
status: accepted
date: 2026-11-19
requires: [ADR-0159, ADR-0160, ADR-0210, ADR-0254]
supersedes: []
superseded-by:
issues: ["#0110", "#0117"]
---

# ADR-0255: 計算のステージは、ソースの行と列の番号を除いた実行ファイルの正規形を deps にする

## 背景

計算のステージは、`build` が `bin/` に作る実行ファイルを deps にする（ADR-0159）。計算を変えない編集でも、ソースの行の位置が
変わると実行ファイルのバイトが変わり、その実行ファイルを呼ぶステージがすべて古くなる（#0110）。ADR-0144 のリファクタリングや
注釈の修正のたびに、出力がビット一致するのにステージを再実行するか、古いまま残すことになる。

main の 92b0751ef で crates/lattice/src/moc.rs の先頭に注釈の行を 1 つ足すと、`bin/lattice` の md5 が c3b3bba1 から
69093522 に変わった（xs、core、chain は同じ。2026-11-19）。大きさは同じで、違う 310 バイトの内訳は次のとおり（節は
`readelf -S` で引き、値を前後の 8 バイトと合わせて読んだ。`.text` とほかの節は同じ。記号は ADR-0160 で落としてあり、
デバッグ情報の節は無い）:

- `.data.rel.ro` の 277 か所: パニックの位置 `core::panic::Location`（配列の境界の検査など）の行。並びは（ファイル名の ptr、
  長さ、行 u32、列 u32）で、ファイル名は `Location::file_as_c_str` のために NUL で終わる（rust-src の
  library/core/src/panic/location.rs）
- `.data.rel.ro` の 22 か所: `tracing` のイベントとスパンの `Metadata` の `line: Option<u32>`
- `.rodata` の 11 バイト: `tracing` のイベントの名前の行の数字。名前を指定しないイベントは
  `"event <ファイル>:<行>"`、`tracing::enabled!` は `"enabled <ファイル>:<行>"`（tracing 0.1.44 の src/macros.rs）

同じファイルの 2 つの関数（`returned` と `add_current`）を入れ替えると、`.data.rel.ro` の 4 バイト（Location の行）だけが
変わった。

crates の追跡しているすべての .rs の先頭に注釈の行を足すと、xs、lattice、chain は上の 3 か所だけが違ったが、core は
`.text` の 53 バイトと `.rodata.cst32` の 398 バイトも違った。.rs を 1 つずつ足して、原因は crates/core/src/nem.rs の
474 行のイベントと分かった。名前 `"event crates/core/src/nem.rs:474"` がちょうど 32 バイトで、リンカがまとめる 32 バイトの
定数の節に入り、その節の並びが中身で変わったので、その節を指す 34 の関数のコードの番地が変わった。行の番号を後から 0 に
しても、並びはリンクの時に決まっているので戻らない。

判断の軸: 正確さ（計算に使う実行ファイルを変えない。変わった実行ファイルを同じとみなさない）、実行速度（リリースの最適化を
変えない）、保守の容易さ。エラーを握りつぶさない（パニックの文の位置を消さない）。

## 決定の要因

- 注釈の行、空行、関数の並べ替えでは、ステージが依存するものを変えない
- ステージが実行するのは cargo が作った実行ファイルのままにする。位置の情報を探す推定が外れたときに起きうるのは「変わったのに
  changed にならない」までで、計算の結果は変わらない
- 推定が外れにくく、外れそうなら止まる（構造体のレイアウトは、ツールチェーンや `tracing` の更新で変わりうる）
- 既にある正規形（.py の python-ast は ADR-0210、入力の TOML の toml-values は ADR-0254）と同じ形にする: cmd は元を読み、deps は
  正規形
- Rust は rust-toolchain.toml の stable の 1.99.0（flake.nix は rust-overlay でこの版を使う）

## 検討した案

パニックの位置:

- 案 A: コンパイラで位置を外す（`-Z location-detail=none`）。1.99.0 の stable では「the option `Z` is only accepted on the
  nightly compiler」で止まり、`-C location-detail` は無い。nightly を使う（または `RUSTC_BOOTSTRAP`）と、ステージの実行ファイルを
  作るツールチェーンが安定版でなくなるので採らない。パニックの位置も失う
- 案 B: `--remap-path-prefix`。ファイル名を変えるだけで、行と列は残る
- 案 C: パニックの経路を位置の要らない形に寄せる（添字の代わりに `get` とエラーを返す）。lattice だけで Location が 4,087 ある
  （「確認」の (a)）。添字の検査はコンパイラが入れるので、本番のコードのすべての添字を書き換えることになり、読みにくく、
  最適化も変わりうる
- 案 D: `build` のステージが、リンクした実行ファイルの番号を 0 にしたコピーを `bin/` に置き、ステージはそれを実行する。
  計算に使う実行ファイルのバイトを、レイアウトの推定で書き換えることになる。推定が外れて別のバイトを 0 にすると、計算の結果が
  黙って変わりうる。今のレイアウトで外れないという確かめ（examples/m0 の 1 例と整合の検査）は、rustc か `tracing` の版が変わると
  前提が動く。パニックの文の位置も `0:0` になる。このブランチで一度この形で作り、マージの前に親が退けた（2026-11-19）
- 案 E（採用）: `bin/` は cargo の出力のまま置き、番号を 0 にしたコピーの md5 を正規形のファイルに書く。計算のステージの cmd は
  `bin/` を実行し、deps には正規形を書く。推定が外れても、起きうるのは変わった実行ファイルを同じとみなすことだけで、ステージが
  実行するプログラムは変わらない。パニックの位置も残る。`dvc.lock` の計算のステージの deps は実行ファイルでなく正規形を指すが、
  python-ast と toml-values と同じ性質で、どの実行ファイルで走ったかは `build` の lock の outs（`bin/` の md5）からたどれる
- 案 F: 番号の場所をレイアウトから探さず、行を 1 つずつずらしたソースでもう一度ビルドし、違ったバイトを消す。レイアウトに
  依らないが、ビルドが倍になり、イベントの名前の桁が変わる位置（9、99、999 行のイベント）では 2 つのビルドの並びが揃わない
- 案 G: 何もしない（#0110 の完了の条件の 2 つ目）。リファクタリングと注釈の修正のたびに、出力が同じステージを再実行するか
  古いまま残す

`tracing` の `Metadata` の行:

- 案 H（採用）: 案 E と同じコピーで 0 にする（正規形にだけ効き、実行するファイルとログは変わらない）
- 案 I: `tracing` の機能で行を出さない。`release_max_level_*` はイベントごと消し、行だけを外す機能は無い（0.1.44 のマクロは
  `line: Some(line!())` を常に書く）

`tracing` のイベントの名前:

- 案 J（採用）: イベントのマクロを xlc-log のマクロ（`xlc_log::info!` など）に替え、名前を `"event <ファイル>"`（行なし）にして
  `tracing::event!` に渡す。名前の長さが行で変わらないので、定数の節の並びも変わらない。clippy.toml の `disallowed-macros` で
  `tracing::trace`〜`tracing::error` を止める（`name:` を付けて `tracing::info!` を呼ぶ形は、clippy が展開の中の呼び出しも
  止め、展開の中の `#[allow]` も効かなかったので、`tracing::event!` を呼ぶ）。ログの出力（`tracing_subscriber::fmt` の既定）は
  名前を出さないので、ログは変わらない
- 案 K: 案 E と同じコピーで数字を 0 にする。上の nem.rs の 474 行のように、定数の節の並びが中身で決まると戻らない。桁が変わる
  ときも、文字列の長さが変わって後ろの並びがずれる
- 案 L: 呼ぶ場所ごとに `name:` を書く。113 の呼び出しに同じ引数を足し、書き忘れを機械で止められない

`tracing::enabled!` は `name:` を取らない（0.1.44）。lattice の 12 か所の名前は 37〜44 バイトで、定数の節の大きさ（4、8、16、
32 バイト）に合わないので、数字だけをコピーで 0 にする（案 K を `enabled!` に限って使う）。自前の `enabled!` は `tracing` の隠した
マクロ（`callsite2!`）を使うことになるので作らない。

正規形の置き場（ADR-0210 の案 G、ADR-0254 と同じ考え方）:

- 採用: リポジトリの直下の `dvc.yaml` の `binary-forms`（foreach。xs、lattice、core、chain）が、`bin/<名前>` から
  `binary-forms/<名前>.md5`（コピーの md5 の 16 進を 1 行）を書く（`cache: false`、git で追跡）。worktree で `dvc checkout`
  しなくても揃う。コピーそのもの（数 MB）は git で追跡せず、md5 だけで比べられる。`build` から分けたので、正規形を書く
  スクリプトの変更では cargo のビルドが走らない
- `build` の outs に正規形を足す形は、正規形のスクリプトの変更（docstring を除く）で cargo の `clean` とビルドが走る

## 決定と理由

案 E、H、J を採る。

- scripts/build-binaries.sh は、ビルドした実行ファイルを `bin/<名前>` にそのまま置く（ADR-0159、ADR-0160 のとおり）。
  `binary-forms` が scripts/normalize_binary.py で正規形を書く。コピーは Location の行と列、`Metadata` の行、`enabled!` の
  名前の行の数字を 0 にし、ほかのバイトを変えない。ほかのバイトが 1 つでも違えば正規形も違う
- ADR-0159 の deps の規則（ステージの入力と、ステージが呼ぶ実行ファイル）のうち、実行ファイルを並べる部分を、その正規形に
  置き換える（ADR-0254 が入力の TOML を置き換えたのと同じ形）。例のステージの deps は `bin/<名前>` の代わりに
  `binary-forms/<名前>.md5` を書き、cmd は `bin/<名前>` を実行する。`dvc.yaml` を
  書くツール（`examples/core/rods_cases.py`、`examples/vera/p3/p3_inputs.py`、`p5_inputs.py`、`p6_inputs.py`）も正規形の
  パスを書く。scripts/test_normalize_binary.py の `RepositoryDepsTest`（pre-commit の hook `binary-forms` も実行する）が、
  生の `bin/` の deps、`binary-forms` が書かない正規形、cmd が実行する実行ファイルと deps の正規形の食い違い、リポジトリの直下の
  `dvc.lock` の `build` の outs と食い違う正規形を止める
- 構造体のレイアウトは `repr(Rust)` で保証されないので、値どうしの整合で見つける。Location は「NUL で終わる `.rs` の文字列を
  指す &str の後に、0 でない行と列」。`Metadata` は「名前 `"event <ファイル>"` の &str の直前に Some(行)、名前から決まった
  距離に同じファイルの &str」で、距離はイベントから決めてスパンにも使う。Metadata のファイル名を先に決め、Location はそれ以外から
  探す。Location が 1 つも無い、名前の直前が Some でない、距離がイベントごとに違う、1 つの &str が Location と Metadata の
  ファイル名の両方に当たる、のどれかでエラーで止める（ツールチェーンか `tracing` のレイアウトが変わった。黙って消し残すと
  注釈の行でステージが古くなるのに気付けず、別のバイトを消すと変わった実行ファイルを同じとみなす）。名前を指定しないイベント
  （`"event crates/…:<行>"`）が残っていてもエラーで止める

守れる範囲: 注釈と空行の追加と削除、行の中の位置（列）の変化、ファイルの中の関数の並べ替え（1 組で確かめた）。守れないもの:
コードの変更（別の関数の変更を含む）、ファイルの名前と、関数のファイルの間の移動（Location のファイル名が変わる）、
lattice の `enabled!` の 12 か所の行の桁の変化（文字列の長さが変わる）。守れないときは正規形が変わり、ステージが
古くなる（安全な側）。

### 移行

deps のパスを替えると、実行ファイルを deps に持つステージの lock が食い違う。ADR-0254 の移行（案 M1）と同じく、例の dvc.lock の
`bin/<名前>` の deps の項目を、lock が記録した版（DVC のキャッシュに md5 で残っているもの）の正規形の項目（パス、md5、大きさ）に
書き換える。キャッシュに無い版と、正規形を求められない版の項目は書き換えない。xlc-log の前の版はイベントの名前に行を含むので、
移行の一度だけ、その数字も 0 にする版の normalize_binary.py（b6f975eb2 の `unnamed_events`）で求めた。移行の後は使わないので、
8bcd3fecf で外した（CLAUDE.md の方針の、古い形式の読み込みを残さない）。
書き換えるのは実行ファイルの deps の項目だけで、outs、cmd、ほかの deps は変えない。

### 帰結

- xlc-log の追加とマクロの置き換えで `bin/` は一度変わり、計算のステージは正規形の差で changed になる（ADR-0130 のとおり
  再実行しない）
- 例のステージに Rust の実行ファイルを足すときは、cmd は `bin/<名前>`、deps は `binary-forms/<名前>.md5` を書く。`build` を
  実行したら `binary-forms` も実行し、`binary-forms/` とリポジトリの直下の `dvc.lock` をコミットする（/dvc-relock）
- ステージのパニックの文は、`bin/` の行と列をそのまま出す
- 使わないモジュールで実行ファイルが変わる件（#0117）は、別の原因。xlc-format に #0117 のモジュールを足すと、番号を 0 にした
  コピーでも 4 つの実行ファイルが変わり、大きさも増える（xs +448、lattice +832、core +2,752、chain +256 バイト。「確認」の (d)）ので、行の番号では
  なく、コードかデータが実際に増えている。この ADR では扱わない
- ツールチェーンか `tracing` を上げてレイアウトが変わると、`binary-forms` が normalize_binary.py のエラーで止まる。そのときは
  整合の条件を新しいレイアウトに合わせ、「確認」の (c) をやり直す

### 確認

(a)〜(e) は main の 92b0751ef からのブランチ（fix/0110-line-positions-in-binaries）、(f) と (g) は main の 9f1a69e72 に rebase した
同じブランチ（1fd31e8f6）で、default の devShell の `scripts/build-binaries.sh`、2026-11-19。(h) はブランチの先端のテスト。

- (a) 0 にした数（xlc-log の前の main の実行ファイル）: xs の Location 1,255、`Metadata` 16、lattice 4,087 と 142、core 2,525 と
  33、chain 1,024 と 3。lattice の `enabled!` は 12
- (b) 「背景」の編集（moc.rs の先頭の 1 行、関数の入れ替え）と #0117 のモジュールのそれぞれの前後で、コピーの md5 は、moc.rs の
  1 行と入れ替えで 4 つとも同じ、#0117 のモジュールで 4 つとも違った（xlc-log の前の main の実行ファイルに、イベントの名前の
  数字も 0 にする版の normalize_binary.py を当てた）
- (c) xlc-log の後、crates の追跡しているすべての .rs と crates/log/src/lib.rs の先頭に注釈の行を 1 つ足してビルドすると、
  コピーの 4 つの md5 が足す前と同じ（xs b5c605b7、lattice 98a772e0、core e6fc72e3、chain 16a674a2。今の `binary-forms/` の値）
- (d) #0117 のモジュールの大きさ: 足す前と後の実行ファイルの大きさは xs 3,072,016 と 3,072,464、lattice 11,234,944 と
  11,235,776、core 7,443,832 と 7,446,584、chain 1,953,016 と 1,953,272 バイト
- (e) 定数を 1 つ変えた実行ファイル（crates/lattice/src/moc.rs の `POLAR_SIN` の 0.166_648 を 0.166_649 に）の lattice の正規形は
  2d50279f で、変える前（98a772e0）と違う
- (f) 計算が変わらないこと: examples/m0 の xs、lattice、core を、main の `bin/`（xs 28e28365、lattice c3b3bba1、core 22e10bb8）と
  このブランチの `bin/`（cargo の出力のまま。xs 66934bbc、lattice 9a4b0ad3、core c7002c64）で流し、出力の md5 が同じ
  （wims69.msgpack 313819dc、pincell-fg.msgpack 09abc432、lattice.out 53445db5、core.out b14708c7。前の 3 つは
  examples/m0/dvc.lock のとおり）。ログは時間の値と、核種の並列の順だけが違う
- (g) 移行: 例の lock の `bin/` の項目は 1,048 で、どれも 4 つの版（xs f98234e2、lattice b54520ab、core d6b0932c、chain 4a9d4768）の
  どれかを記録し、4 つともキャッシュにあって、すべてを書き換えた。DVC の API で解決したステージと lock の deps を静的に照らすと、
  changed は移行の前も後も同じ 1,520 本（例とリポジトリの直下の 1,904 本のうち）。実行ファイルの差だけで changed だった 1,046 本は、
  どれも changed のまま（記録した版はイベントの名前に行を含み、今の版と正規形が違う。位置だけの差で up to date に戻ったものは
  無い）。examples/core の `compare` は deps の rods_cases.py の正規形が変わったので、上流の lock の出力で再実行し、compare.out は
  同じ（8fa7423d）
- (h) テスト: scripts/test_normalize_binary.py（小さな ELF で、0 にするバイトとそれ以外を変えないこと、ほかのバイトを 1 つ
  変えると正規形が変わるかエラーで止まること、行をずらした同じ構造から同じ正規形になること、節の端の Location と名前、偶然 NUL
  で終わる Metadata のファイル名（今のレイアウトの後ろのフィールドでは Location に当たらないことと、両方に当たると止まること）、
  エラーの条件、例の deps と lock の正規形）と、crates/log/src/tests.rs（水準ごとのマクロの名前、target、フィールドと書式の形）
