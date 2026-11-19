---
status: accepted
date: 2026-11-19
requires: [ADR-0159, ADR-0160]
supersedes: []
superseded-by:
issues: ["#0110", "#0117"]
---

# ADR-0255: ステージの実行ファイルから、ソースの行と列の番号を除く

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

判断の軸: 正確さ（計算の結果をビットで変えない。誤ったバイトを消さない）、実行速度（リリースの最適化を変えない）、保守の
容易さ。エラーを握りつぶさない（パニックの文と、場所を引く方法を残す）。

## 決定の要因

- 注釈の行、空行、関数の並べ替えでは、ステージが依存するものを変えない
- ステージが呼ぶ実行ファイルと、`dvc.lock` が記録する実行ファイルの md5 が同じものを指す（ADR-0159 の判断の軸）
- 位置の情報を除く処理が、構造体のレイアウトの変化（ツールチェーンや `tracing` の更新）で黙って外れない
- Rust は rust-toolchain.toml の stable の 1.99.0（flake.nix は rust-overlay でこの版を使う）

## 検討した案

パニックの位置:

- 案 A: コンパイラで位置を外す（`-Z location-detail=none`）。1.99.0 の stable では「the option `Z` is only accepted on the
  nightly compiler」で止まり、`-C location-detail` は無い。nightly を使う（または `RUSTC_BOOTSTRAP`）と、ステージの実行ファイルを
  作るツールチェーンが安定版でなくなるので採らない
- 案 B: `--remap-path-prefix`。ファイル名を変えるだけで、行と列は残る
- 案 C: パニックの経路を位置の要らない形に寄せる（添字の代わりに `get` とエラーを返す）。lattice だけで Location が 4,087 ある
  （「確認」の (a)）。添字の検査はコンパイラが入れるので、本番のコードのすべての添字を書き換えることになり、読みにくく、
  最適化も変わりうる
- 案 D（採用）: `build` のステージが、リンクした実行ファイルの Location の行と列を 0 にしたコピーを `bin/` に置く
  （scripts/normalize_binary.py）。計算のコードは Location を読まない（読むのはパニックの文だけ）ので、計算は変わらない
- 案 E: 実行ファイルは番号を残したまま `bin/` に置き、番号を消したコピーの md5 を別のファイルに書いて、計算のステージの deps を
  それにする。パニックの文に行が残るが、例の 20 の dvc.yaml の 236 行の `bin/` の deps を書き換え、`dvc.lock` が計算に使った
  実行ファイルでなく別のファイルを指す（ADR-0159 の案 D が退けた形と同じく、`dvc checkout` で出力を作った実行ファイルを
  戻せない）
- 案 F: 番号の場所をレイアウトから探さず、行を 1 つずつずらしたソースでもう一度ビルドし、違ったバイトを消す。レイアウトに
  依らないが、ビルドが倍になり、イベントの名前の桁が変わる位置（9、99、999 行のイベント）では 2 つのビルドの並びが揃わない
- 案 G: 何もしない（#0110 の完了の条件の 2 つ目）。リファクタリングと注釈の修正のたびに、出力が同じステージを再実行するか
  古いまま残す

`tracing` の `Metadata` の行:

- 案 H（採用）: 案 D と同じコピーで 0 にする。ログの出力（`tracing_subscriber::fmt` の既定）は行を出さないので、ログは変わらない
- 案 I: `tracing` の機能で行を出さない。`release_max_level_*` はイベントごと消し、行だけを外す機能は無い（0.1.44 のマクロは
  `line: Some(line!())` を常に書く）

`tracing` のイベントの名前:

- 案 J（採用）: イベントのマクロを xlc-log のマクロ（`xlc_log::info!` など）に替え、名前を `"event <ファイル>"`（行なし）にして
  `tracing::event!` に渡す。名前の長さが行で変わらないので、定数の節の並びも変わらない。clippy.toml の `disallowed-macros` で
  `tracing::trace`〜`tracing::error` を止める（`name:` を付けて `tracing::info!` を呼ぶ形は、clippy が展開の中の呼び出しも
  止め、展開の中の `#[allow]` も効かなかったので、`tracing::event!` を呼ぶ）
- 案 K: 案 D と同じコピーで数字を 0 にする。上の nem.rs の 474 行のように、定数の節の並びが中身で決まると戻らない。桁が変わる
  ときも、文字列の長さが変わって後ろの並びがずれる
- 案 L: 呼ぶ場所ごとに `name:` を書く。113 の呼び出しに同じ引数を足し、書き忘れを機械で止められない

`tracing::enabled!` は `name:` を取らない（0.1.44）。lattice の 12 か所の名前は 37〜44 バイトで、定数の節の大きさ（4、8、16、
32 バイト）に合わないので、数字だけをコピーで 0 にする（案 K を `enabled!` に限って使う）。自前の `enabled!` は `tracing` の隠した
マクロ（`callsite2!`）を使うことになるので作らない。

## 決定と理由

案 D、H、J を採る。

- scripts/build-binaries.sh は、ビルドした `target/dvc/release/<名前>` を scripts/normalize_binary.py に通して `bin/<名前>` に
  書く。コピーは、Location の行と列、`Metadata` の行、`enabled!` の名前の行の数字を 0 にし、ほかのバイトを変えない
- 構造体のレイアウトは `repr(Rust)` で保証されないので、値どうしの整合で見つける。Location は「NUL で終わる `.rs` の文字列を
  指す &str の後に、0 でない行と列」。`Metadata` は「名前 `"event <ファイル>"` の &str の直前に Some(行)、名前から決まった
  距離に同じファイルの &str」で、距離はイベントから決めてスパンにも使う。Metadata のファイル名を先に決め、Location はそれ以外から探す。
  Location が 1 つも無い、名前の直前が Some でない、距離がイベントごとに違う、1 つの &str が Location と Metadata の
  ファイル名の両方に当たる、のどれかでエラーで止める（ツールチェーンか `tracing` のレイアウトが変わった。黙って消し残すと
  注釈の行でステージが古くなるのに気付けず、別のバイトを消すと計算が変わりうる）。名前を指定しないイベント
  （`"event crates/…:<行>"`）が残っていてもエラーで止める
- 計算のステージの deps と dvc.yaml は変えない。`bin/` が計算に使う実行ファイルそのもので、`dvc.lock` の md5 もそれを指す

失うもの: ステージのパニックの文の位置が `crates/lattice/src/moc.rs:0:0` の形になる（ファイル名とパニックの文は残る）。行は、
`bin/` を作ったコミットで scripts/build-binaries.sh を流した `target/dvc/release/<名前>`（番号を残した同じビルド）で、同じ
入力を流すと出る。計算のコードは同じなので、同じ場所で止まる。`bin/` をキャッシュから戻した worktree には `target/dvc` が
無いか別のコミットのものなので、そのコミットでビルドし直す。

守れる範囲: 注釈と空行の追加と削除、行の中の位置（列）の変化、ファイルの中の関数の並べ替え（1 組で確かめた）。守れないもの:
コードの変更（別の関数の変更を含む）、ファイルの名前と、関数のファイルの間の移動（Location のファイル名が変わる）、
lattice の `enabled!` の 12 か所の行の桁の変化（文字列の長さが変わる）。守れないときは実行ファイルが変わり、ステージが
古くなる（安全な側）。

### 帰結

- `build` の deps に scripts/normalize_binary.py の正規形（ADR-0210 の python-ast）を足す。xlc-log の追加とマクロの置き換えで今の `bin/` は一度変わり、
  計算のステージは `bin/` の差で changed になる（ADR-0130 のとおり再実行しない）
- 使わないモジュールで実行ファイルが変わる件（#0117）は、別の原因。xlc-format に #0117 のモジュールを足すと、コピーの後も
  3 つの実行ファイルの大きさが変わる（xs +448、lattice +832、core +2,752 バイト。「確認」の (d)）ので、行の番号ではなく、
  コードかデータが実際に増えている。この ADR では扱わない
- ツールチェーンか `tracing` を上げてレイアウトが変わると、`build` が normalize_binary.py のエラーで止まる。そのときは
  整合の条件を新しいレイアウトに合わせ、「確認」の (c) をやり直す

### 確認

いずれも main の 92b0751ef からのブランチ（fix/0110-line-positions-in-binaries）、default の devShell の
`scripts/build-binaries.sh`、2026-11-19。

- (a) コピーの数（xlc-log の前の main の実行ファイル）: xs の Location 1,255、`Metadata` 16、lattice 4,087 と 142、core 2,525 と
  33、chain 1,024 と 3。lattice の `enabled!` は 12
- (b) 「背景」の編集（moc.rs の先頭の 1 行、関数の入れ替え）と #0117 のモジュールのそれぞれの前後で、コピーの md5 は、moc.rs の
  1 行と入れ替えで 4 つとも同じ、#0117 のモジュールで 4 つとも違った（xlc-log の前の main の実行ファイルに、イベントの名前の
  数字も 0 にする版の normalize_binary.py を当てた）
- (c) xlc-log の後、crates の追跡しているすべての .rs と crates/log/src/lib.rs の先頭に注釈の行を 1 つ足してビルドすると、
  `bin/` の 4 つの md5 が足す前と同じ（xs b5c605b7、lattice 98a772e0、core e6fc72e3、chain 16a674a2）。節の端と
  Metadata のファイル名を直した normalize_binary.py（549d77a63）を同じビルドに当てても、4 つの md5 は同じ
- (d) #0117 のモジュールの大きさ: 足す前と後の実行ファイルの大きさは xs 3,072,016 と 3,072,464、lattice 11,234,944 と
  11,235,776、core 7,443,832 と 7,446,584 バイト
- (e) 計算が変わらないこと: examples/m0 の xs、lattice、core を、main の `bin/`（xs 28e28365、lattice c3b3bba1、core 22e10bb8）と
  (c) の `bin/` で流し、出力の md5 が同じ（wims69.msgpack 313819dc、pincell-fg.msgpack 09abc432、lattice.out 53445db5、
  core.out b14708c7。前の 3 つは examples/m0/dvc.lock のとおり）。ログは時間の値と、核種の並列の順だけが違う
- (f) テスト: scripts/test_normalize_binary.py（小さな ELF で、0 にするバイトとそれ以外を変えないこと、行をずらした同じ
  構造から同じコピーになること、節の端の Location と名前、偶然 NUL で終わる Metadata のファイル名、エラーの条件）と、
  crates/log/src/tests.rs（水準ごとのマクロの名前、target、フィールドと書式の形）
