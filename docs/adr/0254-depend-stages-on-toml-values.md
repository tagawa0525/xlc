---
status: proposed
date: 2026-11-19
requires: [ADR-0130, ADR-0159, ADR-0210, ADR-0212]
supersedes: []
superseded-by:
issues: ["#0145", "#0295"]
---

# ADR-0254: DVC のステージは、入力の TOML を tomllib で読んだ値の正規形で deps にする

## 背景

ADR-0159 は Rust のステージの deps をビルドした実行ファイルに、ADR-0210 は Python のファイルを docstring とコメントを除いた
正規形にし、計算に関わらない変更で計算のステージが走らないようにした。入力の TOML は、ファイルそのものを deps にしたまま
残っていた。例のステージの deps にある TOML は 436 ファイル（ステージの出力の `xs-range.toml` の 4 つを除く）で、それを deps に
持つステージは 1,091 ある（2026-11-19、main 92b0751ef の dvc.yaml を DVC の API で解決して数えた）。

- コメントだけの修正で、ライブラリ、NJOY、OpenMC（数時間）のステージが changed になった: #0141 の履歴の書き換え、4aa33c0c9
  の語の修正（p2 の `openmc-index` が changed のまま）、#0291 の語の修正
- #0291 では、DVC のステージの入力だからという理由で 555 のファイルを語の検査から外した。ユーザーはこれを悪手とし、注釈の
  変更がステージに及ばない仕組みがあるはずだと問うた（ユーザー 2026-11-19）。仕組みは Rust（ADR-0159）と Python（ADR-0210）に
  しか無く、TOML に無かった
- ADR-0212 は、コメントだけの差の TOML を `dvc commit` で付け直してよい条件（(b)）を決め、TOML の正規形（案 B5）は決めずに #0145 に
  残した。付け直すたびに、4 つの条件を確かめる手順が人に残る
- 比較のステージの JSON は、読んだ入力の md5 を `inputs` に書く。コメントだけの修正でも JSON の md5 が古くなり、#0291 では 4 本の
  再実行を落とした（#0295）

判断の軸（ADR-0210 と同じ）:

- ステージが走るかどうかを、ステージが実際に読む値で決めること
- 正確さを落とさないこと。値の変更（数、文字列、型、表のキーの順）は必ず下流を changed にする。正規形を作れない入力は
  エラーで止め、元の md5 に黙って戻さない
- 人の判断を手順に残さないこと

## 検討した案

仕組み:

- 案 B5（採用）: TOML の deps を、値の正規形を書くステージの出力に替える（ADR-0210 の延長）。判断は正規形のバイトの同一で
  機械的に決まる
- 案 B2: ADR-0212 の (b) の条件で、コメントだけの修正のたびに `dvc commit` で付け直す。条件の 2（親で changed でない）と 4（親の版の
  バイトが出力に無い）を毎回確かめる手順が残る。比較のステージの JSON は入力の md5 を出力に書くので、条件 4 で付け直せない
- 案 P: TOML を、キーを挙げない `params`（`- file.toml:`）でファイル全体の値の deps にする。DVC は値を Python の `==` で比べるので、
  `1` と `1.0`、`true` と `1` を区別せず、キーの順も比べない。dvc.lock にすべての値が入る（p2 の xs.toml は 80 KB で、545 ステージが
  読む）。DVC は params の TOML を tomlkit で読み、表の見出しの数の 2 乗の時間がかかる（#0268）
- 案 N: 何もしない（ADR-0130 のとおり、結果が要るときに再実行する）。#0291 のように、語を直すために入力を検査から外す対処が残る

正規形の中身:

- 案 J（採用）: 標準ライブラリの `tomllib` で読んだ値を、`json.dumps`（キーの順を保つ、`ensure_ascii=False`、`indent=1`、
  `allow_nan=False`）で書く。コメント、空白、改行、引用符と文字列の書き方（基本と literal、エスケープ）、数の書き方（`1e10` と
  `1.0e10`、`0.10` と `0.1`、`1_000` と `1000`、16 進）、表の書き方（見出し、インライン、ドットのキー）は値に残らない。整数と
  浮動小数、真偽値と数、文字列と数は JSON でも別の書き方になる。浮動小数は Python の `repr`（往復で同じ値になる最短の 10 進）で、
  値の同じ 2 つの書き方は同じ文字列になる。`tomllib` と `json` の書き方は Python の版で変わらない（ADR-0210 の `ast.unparse` と
  違う）
- 案 J-sort: 案 J に `sort_keys=True`。キーの順の差を落とす。核種の並びの順で出力を書くプログラムがあるので（ADR-0212 の (b) の
  条件 3 が順を比べた理由）、順は値に数える
- 案 W: 値を TOML に書き直す（`tomli_w` か tomlkit）。標準ライブラリに無い依存を足し、JSON と比べた読みやすさの差は小さい

正規形の置き場とステージの形は、ADR-0210 の案 G と案 F と同じ理由で同じにする（git で追跡し、worktree で `dvc checkout` しなくても
揃う。並びは 1 か所で、変えたファイルのステージだけが走る）。

比較の JSON の `inputs`（#0295）:

- 案 I1（採用）: TOML の入力は、正規形の md5 を書く（`scripts/normalize_toml.py` の `input_md5`）。ステージの deps が記録する md5 と
  同じ値になる。検査は、追跡している JSON の `inputs` が、それを書いたステージの dvc.lock の deps（TOML は正規形）と同じ md5 を
  持つことを確かめる
- 案 I2: `inputs` の md5 を、今のファイルの md5 と比べる検査を置く（#0295 の書いた形）。上流の値を変えた後、ステージを
  再実行するまで（ADR-0130 のとおり changed のまま残してよい間）検査が止まる。追跡しない出力（テーブルの msgpack）は worktree に
  戻すまで比べられない。JSON と lock が一致すれば、JSON と今のファイルの食い違いはステージの changed と同じで、`dvc status` が示す
- 案 I3: `inputs` から TOML を外す。読んだ入力の記録が欠ける

移行（「移行」の項）:

- 案 M1（採用）: dvc.lock の TOML の deps の項目を、lock が記録した版の正規形の項目に書き換える
- 案 M2: ADR-0159 と ADR-0210 のとおり、分岐点で changed でなかったステージを `dvc commit` で付け直す。出力を worktree に
  戻す必要があり（`dvc commit` は出力の hash を取り直す）、実行ファイルの差で changed のステージ（2026-11-19 の main で、静的に
  数えて 1,047。下の「確認」）は付け直せないので、新しい deps のパスの差が changed の理由に加わる

## 決定と理由

案 B5、J、I1、M1 にする。

- `scripts/normalize_toml.py <元の .toml> <正規形>` が正規形を書く。構文の誤りは `TOMLDecodeError`、TOML の日時は `TypeError`、
  `nan` と `inf` は `ValueError` で止まる。日時と非有限の数は 436 のファイルのどれにも無い（2026-11-19 に `tomllib` で読んで
  確かめた）。日時を文字列で書くと、同じ文字列の値と区別できない。`tomllib`（TOML 1.0）が読めない入力（TOML 1.1 だけの書き方）も
  止まる
- リポジトリの直下の `dvc.yaml` の `toml-values`（foreach）が、どれかの例のステージの deps にある TOML のうちステージの出力でない
  もののすべてについて、`toml-values/<リポジトリの直下からのパス>.json` に正規形を書く（`cache: false`、git で追跡）。例の
  `dvc.yaml` の deps は TOML の代わりにその正規形を書き、cmd は元の TOML を読む。参照解、ライブラリ、lattice、core、比較のどの
  ステージも同じ扱いにし、例外を作らない
- ステージの出力の TOML（`xs-range.toml`）は、そのまま deps にする。バイトを書くのはステージのプログラムで、プログラムが
  変わらなければ同じバイトになる
- `params`（核種ごとのステージの `xs.toml: [nuclide.X]` など）と `vars` は変えない。DVC が挙げたキーの値を比べるので、コメントには
  よらない
- 入力が別のファイルを読む場合（`library`、`chain`、`[history] file`、巨視的断面積の TOML）は、これまでどおり読むファイルを
  それぞれ deps に並べ（ADR-0159）、TOML ならその正規形を並べる。正規形は参照をたどらない
- `dvc.yaml` を書くツール（`examples/core/rods_cases.py`、`examples/vera/p3/p3_inputs.py`、`p5_inputs.py`、`p6_inputs.py`）も正規形の
  パスを書く
- 比較のスクリプトは、TOML の入力の md5 を `input_md5` で書く（JSON の `inputs` に TOML を書く 3 つ:
  `compare_branch_cross_terms.py`、`compare_grid_interpolation.py`、`compare_p6.py`）。`examples/record_provenance.py` の照合も、
  TOML の入力を正規形の deps と比べる

ADR-0159 の deps の規則（ステージの入力と、ステージが呼ぶ実行ファイル）のうち、入力の TOML を並べる部分を、その正規形に
置き換える。

入力のバイトが出力に入るほかの経路（ADR-0212 の (b) の条件 4 が挙げたもの）:

- 入力をコピーするスクリプト: ステージの cmd が呼ぶ Python（python-ast の並びの 38 ファイル）は、入力の TOML を `tomllib` で読み、
  コピーしない（2026-11-19 に grep した）
- lattice のチェックポイントの `Fingerprint`（ADR-0196）: 入力のテキストのバイトを含むが、チェックポイントはステージの出力でなく、
  計算を終えると消える。コメントを直すと、止めた計算を途中から続けられず、始めからになる（安全な側）。Rust は変えない

### 移行

deps のパスを替えると、TOML を deps に持つステージの lock が食い違う。`dvc commit` の代わりに、例の dvc.lock の deps の項目を
書き換える（案 M1）:

- ステージごとに、ステージの出力でない TOML の deps の項目を、正規形のパスと、lock が記録した版（md5 の同じ blob を git の履歴と、
  #0141 の書き換えの前の控え `~/.local/share/history-rewrite/xlc-2026-10-21/original.git` から探す）の正規形の md5 と大きさに
  替える。版が見つからない項目は替えず、そのステージは changed のまま残す
- lock が記録した TOML が今のファイルと同じなら、書き換えた項目は `dvc commit` が書くものと同じになる。違えば、記録した版の
  正規形と今の正規形が比べられる。値の差なら changed のまま、コメントだけの差なら changed でなくなる（ADR-0210 の移行の 2 つ目と
  同じで、記録した deps と今の deps の値が同じなので、ADR-0130 の案 L が退けた形に当たらない）
- 書き換えるのは TOML の deps の項目だけで、outs、cmd、ほかの deps（実行ファイル、上流の出力）は変えない。実行ファイルの差で
  changed だったステージは、その理由で changed のまま残る

比較のスクリプトを変えたステージ（JSON の md5 が正規形の md5 に変わる）と、`rods_cases.py` を deps に持つ examples/core の `compare` は、
数秒から数分の比較なので再実行する。

### 帰結

### 確認
