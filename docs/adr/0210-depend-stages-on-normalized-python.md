---
status: accepted
date: 2026-10-28
requires: [ADR-0130, ADR-0136, ADR-0159, ADR-0160]
supersedes: []
superseded-by:
amends: [ADR-0159]
issues: ["#0145"]
---

# ADR-0210: DVC のステージは、Python のファイルを docstring とコメントを除いた正規形で deps にする

## 背景

ADR-0159 は Rust のステージの deps をビルドした実行ファイルにし（ADR-0160 で strip）、出力を変えないコードの変更では計算の
ステージが走らないようにした。Python のファイル（比較のスクリプト、OpenMC の入力とライブラリを作るスクリプト）は、ソースその
ものを deps にしていたので、docstring とコメントだけの変更でも md5 が変わり、下流のステージが changed になる（#0145）。
13 の `dvc.yaml` が 33 の .py を deps に持っていた。

- 履歴の書き換え（#0141）と語の直し（ed95d3d44）が docstring とコメントだけを直し、examples/depletion と p2 の OpenMC の
  ライブラリと燃焼（数時間）、examples/vera の OpenMC のライブラリが changed になった
- そのたびに、出力を作り直さず `dvc commit` で lock を付け直した（3517e4666、acfe57b01、c09ed7164 の 3 度）。ADR-0159 の案 B、
  ADR-0130 の案 L、ADR-0143 が退けた形を、限り方を決めずにその場の判断で繰り返していた
- 付け直す前に、compare-1C を目標に流した `dvc repro` が OpenMC の燃焼を始め、止めた（Slurm 454。#0145 の経緯）

判断の軸:

- ステージが走るかどうかを、ステージが実際に行う計算で決めること（ADR-0159 と同じ）。docstring とコメントは計算に関わらない
- 正確さを落とさないこと。コードの変更（式、名前、定数）は必ず下流を changed にする。正規形を作れないファイル（構文の誤り）は
  エラーで止め、元の md5 に黙って戻さない
- 人の判断を手順に残さないこと（出力を変えないかの判断を、マージのたびに人が下さない）

## 検討した案

仕組み（#0145 の「参照解の deps の .py の文書だけの変更への案」）:

- 案 a: 文書だけを直したマージの前に `dvc commit` で付け直す手順を /dvc-relock に足す。3 度の手の付け直しを手順にするだけで、
  「文書だけの差か」の判断と確かめ方が人に残る（ADR-0159 の案 B と同じ理由で退ける）
- 案 b（採用）: .py の deps を、docstring とコメントを除いた正規形を書くステージの出力に替える。判断は正規形のバイトの同一で
  機械的に決まる
- 案 c: CLAUDE.md の「`dvc repro` を投げる前に走るステージを数える」規則だけで防ぐ。事故の費用は防ぐが、changed のままの
  ステージが残り、付け直しを繰り返す

正規形の作り方:

- 案 U（採用）: 標準ライブラリの `ast` で読み、モジュール、クラス、関数（async を含む）の本体の先頭にある文字列の定数の式
  （Python が `__doc__` にするもの）を落とし、`ast.unparse` で書き出す。コメント、空行、改行と括弧の付け方、引用符は AST に
  残らない。出力は Python のソースとして読めるので、git の差分でコードの変更を読める
- 案 D: `ast.dump`（`include_attributes=False` で行番号を含めない）。AST をそのまま書くので、unparse の往復に頼らない。差分は
  読みにくい。Python 3.13 で既定の出力（空のリストと None の省略）が変わったように、版による変化は unparse より大きいと見て
  いる（版の間で比べてはいない）
- 案 K: `tokenize` でコメントだけを落とす。docstring の変更（ed95d3d44 と、ADR-0183 に合わせた 3baba57c0 は docstring を
  直した）を拾えない

正規形の置き場:

- 案 G（採用）: outs を `cache: false` にして git で追跡する（`python-ast/<リポジトリの直下からのパス>.ast`）。worktree で
  `dvc checkout` しなくても揃い、`scripts/worktree-restore-outputs.sh` を変えずに済む。commit の hash を含まないので、
  `commit-refs remap` に関わらない。正規形の差分がレビューに出る
- 案 C: DVC のキャッシュに置く（`.gitignore`）。新しい worktree で `dvc checkout` するまで正規形が無く、例のステージの deps が
  消えた状態に見える

ステージの形:

- 案 F（採用）: リポジトリの直下の `dvc.yaml` に、.py ごとの foreach のステージ `python-ast` を置く。並びは 1 か所で、
  変えた .py のステージだけが走る
- 案 S: 1 つのステージで全体をディレクトリの出力に書く。cmd と deps に並びを 2 度書くことになる

## 決定と理由

案 b、U、G、F にする。

- `scripts/normalize_python.py <元の .py> <正規形>` が正規形を書く。docstring でない文字列の式（f-string、2 つ目以降の文の
  文字列）と代入の右辺は残す。本体が docstring だけのクラスと関数は `pass` にする。構文の誤りは `SyntaxError` で止まる
- リポジトリの直下の `dvc.yaml` の `python-ast`（foreach）が、どれかの例の deps にある .py のすべてについて正規形を書く。
  例の `dvc.yaml` の deps は .py の代わりにその正規形を書き、cmd は元の .py を呼ぶ。参照解、ライブラリ、lattice、比較の
  どのステージも同じ扱いにし、例外を作らない。.py が別の .py を import する関係は、これまでどおり各ファイルを deps に並べる
  （正規形もファイルごと）
- `dvc.yaml` を書くツール（`examples/core/rods_cases.py`、`examples/vera/p3/p3_inputs.py`、`examples/vera/p5/p5_inputs.py`）
  も正規形のパスを書く

ADR-0159 の deps の規則（「ステージの入力と、ステージが呼ぶ実行ファイル」）のうち、Python のファイルを入力に並べる部分を、
その正規形に置き換える（`amends`）。ADR-0160 の strip と、ADR-0136 の、OpenMC の実行のステージにスクリプトを deps にしない
決定は変えない。

docstring を落とすと `__doc__` が変わる。deps の 33 のファイルのうち 18 が `__doc__` を使うが、どれも引数の誤りのときの使い方の
表示（`SystemExit(__doc__)` と、`model_run_extract` と `run_or_extract` に渡して `SystemExit` にするもの）で、ステージの出力に
入らない。

`ast.unparse` は Python の版で書き方が変わりうるので、正規形のバイトは Python の版に依る。default と references の devShell は
同じ nixpkgs の `python3`（3.14.7。`flake.lock`）で、`flake.lock` を変えるまで版は変わらない。`python-ast` の deps に
`flake.lock` を入れない: 入れると、Python の版を変えるたびに正規形を書き直し、書き方が変われば参照解までを実行する（この ADR が
防ぐ再実行）。版が変わったあとは、.py を変えたときにその正規形だけが書き直される。

### 移行

deps のパスを替えると、.py を deps に持つすべてのステージの lock が食い違う。再実行せず `dvc commit` で付け直す。対象は次の 2 つに
限る（ADR-0159 の移行と同じ限り方に、2 つ目を足す）:

- 分岐点（main 84517a7f5）の `dvc status` で changed でなかったステージ。lock の deps の .py は今のファイルと同じなので、
  正規形も同じである
- 分岐点で changed だった原因が .py の deps だけで、lock が記録した版（md5 が同じ blob を git の履歴と、#0141 の書き換えの前の
  控え `~/.local/share/history-rewrite/xlc-2026-10-21/original.git` から探す）と今の版の正規形が一致するステージ。コードが同じ
  なので、ADR-0130 の案 L と ADR-0143 が `dvc commit` を退けた理由（記録した deps のコードで実行していない出力を lock が指す）に
  当たらない

分岐点で別の理由（実行ファイル、入力の TOML、上流の出力、出力の欠け）で changed だったステージと、.py の正規形が lock の版と違う
ステージは付け直さず、changed のまま残す。

## 確認

テスト（`scripts/test_normalize_python.py`、`scripts/run-python-tests.sh` が実行する）:

- `InvariantTest`: docstring の変更、コメントの変更、空行と改行と引用符の変更で正規形が変わらない
- `CodeChangeTest`: 式、名前、定数（112.5 → 112.47）の変更で正規形が変わる
- `DocstringTest`: モジュール、クラス、メソッド、ネストした関数、async の関数の docstring を落とし、f-string と 2 つ目の文の
  文字列と代入の右辺を残す。本体が docstring だけのクラスと関数も、Python として読める正規形になる
- `FailureTest`: 構文の誤りは `SyntaxError` で止まる
- `RepositoryDepsTest`: 例の `dvc.yaml` の deps に .py が無いこと、deps の正規形がすべて `python-ast` の並びにあること、並びの
  正規形がどれかの例の deps にあること

2026-10-28、このブランチ（main 84517a7f5 の上）で測った:

- 33 のファイルのすべてで、正規形を読み直した AST が、docstring を落とした元の AST と一致した（`ast.dump` で比べた。unparse の
  往復で区別が失われていない）
- `nix develop -c dvc repro -f -s dvc.yaml:python-ast`（33 ステージ）は 3.2 秒
- ed95d3d44 の直前（586f2c920）と今の examples/openmc_input.py は、正規形が一致した
- 仕組みの試し: examples/openmc_input.py の docstring に 1 文字足して `python-ast` を実行すると、正規形は変わらず、`dvc status`
  の changed は 468 のまま増えなかった。コードの 1 行（`required - table.keys()` を `set(required) - table.keys()`）を変えると、
  正規形が変わり、changed は 730（262 増。p2 の 91、examples/depletion の 69、vera の 68、endfb70/openmc の 34）になった。
  どちらも試したあとに戻した
- 移行: 分岐点の `dvc status`（worktree で `scripts/worktree-restore-outputs.sh` の後）の changed は 687 項目。替えた後に新たに
  changed になった 110 ステージ（理由は .py と正規形の deps だけ）と、正規形が lock の版と一致した 219 ステージ（p2 の
  openmc-neutron 90 と openmc-thermal、examples/vera の OpenMC のライブラリ 68、endfb70/openmc のライブラリと openmc-model 34、
  p3 と endfb70/p3 の compare 6、p5 の比較 20）を `dvc commit` した。`dvc.lock` の差は deps だけで、outs と cmd は変わらなかった。
  付け直しの後の changed は 468 項目で、分岐点の一覧から 219 を除いたものと一致し、各ステージの .py と正規形以外の理由は
  分岐点と同じだった。正規形が lock の版と違ったのは examples/core/rod_table.py と rods_cases.py（書き換えの前の版とコードが
  違う）で、その 10 ステージ（core の table と compare、p5 の hzp-table の 8）は changed のまま残した

## 帰結

- docstring とコメントだけの変更では、`python-ast` のそのファイルのステージだけが走り（数秒）、正規形が変わらなければ計算の
  ステージは走らない。`dvc commit` での付け直しは要らない
- 正規形はコードの変更をすべて拾うので、出力を変えないコードの整理（名前の変更、関数の分割）は、これまでどおり下流を changed に
  する（ADR-0130 のとおり、結果が要るときに再実行する）
- `dvc status examples/<問題>/dvc.yaml` は上流の `python-ast` を示さない（ADR-0159 の帰結の `build` と同じ）。.py を変えたかは
  `dvc status`（引数なし）か `--with-deps` で見る。`dvc repro` は上流の `python-ast` を先に実行する。rebase のあとの手順
  （/dvc-relock）は、`build` に続けて `python-ast` を実行する
- 例のステージに .py を deps として足すときは、`python-ast` の foreach にも足す。`RepositoryDepsTest` が食い違いを止める
- `python-ast/` の正規形は git で追跡する生成物で、手で直さない。.py のコードを変えたブランチは、`python-ast` を実行して正規形と
  リポジトリの直下の `dvc.lock` をコミットする（実行しないまま残すと、`build` と同じく、次の `dvc repro` が書き直す）
- Python の版を変える（`flake.lock` の更新）と、次にそのファイルを変えたときに、書き方の差で正規形が変わりうる。版の間で
  `ast.unparse` の出力がどれだけ変わるかは確かめていない
