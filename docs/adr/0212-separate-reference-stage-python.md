---
status: proposed
date: 2026-10-28
requires: [ADR-0124, ADR-0130, ADR-0159, ADR-0210]
supersedes: []
superseded-by:
amends: [ADR-0130]
issues: ["#0183"]
---

# ADR-0212: 参照のステージの Python を検証用のパッケージから分け、deps の付け直しは値の同じ TOML に限る

## 背景

NJOY と OpenMC を使うステージ（参照のステージ）は、`references` の devShell の store のパスを書いた `references-shell.txt` を
ツールの deps にする（ADR-0130 の案 D）。flake.nix の `pythonPackages`（検証用のツールの Python のパッケージ）は default の
devShell の Python と、`references` の devShell の Python（`referencesPython`）の両方が使っていたので、検証用のスクリプトの
ために足したパッケージも `references` の devShell の store のパスを変える。main の 6634685a6（pyyaml を足した。#0182）で、
`references` の devShell は `qdaa3j64…` から `jrv2kvr1…` に、`referencesPython` は `zwhli51v…` から `wxghvalz…` に変わり、
`openmc`（`r7librmf…`）と `njoy2016`（`jnvmgygh…`）は変わらなかった（`nix eval` で 6634685a6 の親と 6634685a6 を比べた。
2026-10-28）。追跡している 4 つの `references-shell.txt` は前のパスのままで、`references-shell` のステージは `always_changed`
なので、目標の上流にそれを含む `dvc repro` が参照解（数時間）までを再実行する状態だった（#0183）。

もう 1 つ、結果を変えない deps の変更で重いステージが changed になったとき、出力を作り直さずに `dvc commit` で lock を
付け直してよい条件が決まっていない。ADR-0130 の案 L、ADR-0143 の帰結、ADR-0159 の案 B は、`dvc commit` を「記録した deps の
コードで実行していない出力を lock が指す」として退けた。それでも 3517e4666、acfe57b01、c09ed7164 の 3 度、その場の判断で
付け直した（#0145）。.py の docstring とコメントは ADR-0210 の正規形で付け直しが要らなくなったが、TOML のコメント
（ed95d3d44 で p2 の `openmc-index` が changed のまま）と、実行ファイルのバイトだけが変わる refactor（ADR-0211 のブランチは、
ライブラリの作り直しを避けるために `bin/xs` を変える refactor を戻した）が残る。

判断の軸:

- ステージが走るかどうかを、ステージの計算を決めるもので決めること（ADR-0159、ADR-0210 と同じ）
- 正確さを落とさないこと。ステージの環境が変わったら、結果が同じと機械で示せない限り、changed にする
- 人の判断を手順に残さないこと。付け直してよいかは、機械で確かめられる条件で決める
- 黙ったフォールバックを足さず、食い違いはエラーで止めること

## 検討した案

参照のステージのツールの deps（決定の (a)）:

- 案 A1: 今のまま（`references` の devShell の store のパスを deps にし、Python のパッケージの並びを default と共有する）。
  検証用のツールのためにパッケージを足すたびに、参照解が changed になる（#0183）
- 案 A2: deps を、devShell でなくツールの store のパス（`njoy2016`、`openmc`、`referencesPython`）だけにする。6634685a6 で
  変わったのは `referencesPython` そのものなので、並びを共有したままでは今回の種類を防げない。並びを分けた後（案 A3）は、
  devShell のパスとの差は `mkShell` の stdenv（PATH の bash や coreutils）と `shellHook`（`unset PYTHONPATH`。どの
  Python のパッケージを読むかを決める）だけで、どちらもステージの環境の一部である。外すと、それらの変更で参照解が走らない
- 案 A3（採用）: 参照のステージの Python のパッケージの並び（`referencesPythonPackages`）を、検証用のツールとテストの並び
  （`pythonPackages`）から分け、deps は devShell の store のパスのままにする。検証用のツールだけが使うパッケージは
  `pythonPackages` にだけ足し、`references` の devShell を変えない。テスト（OpenMC の Python API と pyyaml の両方を使う）は、
  新しい devShell `python-tests` で実行する
- 案 A4: deps を変えず、devShell が変わったら (b) の条件で付け直す。パッケージを足すと、参照のスクリプトが読むパッケージの
  任意の import（`try: import …`）が変わりうるので、環境の差が結果を変えないことを機械で示せない

結果を変えない deps の変更での付け直し（決定の (b)）:

- 案 B1: 付け直さない（ADR-0130 のまま、結果が要るときに再実行する）。TOML のコメントだけの直しでも、ライブラリや OpenMC
  （数時間）を再実行するか、changed のまま残す。changed のまま残すと、目標の上流に入った `dvc repro` が重いステージを実行する
  （Slurm 454）
- 案 B2（採用）: 機械で確かめられる条件をすべて満たすときだけ `dvc commit` で付け直す。条件は「決定と理由」の (b)
- 案 B3: 人が結果を変えないと判断した変更は付け直す（3517e4666、acfe57b01、c09ed7164 の運用）。判断が人に残り、
  ADR-0159 の案 B と同じ理由で退ける
- 案 B4: 実行ファイルのバイトが変わっても、代表の入力で出力が同じなら付け直す。1 つの入力で同じでも、ほかの入力で同じとは
  言えない。重いステージでは、その入力で流すこと自体が再実行になる。実行ファイルが行の位置で変わる件は #0110 で扱う
- 案 B5: TOML の deps も、.py と同じく正規形を書くステージの出力に替える（ADR-0210 の延長）。付け直しが要らなくなる構造の
  解決だが、TOML を deps に持つすべてのステージの `dvc.yaml` を変える。コメントだけの直しはこれまで 2 度（#0141 の書き換えと
  ed95d3d44）で、この ADR では決めない（#0145 に残す）

## 決定と理由

案 A3 と案 B2 にする。

### (a) 参照のステージの Python の並びを分ける

flake.nix の Python のパッケージの並びを 2 つにする:

- `referencesPythonPackages`: 参照のステージのスクリプト（`references` の devShell で走るもの）が使うパッケージ（numpy、h5py、
  msgpack、lxml、scipy、pandas、uncertainties）。`referencesPython` はこれに OpenMC の Python API を足す
- `pythonPackages`: `referencesPythonPackages` に、検証用のツールだけが使うパッケージ（pyyaml）を足したもの。default の devShell の
  Python が使う
- `testsPython`: `pythonPackages` に OpenMC の Python API を足したもの。devShell `python-tests` に置き、
  `scripts/run-python-tests.sh` がその中でテストを実行する

参照のステージの deps は、これまでどおり `references` の devShell の store のパス（`references-shell.txt`）にする（案 D は
変えない）。この変更で `references` の devShell は pyyaml を足す前と同じ derivation になり、store のパスは 4 つの
`references-shell.txt` と同じ `qdaa3j64…` に戻る。参照のステージの lock は、記録した環境で作った出力をそのまま指すので、
付け直しは要らない。

flake.nix か flake.lock を変えるコミットは、pre-commit の hook `references-shell`（`scripts/check-references-shell.sh`）が、
`references` の devShell の store のパスと、追跡しているすべての `references-shell.txt` を比べ、食い違えば止める。
`references-shell` のステージは `always_changed` で、`dvc repro` が走るまで `references-shell.txt` を書き直さないので、
hook がなければ、devShell の変更が `dvc status` に出ないまま main に入る（#0183 が見つからなかった理由）。
食い違ったら、次のどちらかで直す:

- 足したパッケージを参照のステージが使わないなら、`referencesPythonPackages` でなく `pythonPackages` に足す
- 参照のステージのツール（NJOY、OpenMC、`referencesPython`）を変えたなら、`references-shell.txt` のある例ごとに
  `dvc repro examples/<問題>/dvc.yaml:references-shell`（評価だけで数秒）を実行してコミットする。参照のステージは changed に
  なり、結果が要るときに再実行する（ADR-0130）。(b) の付け直しの対象にしない

### (b) 付け直してよい条件

出力を作り直さずに `dvc commit <yaml>:<ステージ>` で lock を付け直してよいのは、次をすべて満たすステージだけである:

1. 再実行が数分を超えるステージ（ライブラリ、NJOY、OpenMC、xs）。数分以内のステージは再実行する
2. deps を変えたコミットの親で、そのステージが changed でない。lock の deps は、その親のファイルを記録している（ADR-0159 と
   ADR-0210 の移行と同じ限り方）。前から changed だったステージの別の理由を取り込まない
3. 変わった deps がすべて TOML のファイルで、親の版と今の版を `tomllib` で読んだ値が、キーの順序を含めて一致する（両方を
   `json.dumps` した文字列が同じ）。コメント、空白、引用符の書き方だけの差にあたる。順序を比べるのは、核種の並びの順で
   出力を書くスクリプトがあるため
4. 変わった deps の lock の md5 が、そのステージの出力のどのファイルにも現れない（`grep -rqaF <md5> <outs>` が見つけない）。
   入力のファイルの md5 を出力に書くスクリプト（比較のスクリプトの `inputs`）では、バイトの差が出力に入る

付け直したら、`scripts/compare-dvc-lock.py <親> HEAD` で `cmd`、`outs`、`params` が変わっていないことを確かめ、コミットの本文に
付け直したステージと、2〜4 を確かめたコマンドを書く。

次の変更は付け直さず、再実行するか changed のまま残す（ADR-0130）:

- 実行ファイル（`bin/`）のバイトの差。同じ計算かを機械で示せない（ADR-0143 の帰結、ADR-0159 の案 B のまま）
- `references-shell.txt` の差。(a) の後は、参照のステージの環境そのものの差である
- 上流の出力、データ、TOML の値の差
- .py の差は ADR-0210 の正規形が扱う（正規形が変われば、コードが変わっている）

ADR-0130 の「rebase やあとのコミットで deps が変わると…そのままでよい」と、案 L（`dvc commit` での取り込み）を退けた部分の
うち、(b) の条件を満たす付け直しだけを許す形に置き換える。条件を満たす変更では、lock が記録した deps と今の deps から読む値が
同じなので、案 L が退けた「記録した deps のコードで実行していない出力を lock が指す」に当たらない。ADR-0130 のほかの部分と、
ADR-0143 の帰結（実行ファイルが一致しないときに `dvc commit` しない）は変わらない。

## 確認

## 帰結

- 検証用のツールやテストのために Python のパッケージを足しても、参照解は changed にならない。参照のステージのスクリプトが
  新しいパッケージを使うなら、`referencesPythonPackages` に足し、hook が `references-shell.txt` の書き直しを求める
- `references` の devShell の Python は、default の devShell の Python の上位集合でなくなる。default のパッケージを使う
  スクリプトを `nix develop .#references` で実行しない（テストは `python-tests`）
- /dvc-relock は、`references-shell` を対象から外すだけでなく、hook と同じ比較で `references-shell.txt` の中身を確かめる
- TOML のコメントの直しでは、(b) の条件で付け直せる。p2 の `openmc-index`（ed95d3d44 で changed のまま）への適用は別の作業にする
  （#0145）。TOML の正規形（案 B5）は、コメントだけの直しが繰り返すときに #0145 で決める
- 実行ファイルのバイトだけが変わる refactor は、これまでどおり計算のステージを changed にする（#0110）
