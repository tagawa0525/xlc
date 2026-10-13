---
status: accepted
date: 2026-10-14
requires: [ADR-0123, ADR-0130, ADR-0152, ADR-0159]
supersedes: []
superseded-by:
issues: ["#0018"]
---

# ADR-0136: OpenMC の段階を model.xml までと model.xml からに分ける

## 背景

ADR-0130 は OpenMC の参照解を 1 つの段階（`openmc_model.py run`）にした。その deps はスクリプト、共通の部品、入力、
ライブラリ全体で、スクリプトの 1 文字を直しても、モデルが同じなら結果が同じ OpenMC の計算（1B で約 18 分、2B〜2K は
1 問題 30 分ほど）が再実行する。ADR-0130 の「帰結」は、この再実行を避けるために、分岐点から復元する運用を置いた。
根本の直しは、段階を「入力とスクリプト → model.xml」と「model.xml とライブラリ → statepoint」に分けることだった（#0018）。

ADR-0159 は、Rust の計算の段階の deps をビルドした実行ファイルにして、出力を変えないコードの変更で段階が古くならない
ようにした。OpenMC の段階の「実際に呼ぶプログラム」は、OpenMC の実行ファイルと、それに渡す model.xml、ライブラリで、
Python のスクリプトは model.xml を作るツールにすぎない。同じ考えを OpenMC の段階に当てはめる。実際に、標準のタリー
（ADR-0152）の追加とスクリプトの整理で、1B の `openmc-1B`（約 18 分）が古くなった。結果を変えない変更を含む。

分けるときに、model.xml が次の 2 つを含まないことが要る。含むと、モデルが同じでも model.xml が変わり、分けた意味がない。

- 断面積のライブラリのパス: `Materials.cross_sections` は、相対パスも Python の作業ディレクトリから絶対パスに直して書くので、
  worktree のパスが入る（ADR-0130 の「確認」で outs の model.xml に絶対パスが入ると書いた）
- スレッド数: `openmc.run` の引数で、model.xml には入らない

ライブラリの段階（`openmc_library.py`）は、入力を読む関数を `openmc_common.py` から取っていた。`openmc_common.py` は
モデルの部品（物質、tally）も持つので、モデルの部品を直すと、核種ごとの 63 の NJOY の段階が再実行する（逐次で約 30 分）。

## 検討した案

- 案 A: 分けない。分岐点からの復元を続ける。復元は手作業で、`dvc.lock` が HEAD と合わない状態が常に残る
- 案 B: model.xml にライブラリの相対パスを書く。OpenMC が作業ディレクトリから絶対パスに直す（`openmc.config` の
  `resolve_paths`）ので、書く側では防げない。書いたあとに文字列を置き換える手もあるが、OpenMC の出力の整形に頼る
- 案 C: model.xml に `cross_sections` を書かず、`run` が環境変数 `OPENMC_CROSS_SECTIONS` で渡す。model.xml は
  worktree のパスを含まない
- 案 D: ライブラリの段階も、NJOY の入力を作る段階と NJOY を実行する段階に分ける。NJOY の入力は `openmc.data.njoy` の
  `make_ace` の中で組み立てられ、取り出す API がない。取り出すには `njoy.run` を差し替えるか、入力の組み立てをコピーすることに
  なり、OpenMC の版が変わると食い違いを黙って残す。`openmc_library.py` は、入力の項目を `from_njoy` に渡すだけの薄い
  スクリプトで、直すことも少ないので、案 D は採らない。代わりに、入力を読む部分を別のファイルにして、ライブラリの段階の
  deps から、モデルの部品を外す

## 決定と理由

案 C を採る。

- `openmc_model.py model <入力> <モデルのディレクトリ>` が `model.xml` を書く。物質、形状、設定（粒子数、世代数、乱数のシード）、
  tally を含み、ライブラリのパスもスレッド数も含まない。段階の deps はスクリプト、共通の部品、入力
- `openmc_model.py run <入力> <モデルのディレクトリ> <作業ディレクトリ>` が model.xml を解く。ライブラリは入力の
  `cross-sections` から環境変数で渡し、スレッド数は入力の `threads` から `openmc.run` に渡す。段階の deps は model.xml、
  ライブラリのディレクトリ、入力の params（`cross-sections`、`threads`）。スクリプトは deps にしない。このスクリプトの変更が
  結果を変えるのは、モデルを変えるときだけで、そのときは model.xml が変わる
- 入力を読む関数（`check_keys`、`read_input`、`input_path`）を `examples/openmc_input.py` に移す。ライブラリの段階の
  deps は `openmc_input.py` だけにし、`openmc_common.py` の変更で NJOY が再実行しない

threads を params に置くのは、スレッド数が結果を変えないことを確かめていないため（#0060 は xlc の出力の話で、OpenMC は別）。
変えると再実行するが、結果の取りこぼしより時間の損のほうが軽い。

## 帰結

- スクリプトの直しでモデルが変わらなければ、`openmc-<問題>-model` は再実行しても model.xml の md5 が同じで、OpenMC の
  段階は走らない
- 分岐点からの復元（ADR-0130 の「帰結」）は、この ADR のあとに実行した段階では要らない。rebase で deps が変わっても、
  model.xml の md5 が同じなら OpenMC は走らない
- 他の問題（1C〜2K、MOX など）の段階は、`openmc-1B-model` と `openmc-1B` と同じ形で足す
- 実行の段階は OpenMC の `run`（`--restart` も）と statepoint の出力を変えない。外部の計算機の自己完結な実行
  （ADR-0155、ADR-0168。`openmc_resume.py`、問題 3 と 5-2D のスクリプト）は、model.xml を VM の側で作る別の経路で、この
  ADR では変えない。DVC のステージがない（手で流して extract だけを段階にした）問題も対象外
- 既存の `openmc-1B` の出力は、model.xml を deps に結び付けるだけで再実行しない（`dvc commit`）。ただし model.xml が、
  その出力を作ったときのモデルと違えば結び付けない（#0104 の、model.xml を比べる案の前提でもある）

### 確認

- `PYTHONHASHSEED` を変えて 3 回 model.xml を書き、md5 が一致した。`/home` と `cross_sections` を含む行はない
  （1B、2026-10-14。md5 8671ba78）。入力の `threads` と `cross-sections` を変えても model.xml は同じバイトになる
  （examples/vera/test_openmc_model.py）
- 少ない粒子数（2000 粒子、4 世代）で `model`、`run`、`extract` を通して実行し、statepoint から k と標準のタリーを抜けた
- 既存の `openmc-1B` の出力（2026-10-11）を作った model.xml は、新しい model.xml と、標準のタリー（ADR-0152）と
  タリー出力の設定の有無だけが違う。そのためこの出力は新しい model.xml に結び付けず、実行の段階は古いまま残し、作り直し
  （1B の参照解を標準のタリーで再実行する）は別に行う
- 1B の model.xml と statepoint を作り直したときの k が、ADR-0130 の測定（examples/vera/verification.md の 1B）と
  統計誤差の中で一致するかは、その作り直しで書く
