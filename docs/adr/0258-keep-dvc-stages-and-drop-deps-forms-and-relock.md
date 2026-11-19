---
type: ADR
status: proposed
date: 2026-11-19
requires: [ADR-0122, ADR-0124, ADR-0126, ADR-0130, ADR-0159, ADR-0160, ADR-0163, ADR-0176, ADR-0193, ADR-0222, ADR-0223, ADR-0257]
supersedes: []
superseded-by:
issues: []
---

# ADR-0258: DVC のステージと dvc.lock は残し、deps の正規形と relock をやめ、上流を飛ばす実行は上流の出力が dvc.lock のとおりかで決める

## 背景

計算のステージは ADR-0123（2026-10-14）で DVC に載せ、ADR-0159 が deps をビルドした実行ファイルにした。DVC は検討の段階で
親（当時のセッション）が提案して入れたもので、ユーザーは当時から入れるのが早すぎないかと見ていた（ユーザー 日付不明）。

DVC は deps の md5 で出力の古さを判定する。deps にはコード（`bin/` の実行ファイル。ADR-0159）と Python のツール
（ADR-0210）が入るので、xs、lattice、core と比較のツールを直すと、それを deps に持つステージがすべて changed になる。
`dvc repro <目標>` は目標の上流の changed を実行する（ADR-0184）ので、直した後の最初の計算が上流の作り直し
（核種ごとのライブラリで約 2 時間。#0150）を背負う。そのため、直すべきものを直さない判断が繰り返し出ている:

- #0149（閉じた）: ADR-0186 のレビューの 2 周目の指摘を、xs を変えるとライブラリ（約 2 時間）が changed になるので直さなかった
- #0150: 同じブランチで、作り直しを避けるため 2 周目の指摘を issue に残した
- #0241: lattice の CMFD の複製を共通のクレートに移すのを、15 の dvc.lock の約 380 行が changed になるので後回しにした
- #0174: p2 の入力に `few-group` を書くかを、入力を変えると changed になるので ADR-0196 で決めなかった
- #0106: 実行中の statepoint を間引く修正を、`openmc_resume.py` の正規形が変わって `openmc-5A-2D-extract` が changed になるので、
  その extract を再実行できるブランチまで待つ
- 2026-11-19 の issue の片付けで、ツールの誤り #0104、#0106、#0236 を、直すと extract と compare のステージが changed になる
  という理由で残し、xs と lattice のコードを触る修正も同じ理由で範囲から外した（親の依頼による。文書に残るのは #0106 の本文の
  条件だけ）

ユーザーは、再実行の引き金になるから問題を放置するのは本末転倒で、その判断が出てくるなら DVC を消すべきだと指摘した
（ユーザー 2026-11-19）。親も、ツールが直すべきものを直さない向きに判断を曲げていると見ている。

結果を変えない変更で changed にしないための仕組みも増え続けた: Python の正規形（ADR-0210）、TOML の値の正規形（ADR-0254）、
実行ファイルの正規形（ADR-0255）、`dvc commit` で付け直す条件（ADR-0212 の (b)）、上流を飛ばす `-s` の条件（ADR-0184）、
rebase のあとの relock（スキル dvc-relock）、生の `dvc checkout` を止める hook（ADR-0223）、Slurm に投げる `dvc repro` の
走るステージを数えるラッパーと hook（ADR-0257）。これらは結果を変えない変更だけを見えなくするもので、結果を変える修正
（誤りを直すこと）の費用は減らさない。さらに DVC の操作そのものが遅い: p2 の 545 核種の作り直し（Slurm 887）の 2 時間 3 分の
うち 16 コアの計算は約 27 分で、`dvc commit` が 1 回約 100 秒で 64 回、全体の `dvc status` は手元の 15 分に収まらない（#0268）。

根本の原因は、出力の古さをコードの md5 で機械が判定し、その判定が実行の順序（上流の changed を先に実行する）と検査の
対象（relock、正規形の食い違い）に結びついていることである。判定が実行の順序か検査に結びついていれば、どの駆動のツールでも、
誤りを直すことと重い出力を作り直すことが結びつく。結びつきを切るのは、判定をやめることか、判定を残して実行の順序と検査から
外すことのどちらかである（検討した案の案 B と案 C2）。

判断の軸（(c) を満たさない形は採らない）:

- (a) 記録の数値がどのコードと入力から出たかを辿れること（ADR-0122。記録は数値を出したコミットと日付を引く）
- (b) 重い出力（ライブラリ、参照解、燃焼の履歴、テーブル）を作り直さずに使い回せること
- (c) 直すべきものを直すことが、手続きの費用で妨げられないこと

## 棚卸し

main d17396374（2026-11-19）で数えた。数え方は「## 補足」にある。

- `dvc.yaml` は 24（リポジトリの直下の 1 と examples の 23）。ステージの定義は 354（うち foreach が 30）、`dvc.lock` の
  ステージは展開して 1,897。lock の cmd で分けると xs 745、正規形 476、参照解（cmd が references の devShell で包まれた NJOY と
  OpenMC）281、lattice 207、core 98、比較 75、ほか 15
- `dvc.lock` の outs は 2,326。git で追跡しないものが 1,768（すべてキャッシュにある）、`cache: false` で git が追跡するものが
  558（正規形の 480 = python-ast 40、toml-values 436、binary-forms 4。`references-shell.txt` 4、比較の JSON 63、ほかの JSON 7、
  `xs-range.toml` 4）。追跡しない outs の大きさの和は 9.5 GB
- `.dvc` のファイル（ADR-0222）は 765（examples/depletion 749、examples/vera 12、examples/mox 2、examples/branches 2）。
  outs は 765（ディレクトリ 8）、大きさの和は 9.5 GB。拡張子は msgpack 408、log 263、out 91、ほか 3
- キャッシュ `references/dvc-cache`（reflink。ADR-0123、ADR-0126）は 36 GB で、うち run-cache は 37 MB
- 比較の JSON（examples の下の追跡する `compare-*.json`）は 101 で、71 が読んだ入力の md5（`inputs`）を持つ。lock の outs は 63 で、
  ほかの 38 はステージの出力でない。検証の記録は 14 本で、表の目印は 222
- ADR: 番号のある 64 本が `dvc` の語を含む（大小を区別しなければ 71 本）。決定が DVC のステージの仕組みそのものなのは
  ADR-0184、ADR-0210、ADR-0223、ADR-0254、ADR-0255、ADR-0257 の 6 本（ADR-0123 は既に superseded）。決定の一部が DVC の
  ステージに拠るのは ADR-0130、ADR-0136、ADR-0159、ADR-0160、ADR-0163、ADR-0176、ADR-0193、ADR-0212、ADR-0222 の 9 本。ほかは、確認や
  帰結に `dvc repro` で流した経緯を書くだけで、凍結した記録として変えない
- script、hook、テスト、設定: `dvc` の語（大小を区別しない）を含む追跡したファイルは、`scripts/` 21、`.claude/` 11（hook 5、
  `settings.json`、スキル 5）、`crates/` 4（コメントと ignore の理由）、`flake.nix`、`.pre-commit-config.yaml`。`examples/` では
  検証の記録の外の 74（`dvc.yaml` 23 を含む yaml 24、py 21、入力の TOML 21、sh 6、`.gitignore` 2）と検証の記録 12。pre-commit
  の hook のうち `dvc-lock`、`python-ast`、`toml-values`、`binary-forms`、`references-shell` の 5 つと、`verification-tables` の
  lock の md5 の照合が dvc.lock に拠る
- 文書: スキル dvc-relock（259 行）の全体と、compute、topic-start、openmc-reference、verification-record の DVC の項。
  CLAUDE.md の 3 行、docs/environment.md の 25 行。開いている issue のうち 41 が `dvc` の語を含む
- crates/docs の検査: `record_tables_match_the_compare_outputs` は `examples/verification_tables.py --check` を呼び、JSON の
  md5 を同じ場所の `dvc.lock` と照らす。`linked_worktrees_share_references` はキャッシュの共有を確かめ、ステージに拠らない
- 履歴: DVC を入れた 2026-10-11 からの main のマージでないコミット 5,792 のうち 373 が件名か本文に relock を含み、428 が
  `dvc.lock` を変える
- 走っているもの（2026-11-19 20 時）: Slurm 961（fix/0131-endfb70-keno-offset）と 962（feat/0293-energy-normalization-check）。
  マージしていないブランチのうち、fix/0131-endfb70-keno-offset、feat/0293-energy-normalization-check、
  fix/0104-0106-0236-tool-defects が `dvc.lock` か `.dvc` を変える。docs/branches.md の残すブランチ 5 本はどれも DVC を入れる
  前のもので、`dvc.yaml` を持たない

## 検討した案

- 案 A: DVC をすべてやめる。重い出力は `references/` などの追跡しない置き場に置き、出どころは記録と、出力に添える小さな
  出どころのファイル（コミット、入力の md5、コマンド）で辿る
  - (a): 出どころのファイルと記録のコミットで辿れる
  - (b): 置き場の 2 つの形で分かれる。A1（同じパスに上書き）は前の版を失い、記録が引く版を戻せない（ADR-0222 が案 C、D を
    退けた理由と同じ）。A2（中身の md5 で置く自前の置き場と、worktree への reflink の戻し）は (b) を満たすが、DVC のキャッシュと
    `dvc add`、`dvc checkout`、`dvc gc` を作り直すことになり、ADR-0123 が案 G で退けた保守を自分で持つ
  - (c): 満たす
  - 要らなくなるもの: 案 B の一覧に加えて、`flake.nix` の dvc のラッパー（#0101）と `.dvc/config`。`.dvc` の 765 ファイル、
    `scripts/worktree-restore-outputs.sh` と `scripts/dvc-checkout.py` は、A2 では自前の置き場の形と戻しに替わる
  - 移行: 案 B の移行に加えて、キャッシュの 36 GB と `.dvc` を新しい置き場の形に移す（A2）か、記録が引く前の版を諦める（A1）
- 案 B: ステージと lock（`dvc.yaml`、`dvc.lock`、正規形、relock）をやめ、`dvc add` と `.dvc` のファイルによる保管
  （キャッシュ、worktree への戻し、前の版が残ること。ADR-0222、ADR-0126）だけを残す。実行の手順は、例ごとの plain なコマンドの
  一覧と、それを実行するテストできる Python の runner に移す
  - (a): 出力の `.dvc` に、それを作ったコミット、コマンド、入力と実行ファイルの md5 を添える（補足の「案 B の設計」の「出どころ」）。記録は
    今までどおりコミットと日付を引く。今の `dvc.lock` の deps は、出力が実際に何から作られたかの記録なので、移すときに `.dvc` に
    そのまま移せる
  - (b): 今の `.dvc` と同じ仕組みで、キャッシュ、reflink、worktree での戻し、前の版の保持はそのまま残る
  - (c): 古さを機械が判定しないので、コードを直しても何も changed にならず、直すことに作り直しの費用が付かない。作り直すのは、
    その出力を読む作業が来たときに、その作業が決める（CLAUDE.md の後回しの規則と同じ）
  - 失うもの: 何が古いかの自動の一覧（`dvc status`）、目標の上流を自動で実行すること、run-cache、`dvc dag`。古い出力を今の
    コードの結果と取り違える危険は、出どころのコミットを記録に書く手順（record_provenance.py）で抑える
  - 新しく作るもの: 例ごとの手順の一覧（`steps.toml` 23）と、それを書くスクリプトの付け替え、runner（指定した手順だけを実行し、
    上流を自動で実行しない。汚れた作業ツリーで止まり、`bin/` を作り、入力の md5 を前後で比べ、worktree ごとの排他を持つ）、
    `.dvc` の `meta` と出どころのファイルの書き手、`record_provenance.py`、`verification_tables.py`、JSON の入力の照合の付け替え、
    lock の 1,768 の outs を `.dvc` に移すスクリプト。runner の働きのうち、指定した手順だけを上流なしに実行すること、入力の md5 を
    出どころに書くこと、同じ worktree の 2 つ目を止めることは、`dvc repro -s`、dvc.lock の deps、dvc のリポジトリのロックが今も持つ
  - 要らなくなるもの: 補足の「案 B の設計」の「要らなくなるもの」
  - 移行: 補足の「案 B の設計」の「移行の順と見積もり」（4 段）
- 案 C: DVC のステージを残し、changed を「再実行の義務」として扱う規則だけをやめる。ADR-0130 は既に「changed は結果が要るとき
  に再実行する」と決めているので、実際に変えるのは、`-s` を ADR-0184 の 2 つの例外に限らず使ってよいこと、rebase のあとの
  relock をやめること、changed が残っていても検査で止めないことになる
  - (a): 満たす（dvc.lock がつながる限り）
  - (b): 満たす
  - (c): 満たさない。誤りを直したコミットは、正規形を書き直して lock をコミットする手続き（直下の `dvc.lock` と `python-ast/`
    などの書き直し）を毎回伴い、rebase で lock が衝突する（どの hook がどの修正を止めるかは案 C2 の (c)。この手続きは正規形を
    残すことから来る）。changed が積もる（examples/vera の yaml だけで 83。ADR-0257 の確認）と
    `dvc status` の一覧は意味を失い、それでも DVC の遅さ（#0268）と数えるラッパー（ADR-0257）は残る。`-s` を広げると、
    ADR-0184 が防いだ、上流の outs を戻していない作業ツリーで下流を作る誤りが戻る
  - 要らなくなるもの: スキル dvc-relock の前半と ADR-0184 の `-s` の条件の一部だけ。正規形、hook、数える script は残る
- 案 C2: DVC のステージと dvc.lock は残し、changed の判定を手続きと検査から外す。やめるのは 3 つで、(1) 正規形の 3 種
  （python-ast、toml-values、binary-forms。ADR-0210、ADR-0254、ADR-0255）と、それを守る hook、検査、追跡する 480 ファイル、
  (2) changed を直す手続き（スキル dvc-relock の前半の 2〜6、付け直しを ADR の移行に限る決まり。ADR-0130、ADR-0212 の (b)）、
  (3) `-s` を ADR-0184 の 2 つの例外に限る決まり。ステージの deps は元のファイル（`bin/`、.py、TOML）に戻す（コードを deps から
  外す形は案 D で、(a) が弱まる）。直下の `dvc.yaml` は、正規形のステージが無くなると `build` だけになるので、`build` もやめ、
  `bin/` は実行の前に `scripts/build-binaries.sh` で作る（案 B と同じ）
  - (a): 満たす。lock の deps は、出力を作ったファイルそのもの（実行ファイルを含む）の md5 で、正規形より細かい。上流を飛ばして
    下流を作っても、上流の outs が lock のとおりなら、lock の上流と下流の記録はつながったまま残り、上流の古さは上流の lock の項の
    deps から読める（ADR-0184 の決定の、`-s` の 2 つ目の例外の後の段落。そこで下流を数秒のものに限る理由は、上流を作り直した
    ときに下流も数秒で作り直せることで、つながりの条件ではない）。`record_provenance.py` と `verification_tables.py` は、今のまま
    dvc.lock の履歴と outs を読む
  - (b): 満たす（キャッシュと `.dvc` は今のまま）
  - (c): 満たす。今の main で誤りを 1 つ直してマージするまでの手続きは、直す対象で分かれる（main d17396374 のこの worktree で
    2026-11-19 に変更を当てて hook の検査を実行し、元に戻して確かめた）:
    - .py: ステージの deps の .py の中身を変えると（examples/vera/compare_k.py の `return values[0]` を `return float(values[0])`）、
      pre-commit の hook `python-ast`（`test_normalize_python.py RepositoryDepsTest`）が
      `test_forms_of_python_files_changed_since_python_ast_ran_are_current` で止める。通すには `dvc repro -s dvc.yaml:python-ast` を
      実行し、直下の `dvc.lock` と `python-ast/` をコミットする
    - TOML: 入力の値を変えると（examples/vera/p1/pincell-1C.toml の `azimuthal = 64` を 32）、hook `toml-values` が
      `test_every_form_matches_its_toml_file` で止める。`toml-values` を実行して同じくコミットする
    - Rust: hook は止めない（`binary-forms` の検査は、直下の `dvc.lock` の `build` と `binary-forms` の項と追跡した正規形を
      比べるだけで、`bin/` を見ない）。手続きは /dvc-relock の前半の 2（`build` と `binary-forms` を実行し、直下の `dvc.lock` と
      `binary-forms/` をコミットする）で、2026-10-28 からの main のマージ 184 のうち `crates/` を変えた 82 の 66 が直下の
      `dvc.lock` も変えた
    - どれも、rebase のあとに /dvc-relock の前半（出力を戻す、正規形を書き直す、changed を数える、再実行する数秒のステージを
      決める、lock をコミットする、lock が読めるかを確かめる）が付く

    案 C2 では正規形と relock が無くなるので、これらの手続きは無くなる。残る検査のうち dvc.lock を読むのは、`dvc-lock`（dvc が
    lock を読めるか）、`verification-tables`（表と JSON、JSON と lock の outs の md5）、JSON の `inputs` と lock の deps の照合
    （#0295）で、どれも追跡した JSON と lock を比べる。修正は JSON も lock も変えないので、これらは止めない。`references-shell` は
    flake.nix、flake.lock、docs/patches/ の変更でだけ `references-shell.txt` の書き直しを求め、コードの修正では動かない
  - 失うもの: 結果を変えない変更（注釈、docstring、TOML の書き方、実行ファイルに入るソースの行の番号）で changed にならないこと。
    ただし今も、例の lock のステージ 1,416 のうち `binary-forms/` を deps に持つ 1,048 のすべてが、今の正規形と違う md5 を持つ
    （main d17396374）。正規形を入れた移行のコミット 37e3044cb の本文でも、静的に数えた changed のステージは移行の前後とも 1,520 で、
    実行ファイルだけで changed の 1,046 は 1 つも changed でなくならなかった。changed の一覧は、今も作業の一覧として読めない
  - 要らなくなるもの: 正規形の 480 ファイルと直下の `dvc.yaml`、`dvc.lock`（`build` と正規形の 481 の項）。pre-commit の
    `python-ast`、`toml-values`、`binary-forms`。`scripts/normalize_binary.py`、`normalize_python.py`、`normalize_toml.py` と、
    それぞれの `test_*.py`（script とテストの 6 本で 1,342 行。`test_normalize_toml.py` の `JsonInputsTest` は別のテストに移す）と、移行を確かめた後の
    `compare-dvc-lock.py` とそのテスト。スキル dvc-relock の前半の 2〜6（前半の 1 の出力の戻し、「マージが止まったら」、後半の
    出力の戻しは残る）。ADR は 3 本を置き換え（ADR-0210、ADR-0254、ADR-0255）、7 本の一部を変える（ADR-0130 の rebase のあとの
    relock、ADR-0159、ADR-0160、ADR-0163 の `bin/` を `build` の出力にする部分、ADR-0184 の `-s` の条件、ADR-0212 の (b)、
    ADR-0257 のラッパーと hook の範囲）
  - 残るもの: 例の `dvc.yaml` と `dvc.lock` の 23 ずつ、pre-commit の `dvc-lock`、`references-shell`、`verification-tables`、
    `scripts/dvc-repro.py`、`check-skip-upstream.py`、`count-changed-stages.py`、`dvc-checkout.py`、`check-dvc-lock.sh`、
    `check-references-shell.sh`、`.claude/hooks` の 2 つの hook、`examples/stage_runner.py`（`dvc status` と `dvc commit` を含む）、
    DVC の遅さ（#0268）
  - 新しく作るもの: `scripts/dvc-repro.py` の `-s` の経路と、hook の範囲を手元の `dvc repro` に広げること（決定と理由の「実行」）。
    lock の deps を元のファイルに付け替える 1 回のスクリプト（決定と理由の「移行の順と見積もり」）
  - 移行: 決定と理由の「移行の順と見積もり」（段 0 の測りと 3 段）
- 案 D: DVC のステージを残し、deps からコード（`bin/` の正規形、python-ast）を外して、入力の TOML と上流の出力だけにする。
  コードの版は出どころに書く
  - (a): 弱まる。lock は最新と示しながら、出力は前のコードのものでありうる（ADR-0159 の案 C と #0150 の 1 つ目の案が残した、
    人が動くまで DVC が古さを知らない形）
  - (b): 満たす
  - (c): コードの修正では満たすが、入力の修正（#0291 の語の修正、#0174 の入力の追加）では満たさない。TOML の正規形
    （ADR-0254）、relock、DVC の遅さ（#0268）は残る
- 案 E: ほかの駆動のツール（Snakemake、make、Nix の derivation）に替える。ADR-0123 の案 C、E、F で退けた理由に加え、
  コードを deps に持てば (c) は同じく満たさず、持たなければ案 D と同じになる

案 B と案 C2 の比較（数は上の各案の項と棚卸しから）:

|                     | 案 B                                                                                                                                                                         | 案 C2                                                                                                                                                                                          |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 失うもの            | `dvc status` の古さの一覧、目標の上流を自動で実行すること、run-cache、`dvc dag`                                                                                              | 結果を変えない変更で changed にならないこと（今も Rust の 1,048 のステージのすべてが changed で、効いていない）                                                                                |
| 要らなくなるもの    | `dvc.yaml` と `dvc.lock` の 24 ずつ、正規形 480、script 9 とテスト 8、Claude Code の hook 2、pre-commit の hook 5、スキル dvc-relock の全体、ADR の置き換え 6 と一部の変更 9 | 直下の `dvc.yaml` と `dvc.lock`、正規形 480、script 3 とテスト 3（と移行の後の `compare-dvc-lock.py`）、pre-commit の hook 3、スキル dvc-relock の前半の 2〜6、ADR の置き換え 3 と一部の変更 7 |
| 新しく作るもの      | 手順の一覧 23 とそれを書くスクリプトの付け替え、runner、出どころの書き手、記録と JSON の照合のツールの付け替え、1,768 の outs の移し替え                                     | `dvc-repro.py` の `-s` の経路、hook の範囲の拡張、lock の deps の付け替え（1 回）                                                                                                              |
| DVC の遅さ（#0268） | ステージの操作は無くなる（代わりの `dvc add` の時間は測っていない）                                                                                                          | 残る（段 0 で測る）                                                                                                                                                                            |
| 移行                | 4 段                                                                                                                                                                         | 段 0 の測りと 3 段                                                                                                                                                                             |

## 決定と理由

案 C2 を勧める（proposed。決定はユーザーが読んでから）。ただし、DVC の固定費（#0268）を DVC のまま計算の時間より短くできるかを、
移行の前に測る（下の「移行の順と見積もり」の段 0）。できなければ、案 B と比べ直す。

- (c) を妨げていたのは、changed の判定そのものでなく、判定に付いた手続き（正規形の書き直しと relock）と、目標の上流を先に実行する
  順序である。案 C を (c) で退けた理由（誤りを直したコミットが、正規形を書き直して lock をコミットする手続きを伴う）は、正規形を
  残す前提に拠る。案 C2 の (c) に数えたとおり、正規形と relock をやめると、修正のコミットとマージに付く DVC の手続きは無くなる。
  上流を先に実行する順序は、`-s` の経路で外す
- 案 B と案 C2 は、どちらも (a)、(b)、(c) を満たす。案 B が新しく作る runner と出どころの書き手の働きのうち、指定した手順だけを
  上流なしに実行すること、入力の md5 を出どころに書くこと、同じ worktree の 2 つ目を止めることは、`dvc repro -s`、dvc.lock の deps、
  dvc のリポジトリのロックが今も持つ。事実上の標準のツールと手書きのうち保守の小さい方を選ぶ（グローバルの CLAUDE.md の設計
  指針）と、作り直さない案 C2 になる。案 C2 が新しく作るのは `-s` の経路と 1 回の移行のスクリプトで、1,768 の outs の移し替えも、
  `record_provenance.py` と `verification_tables.py` の付け替えも無い。移行は 3 段で、案 B の 4 段より短い（上の表）
- 案 C2 が案 B に劣る点は 2 つある。1 つは DVC の遅さ（#0268）が残ること。Slurm 887（p2 の 545 核種、2026-11-15）の 2 時間 3 分の
  うち、16 コアの計算は約 27 分で、`dvc commit` が約 1 時間 40 分、上流の確認が 547 秒だった。これは判断の軸の実行速度に当たり、
  修正の手続きでなく作り直しのジョブの時間に効く。#0268 は直す候補（commit を最後の 1 回にする、投入と回収から切り離す、vars の
  TOML の表の見出しを減らす）を挙げているが、確かめていない。そのため段 0 で測る。もう 1 つは changed が見え続けること。修正の
  たびに `dvc status` の changed が増え、それを理由に直さない判断が戻る危険がある。今も 1,520 のステージが changed で（37e3044cb）、
  正規形はそれを減らしていないので、案 C2 で見え方が悪くなるわけではない。changed は再実行の義務でないこと（ADR-0130）と、changed を
  修正の後回しの理由にしないことを、スキル compute の DVC の項に書く。changed を見る検査と hook は足さない
- rebase のときの衝突は、案 C2 と案 B で大きく違わない見込みである。2026-10-28 からの main のマージ 184 のうち、例の `dvc.lock` を
  変えたのは 28（examples/vera は 5）、`.dvc` を変えたのは 47。DVC を入れてからの、`dvc.lock` を変えたマージでないコミット 428 の
  うち 220 は直下の `dvc.lock` だけを変えたもので、案 C2 も案 B も直下の lock を無くす。relock のコミットの本文が衝突に触れたのは
  2 つ（4b7dbed36 は直下の lock、4c770bc5a）。同じ例の lock の別々のステージを 2 本のブランチが実行したときに git が混ぜられるか
  （ステージの項が離れていれば混ぜる）は確かめていない。案 B では出力ごとの `.dvc` に分かれるので、同じ出力を作り直したときだけ
  衝突する

### deps

- 例のステージの deps は元のファイルにする。実行ファイルは `../../bin/<名前>`、Python は .py、入力の TOML はそのファイル
  （ADR-0159 の deps の形に戻し、ADR-0210、ADR-0254、ADR-0255 を置き換える）
- コードを deps から外す形（案 D）は採らない。lock が最新と示しながら、出力が前のコードのものでありうる。案 C2 では changed が
  手続きと検査に結びつかないので、コードを deps に持つ費用は changed が増えることだけである
- 直下の `dvc.yaml` の `build` をやめ、`bin/` は実行の前に `scripts/build-binaries.sh` で HEAD のソースから作る。使った実行
  ファイルの md5 は、ステージの deps として lock に残る。ADR-0159 の案 P1 と、案 T1 のうち ADR-0163 が置き換えなかった部分、
  ADR-0160 の strip、ADR-0163 の 1 回の cargo は残し、同じソースから worktree によらず同じ実行ファイルを作る（lock の md5 を
  worktree の間で比べられるようにする）

### 実行

- `scripts/dvc-repro.py` に `-s` の経路を足す。`-s` のときは、(1) 追跡しているファイルに変更のある作業ツリーで止まり（lock を
  コミットするコミットで出し直せることを保つ。ADR-0122。dvc が書く `dvc.lock` の差と、出力のパスを足す `.gitignore` の行の差は
  除く。案 B の runner と同じ）、(2) `scripts/build-binaries.sh` を実行し、(3)
  `scripts/check-skip-upstream.py` で目標の上流が `changed outs` を持たないことを確かめ、(4) `dvc repro -s <目標>` を実行する。
  `-s` の無いときも (1) と (2) を行う
- `check-skip-upstream.py` は、上流の outs が lock のとおりかで判定する（消えた出力は `deleted`、中身が違えば `modified` の
  `changed outs` に出る）。正規形の項（`RUN_FIRST`）は消す。ADR-0184 が防いだ、上流の outs を戻していない作業ツリーで下流を作る
  誤りは、この判定で止まる（ADR-0184 で `-s` の条件にしていた判定と同じ）。`changed command` と `changed checksum` でも止めるかは、
  その種類が lock のつながりを壊すかを段 2 で確かめて決める
- 目標が数秒で終わるという条件はやめる。重いステージも、上流の outs が lock のとおりなら、`-s` で 1 本だけ実行してよい。記録には、
  ADR-0184 のとおり、上流の出どころ（上流を実行したコミットと分岐点）を下流と別に書く
- hook `block-raw-dvc-repro-in-slurm.py` の範囲を、Slurm の外の `dvc repro` にも広げる（手元の `-s` も (1)〜(3) を通すため）。
  `--max-stages` と参照解の許可（ADR-0257）は残す
- 核種ごとの xs（`examples/stage_runner.py`、ADR-0193）は今のまま、changed のものを実行する。xs のコードを直すとすべてが changed に
  なり、作り直すかは、ライブラリを読む作業が決める

### 移行の順と見積もり

段ごとに 1 本のブランチにし、各段の後で main は run-checks.sh と `cargo test` を通る。重いステージの再実行は、どの段にも無い。
時間は担当の作業の目安で、測っていない。

0. DVC の固定費を測る。小さい入力（核種の数を変えた xs-nuclide の組）で、`stage_runner.py` の `dvc commit` を最後の 1 回にした
   形と今の形の、上流の確認、commit、投入の前の処理の時間を、核種の数を変えて測る（時間を記録する測定なので Slurm の timed）。
   p2 の 545 核種で固定費が計算（約 27 分）を下回る見込みが立てば段 1 に進み、直す形は #0268 で実装する。立たなければ、案 B の
   段 1 の `dvc add` の時間の測りと並べて比べ直す。目安は担当 1 本
1. 正規形をやめ、deps を元のファイルに戻す。例の `dvc.yaml` の deps を元のファイルに書き換え、lock の 1,405 のステージ（正規形を
   deps に持つもの）の deps を付け替える。付け替える md5 は、正規形の元のファイルの md5 を git の履歴から引く: 実行ファイルは
   37e3044cb の親の例の lock の `bin/` の md5（4 つの版。そのコミットの本文）、TOML は e327c56ca の親の lock、.py は直下の
   `dvc.lock` の履歴の `python-ast@<パス>` の項（deps の元のファイルと outs の正規形の md5 の組）のうち、ステージの lock の項を
   入れたコミットのもの。outs と cmd は変えず、`compare-dvc-lock.py` で変わったのが deps だけと確かめる。比較の JSON の
   `inputs` の TOML の md5 は正規形の md5（ADR-0254）なので、JSON を書き換えるか、照合で書いたコミットの木の正規形を読むかを
   この段で決める（案 B の段 2 と同じ問題）。比較のツール 3 本（compare_branch_cross_terms.py、compare_grid_interpolation.py、
   compare_p6.py）は元のファイルの md5 を書く形にする。直下の `dvc.yaml` の正規形のステージ、正規形の 480 ファイル、pre-commit の
   hook 3、normalize の script とテストを消し、`check-skip-upstream.py` の正規形の項を消す。目安は担当 1〜2 本、レビュー 2〜3 周
2. `build` をやめ、`-s` の経路を足す。直下の `dvc.yaml` と `dvc.lock` を消し、`scripts/worktree-restore-outputs.sh` は `bin/` を
   `scripts/build-binaries.sh` で作る。`dvc-repro.py` の `-s` の経路と hook の範囲（決定と理由の「実行」）を TDD で足す。スキル
   dvc-relock の前半の 2〜6 を消し、compute、topic-start、openmc-reference、verification-record の該当の項と、CLAUDE.md の DVC の
   行、docs/environment.md を直す。目安は担当 1 本、レビュー 2 周
3. ADR と issue を今の状態に合わせる。ADR-0210、ADR-0254、ADR-0255 を、この ADR で置き換える（superseded）。一部を変える 7 本は
   「## 変更」で直す。#0150 など、正規形と relock があるために開いている issue を閉じるか書き直す。目安は担当 1 本

段 1 はすべての例の lock を書き換えるので、`dvc.lock` を変えるマージしていないブランチ（今は 3 本）と衝突する。段 1 は、それらが
マージされるか、段 1 の後に rebase して lock の差を付け替えるかを決めてから始める。走っているジョブ（961、962）の worktree は
変えない。

### 帰結

- 誤りの修正は、出力の作り直しも DVC の手続きも伴わずにマージできる。作り直しは、その出力を読む作業か記録を書き直す作業が、
  時間を見積もって決め、上流は `-s` で飛ばせる
- changed は、出力が今のコードと入力より前のものだという情報として残り、`dvc status` で読める。検査と hook はそれを見ない
  （見ると (c) を妨げる形に戻る）
- 結果を変えない変更かを確かめたいときは、機械の判定でなく、小さな入力で新旧の実行ファイルの出力を比べる。
  feat/0150-relock-unchanged-results の worktree に書きかけの `examples/relock-check/xs.toml`（コミットなし。核種の処理の種類ごとの
  7 核種を 600 K で処理する組）は、この比べに使える。その注釈が引く「ADR-0258」は、この ADR とは別の、書かれていない決定を指す
- DVC の遅さ（#0268）は案 C2 では残り、段 0 の結果で直す形を決める
- 正規形と relock に拠る開いた issue（#0150 など）は段 3 で閉じるか書き直す。この ADR は issue を起こさず、`issues` は決めてから結ぶ

### 確認

段ごとに、次を通してからマージする。

- 段 0: 核種の数ごとの固定費の値（commit、上流の確認、投入の前の処理）と、545 核種への見込み
- 段 1: `compare-dvc-lock.py` で、移す前後の lock の違いが deps（と、段 1 で決めた JSON の扱いの分）だけであること。changed の
  ステージを移す前後で静的に数え、増えたものが、lock の後に注釈などだけを変えた .py と TOML を deps に持つステージであること。
  .py、TOML、Rust の 1 行の変更で、pre-commit の hook がどれも止めないこと（案 C2 の (c) と同じ手順）。`verification_tables.py --check`、
  `cargo test`、`scripts/run-python-tests.sh`、新しい worktree での `scripts/worktree-restore-outputs.sh`（MISSING なし）
- 段 2: `dvc-repro.py` の `-s` の経路のテスト（汚れた作業ツリーで止まる、`bin/` を先に作る、上流の `changed outs` で止まる、
  上流が `changed deps` だけなら実行する）、hook のテスト（Slurm の外の生の `dvc repro` も止める）。消した名前（正規形、relock、
  `dvc.yaml:build`）を `git grep` して、残る箇所が凍結した ADR と検証の記録の経緯、閉じた issue だけであること
- 段 3: `design-records check`、docs の整合性の検査（ADR と issue の両方向）

## 補足

`requires` には、案 C2 が残して拠る決定（記録の出どころ、references の固定、キャッシュの共有、参照解のステージ、実行ファイルの
ビルドの決まり、表の生成と lock の照合、核種の組の同時実行、`dvc add` の保管、生の `dvc checkout` を止める hook、走るステージを
数えるラッパー）を挙げた。一部を変える ADR と置き換える ADR は、決めたときに `supersedes` と「## 変更」で
結ぶ（proposed のうちに `supersedes` を書くと、置き換えられる側の `superseded-by` が要る）。

棚卸しの数え方（main d17396374 の worktree の直下。2026-11-19）:

- `dvc.yaml`: `git ls-files '*dvc.yaml'`。ステージの定義と foreach は各 `dvc.yaml` の `stages` を YAML で読み、lock のステージは
  各 `dvc.lock` の `stages` の数。cmd の分け方は、`#references` を含むものを参照解、`bin/xs` などを含むものをその実行ファイル、
  `normalize` を含むものを正規形とした
- lock の outs: 各 `dvc.lock` の outs を、パスが `git ls-files` にあるかと、`references/dvc-cache/files/md5/<md5 の先頭 2 文字>/<残り>`
  があるかで分けた
- `.dvc`: `git ls-files '*.dvc'`（`find . -name '*.dvc'` の 766 は `.dvc` のディレクトリを含む）。大きさは各ファイルの `size` の和
- キャッシュ: `du -sh references/dvc-cache references/dvc-cache/runs`
- 比較の JSON: `git ls-files 'examples/**/compare-*.json'` と、そのうち `"inputs"` を含むもの
- ADR: `grep -l dvc docs/adr/[0-9]*.md`（`-i` で 71）。決定の分け方は各 ADR の「決定と理由」を読んで決めた
- ファイル: `git ls-files 'scripts/*' '.claude/*' '*.rs' flake.nix .pre-commit-config.yaml | xargs grep -li dvc` と、
  `git ls-files examples | xargs grep -li dvc` を拡張子と検証の記録で分けたもの
- 開いた issue: `grep -lis dvc docs/issues/[0-9]*.md` のうち `status: open` のもの
- 履歴: `git log --no-merges -i --grep=relock main` と `git log --no-merges main -- '*dvc.lock'` の行の数
- ブランチ: マージしていないブランチごとに `git diff --name-only $(git merge-base main <ブランチ>) <ブランチ>` が `dvc.lock`、
  `dvc.yaml`、`.dvc` を含むか
- 修正の手続き（案 C2 の (c)）: この worktree（main d17396374 と同じ木）で、examples/vera/compare_k.py の 22 行目と
  examples/vera/p1/pincell-1C.toml の `azimuthal` を書き換え、`python3 scripts/test_normalize_python.py RepositoryDepsTest` と
  `python3 scripts/test_normalize_toml.py RepositoryDepsTest JsonInputsTest` を実行して失敗を見てから、`git checkout --` で戻した。
  変更の前はどちらも、`python3 scripts/test_normalize_binary.py RepositoryDepsTest` とともに通る。Rust は
  `test_normalize_binary.py` の `RepositoryDepsTest` を読んで、`bin/` を読む検査が無いことを確かめた
- マージの数: 2026-10-28 からの `git log --merges --first-parent main` のマージごとに `git diff --name-only <マージ>^1 <マージ>` が
  `crates/`、直下の `dvc.lock`、`examples/**/dvc.lock`、`*.dvc` を含むかを数えた
- lock を変えたコミット: `git log --no-merges main -- '*dvc.lock'`（428）のうち、`git diff-tree` の lock が直下の `dvc.lock`
  だけのもの（220）と、`examples/` の lock を含むもの（208）。relock の本文の衝突は、`git log --no-merges -i --grep=relock` の
  本文に「衝突」か conflict を含むもの
- 正規形と changed: 例の各 `dvc.lock` のステージのうち deps に `python-ast/`、`toml-values/`、`binary-forms/` を持つもの（1,405。
  `binary-forms/` 1,048、`python-ast/` 384、`toml-values/` 355）と、その deps の md5 が今の正規形のファイルの md5 と違うもの
  （`binary-forms/` 1,048、`toml-values/` 3）を YAML で読んで数えた。changed の 1,520 は 37e3044cb の本文の値で、数え直していない
- 正規形の script とテストの行の数: `wc -l scripts/normalize_*.py scripts/test_normalize_*.py`

### 案 B の設計（勧めない案）

前の版（fb2b5b72f〜8bbbe47f0、2026-11-19）が勧めた案 B の設計で、レビューを 3 周受けた。段 0 で案 B と比べ直すときに使う。

#### 実行

- 例ごとに、実行の手順を `examples/<問題>/steps.toml`（名前、作業ディレクトリからのコマンド、入力、出力）に書く。今の
  `dvc.yaml` の foreach と vars は、展開した行にする（examples/vera/p5/p5_inputs.py など `dvc.yaml` を書くスクリプトは、
  この一覧を書く）。tomllib で読むので、#0268 の tomlkit の遅さは無い
- runner（`examples/run_steps.py <steps.toml> <手順か組>...`）は、指定した手順だけを、指定した順に実行する。上流を自動で
  実行しない。古さも判定しない（入力の一覧は出どころの md5 を取るためだけに使う）
- 追跡しているファイルに変更がある作業ツリーでは止まる（出どころのコミットで出し直せることを保つため。ADR-0122）。runner 自身が
  書く追跡したファイル（一覧の手順の出力そのもの、その `.dvc` と出どころのファイル、`stage-cpu-seconds.json`）と、`.gitignore` の
  差のうち一覧の出力のパスを足す行は除く（`.gitignore` はファイルでなく行の差で見る。ほかの行の差があれば止まる）。前の手順の出力をコミットせずに次の手順の入力にしたときは、その入力の md5 と、それを書いた手順の出どころで辿る
- 始めに `scripts/build-binaries.sh` を実行し、`bin/` を HEAD のソースから作る。`bin/` は追跡しないので、作業ツリーの検査では
  古い実行ファイルを見分けられない（今は `dvc repro` が `build` のステージで作り直している）。実行の間は worktree ごとの排他の
  ロックを持ち、同じ worktree で 2 つめの runner が始まれば止まる（走っている runner の `bin/` を作り直さないため。今は dvc の
  リポジトリのロックが同じことをしている）。実行ファイルも前後の md5 の照合の入力に入れる
- 実行の前後で入力の md5 を取り、食い違えば止まる（ADR-0193 の案 D1 と同じ）。成功したら、git で追跡しない出力は `dvc add` して
  `.dvc` の `meta` に出どころを書き、git で追跡する出力（比較の JSON など、今の `cache: false`）は隣の出どころのファイル
  （`<出力>.provenance.json`）に同じ中身を書く。`dvc add` は git が追跡するファイルを受けない
- 核種ごとの xs の組は、ADR-0193 の同時実行（前回の CPU 時間の長い順、予算、SIGTERM の扱い、CPU 時間の記録）を runner に移す。
  `dvc status` と `dvc commit` の部分は要らなくなる
- `bin/` は `scripts/build-binaries.sh` で作るビルドの生成物のまま残し、DVC の出力にしない。同じソースから worktree によらず同じ
  実行ファイルを作る決まり（ADR-0159 の案 P1 と、案 T1 のうち ADR-0163 が置き換えなかった `target/dvc` と毎回の `cargo clean -p`、
  ADR-0160 の strip、ADR-0163 の 3 つのパッケージだけの 1 回の cargo）は残し、出どころの実行ファイルの md5 を worktree の間で
  比べられるようにする

#### 出どころ

- `.dvc` の `meta`（DVC が読まない欄）に、実行したコミット、日付、コマンド、入力ごとの（パス、md5）、実行ファイルの md5、
  参照解の手順では references の devShell の store のパス（今の `references-shell.txt` の中身。ADR-0212）、Slurm のジョブの番号を書く。
  md5 はファイルそのものの md5 で、正規形を使わない
- 段 2 で lock から移した `meta` は、lock の deps をそのまま持つので、入力の md5 が正規形（`toml-values/`、`binary-forms/`、
  `python-ast/` のパス）のものを含む。`meta` に出どころの種類（runner か lock か）と lock を入れたコミットを書いて分け、lock から
  移した入力は、そのコミットの木の正規形のファイルで読む（段 3 で正規形を消しても git の履歴に残る）
- 記録は今までどおり、数値を出したコミットと日付を書く（ADR-0122）。`examples/record_provenance.py` は `dvc.lock` の履歴の
  代わりに `.dvc` の `meta` と出どころのファイルを読む
- 比較の JSON は git で追跡したまま、読んだ入力の md5（`inputs`）を持つ。`verification_tables.py --check` は表と JSON の一致を
  確かめ続ける。JSON の md5 を `dvc.lock` と照らす部分は、照らす相手が無くなるので外す。JSON の `inputs` のうち `.dvc` で追う
  出力の md5 は、今の出力と一致するかでなく、キャッシュにあって戻せるかを確かめる（記録が引く版を出し直せることの検査で、古さの
  検査にしない）。git で追跡する入力（入力の TOML など）はコミットから戻り、`references/` の入力（AWS の結果など。ADR-0174）は
  ADR-0124 の固定に任せるので、この検査に入れない。今の JSON が TOML の入力に持つ md5 は値の正規形の md5（ADR-0254。
  `scripts/normalize_toml.py` の `input_md5`）で、書いたコミットの木の `toml-values/` で読む。段 3 からは、比較のツールも入力の
  ファイルそのものの md5 を書く
- ステージの出力でない比較の JSON（38）とその入力の出どころは、今と同じく記録の「計算」と JSON の `inputs` で辿り、
  移行で変えない

#### 移行の順と見積もり

段ごとに 1 本のブランチにし、各段の後で main は run-checks.sh と `cargo test` を通る。計算（重いステージの再実行）は
どの段にも無い。時間は担当の作業の目安で、測っていない。

1. runner と出どころの書き手を足す（足すだけで、何も消さない）。TDD で、手順の実行、汚れた作業ツリーで止まること、入力の
   md5 の前後の照合、`.dvc` の `meta`、核種の組の同時実行を確かめる。DVC 3.67.1 が `meta` のある `.dvc` を読み、`dvc add` の
   やり直しで `meta` をどう扱うか（runner は毎回書き直すので、残らなくてよい）と、約 2,500 の `.dvc` での `dvc checkout`
   と `dvc add` の時間を、使い捨ての複製で測る。遅ければ、案 A2 との比べ直しをこの段で決める。目安は担当 1 本、レビュー 2 周
2. 例の `dvc.lock` の outs を `.dvc` に移し、例の `dvc.yaml` と `dvc.lock` を消す。追跡しない outs の 1,768 は、同じパスと
   md5 の `.dvc` にし、`meta` に lock の cmd、deps、その md5 を lock に入れたコミット（record_provenance.py の求め方）を移す。
   md5 が同じなのでキャッシュの中身はそのまま使え、再実行は要らない。changed のステージの出力も、lock の deps が実際に
   作った入力なので、そのまま移す。git で追跡する outs のうち正規形と `references-shell.txt` を除く 74（比較の JSON 63、ほかの
   JSON 7、`xs-range.toml` 4）は、lock の同じ欄を出どころのファイルに移す。lock に hash の無い参照解（未実行）は手順の一覧に
   だけ載る。同じブランチで
   `verification_tables.py`、`record_provenance.py`、#0295 の JSON の入力の照合（`scripts/test_normalize_toml.py` の
   `JsonInputsTest`）を `.dvc` に向け、`dvc.yaml` を書くスクリプト（p3、p5、p6 の inputs、`rods_cases.py`）を steps.toml を
   書く形にし、スキル compute、topic-start、openmc-reference、verification-record の該当の項を直す。目安は担当 1〜2 本、
   レビュー 2〜3 周
3. 直下の `dvc.yaml` と `dvc.lock`（`build`、正規形 3 種）と、正規形を守る仕組みを消す（下の一覧）。`normalize_toml.py` の
   `input_md5` を読む比較のツール（compare_branch_cross_terms.py、compare_grid_interpolation.py、compare_p6.py）は、
   同じ段で入力のファイルそのものの md5 を書く形にする。スキル dvc-relock を消し、
   CLAUDE.md の 3 行と docs/environment.md の DVC の項を直す。段 2 の後に例の `dvc.yaml` が無いと正規形の検査が空で通るか
   落ちるかは確かめていない。落ちるなら段 2 と 3 を 1 本にする。目安は担当 1 本、レビュー 2 周
4. ADR と issue を今の状態に合わせる: 決定がステージの仕組みそのものの 6 本を、案 B を決めた ADR で置き換える（superseded）。一部が
   拠る 9 本は「## 変更」で直す。#0150、#0268 など、ステージがあるために開いている issue を閉じるか書き直す。目安は担当 1 本

段 2 は、`dvc.lock` を変えるマージしていないブランチ（今は 3 本）と衝突する。段 2 は、それらがマージされるか、段 2 の後に
rebase して lock の差を `.dvc` に移すかを決めてから始める。走っているジョブ（961、962）の worktree は変えない。

#### 要らなくなるもの

- ADR（置き換え）: ADR-0184、ADR-0210、ADR-0223、ADR-0254、ADR-0255、ADR-0257
- ADR（一部の変更）: ADR-0159（DVC で駆動し `bin/` を deps にする部分。案 P1 と T1 のビルドの決まりは残す）、ADR-0193（`dvc status`
  と `dvc commit` の部分。同時実行の順、予算、CPU 時間の記録、失敗の扱いは残す）、ADR-0130（参照解をステージで作る部分。
  references の devShell で包むことと入力の TOML の決まりは残す）、ADR-0136（model.xml で分けることは手順の分け方として残る）、
  ADR-0176 の 4 と 6（lock の md5 の照合）、ADR-0212 の (b)（付け直し）、ADR-0222（対象を記録が引く生成物から重い出力のすべてに
  広げる）、ADR-0160 と ADR-0163（`bin/` が DVC の出力でなくなる）
- スキル: dvc-relock の全体と、compute、topic-start、openmc-reference、verification-record の DVC のステージの項
- hook: `.claude/hooks/block-raw-dvc-checkout.py`、`block-raw-dvc-repro-in-slurm.py`、両者が共有する `shell_words.py` と、それぞれの
  テスト（`.claude/settings.json` の 2 つの項）。pre-commit の `dvc-lock`、`python-ast`、`toml-values`、`binary-forms`、`references-shell`
- script とテスト: `scripts/check-dvc-lock.sh`、`check-skip-upstream.py`、`compare-dvc-lock.py`、`count-changed-stages.py`、
  `dvc-repro.py`、`normalize_binary.py`、`normalize_python.py`、`normalize_toml.py`、`check-references-shell.sh`（中身は出どころに
  移る）と、それぞれの `test_*.py`（8 本。`test_normalize_toml.py` の `JsonInputsTest` は段 2 で移す）。`examples/stage_runner.py`
  とそのテストは runner に移す
- 追跡したファイル: `dvc.yaml` 24、`dvc.lock` 24、正規形 480（`python-ast/`、`toml-values/`、`binary-forms/`）、
  `references-shell.txt` 4
- 残るもの: `scripts/worktree-restore-outputs.sh` と `scripts/dvc-checkout.py`（引数なしの `dvc checkout` は `.dvc` だけを戻し、
  未実行の参照解による 255 は無くなるので、MISSING の一覧だけを残す）、`flake.nix` の dvc、`.dvc/config`、キャッシュ、
  `linked_worktrees_share_references`

#### 案 B の帰結

- 誤りの修正は、出力の作り直しを伴わずにマージできる。作り直しは、その出力を読む作業か記録を書き直す作業が、時間を見積もって
  決める。同じ出力を前のコードの結果のまま読むときは、記録に出どころのコミット（`meta`）を書く
- 結果を変えない変更かを確かめたいときは、機械の判定でなく、小さな入力で新旧の実行ファイルの出力を比べる。
  feat/0150-relock-unchanged-results の worktree に書きかけの `examples/relock-check/xs.toml`（コミットなし。核種の処理の種類ごとの
  7 核種を 600 K で処理する組）は、DVC と無関係にこの比べに使える。その注釈が引く「ADR-0258」は、この ADR とは別の、書かれて
  いない決定を指す
- 古い出力の自動の一覧は無くなる。必要になれば、`meta` の入力の md5 と今の入力を比べて示す読み出しだけのツールを足せるが、
  検査と hook には入れない（入れると (c) を妨げる形に戻る）
- ステージに拠る開いた issue（#0150、#0268、#0211、#0269 など）は段 4 で閉じるか書き直す。この ADR は issue を起こさず、
  `issues` は決めてから結ぶ

#### 案 B の確認

段ごとに、次を通してからマージする。

- 段 1: runner のテスト（手順の実行、汚れた作業ツリーで止まり runner 自身が書くファイルでは止まらない、始めに `bin/` を作る、
  同じ worktree の 2 つめの runner が止まる、入力の md5 の食い違いで止まる、`meta` と出どころのファイルの中身、git で追跡する出力を `dvc add` しない、組の同時実行と
  SIGTERM）。DVC 3.67.1 が `meta` のある `.dvc` を `dvc checkout` で戻すこと。約 2,500 の `.dvc` での `dvc checkout` と `dvc add` の時間
- 段 2: 移す前の 2,326 の outs のうち追跡しない 1,768 のすべてが、同じパスと md5 の `.dvc` にあること（移す前の lock と比べる
  スクリプト）。新しい worktree で `scripts/worktree-restore-outputs.sh` が MISSING なしで終わること。記録が引くステージの
  `record_provenance.py` の出力（コミットと md5）が、移す前と同じこと。`verification_tables.py --check`、`cargo test`、
  `scripts/run-python-tests.sh`。キャッシュの大きさが増えないこと（中身は同じ md5）
- 段 3: run-checks.sh、`cargo test`（`linked_worktrees_share_references` を含む）、pre-commit の全部の hook。消した script、hook、
  正規形の名前と `dvc repro`、`dvc.lock` を `git grep` して、残る箇所が保管（`.dvc`、キャッシュ、戻し）、凍結した ADR と検証の記録の
  経緯、閉じた issue だけであること（入力の TOML と例のスクリプトの注釈は、この段で直す）
- 段 4: `design-records check`、docs の整合性の検査（ADR と issue の両方向）
