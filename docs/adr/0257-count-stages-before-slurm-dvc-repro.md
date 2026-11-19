---
type: ADR
status: accepted
date: 2026-11-19
requires: [ADR-0184, ADR-0212, ADR-0223]
supersedes: []
superseded-by:
issues: []
---

# ADR-0257: Slurm に投げる dvc repro は、走りうるステージを数えて上限と参照解の許可を確かめるラッパーを通し、生の dvc repro を hook で止める

## 背景

スキル compute は、`dvc repro` を投げる前に走りうるステージを数えて見積もり、走るステージに参照解（数時間の OpenMC）が
入っていれば投げずに原因を調べる、と文で定めていた（数え方は `scripts/count-changed-stages.py` と、目標の上流のうち changed の
下流にあるもの）。守られなかった例:

- 2026-10-22: 文書だけを直した deps で changed になった 5 時間の OpenMC を目標の上流に入れて投げ、12 分で止めた（Slurm 454。スキル compute の記載）
- 2026-11-19: `examples/vera/endfb70/dvc.yaml` の全体を目標にした `dvc repro` が、`xs-range` と 53 の `xs-nuclide` から走り始め、
  7 分 45 秒で取り消した（Slurm 960。ログ `~/github/slurm-logs/openmc-endfb70-1B-1D-2B-960.out` の `Running stage` と `TIME:` の行）

数える手順は投げる前に手元で行うもので、dvc の呼び出しが遅い（#0268。全体の `count-changed-stages.py` が手元の 15 分に
収まらなかった例がある）ので、省かれやすい。

判断の軸:

- 数えずに投げることが、覚えていなくてもできないこと（止めたときに使うものが分かること）
- 数えるのに投げる側を待たせないこと
- 上流を飛ばす数秒の `dvc repro -s`（/dvc-relock、ADR-0184）と、手元の数秒のステージを止めないこと

## 検討した案

- 案 A（採用）: `scripts/dvc-repro.py` が、ジョブの中で走りうるステージを数えて一覧を出し、呼ぶ側が書いた上限
  （`--max-stages`）を超えるか、許していない参照解のステージが入っていれば、何も実行せずに終了コード 1 で止まる。通れば
  `dvc repro` を実行する。xlc の `.claude/settings.json` の hook が、`slurm-run.sh` に生の `dvc repro` を渡すコマンドを止める
- 案 B: スキルの文で守る（今まで）。上の背景のとおり守られなかった
- 案 C: hook が投げる前に数える。hook は Bash のコマンドの前に 10 秒の timeout で走り、数えるための dvc の呼び出しは
  1 つの yaml で 1 分を超える（下の確認）ので、hook の中では数えられない
- 案 D: ラッパーが先に正規形のステージ（`python-ast`、`toml-values`、`build` と `binary-forms`）を実行してから数える。docstring
  だけの .py の変更などで正規形のバイトが変わらない場合に下流まで数える分（下の決定の「上限」）が減るが、`build`（数分の cargo）が
  `bin/` を書き換えてから止まることになり、「何も実行せずに止まる」が崩れる。正規形の作り直しは /dvc-relock の前半が手元で
  行う
- 案 E: 参照解のステージを、`references-shell.txt` を deps に持つステージ（`dvc dag` で `references-shell` の子）として見分ける。
  `dvc dag` だけで決まるが、拠るのは deps の書き方の取り決め（ADR-0212）で、守らせる検査が無い。cmd が references の devShell で
  包まれていること（NJOY と OpenMC を実行する条件そのもの）で見分ける方が直接である
- 案 F: 核種ごとの xs の `examples/stage_runner.py` にも同じ上限を足す。runner の目標は組（`xs-nuclide`）そのもので、走る核種は
  投げる側が作り直すと決めたものである。組の外の上流は直下の dvc.yaml の正規形、`build`、`xs-range` だけで参照解を含まない
  （`dvc dag --dot examples/vera/dvc.yaml` で `xs-nuclide@U-235` の上流をたどった。2026-11-19）。runner は今のままにし、runner に
  続けて流す下流（`dvc repro <yaml>:xs`）をラッパーにする

## 決定と理由

案 A にする。

- 使い方: `scripts/dvc-repro.py --max-stages <N> [--allow-reference <yaml>:<ステージか foreach の組>]... [--count-only]
  [--no-run-cache] (-P | <dvc repro の目標>...)`（default の devShell で、リポジトリの直下から。目標は cwd からの相対で、ステージの名前は
  直下からのパスなので、ほかの場所では止まる）。`-s` と、ほかの dvc repro のオプションは、走るステージの
  数え方を変えうるので受けない（`-s` は手元の数秒の再実行に使う）
- 範囲: `dvc repro --dry -f <目標>` が並べるステージ（目標とその上流。`-f` で changed によらずすべてを並べる）と、その cmd
- 走りうるステージ: 範囲のうち、`dvc status --json --with-deps <目標>` が changed と示すものと、`dvc dag --dot <目標>`（目標ごと）
  の辺でその下流にあるもの。dvc status と `dvc repro --dry` は上流が出力を書き換えた後に走る下流を示さないので足す。上流が
  同じ出力を書けば下流は飛ぶので、この数は上限である
- references-shell: always_changed でいつも changed だが、書き出す `references-shell.txt` が references の devShell と同じなら
  参照解は飛ぶ。`scripts/check-references-shell.sh`（ADR-0212 の hook と同じ検査）が通るときだけ、その下流を足さない。通らなければ
  下流の参照解を数え、許していなければ止まる
- 参照解のステージ: cmd が `nix develop <flake>#references` で包まれたもの（NJOY と OpenMC のステージ。CLAUDE.md の「計算」）。
  抽出のような短いものも含め、`--allow-reference` に挙げたものだけを実行してよい。挙げたものが範囲に無ければ（書き違い）止まる
- 出力: ステージごとに `RUN:`（changed の種類か、どの changed の下流か。参照解には `[reference]`）、`COUNT:`、dvc の呼び出しごとの
  `TIME:`。止まるときは `STOP:` と `NEXT:`。`--count-only` は数えて終わる（投げる前に手元で数を見るとき）
- hook: `.claude/hooks/block-raw-dvc-repro-in-slurm.py`。Bash のコマンドのうち `slurm-run.sh` の呼び出しの最後の引数（ジョブの
  コマンド）を ADR-0223 と同じ読み方で読み、dvc の repro を直接呼べば終了コード 2 で止め、ラッパーの使い方を示す。`slurm-run.sh` の
  外の `dvc repro`（手元の `-s` と数秒のステージ）は止めない。stdin に `slurm-run` の文字列が無ければ python3 を起動しない
- コマンドの読み方は ADR-0223 の hook から `.claude/hooks/shell_words.py` に移し、2 つの hook が共有する。ジョブのコマンドを
  最後の語で見分けるので、`slurm-run.sh … '<コマンド>' 2>&1` の fd の番号を語から外す

理由:

- 案 B は守られなかった。止めるには、投げるコマンドを機械で判定する hook と、数えて止まるスクリプトが要る
- 数えるのをジョブの中にすれば、#0268 の遅い dvc の呼び出しで投げる側を待たせない（案 C は hook の時間に収まらない）
- 上限を呼ぶ側に書かせると、投げる前の見積もりの本数がコマンドに残り、外れたときに走る前に止まる

### 帰結

- Slurm 454 と 960 の形（目標の上流の参照解、目標の全体のライブラリの作り直し）は、上限か参照解の許可で、走る前に止まる
- 数は上限なので、正規形（python-ast など）が古いまま数えると、正規形のバイトの変わらない変更でも下流を数えて止まる。その
  ときは /dvc-relock の前半のとおり正規形を手元で作り直してから数え直す
- ジョブのコマンドの中の生の `dvc repro`（`a && dvc repro …` を含む）は hook が止める。`slurm-run.sh` を通さずに Slurm に渡す形
  （直接の sbatch は claude-hooks が止める）、シェルにスクリプトとして渡す `bash ~/.claude/scripts/slurm-run.sh …`、ADR-0223 の帰結に
  挙げた読まない書き方は止め損なう
- 数えるのに dvc を 3 種類（dag は目標ごと）呼ぶので、ジョブの始めに 1 分ほどかかる（下の確認の examples/vera/dvc.yaml で 57 秒）

### 確認

- `scripts/test_dvc_repro.py`: 参照解のステージを含む目標、上限を超える目標で止まって dvc repro を呼ばず、通る目標で
  dvc repro を呼ぶこと（PATH に置いた偽の dvc で）。下流の数え、references-shell の扱い、組の許可、`-s` とリポジトリの直下の外を
  受けないこと、ログで `STOP:` が一覧の後に来ること
- `.claude/hooks/test_block_raw_dvc_repro_in_slurm.py`: `slurm-run.sh` に生の dvc repro を渡すコマンドを止め、ラッパー、
  stage_runner、手元の `dvc repro` と `dvc repro -s` を止めないこと。`slurm-run.sh … 2>&1 | tail` のように fd の番号が後に続いても
  止めること。`.claude/hooks/test_block_raw_dvc_checkout.py` は `2>/dev/null dvc checkout` を止めること
- 2026-11-19、xlc 4d6af2c0f（分岐点は main 3eaeb4dcd。出力は `scripts/worktree-restore-outputs.sh` で dvc.lock のとおりに戻した）で
  `scripts/dvc-repro.py --max-stages 100 --count-only examples/vera/dvc.yaml`: 走りうるステージは 83（xs-nuclide 60、keno-assembly 12、
  keno-pincell 3、xs-range、xs、lattice-1B、compare-1B、references-shell、参照解 3）で、どれも dvc status の changed（下流から足した
  ものは無い）。参照解の openmc-1B-model、openmc-1B、openmc-1B-extract を許していないので、何も実行せずに終了コード 1 で止まった。
  changed の deps は xs と lattice の正規形（`binary-forms/xs.md5` が 62 ステージ、`binary-forms/lattice.md5` が 16）、openmc-1B-model と
  openmc-1B-extract は python-ast の 5 つと入力の正規形、openmc-1B は model.xml、ライブラリ、入力と cmd（同じ日に `dvc status --json` で数えた）。dvc の呼び出しは repro --dry 7.2 秒、dag 9.0 秒、status 41.1 秒（手元。references-shell.txt は devShell と同じ）
