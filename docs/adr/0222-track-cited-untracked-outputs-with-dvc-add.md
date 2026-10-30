---
status: accepted
date: 2026-10-31
requires: [ADR-0126, ADR-0130]
supersedes: []
superseded-by:
issues: ["#0207", "#0209"]
---

# ADR-0222: 記録と ADR が引く追跡しない生成物は dvc add で DVC に追う

## 背景

検証の記録（examples/*/verification*.md）と ADR は、DVC のステージの出力でない生成物（燃焼の履歴、分岐のテーブル、core の
履歴、OpenMC の参照解のピン出力など。`.gitignore` が無視する `*.msgpack`）を md5 か名前で引く。生成物は作ったブランチの
worktree にだけあり、/git-merge の後片付け（`git worktree remove`）が無視されたファイルごと消す。DVC の出力はキャッシュ
（`references/dvc-cache`。main と共有。ADR-0126）から `scripts/worktree-restore-outputs.sh`（`dvc checkout`。ADR-0130）で戻るが、
それ以外に戻す経路は `examples/vera/p2/openmc-2K.msgpack` を main からコピーする 1 件しかなかった。

経緯は #0207 に 3 件ある。micro のテーブル（`table-nominal.msgpack` md5 fde256d9 など）が後片付けで消えて main にも無い、
置き換えた組成の履歴を作り直し、テーブルを Slurm 618（20 分 9 秒）で解き直した、`histories` の形が変わって親が main の
`examples/depletion/history/` に新しい版を `cp` で上書きし、記録が引く前の形の `table-residual-nominal-histories.msgpack`
（md5 3d73faee）が残らなかった。手順の文で補う案は /retro で試して捨てた（#0207 の「経緯」）。

判断の軸:

- 手順を覚えていなくても、記録が引く生成物がマージの後に残るか戻せること
- 形が変わって置き換えても、前の版が残ること
- 生成物を作り直す計算を増やさないこと（CLAUDE.md の「計算」）

## 検討した案

- 案 A（採用）: 引く生成物を、置き場のまま `dvc add` で追う。ファイルごとに `.dvc` ファイル（md5 と大きさ）を git で追跡し、
  中身はキャッシュに置く
- 案 B: 既存の dvc.yaml（examples/depletion/dvc.yaml など）にステージを足し、生成物をその出力にする。記録は、分岐の
  テーブルをステージにしない理由を「履歴（p2 の `lattice@2A`）を deps にすると、このブランチで lattice の実行ファイルが
  変わったために履歴のステージ（約 1 時間 41 分）が changed になり、目標を指定した `dvc repro` がそれを再実行するため（ADR-0184）」
  と書き（examples/depletion/verification.md の「燃焼度の軸を持つテーブル（2A）」）、「履歴の効き（2A）」、「ミクロ燃焼の
  線形の補正（2A）」、「履歴の残り R（2A の減速材の密度）」も同じ理由でステージにしていない。生成物の多くは、ブランチの
  `target/release` の実行ファイルで手で流した試し（examples/depletion/branches、history、micro には dvc.yaml が無い）で、
  ステージにすると deps の実行ファイルとコマンドをステージの形に直すことになり、ステージの再実行が要る
- 案 C: `worktree-restore-outputs.sh` の逆向きのスクリプトで、ブランチの記録が引く無視されたファイルを、マージの前に
  main にコピーする。/retro の 25 回目にシェルの手順として試し、ローカルの /code-review 2 周で cwd、glob、`git diff` の行の
  拾い方、DVC の出力の除き方、main のファイルの上書き、コピーの段が /git-merge に無いことの指摘が出て捨てた（#0207 の
  「経緯」）。同じ名前のファイルを上書きするので、前の版は残らない
- 案 D: マージの前に、記録が引く無視されたファイルを持つ worktree の削除を止める。`post-merge-cleanup.sh` はすべての
  リポジトリで共有するので xlc の側の検査として置くことになる。止めた後の置き場は手で決めることになり、前の版も残らない

## 決定と理由

案 A にする。ユーザーが「8は追ったらいいと思います」（ユーザー 2026-10-31。「8」は判断待ちの一覧の番号で、#0207 の修正方法の
`dvc add` で追う案を指す）と決めた。

- 対象: 記録か ADR が md5 で引くか、記録のコマンドの入力や入力の表で名前で引く、git で追跡しない生成物で、DVC のステージの
  出力でないもの。置き場は記録が引くパスのまま。コマンドの出力先（`> x.out 2> x.log`）としてだけ名前が出るログと k の出力は
  追わない。`references/` の下（OpenMC の作業ディレクトリ、AWS の結果。#0018）は main の worktree にだけあり、後片付けで
  消えないので対象にしない
- 時機: 記録に md5 か名前を書く前に `dvc add <ファイル>` し、`.dvc` ファイルを同じブランチでコミットする。同じ名前で作り
  直したら `dvc add` をやり直す
- 戻し方: 引数なしの `dvc checkout` は dvc.yaml の出力と `.dvc` ファイルの対象の両方を戻すので、
  `scripts/worktree-restore-outputs.sh` はそれだけを実行する（`openmc-2K.msgpack` のコピーは、p2 の参照解を `.dvc` で追う
  ので要らない）

理由:

- 案 B と違い、`.dvc` ファイルは deps もコマンドも持たないので、実行ファイルや上流の履歴が変わっても changed にならず、
  再実行を招かない。生成物を作ったコマンドと実行ファイルは、今までどおり記録の「計算」に書く
- キャッシュは中身の md5 で置くので、置き換えても前の版はキャッシュに残り、前の版を指す `.dvc` は git の履歴に残る
  （`dvc gc --all-commits` も消さない。ADR-0123 の帰結）。案 C と D は、同じパスの上書きで前の版を失う
- キャッシュは worktree の間で共有し（ADR-0126）、reflink なので `dvc add` と `dvc checkout` は中身をコピーしない。費用は
  ファイルごとの `.dvc`（数行）と、記録に書く前の `dvc add` の 1 手である
- 戻す経路が dvc.yaml の出力と同じ `dvc checkout` になり、生成物ごとの特別なコピーが要らない

### 帰結

- 追う生成物（2026-10-31。main の作業ツリーから `cp --reflink` でコピーし、md5 が同じことを確かめて `dvc add` した）。
  43 件、543,697,755 バイト:
  - examples/depletion/history/（13 件）: core の履歴 `core-2A-{plain-nominal,plain-switch,residual-branch,residual-nominal,residual-switch,sigma-branch}.msgpack`、
    格子の履歴 `lattice-residual-density.msgpack`（a30302d3）、`lattice-switch-density.msgpack`（1840f666。どちらも 55 MB）、
    テーブル `table-residual-{density-branch,density-reverse,nominal-histories,nominal}.msgpack`、`table-switch-density-branch.msgpack`
    （`table-residual-nominal-histories.msgpack` は今の形の 222fabb0）
  - examples/depletion/branches/（12 件）: `core-2A-coarse.msgpack`、`core-2A-steps.msgpack`、履歴 `lattice-2A-6272e937.msgpack`（55 MB）、
    テーブル `table-2A-{all,cross,cross-rods,fine,micro-boron,micro-boron-rods,micro-boron-rod-axis}.msgpack`、`table-2O-all.msgpack`、
    `table-reflector-I-boron.msgpack`
  - examples/depletion/issue-0133/ufg-free-gas-assemblies/（2 件）: `lattice-2A-fg.msgpack`（32645292、55 MB）、`lattice-2O-fg.msgpack`
    （d7692136、130 MB）。記録はこの置き場で読むコマンドを書くが、main の `trials/ufg-free-gas-assemblies/` にだけ残っていた
  - examples/branches/（2 件）: 多群ライブラリ `branches.msgpack`（188 MB）、テーブル `table-2B.msgpack`（7d4dfed5）
  - examples/vera/p2/（12 件）: OpenMC の参照解のピン出力 `openmc-{2A,2B,2C,2D,2E,2E-empty,2F,2G,2H,2I,2K,2Q}.msgpack`（ADR-0121）
  - examples/mox/（2 件）: OpenMC の燃料の実効断面積 `openmc-mox1.msgpack`、`openmc-mox2.msgpack`（ADR-0121。作業ディレクトリに残っていない）
- 追わなかったもの: 記録のコマンドの出力先としてだけ名前が出る `.log` と `.out`（branches、history、examples/branches）、
  記録が引かない `table-reflector-I-boron-fuel.msgpack`、`references/openmc/` の下の OpenMC の結果（記録が md5 で引く
  `aws/m4-assembly/openmc-2A/openmc-2A.msgpack` など）
- 失われたもの（記録が md5 で引くが、`/home/tagawa/github` の下の `*.msgpack` と DVC のキャッシュを md5 で探して 0 件。
  2026-10-31）。作り直さず、それを読む作業のときに、記録のコミットと入力で再実行する（CLAUDE.md の「計算」の後回し）:
  - micro（「ミクロ燃焼の線形の補正（2A）」）: テーブル nominal fde256d9、置き換えた組成 da7265c2、ε = 0.05 f072e501、
    置き換えた組成の履歴 `shifted-2A.msgpack` 456654d5
  - history: 前の形の `table-residual-nominal-histories.msgpack` 3d73faee
  - branches: `table-2A.msgpack` 3ac7a83e、`table-2F-all.msgpack` 78ad62e3、`table-2O-guide.msgpack` 7f323f0e、
    `table-2A-boron.msgpack` 5f443d2e、`table-2A-micro.msgpack` be25725e と 646d0ec2、制御棒を入れた組 70b204f3、core の履歴
    b46d7ecb、f78e0822、8e547efb、3a9232a1、e545a204、0b4a44cd、9ac68868、72565c38
  - issue-0129 ほかの燃焼の試し: `2A-step-start.msgpack` 355ddae2、`2O-step-start.msgpack` 72b2b84c、
    `2A-10-every-transport.msgpack` df58a8b8、`2A-10-step-start.msgpack` 803bce96
- main の作業ツリーには、コピーした元の生成物が無視されたファイルとして残る（md5 は `.dvc` と同じ）。`trials/ufg-free-gas-assemblies/`
  の元も消していない
- `.dvc` ファイルは `dvc status` に出る（ファイルが無いか変わったとき）。`scripts/count-changed-stages.py` は、コロンを持たない
  `.dvc` のキーを直下の dvc.yaml のステージと分けて `*.dvc` の行に数える
- `.dvc` を足すか変えたブランチをマージしたあとは、main で `scripts/dvc-checkout.py <その .dvc>` で置き場をキャッシュに合わせる
  （スキル dvc-relock の後半）。ブランチで作り直して `dvc add` をやり直していない生成物があると、`dvc checkout` は
  「Can't remove the following unsaved files」で全体が止まる（ほかの出力も戻らない。スキル dvc-relock の前半 1）
- `dvc repro` は、目標にした `.dvc` の対象が変わっていると、`.dvc` の md5 を今のファイルに書き直す（コミットはしない。`-P` と上流の扱いは確かめていない。
  2026-10-31 に `dvc repro examples/mox/openmc-mox2.msgpack.dvc` で確かめた）。前の版はキャッシュに残る
- 規則の 1 行は CLAUDE.md の「計算」にある。新しい生成物を追うのは記録を書く担当の 1 手に頼り、記録の md5 と `.dvc` を
  突き合わせる検査は無い（#0209）

### 確認

- `scripts/test_count_changed_stages.py` の `test_dvc_files_are_counted_apart_from_the_root_yaml`
- 新しい worktree（`git worktree add --detach` に `references` の symlink を張ったもの）で `scripts/worktree-restore-outputs.sh` を
  実行し、終了コード 0（未実行のステージの出力 7 件は `UNRUN:` に数えるだけ）で、追った 43 件の md5 が main の作業ツリーと
  一致した（2026-10-31）
