---
status: accepted
date: 2026-10-23
requires: [ADR-0067, ADR-0130, ADR-0138, ADR-0159, ADR-0184, ADR-0189]
supersedes: []
superseded-by:
issues: ["#0157"]
---

# ADR-0193: 核種ごとのステージを、DVC の外の runner が前回の CPU 時間の長い順に同時に実行し、dvc commit で lock に取り込む

## 背景

[ADR-0189](0189-fill-cores-within-nuclide-stages.md) は核種ごとのステージ（[ADR-0138](0138-split-xs-into-nuclide-stages.md)）を保ち、
核種の中の並列でコアを埋めた（案 D）。`dvc repro` はステージを 1 つずつ実行するので、examples/depletion の 65 核種で核種のプロセスの和は
1,154 秒で、CPU 時間の和の 1/16（710 秒）の 1.6 倍かかる。changed の核種のステージを、前回の CPU 時間の長い順に、スレッドの数
（`xs --threads`）を決めて同時に実行する試行（ADR-0189 の案 C）は、同じジョブ（Slurm 479）で 720 秒だった。DVC と `dvc commit` の時間は
含まない。ユーザーは案 C を採ると決めた（2026-10-23「速さは正義だと思います。そこまで複雑でもないですし」。元の案もユーザー、2026-10-23）。

決めることは、DVC の外で実行するときに `dvc.lock` が出力を作った deps を指し続けること（[ADR-0159](0159-depend-stages-on-built-binaries.md) の
案 B が退けた「lock が、記録した deps のコードで実行していない出力を指す」を起こさないこと）、前回の CPU 時間の置き場、スレッドの
割り当て、runner の形、失敗の扱いである。

判断の軸:

- `dvc.lock` が、出力を作った cmd と deps を指すこと（ADR-0122 の「hash から出し直せる」）
- 出力が `dvc repro` と同じこと（ADR-0189 の確認で、`xs nuclide` の出力はスレッドの数によらない）
- 65 核種の経過時間（DVC と `dvc commit` を含む）
- 保守の層（例ごとの dvc.yaml の変更、DVC の内部への依存、手で書く写し）を増やさないこと

## 検討した案

前回の CPU 時間の置き場:

- 案 R1（採用）: dvc.yaml の隣の `stage-cpu-seconds.json`（ステージの名前 → CPU 秒）を git で追跡する。runner がステージを
  `dvc commit` するたびに更新し、`dvc.lock` と一緒にコミットする。どの worktree と Slurm のジョブからも同じ値が見え、値を出した実行が
  コミットで分かる
- 案 R2: Slurm のログか DVC のキャッシュの隣など、git の外に置く。ログは消えうる、worktree ごとに違いうる、読み取りの形が決まらない
- 案 R3: 記録を持たず、ENDF の大きさ × 温度の数から見積もる。CPU 時間とほとんど相関しない（ADR-0189。対数の相関係数 0.12）

スレッドの割り当て:

- 案 T1（採用）: 実行するステージの CPU 時間の和を C、予算を B とし、各ステージのスレッドの数を ceil(c / (C / B)) を 1〜B に抑えたものにする
  （#0157 と ADR-0189 の試行の式）。長い順に並べ、空いたスレッドに収まるものを前から順に投入する（収まらないものを飛ばして、後ろの
  小さいものを先に入れる。試行と同じ）。収まらないものは走っているものが終わるのを待つ（空きの少ないスレッドで始めると、そのステージが
  最後まで遅い）。予算より多くを求めるステージは、何も走っていないときに予算の全部で投入する
- 案 T2: 全ステージを 1 スレッドにする。最長の核種の CPU 時間（Am-242m 2,278 秒）が経過時間の下限になる（ADR-0189）

runner が cmd と deps を知る方法:

- 案 C1（採用）: DVC の CLI に聞く。changed のステージは `dvc status --json --with-deps <yaml>:<foreach の名前>`、cmd は
  `dvc repro --dry -f -s <ステージ>...` の出力（`Running stage '<名前>':` と `> <cmd>`）、記録した deps は `dvc.lock` から読む
- 案 C2: dvc.yaml のテンプレート（`foreach`、`vars`、`${item.endf}`）を runner が展開する。DVC の解決を作り直すことになり、default の devShell の
  Python には YAML の読み手がない
- 案 C3: DVC の Python API（`dvc.repo.Repo`）を使う。default の devShell は、DVC の依存が検証用の Python より先に読まれないように、
  DVC の Python を PATH の python3 から分けている（flake.nix の dvc の注釈、#0101）。runner のためにそれを戻すことになる

`--threads` の置き場:

- 案 A1（採用）: runner が実行するときだけ、cmd の実行ファイルの直後に `--threads <数>` を足す。dvc.yaml の cmd は変えない
- 案 A2: dvc.yaml の cmd に `--threads` を書く。`dvc.lock` が cmd を記録するので、スレッドの数が変わるたびにステージが changed になる

実行中に deps が変わったとき:

- 案 D1（採用）: runner の始めに、見るべきパス（dvc.yaml、`dvc status` が挙げた deps、`dvc.lock` が記録した deps と params のファイル）の
  指紋（ファイルごとの大きさ、mtime、inode。ディレクトリは中のファイルすべて）を取る。`dvc commit` の直前に指紋を比べ、変わっていれば
  commit せずに止まる。commit の直後に、`dvc.lock` が記録したそのステージの deps がすべて見たパスに入り、指紋が変わっていないことを
  確かめる（commit が hash を取った間も deps が始めのままだったことの確認）。dvc.lock の cmd が実行した cmd と同じことも確かめる（dvc.yaml の
  vars のファイルは指紋で見ないので、cmd はこちらで見る）。成り立たなければ dvc.lock を commit の前に戻し、`dvc commit` が足した run-cache の
  項目（DVC 3.67.1 の `Stage.save` が保存する。新しい deps を、それで作っていない出力に結ぶ）を消して止まる。消すのは commit の間に増えた
  項目のうち commit した出力の md5 を持つものだけ（キャッシュは worktree で共有し、ほかのセッションの項目を消さない）
- 案 D2: commit の前後で md5 を計算して比べる。ステージごとに熱中性子散乱のディレクトリ（225 MB）などを読み直す。DVC 自身も、
  hash を計算し直すかを大きさ、mtime、inode で決めている（state）ので、指紋で足りる
- 案 D3: 確かめない。runner の外で同じ worktree の deps を変えたとき（ビルドし直した、xs.toml を直した）、`dvc.lock` が実行していない deps を指す

失敗の扱い:

- 案 F1（採用）: ステージのプロセスが失敗したら、新しいステージは投入せず、走っているものは終わらせて commit し、失敗したステージを
  並べて終了コード 1 で終わる。`dvc commit` できないステージ（出力が無い、など。まとめた commit が失敗したら 1 つずつ commit し直して
  見分ける）も同じくそのステージの失敗にする。deps の変化など runner の側の誤りは、走っているプロセスを止めて終了コード 1 で終わる
- 案 F2: 1 つ失敗したら全部を止める。ほかの核種は失敗したステージと独立で、終わっていれば正しい出力なので、捨てる理由がない
- 案 F3: 失敗しても残りを全部投入する。`dvc repro` の既定（最初の失敗で止まる）と違い、同じ原因の失敗を繰り返しうる

## 決定と理由

案 R1、T1、C1、A1、D1、F1 にする。runner は `examples/stage_runner.py`（標準ライブラリだけ）で、default の devShell で実行する
（dvc は references の devShell では実行しない。ADR-0130）。

```sh
python3 examples/stage_runner.py examples/<問題>/dvc.yaml xs-nuclide [--threads <予算>]
```

- 対象は foreach のステージの組 1 つ（`<yaml>:<名前>`）。`dvc status --with-deps` が組の外（上流の build、xs-range）を changed と示せば、先に
  `dvc repro <それらのステージ>` で実行する（`dvc repro <yaml>:<名前>` と同じく上流を飛ばさない。ADR-0184）。`dvc status` は上流が出力を
  書き換えた後に changed になる下流を示さない（build の後の xs-range）ので、組の外が changed でなくなるまで繰り返す。同じ上流が続けて
  changed なら（always_changed など）止まる
- deps も cmd も変わらず、出力だけが lock と違うステージ（出力を戻していない worktree など）があれば止まる。`dvc repro` は run-cache から
  戻すものを、runner は作り直さない
- 予算の既定は、このプロセスに許された CPU の中の物理コアの数（ADR-0067）。`--threads` で変える。xs の既定（ノードの物理コアの数を許された
  CPU の数で抑えたもの。ADR-0189）は、Slurm の `-c 16` で 8 コアに 16 を返すので、使わない（xs の側は #0158）。
  同時に走るプロセスの数は予算以下になる（各ステージに 1 スレッド以上）。メモリは見ない: 1 核種のプロセスは最大 998 MB（ADR-0189）で、
  16 本でも r995 の 91 GB に収まる
- 記録の無いステージ（初めての核種）の CPU 時間は、同じ組の記録の最大値とみなす（長いものとして先に、多めのスレッドで流す）。組の記録が
  1 つも無ければ、すべて同じ値とみなし、予算を等分する。記録は 0.1 秒の刻みで、0 のステージ（刻みより短い）は 0 として扱い、すべて 0 なら
  予算を等分する。負の値、数でない値は誤りとして止まる
- cmd の実行ファイルの名前が `xs` でなければ止まる（`--threads` を受けるのは xs だけ）。dvc.yaml に `wdir` があれば止まる（cmd を
  dvc.yaml のディレクトリで実行し、`dvc.lock` のパスをそこから読む）
- CPU 時間は `wait4` の rusage（ユーザーとシステム。シェルが待った子を含む）
- リポジトリの直下の dvc.yaml の組は受けない（DVC はそのステージを `dvc.yaml:` を付けずに示すので、組と上流を見分けられない）
- runner が SIGTERM を受けたら、走っているステージ（それぞれ自分のセッションに置く）を止め、待っている dvc は SIGTERM で止めて終わるのを待つ
  （SIGKILL すると `.dvc/tmp/rwlock` と dvc.lock を書きかけで残す）
- 終わったステージは、空いたスレッドに次を投入してから、まとめて `dvc commit -f` する。commit できないステージが出たら、その後は投入
  しないが、同じ回に commit の前に投入したステージは走り切る
- ステージを投入する直前に、その出力（dvc.lock の outs と、`dvc status` の changed outs）を消す（`dvc repro` と同じ。投入しない
  ステージの出力は残す）。xs が出力を書かずに
  0 で終わったとき、古い出力を新しい deps で commit しない（commit が出力の無さで失敗し、そのステージの失敗になる）
- 記録のうち、今の組（`dvc stage list`）に無いステージの値は消す（記録の無いステージの見積もりに、消した核種の値を使わない）
- 後始末で、SIGTERM の後 10 秒で終わらないステージはセッションごと SIGKILL する。commit、確認、dvc.lock の戻し、記録の書き込みの間は SIGTERM と SIGINT を保留する

ADR-0159 の案 B（出力を変えないとみなした変更のあとに `dvc commit` で deps だけを取り込む）は、記録する deps で実行していない出力を
lock が指すので退けた。runner は、今の deps（今の `bin/xs`、入力、ENDF）で、dvc.yaml の cmd と同じ計算（`--threads` だけが違い、出力は
スレッドの数によらない。ADR-0189）を実行し、その deps を commit する。実行の間に deps が変わらなかったことを案 D1 で確かめるので、案 B の
理由には当たらない。

## 確認

- examples/depletion の 65 核種（examples/depletion/verification.md の「燃焼の多群ライブラリ」。xlc c92e8cac、分岐点 2b6e1201、2026-10-23、Slurm の
  排他のジョブ 485）: main の lock で bin/xs が変わって changed だった xs-range と 65 核種を runner が実行して commit し、続けて
  `dvc repro examples/depletion/dvc.yaml:xs` を実行した。67 のステージの出力の md5 は main の dvc.lock（Slurm 473）と同じで、lock の差は
  deps の bin/xs だけ。runner の始めから最後の `dvc commit` まで 827 秒（上流 54 秒、核種の最初の投入から最後の回収まで 753 秒、最後の commit 20 秒）。
  CPU 時間の和の 1/16 は 735 秒で、ADR-0189 の案 D（核種のプロセスの和 1,154 秒と、DVC のステージの間）より短い
- 単体テスト（examples/test_stage_runner.py。scripts/run-python-tests.sh が実行する）: スレッドの割り当てと順序、記録の無いステージの扱い、空きへの
  投入と待ち、`dvc status`、`dvc repro --dry`、dvc.lock（折り返した cmd、deps、params、outs）の読み方、指紋（書き直し、追加、symlink の先、読めない
  ディレクトリ）、commit の前後の確認と dvc.lock と run-cache の戻し、commit できないステージ、上流の繰り返し、失敗と SIGTERM の扱い、偽の xs と dvc での実行
- 一時の DVC のリポジトリ（偽の xs）で、上流の実行、commit、実行の途中で deps を変えたときに止まること、出力だけが違うステージで止まることを確かめた
  （手元、2026-10-23）

## 帰結

- 例ごとの dvc.yaml は変えない。`dvc repro` だけでも従来どおり作れる（直列で遅いだけ）
- `dvc.lock` の cmd には `--threads` が入らない。lock の cmd を `dvc repro` で実行し直しても、出力は同じ（ADR-0189）
- 記録の `stage-cpu-seconds.json` は計算の結果に効かない（並べる順とスレッドの数だけを決める）。DVC の deps にも outs にもしない
- 指紋の比べ方は、内容を変えずに mtime だけが変わったとき（`touch`、`dvc checkout` で戻したとき）も止まる（誤って止まる側に倒す）
- 空いたスレッドに収まらないステージを飛ばして後ろの小さいものを先に入れるので、収まらないステージは後回しになりうる。ceil の和は
  予算を超えうるので、ステージが少ない組では、収まらないステージを待つ間コアが遊ぶ（同じ長さの 3 本を 16 スレッドで 6、6、6 にすると、
  3 本目は 1 本が終わるまで待つ）。スレッドを多く取るステージは CPU 時間の長いもので、並びの先頭にあり、始めにすべての予算が空いている
  ときに投入される（examples/depletion の記録で、Am-242m 4、Eu-151 2、Pu-239 2 と 1 スレッドの 8 本が始めに入る）
- `dvc commit` は投入と同じループで順に実行する。commit の間に終わったステージの回収と次の投入は commit の後になる。65 核種で commit は 38 回、
  1 回 2〜4 秒（最後の 15 核種で 20 秒）で、回収の遅れは 1 スレッドのステージの経過時間 − CPU 時間の和で 172 秒（16 で割って 11 秒）以下
  （Slurm 485）。commit を別のスレッドに分けて縮むのはこの分と最後の commit だけなので、分けない
- `dvc.lock` は、record_provenance.py の `parse_stage`（DVC が書く形を前提にした、標準ライブラリだけの読み手）で読む。DVC の lock の形が
  変われば、読めずに止まる。params のファイルの読み取りは runner が持つ（同じ形の読み手が 2 つ。3 つ目が要るときに parse_stage へまとめる）
