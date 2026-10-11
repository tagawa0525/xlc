---
status: accepted
date: 2026-10-11
supersedes: []
superseded-by:
issues: ["#0018", "#0060"]
requires: []
---

# ADR-0123: 計算の段階は DVC で駆動し、dvc.yaml と dvc.lock をコミットする

## 背景

xlc の計算は、入力のファイルから出力のファイルを作る段階の連なりである（xs → ライブラリ → lattice → core、NJOY →
OpenMC のライブラリ → OpenMC、xlc と OpenMC の比較）。どの段階を実行するか、もう実行してあるかは、人と run shell
（examples/m1/run.sh など）が判断していた。何が何から作られたかは、記録に添えたコミットの hash と手順から辿るしかなく、
手元の生成物（ライブラリ、OpenMC の結果）がどのコードと入力から作られたものかは、どこにも残っていない。

ADR-0122 は、記録の hash を checkout すれば当時の結果を出し直せることを求め、結果を変えない変更（rebase、文書）では
再実行しないと決めた。その手段は #0018 で整える。

判断の軸:

- 実行するかどうかを、入力とコードの内容で決めること（rebase、checkout、文書の変更で走らないこと）
- 何が何から作られたかを、git にコミットできる形で残すこと
- デファクトのツールを使い、保守を自分で持たないこと
- 10 GB 級の生成物と、Slurm で 1 本ずつ流す運用に合うこと

## 検討した案

駆動のツール:

- 案 A（採用）: DVC（nixpkgs の python3Packages.dvc 3.67.1）。段階を `dvc.yaml`（cmd、deps、outs）に書き、各段階の
  deps と outs の md5 を `dvc.lock` に書く。判定は内容の md5 で、`dvc.lock` は YAML なので git の差分で読める。
  出力は内容でアドレスするキャッシュに置き、`dvc checkout` で戻せる。同じ deps で一度実行した段階は、run-cache から
  出力を戻して実行しない。利用者が多く、保守は外にある
- 案 B: xvc（0.7.1-alpha）。作者が 1 人で、2026-10-11 の試験で、TOML のキーへの依存が panic して hang した。出力の
  digest を記録せず、`.xvc/` は追記型の JSON で、ツールなしには読めない
- 案 C: make。mtime で判定するので、checkout と rebase のたびに無駄に走る。何から作ったかの hash を残せない
- 案 D: just。コマンドの短縮であって、依存の状態を持たないので駆動にならない
- 案 E: Nix で段階を derivation にする。環境の固定には使う（flake.nix）が、10 GB の出力を store に置くことと、
  Slurm でジョブを流すことに合わない
- 案 F: Snakemake。Slurm の executor は良いが、判定は mtime が主体で、DSL（Python の上の独自の文法）が大きい
- 案 G: 自前の Rust の駆動。デファクトのツールを使う方針に反し、判定と記録の保守を自分で持つことになる

段階の deps に、その段階のバイナリを表すものとして何を書くか:

- 案 H（採用）: コードのソース。段階のクレートと、パスで依存するクレートの `src` と `Cargo.toml`、ワークスペースの
  `Cargo.toml` と `Cargo.lock`、ツールチェーン（`rust-toolchain.toml`）、コンパイルの設定（`.cargo/config.toml`、
  ADR-0109）、リンクする glibc とリンカを決める `flake.lock`。リポジトリの外の設定（`~/.cargo/config.toml`、環境変数の
  `RUSTFLAGS`）は deps にできないので、変えるときは `dvc repro --force` で再実行する。内容の md5 で判定するので、rebase や文書の変更では走らない。`src` の中のテスト（`#[cfg(test)]` の
  モジュール）を変えても古くなるが、これは許す（テストを `src` の外に出すのは、テストの置き方を変える別の決定になる）
- 案 I: ビルドしたバイナリ（`target/release/xs` など）。結果を変える変更だけで走りそうに見えるが、バイナリには
  worktree の絶対パスが埋まる（`target/release/xs` に `/home/tagawa/github/xlc-feat-dvc-pipeline/outputs/out/lib` が
  あった。2026-10-11）ので、同じコードでも worktree ごとに md5 が変わり、別の worktree では再実行する。`dvc.lock` から
  どのコードかも読めない。同じパスでビルドし直したときにビット一致するかは確かめていない

キャッシュの置き方:

- 案 J（採用）: reflink。`~/github` は btrfs（`df -T`、2026-10-11）で、reflink の写しは書き込むまでブロックを共有するので、
  ワークスペースとキャッシュで容量が倍にならない。写しは互いに独立なので、ツールが出力を同じパスに上書きしてもキャッシュは
  壊れない。reflink のできないファイルシステムでは、コピー（copy）に黙って戻らず、DVC がエラーで止まるよう `copy` を
  並べない
- 案 K: hardlink。容量は倍にならないが、ワークスペースのファイルがキャッシュそのものなので、DVC はそれを読み出し専用に
  する。ツールがその場で書き換えるとキャッシュが壊れ、書き換えるには `dvc unprotect` が要る
- 案 L: copy（DVC の既定は reflink、copy の順）。どのファイルシステムでも動くが、10 GB 級の出力で容量が倍になる

キャッシュの置き場は、DVC の既定（`.dvc/cache`）だと worktree ごとに分かれ、worktree を作るたびに段階を作り直すことに
なる。`.dvc/config` の `cache.dir` を `references/dvc-cache` にし、worktree では `references` を main の `references/`
への symlink にして、キャッシュを共有する。キャッシュは `references/` と同じく、消しても `dvc repro` で作り直せる
作業の置き場である。

大きな出力（OpenMC の statepoint）の扱い:

- 案 M（採用）: statepoint を OpenMC の段階の outs にしてキャッシュし、比べる小さな結果（ピン出力の msgpack など）は、
  statepoint を読む別の段階で抽出する。statepoint は 1 問題に約 100 MB（references/openmc/vera-xs-2B/statepoint.400.h5 が
  109 MB、21 本、2026-10-11）で、reflink なら容量は倍にならない。作り直すと OpenMC の計算に 1 問題 16〜34 分かかる
  （ADR-0121）ので、抽出の仕方を変えても OpenMC を走らせ直さずに済む
- 案 N: 抽出した小さな結果だけを outs にし、statepoint は作業ディレクトリに置いて追跡しない。キャッシュは小さいが、
  抽出を変えるには OpenMC を再実行するか、DVC の外で statepoint を読むことになり、何から作ったかが `dvc.lock` に残らない
- 案 O: statepoint を `cache: false` の outs にする。md5 は `dvc.lock` に残るが、キャッシュにないので `dvc checkout` で
  戻せず、消えたら OpenMC を再実行することになる

## 決定と理由

案 A、H、J、M にする。

- 段階は `examples/<問題>/dvc.yaml` に書く（問題ごとの置き場所は ADR-0116）。段階のコマンドは作業ディレクトリをその
  ディレクトリにして、素の `cargo run --release -p <クレート> -- <入力>` にする。結果を標準出力に出す段階（lattice の
  k∞、core の k_eff）は、標準出力をファイルに受けて outs にする（run-cache や `dvc checkout` で戻したときにも結果が
  残るように）。Slurm や手元のツール（~/.claude のスクリプト）に依存させない。実行するのは
  `dvc repro examples/<問題>/dvc.yaml`（すべては `dvc repro -P`）で、1 分を超えるときは `dvc repro` 全体を Slurm で流す。
  段階は逐次に実行する（Slurm で 1 本ずつ流す今の運用と同じ）
- deps には、段階の入力（入力の TOML と、そこから読むファイル）と、案 H のコードを並べる。入力の TOML は、`params` で
  キーごとに選ばず、ファイル全体を deps にする。入力のどのキーも結果を変えうるので、キーを並べ直す手間と漏れを増やす
  だけである。コメントだけの変更でも走るが、戻せば run-cache から出力を戻す
- `dvc.yaml` と `dvc.lock` をコミットする。記録の hash を checkout すれば、`dvc.lock` が、どの出力がどの入力とコードから
  作られたかを示す。`dvc checkout` でキャッシュから出力を戻し、キャッシュになければ `dvc repro` で作り直す。これが
  ADR-0122 の「hash から出し直せる」の手段になる
- DVC の利用状況の送信（`core.analytics`）と更新の問い合わせ（`core.check_update`）は、`.dvc/config` で切る

内容の md5 で判定するので、ADR-0122 の「結果を変えない変更で再実行しない」が、ツールの判定として成り立つ。`dvc.lock` は
テキストなので、何が変わって走ったかを git の差分で読める。

## 確認

M0（examples/m0/dvc.yaml の xs → lattice → core）で確かめた（fef2957a、2026-10-11）。

- 最初の `dvc repro examples/m0/dvc.yaml`（Slurm）が 3 段階を順に走らせ（U-238 の xs が 55 s、core が 9.4 s）、
  `dvc.lock` を書いた。`wims69.msgpack`、`pincell-fg.msgpack`、`core.out` の md5 は、xs から実行した 599870ee の
  1 回目と一致した（2 回とも Slurm の上で同じスレッド数）
- 2 回目は 3 段階とも「didn't change」で何も走らせなかった（1.6 s）
- `pincell.toml` のピッチを 1.26 から 1.27 cm に変えると、xs は走らず、lattice と core が走った（k∞ 0.731552 →
  0.741536）。戻すと、lattice と core は run-cache から出力を戻して走らず、`dvc.lock` は最初のものとバイト単位で一致した
- `pincell.toml` のコメントだけを変えると lattice は走ったが、出力の `pincell-fg.msgpack` の md5 が変わらないので、
  core は走らなかった
- キャッシュは `references/dvc-cache` にでき、消した `wims69.msgpack` を `dvc checkout` が reflink で戻した（リンク数 1 で、
  hardlink でない）。reflink のできない tmpfs（/dev/shm）で同じ設定にすると、`dvc checkout` は copy に戻らず
  「no more link types left to try out」で止まった

## 帰結

- 段階を DVC で実行しているのは M0 だけである（M0 の run.sh は消した）。ほかの例（M1、VERA、MOX など）と、参照の結果
  （NJOY → OpenMC のライブラリ → OpenMC、比較）を DVC のステージにするのは #0018 の次のブランチで行う。ツール（NJOY、
  OpenMC）の作り方は ADR-0031 のままで、変えるならそのブランチの ADR が置き換える
- キャッシュは、同じ出力の古い版も残すので増えていく。`dvc gc` で消すが、キャッシュは worktree の間で共有するので、
  `--workspace` だけで消すと、ほかの worktree とブランチの出力が消える。消すときは `--all-commits`（どのブランチの
  コミットの `dvc.lock` が指すものも残す）にし、まだコミットしていない出力のある worktree がないことを確かめる
- deps のクレートの並びは、各 `Cargo.toml` のパスの依存を手でコピーしたもの、ENDF のファイルは入力の TOML をコピーしたもので
  ある。片方だけを直すと、古い出力が黙って最新とみなされる。ほかの例に広げるとき（#0018 の次のブランチ）に、
  コピーを作るツールか、食い違いを見つける検査を考える
- `.dvcignore` で `target/` と `references/dvc-cache` をたどらない。`references/` のほかの場所は、ENDF などを deps に
  するので無視しない
- main の examples/m0/ に DVC を入れる前の生成物（`wims69.msgpack`、`pincell-fg.msgpack`）が残っていると、`dvc.lock`
  と合わない出力になる。マージの後に消し、`dvc checkout examples/m0/dvc.yaml` でキャッシュから戻す
- `src` の中のテストを変えた段階は、結果が変わらなくても走る（案 H）
- DVC は出力をバイトの一致（md5）で記録し、数値の許容幅をもたない。再実行の引き金は deps なので、走らせなければ
  出力はキャッシュから戻り、揺らがない。並列の和の順序が実行方法（スレッド数）で変わって出力のビットが揺らぐと、
  効くのは次の 3 つである: 上流を再実行したときに下流の段階まで走る連鎖、`dvc.lock` の意味のない差分（ブランチの
  間の衝突）、DVC の外で実行した出力の「changed outs」。本筋の対策は、xlc の出力をスレッド数によらずビット一致させる
  ことである（和の順序の固定。MOC のスイープは ADR-0064 で固定している）。段階の出力でそれを確かめるのは #0060
- DVC の排他（`.dvc/tmp/lock`）は worktree ごとなので、共有のキャッシュに 2 つの worktree から同時に `dvc repro` や
  `dvc gc` を実行すると、互いに排他されない。重い段階は Slurm の排他で 1 本ずつ流すので重ならないが、手元で実行する
  ときと `dvc gc` は、ほかの worktree で DVC が動いていないときに行う
- 共有のキャッシュは、worktree で `references` を main への symlink にすることに頼っている。symlink がないと、DVC は
  worktree の中に別のキャッシュを黙って作る（2026-10-11 に一時的なリポジトリで確かめた。ENDF を deps にする段階は、deps がなくて止まる）。worktree を作るツールで
  symlink を作るか、置き場が違えば止める仕組みは、まだない
