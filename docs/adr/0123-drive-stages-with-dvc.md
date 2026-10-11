---
status: proposed
date: 2026-10-11
supersedes: []
superseded-by:
issues: ["#0018"]
---

# ADR-0123: 計算の段階は DVC で駆動し、dvc.yaml と dvc.lock をコミットする

## 背景

xlc の計算は、入力のファイルから出力のファイルを作る段階の連なりである（xs → ライブラリ → lattice → core、NJOY →
OpenMC のライブラリ → OpenMC、xlc と OpenMC の比較）。どの段階を走らせるか、もう走らせてあるかは、人と run shell
（examples/m0/run.sh など）が判断していた。何が何から作られたかは、記録に添えたコミットの hash と手順から辿るしかなく、
手元の生成物（ライブラリ、OpenMC の結果）がどのコードと入力から作られたものかは、どこにも残っていない。

ADR-0122 は、記録の hash を checkout すれば当時の結果を出し直せることを求め、結果を変えない変更（rebase、文書）では
流し直さないと決めた。その手段は #0018 で整える。

判断の軸:

- 走らせるかどうかを、入力とコードの内容で決めること（rebase、checkout、文書の変更で走らないこと）
- 何が何から作られたかを、git にコミットできる形で残すこと
- デファクトの道具を使い、保守を自分で持たないこと
- 10 GB 級の生成物と、Slurm で 1 本ずつ流す運用に合うこと

## 検討した案

駆動の道具:

- 案 A（採用）: DVC（nixpkgs の python3Packages.dvc 3.67.1）。段階を `dvc.yaml`（cmd、deps、outs）に書き、各段階の
  deps と outs の md5 を `dvc.lock` に書く。判定は内容の md5 で、`dvc.lock` は YAML なので git の差分で読める。
  出力は内容でアドレスするキャッシュに置き、`dvc checkout` で戻せる。同じ deps で一度走らせた段階は、run-cache から
  出力を戻して走らせない。利用者が多く、保守は外にある
- 案 B: xvc（0.7.1-alpha）。作者が 1 人で、2026-10-11 の試験で、TOML のキーへの依存が panic して hang した。出力の
  digest を記録せず、`.xvc/` は追記型の JSON で、道具なしには読めない
- 案 C: make。mtime で判定するので、checkout と rebase のたびに無駄に走る。何から作ったかの hash を残せない
- 案 D: just。コマンドの短縮であって、依存の状態を持たないので駆動にならない
- 案 E: Nix で段階を derivation にする。環境の固定には使う（flake.nix）が、10 GB の出力を store に置くことと、
  Slurm でジョブを流すことに合わない
- 案 F: Snakemake。Slurm の executor は良いが、判定は mtime が主体で、DSL（Python の上の独自の文法）が大きい
- 案 G: 自前の Rust の駆動。デファクトの道具を使う方針に反し、判定と記録の保守を自分で持つことになる

段階の deps に、その段階のバイナリを表すものとして何を書くか:

- 案 H（採用）: コードのソース。段階のクレートと、パスで依存するクレートの `src` と `Cargo.toml`、ワークスペースの
  `Cargo.toml` と `Cargo.lock`、ツールチェーン（`rust-toolchain.toml`）、コンパイルの設定（`.cargo/config.toml`、
  ADR-0109）。内容の md5 で判定するので、rebase や文書の変更では走らない。`src` の中のテスト（`#[cfg(test)]` の
  モジュール）を変えても古くなるが、これは許す（テストを `src` の外に出すのは、テストの置き方を変える別の決定になる）
- 案 I: ビルドしたバイナリ（`target/release/xs` など）。結果を変える変更だけで走りそうに見えるが、cargo のビルドは
  ビット一致しない（パスや時刻が埋まる）ので、ビルドし直すたびに古くなる。`dvc.lock` からどのコードかも読めない

キャッシュの置き方:

- 案 J（採用）: reflink。`~/github` は btrfs（`df -T`、2026-10-11）で、reflink の写しは書き込むまでブロックを共有するので、
  ワークスペースとキャッシュで容量が倍にならない。写しは互いに独立なので、道具が出力を同じパスに上書きしてもキャッシュは
  壊れない。reflink のできないファイルシステムでは、写し（copy）に黙って戻らず、DVC がエラーで止まるよう `copy` を
  並べない
- 案 K: hardlink。容量は倍にならないが、ワークスペースのファイルがキャッシュそのものなので、DVC はそれを読み出し専用に
  する。道具がその場で書き換えるとキャッシュが壊れ、書き換えるには `dvc unprotect` が要る
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
  抽出を変えるには OpenMC を走らせ直すか、DVC の外で statepoint を読むことになり、何から作ったかが `dvc.lock` に残らない
- 案 O: statepoint を `cache: false` の outs にする。md5 は `dvc.lock` に残るが、キャッシュにないので `dvc checkout` で
  戻せず、消えたら OpenMC を走らせ直すことになる

## 決定と理由

案 A、H、J、M にする。

- 段階は `examples/<問題>/dvc.yaml` に書く（問題ごとの置き場所は ADR-0116）。段階のコマンドは作業ディレクトリをその
  ディレクトリにして、素の `cargo run --release -p <クレート> -- <入力>` にする。Slurm や手元の道具（~/.claude の
  スクリプト）に依存させない。走らせるのは `dvc repro examples/<問題>/dvc.yaml`（すべては `dvc repro -P`）で、
  1 分を超えるときは `dvc repro` 全体を Slurm で流す。段階は逐次に走らせる（Slurm で 1 本ずつ流す今の運用と同じ）
- deps には、段階の入力（入力の TOML と、そこから読むファイル）と、案 H のコードを並べる。入力の TOML は、`params` で
  キーごとに選ばず、ファイル全体を deps にする。入力のどのキーも結果を変えうるので、キーを並べ直す手間と漏れを増やす
  だけである。コメントだけの変更でも走るが、戻せば run-cache から出力を戻す
- `dvc.yaml` と `dvc.lock` をコミットする。記録の hash を checkout すれば、`dvc.lock` が、どの出力がどの入力とコードから
  作られたかを示す。`dvc checkout` でキャッシュから出力を戻し、キャッシュになければ `dvc repro` で作り直す。これが
  ADR-0122 の「hash から出し直せる」の手段になる
- DVC の利用状況の送信（`core.analytics`）は、`.dvc/config` で切る

内容の md5 で判定するので、ADR-0122 の「結果を変えない変更で流し直さない」が、道具の判定として成り立つ。`dvc.lock` は
テキストなので、何が変わって走ったかを git の差分で読める。

## 確認

M0（examples/m0/dvc.yaml の xs → lattice → core）で、次を確かめる（実装の PR で書く）。

- 最初の `dvc repro` が 3 段階を順に走らせ、`dvc.lock` を書く
- 2 回目の `dvc repro` が何も走らせない
- 入力を変えると、その段階と、出力が変わった下流の段階だけが走る

## 結果

- 段階を DVC で走らせるのは、この ADR のブランチでは M0 だけである。ほかの例（M1、VERA、MOX など）と、参照の結果
  （NJOY → OpenMC のライブラリ → OpenMC、比較）を DVC の段階にするのは #0018 の次のブランチで行う。道具（NJOY、
  OpenMC）の作り方は ADR-0031 のままで、変えるならそのブランチの ADR が置き換える
- キャッシュは、同じ出力の古い版も残すので増えていく。`dvc gc` で消す（どの版を残すかは `--workspace`、
  `--all-commits` などで選ぶ）
- `src` の中のテストを変えた段階は、結果が変わらなくても走る（案 H）
