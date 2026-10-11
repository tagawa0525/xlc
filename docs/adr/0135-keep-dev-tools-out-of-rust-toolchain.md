---
status: accepted
date: 2026-10-12
supersedes: []
superseded-by:
issues: []
requires: [ADR-0123]
---

# ADR-0135: 開発のツールは rust-toolchain.toml でなく devShell に足す

## 背景

Claude Code は Rust のファイルの言語サーバーとして、lsp-det 経由で rust-analyzer を PATH から起動する（nixfiles の
modules/home/parts/claude-code.nix が ~/.claude/skills/lsp-det-dogfood に置く lsp-det の dogfood/claude-plugin/.lsp.json）。xlc の devShell には rust-analyzer がなく、LSP の要求は
`lsp-det: the upstream language server exited before answering initialize` で失敗していた（2026-10-12）。
devShell のツールチェーンは `rust-toolchain.toml`（components は rustfmt、clippy、rust-src）から rust-overlay で作る。

`rust-toolchain.toml` は、Rust を使う DVC のステージすべての deps に入っている（[ADR-0123](0123-drive-stages-with-dvc.md)。
ツールチェーンが変われば出力が変わりうるため）。`flake.nix` は deps に入っておらず、`flake.lock` だけが入っている。
DVC のステージは default の devShell で実行する。

## 検討した案

- 案 A: `rust-toolchain.toml` の components に `rust-analyzer` を足す。1 行で済み、rustup でも同じものが入る。
  ただしファイルが変わるので、計算に関係しない変更で Rust の段階がすべて changed になり、再実行するか、実行せずに
  `dvc.lock` を書き換えるかを迫られる。xlc の利用者は Nix の devShell だけで、rustup の利点は使わない
- 案 B: flake でツールチェーンを override し、components に `rust-analyzer` を足す。ブランチで 2 つの形を試して捨てた。
  `rust-toolchain.toml` を `lib.importTOML` で読み直して components に足す形は、rust-overlay が解いた一覧を二重に読む。
  `override (old: …)` で `old.extensions` に足す形はそれを避けるが、どちらも段階をビルドするツールチェーンの derivation を
  変える（版は同じでも store パスが変わる）。その変更は `flake.nix` にあって deps に表れないので、同じ override に後から
  targets などを足すとコンパイラが黙って変わる。防ぐのは文書の規則だけになる
- 案 C: ビルド用と開発用の devShell を分け、DVC のステージと Slurm のジョブはビルド用で実行する。段階のビルドが開発のツールから
  切り離されるが、すべての段階の cmd を `nix develop ..#build -c` で包み直し、direnv の shell と段階の shell を
  使い分ける手間が増える。開発のツールがコンパイラを変えないなら、ビルドする rustc は案 E と同じ。shell の設定の変更
（帰結）から段階を切り離せる利点はあるが、それはこの決定の前からある別の問題で、言語サーバーを足すことには要らない
- 案 D: nixpkgs の rust-analyzer を入れる。手続き型マクロの展開はツールチェーンの sysroot の proc-macro-srv を使うので
  動くが、rust-analyzer の版が nixpkgs の更新で動き、ツールチェーンの版と別に変わる
- 案 E（採用）: ツールチェーンはそのままにし、rust-overlay の同じ版の `rust-analyzer` の component
  （ツールチェーンの `availableComponents`）を default の devShell に横に置く

## 決定と理由

`rust-toolchain.toml` には、ビルドとチェック（run-checks.sh の rustfmt、clippy）に要る components と、sysroot に
置く必要のあるもの（rust-src。rust-analyzer が std のソースを sysroot の lib/rustlib/src から読む）だけを書く。
言語サーバーのような開発のツールは、ツールチェーンを変えずに、同じ版の component として default の devShell に足す。
案 C より案 E を採るのは、開発のツールがコンパイラを変えないなら、shell を分けても段階のビルドは同じで、言語サーバーを
足すためには cmd の包み直しと shell の使い分けが増えるだけのため。

DVC の deps はビルドの入力を表すので、計算に関係しないツールの変更で段階を changed にしない。段階をビルドする
コンパイラの derivation は `rust-toolchain.toml` と `flake.lock` だけで決まる形を保ち、`flake.nix` の文書の規則に頼らない。
rust-analyzer は `rustc --print sysroot` でツールチェーンの sysroot と proc-macro-srv を見つけるので、component を
分けても版はそろう。

## 確認

6b75c94b（分岐点は main の 5f984f72。2026-10-12）の worktree で確かめた。

- devShell の rustc は main の devShell と同じ store パス（`24xl8c07…-rust-minimal-1.99.0`、1.99.0、b940084d7）
- `rust-analyzer --version` は `rust-analyzer 1.99.0 (b940084 2026-09-28)`。sysroot の libexec に
  `rust-analyzer-proc-macro-srv` がある
- `grep flake.nix examples/*/dvc.yaml` は空で、`flake.nix` は DVC の deps にない
- worktree で `nix develop -c claude -p` として起動した Claude Code から、LSP の findReferences を
  crates/numeric/src/lu.rs の `Lu` に投げ、7 ファイルの 23 件の参照が返った（lsp-det 経由。crates/lattice/src/b1.rs、
  cmfd.rs、crates/numeric/src/pins/node.rs など）
- 同じ方法で、crates/numeric/src/lu.rs の `Vec` の goToDefinition が sysroot の
  `lib/rustlib/src/rust/library/alloc/src/vec/mod.rs` を返した（rust-src が sysroot から読まれる）
- 手続き型マクロの展開（proc-macro-srv）は確かめていない。derive した構造体への goToImplementation は構造体の位置を
  返すだけで、展開したかを見分けられなかった

## 帰結

- 開発のツールを足しても `rust-toolchain.toml`、`flake.lock`、コンパイラの store パスは変わらず、DVC のステージは
  changed にならず、cargo のビルドもやり直さない
- 段階は default の devShell で走るので、`flake.nix` の shell の設定（shellHook の環境変数、別の C コンパイラなど）を
  変えると、deps に表れずに段階のビルドが変わりうる。これはこの決定の前からある限界で、この決定はツールチェーンの
  derivation だけを `rust-toolchain.toml` と `flake.lock` で決まる形に保つ
- rust-analyzer は開いたワークスペースで `cargo check` を実行し、既定では同じ `target/` を使う。run-checks.sh の
  clippy やテストとビルドのロックを取り合いうる（確かめていない）
- 走っている Claude Code の LSP のサーバーは起動時の PATH を引き継ぐので、direnv で devShell を読み直してから
  Claude Code を起動し直すと Rust の LSP が使える
