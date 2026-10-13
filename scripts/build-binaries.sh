#!/usr/bin/env bash
# build-binaries.sh — DVC の build の段階（リポジトリの直下の dvc.yaml）。xs、lattice、core の release の実行ファイルを
# 記号を落としてビルドし、bin/ に置く。examples/*/dvc.yaml の Rust の段階は、この実行ファイルを deps にして呼ぶ（ADR-0159）
#
# Usage: scripts/build-binaries.sh（リポジトリの直下で、default の devShell の中から。dvc repro が呼ぶ）
#
# cargo はリンクの環境（NIX_LDFLAGS、ホームの設定のリンカ）を fingerprint に含めず、環境だけが変わってもリンクし直さない。
# そのため、手で使う target/release と分けた target/dvc で、実行ファイルの 3 つのパッケージを毎回消してからビルドする
# （コンパイルとリンクを今の環境でやり直す。依存のクレートは再利用する）。

set -euo pipefail

top=$(git rev-parse --show-toplevel)
cd "$top"

# stdenv は devShell に入ったディレクトリの outputs/out/lib を RUNPATH に足す。足すと同じソースでも worktree ごとに
# 実行ファイルのバイトが変わるので、flake.nix の default の devShell は NIX_NO_SELF_RPATH = 1 で止めている。下の grep も
# 同じ漏れを見つけるが、ビルドの前に、どの shell で実行するべきかを示して止める
if [[ "${NIX_NO_SELF_RPATH:-}" != 1 ]]; then
  echo "ERROR: NIX_NO_SELF_RPATH が 1 でない。flake.nix の default の devShell の中で実行する（ADR-0159）" >&2
  exit 1
fi

# シェルで設定したコンパイルとリンクの設定（手で試すときの RUSTFLAGS など）は、DVC の deps に表れないまま実行ファイルを
# 変えるので、設定されていれば止める
flags='^(RUSTFLAGS|CARGO_ENCODED_RUSTFLAGS|CARGO_BUILD_RUSTFLAGS|CARGO_PROFILE_[A-Z0-9_]+|CARGO_TARGET_[A-Z0-9_]+_(RUSTFLAGS|LINKER))='
if env | grep -qE "$flags"; then
  echo "ERROR: コンパイルかリンクの設定が環境変数にある（$(env | grep -oE "$flags" | tr '\n' ' ')）。外してから実行する（ADR-0159）" >&2
  exit 1
fi

# 実行ファイルの名前（パッケージは xlc-<名前>）。clean、build、コピーする先で同じ並びを使う
names=(xs lattice core)
packages=()
for name in "${names[@]}"; do
  packages+=(-p "xlc-$name")
done

cargo clean --release --target-dir target/dvc "${packages[@]}"
# 記号表（.symtab、.strtab）は、計算に使わない内部の記号の名前（`.llvm.<数>` の添字）が、テストを tests.rs に移すと変わった
# （モジュールの構成で決まると見ているが、確かめていない）。計算する節は同じでも md5 が変わるので、リンクのときに記号を
# 落とす（ADR-0160）。手で使う target/release は記号を残す
# ワークスペースの全体でなく 3 つのパッケージだけを 1 回の cargo でビルドする。全体をビルドすると、実行ファイルが依存しない
# クレートの依存が共有の依存の feature を足し、実行ファイルのバイトが変わる。3 つを別々の cargo にしないのは、3 つの間の
# feature の統合で今の実行ファイルが決まっているため（ADR-0163）
cargo build --release --locked "${packages[@]}" --target-dir target/dvc --config 'profile.release.strip="symbols"'

mkdir -p bin
for name in "${names[@]}"; do
  built="target/dvc/release/$name"
  # worktree のパスや、devShell に入ったディレクトリの outputs/out/lib（入ったのが別の場所なら $top を含まない）が埋まった
  # 実行ファイルは、ほかの worktree と md5 が合わず、段階を再実行することになる
  if grep -qaF -e "$top" -e /outputs/out/lib -- "$built"; then
    echo "ERROR: $built に worktree のパスか outputs/out/lib が埋まっている（ADR-0159）" >&2
    exit 1
  fi
  cp -- "$built" "bin/$name"
done
