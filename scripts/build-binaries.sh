#!/usr/bin/env bash
# build-binaries.sh — DVC の build の段階（リポジトリの直下の dvc.yaml）。xs、lattice、core の release の実行ファイルを
# ビルドして bin/ に置く。examples/*/dvc.yaml の Rust の段階は、この実行ファイルを deps にして呼ぶ（ADR-0159）
#
# Usage: scripts/build-binaries.sh（リポジトリの直下で、default の devShell の中から。dvc repro が呼ぶ）
#
# target/dvc に分けてビルドする。cargo はリンクの環境（NIX_LDFLAGS、ホームの設定のリンカ）を fingerprint に含めないので、
# 手で走らせる target/release を使い回すと、前の環境でリンクした実行ファイルが残りうる。

set -euo pipefail

top=$(git rev-parse --show-toplevel)
cd "$top"

# stdenv は devShell に入ったディレクトリの outputs/out/lib を RUNPATH に足す。足すと同じソースでも worktree ごとに
# 実行ファイルのバイトが変わるので、flake.nix の default の devShell は NIX_NO_SELF_RPATH = 1 で止めている
if [[ "${NIX_NO_SELF_RPATH:-}" != 1 ]]; then
  echo "ERROR: NIX_NO_SELF_RPATH が 1 でない。flake.nix の default の devShell の中で走らせる（ADR-0159）" >&2
  exit 1
fi

cargo build --release --locked --target-dir target/dvc

mkdir -p bin
for name in xs lattice core; do
  built="target/dvc/release/$name"
  # worktree のパスが埋まった実行ファイルは、ほかの worktree と md5 が合わず、段階を走らせ直すことになる
  if grep -qaF -- "$top" "$built"; then
    echo "ERROR: $built に worktree のパス $top が埋まっている（ADR-0159）" >&2
    exit 1
  fi
  cp -- "$built" "bin/$name"
done
