#!/usr/bin/env bash
# ADR-0157 の測定の実行ファイルを作る: コミットごとのソースを、整列の設定ごとに release でビルドし、bin/ に置く。
# 使い方: build.sh <作業ディレクトリ> <コミット>... （例: build.sh ~/.cache/xlc-0157 6b39ad13 ab7c70b3 b501c22e）
# 作業ディレクトリはリポジトリの外に置く。cargo は上のディレクトリの .cargo/config.toml も読むので、リポジトリの中に
# ソースを広げると、このリポジトリの rustflags（ADR-0157 の設定を含む）が足される。
# 各ソースの .cargo/config.toml の rustflags はそのコミットのもの（ADR-0157 より前のコミットは target-cpu=znver5 だけ）で、
# 下の設定はそれに `--config` で足す。
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(git -C "$here" rev-parse --show-toplevel)
work=$1; shift
# 作業ディレクトリの上に .cargo/config.toml か .cargo/config があれば止める（ホームの ~/.cargo/config.toml はどのビルドにも効くので除く）
d=$(realpath -m "$work")
while [ "$d" != / ]; do
  if [ "$d" != "$HOME" ] && { [ -e "$d/.cargo/config.toml" ] || [ -e "$d/.cargo/config" ]; }; then
    echo "$d/.cargo/config(.toml) would add its rustflags; use a work directory outside it" >&2
    exit 1
  fi
  d=$(dirname "$d")
done
declare -A flags=(
  [v0]=''
  [v1]='"-C","llvm-args=-align-all-functions=6"'
  [v2]='"-C","llvm-args=-align-all-functions=6","-C","llvm-args=-align-all-nofallthru-blocks=6"'
  [v3]='"-C","llvm-args=-align-loops=64"'
)
# 環境変数の rustflags は設定ファイルと --config の rustflags を上書きし、設定ごとの差が消えるので止める（scripts/build-binaries.sh と同じ確認）
env_flags='^(RUSTFLAGS|CARGO_ENCODED_RUSTFLAGS|CARGO_BUILD_RUSTFLAGS|CARGO_PROFILE_[A-Z0-9_]+|CARGO_TARGET_[A-Z0-9_]+_(RUSTFLAGS|LINKER))='
# grep -q は最初の一致で抜け、pipefail の下で env の SIGPIPE がパイプラインを失敗にするので、-q を使わない
if env | grep -E "$env_flags" >/dev/null; then
  echo "compile or link settings are set in the environment: $(env | grep -oE "$env_flags" | tr '\n' ' ')" >&2
  exit 1
fi
# devShell に入ったディレクトリが RUNPATH に入らないよう、scripts/build-binaries.sh と同じく default の devShell の中で走らせる
if [[ "${NIX_NO_SELF_RPATH:-}" != 1 ]]; then
  echo "NIX_NO_SELF_RPATH is not 1; run this in the default devShell of flake.nix (ADR-0159)" >&2
  exit 1
fi
mkdir -p "$work" "$here/bin"
for c in "$@"; do
  rm -rf "$work/src-$c"; mkdir -p "$work/src-$c"
  git -C "$repo" archive "$c" | tar -x -C "$work/src-$c"
  # ADR-0157 のあとのコミットは、そのコミットの設定に整列が入っており、v0 が「なし」にならないので止める
  if grep -q 'align' "$work/src-$c/.cargo/config.toml"; then
    echo "$c already sets an alignment in .cargo/config.toml; v0 would not be the build without it" >&2
    exit 1
  fi
  for v in v0 v1 v2 v3; do
    f=${flags[$v]}
    cfg=()
    [ -n "$f" ] && cfg=(--config "target.x86_64-unknown-linux-gnu.rustflags=[$f]")
    # コミットごとに target を分ける（同じ target を使い回すと、パッケージの名前が同じなので作り直さないことがある）
    # cargo はリンクの環境を fingerprint に含めないので、使い回す target でも 3 つのパッケージを消してから作る
    (cd "$work/src-$c" && cargo clean -q --release --target-dir "$work/target-$c-$v" -p xlc-xs -p xlc-lattice -p xlc-core)
    (cd "$work/src-$c" && cargo build -q --release --locked -p xlc-xs -p xlc-lattice -p xlc-core \
      --target-dir "$work/target-$c-$v" "${cfg[@]}")
    for b in xs lattice core; do cp "$work/target-$c-$v/release/$b" "$here/bin/$b-$c-$v"; done
  done
done
