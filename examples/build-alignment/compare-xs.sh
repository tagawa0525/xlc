#!/usr/bin/env bash
# ADR-0144 の判定（xs）: 2 つのコミットの xs を release でビルドし、examples/m0/xs.toml で 2 回ずつ交互に走らせて、出力の
# msgpack のバイト一致と時間を比べる。変更後は 1 CPU に縛っても 1 回走らせる。J, K の表の構築を含むテスト
# （table_matches_direct_integration）の時間も 3 回ずつ測る。
# 使い方: compare-xs.sh build <作業ディレクトリ> <分岐点> <変更後>（default の devShell の中で）
#         compare-xs.sh run（Slurm の排他のジョブで）
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
case "$1" in
build)
  work=$2 base=$3 new=$4
  repo=$(git rev-parse --show-toplevel)
  # 作業ディレクトリの上の .cargo/config.toml は rustflags を足すので、リポジトリの外に置く（ホームの設定は除く）
  d=$(realpath -m "$work")
  while [ "$d" != / ]; do
    if [ "$d" != "$HOME" ] && { [ -e "$d/.cargo/config.toml" ] || [ -e "$d/.cargo/config" ]; }; then
      echo "$d/.cargo/config(.toml) would add its rustflags; use a work directory outside it" >&2
      exit 1
    fi
    d=$(dirname "$d")
  done
  if [[ "${NIX_NO_SELF_RPATH:-}" != 1 ]]; then
    echo "NIX_NO_SELF_RPATH is not 1; run this in the default devShell of flake.nix (ADR-0159)" >&2
    exit 1
  fi
  mkdir -p "$work" bin-compare
  for tag in base new; do
    if [ $tag = base ]; then c=$base; else c=$new; fi
    rm -rf "$work/src-$tag"; mkdir -p "$work/src-$tag"
    git -C "$repo" archive "$c" | tar -x -C "$work/src-$tag"
    (cd "$work/src-$tag" && cargo clean -q --release --locked --target-dir "$work/target-$tag" -p xlc-xs &&
      cargo build -q --release --locked -p xlc-xs --target-dir "$work/target-$tag")
    cp "$work/target-$tag/release/xs" "bin-compare/xs-$tag"
    test=$(cd "$work/src-$tag" && cargo test --release --locked -p xlc-xs --lib --no-run --target-dir "$work/target-$tag" 2>&1 |
      grep -oE "$work/target-$tag/release/deps/xlc_xs-[0-9a-f]+" | head -1)
    cp "$test" "bin-compare/xs-test-$tag"
    echo "$tag $c $(md5sum < "bin-compare/xs-$tag")"
  done
  ;;
run)
  if [ -e ../m0/wims69.msgpack ]; then
    echo "../m0/wims69.msgpack exists; run this in a worktree without the M0 outputs" >&2
    exit 1
  fi
  rm -rf out/compare; mkdir -p out/compare
  trap 'rm -f ../m0/wims69.msgpack' EXIT
  run() { # <tag> <command...>
    local tag=$1; shift
    local s e dt
    s=$(date +%s.%N)
    (cd ../m0 && "$@" xs.toml 2>"$here/out/compare/$tag.log")
    e=$(date +%s.%N)
    dt=$(echo "$e - $s" | bc)
    mv ../m0/wims69.msgpack "out/compare/$tag.msgpack"
    echo "TIME xs $tag $dt"
  }
  run base1 "$here/bin-compare/xs-base"
  run new1 "$here/bin-compare/xs-new"
  run new2 "$here/bin-compare/xs-new"
  run base2 "$here/bin-compare/xs-base"
  run new-1cpu taskset -c 0 "$here/bin-compare/xs-new"
  for t in base2 new1 new2 new-1cpu; do
    if cmp out/compare/base1.msgpack "out/compare/$t.msgpack"; then echo "CMP base1 $t identical"; else echo "CMP base1 $t DIFFER"; fi
  done
  md5sum out/compare/*.msgpack
  for i in 1 2 3; do
    for v in base new; do
      "bin-compare/xs-test-$v" --exact resonance::line_shape::tests::table_matches_direct_integration 2>&1 |
        grep -oE "finished in [0-9.]+s" | sed "s/^/TIME table $v $i /"
    done
  done
  ;;
*) echo "usage: compare-xs.sh build <work> <base> <new> | run" >&2; exit 1 ;;
esac
