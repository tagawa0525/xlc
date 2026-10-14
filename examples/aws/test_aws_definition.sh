#!/usr/bin/env bash
# 実行の定義（examples/<問題>/aws.toml。ADR-0174）を読む部品 def_parse と def_find の検査。flake の check（aws-definition）が走らせる。
#
# 使い方: examples/aws/test_aws_definition.sh [<aws_common.sh のパス> [<aws.toml>...]]
#   aws.toml を渡すと、その定義が読めることを確かめる。リポジトリの直下から走らせれば、定義が挙げる道具のファイルの存在も確かめる
set -euo pipefail
common=${1:-$(dirname "$0")/aws_common.sh}
shift || true
# shellcheck source=aws_common.sh
source "$common"

failures=0
fail() { echo "FAIL: $*" >&2; failures=$((failures + 1)); }

valid='# 実行の定義
tools = [
  "examples/a.py",   # 道具
  "examples/b/c.py",
]
entry = "examples/b/run_remote.sh"
progress = "steps"
library = "examples/b/openmc-library"
data = ["examples/b/chain.xml"]
jobs = '"'"'^openmc-1C-[a-z0-9]+$'"'"'
'

# def_parse は変数に値を入れるので、戻りの値を見るものは現在のシェルで呼ぶ（失敗は exit なので subshell）
(
  def_parse "$valid"
  [[ ${#DEF_TOOLS[@]} == 2 && ${DEF_TOOLS[0]} == examples/a.py && ${DEF_TOOLS[1]} == examples/b/c.py ]] || fail "tools: ${DEF_TOOLS[*]}"
  [[ $DEF_ENTRY == examples/b/run_remote.sh ]] || fail "entry: $DEF_ENTRY"
  [[ $DEF_PROGRESS == steps ]] || fail "progress: $DEF_PROGRESS"
  [[ $DEF_LIBRARY == examples/b/openmc-library ]] || fail "library: $DEF_LIBRARY"
  [[ ${#DEF_DATA[@]} == 1 && ${DEF_DATA[0]} == examples/b/chain.xml ]] || fail "data: ${DEF_DATA[*]}"
  [[ $DEF_JOBS == '^openmc-1C-[a-z0-9]+$' ]] || fail "jobs: $DEF_JOBS"
  [[ ! openmc-1C-rings1 =~ $DEF_JOBS ]] || fail "jobs が 1C-rings1 を受け付けた"
  [[ openmc-1C-rings1x =~ $DEF_JOBS ]] || fail "jobs が openmc-1C-rings1x を受け付けない"
) || fail "有効な定義を読めない"

# data は省略できる（空の配列）
(
  def_parse "$(grep -v '^data' <<< "$valid")"
  [[ ${#DEF_DATA[@]} == 0 ]] || fail "data を省いても空でない"
) || fail "data のない定義を読めない"

# 前の定義の値を持ち越さない
(
  def_parse "$valid"
  def_parse "$(grep -v '^data' <<< "$valid")"
  [[ ${#DEF_DATA[@]} == 0 ]] || fail "前の定義の data を持ち越した"
) || fail "続けて読めない"

rejects() { # <名前> <定義の本文>
  if (def_parse "$2") > /dev/null 2>&1; then fail "$1: 受け付けた"; fi
}
rejects "知らないキー" "${valid}extract = \"examples/x.py\""
rejects "キーの重複" "${valid}entry = \"examples/b/other.sh\""
for key in tools entry progress library jobs; do
  rejects "$key がない" "$(grep -v "^$key" <<< "$valid")"
done
rejects "tools がない（複数行の中だけ消す）" "$(sed '/^tools/,/^]/d' <<< "$valid")"
rejects "progress が不明" "${valid//progress = \"steps\"/progress = \"batches\"}"
rejects "絶対パス" "${valid//examples\/b\/run_remote.sh//tmp/run_remote.sh}"
rejects ".. を含む" "${valid//examples\/b\/run_remote.sh/examples\/..\/run_remote.sh}"
rejects "使えない文字" "${valid//examples\/b\/chain.xml/examples\/b\/chain xml}"
rejects "配列が閉じない" "${valid//\"examples\/b\/c.py\",/\"examples\/b\/c.py\"}"
rejects "tools が空" "$(sed '/^tools/,/^]/c\tools = []' <<< "$valid")"
rejects "文字列でない" "${valid//entry = \"examples\/b\/run_remote.sh\"/entry = 3}"
rejects "jobs が正規表現として読めない" "${valid//\'^openmc-1C-[a-z0-9]+\$\'/\'^openmc-[\'}"
rejects "行を読めない" "${valid}これは定義ではない"

# def_find: 入力から上へ、aws.toml のある最初のディレクトリ（リポジトリの直下から相対）
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
(
  cd "$tmp"
  mkdir -p examples/p/q examples/other
  : > examples/p/aws.toml
  : > examples/p/q/in.toml
  : > examples/other/in.toml
  [[ $(def_find examples/p/q/in.toml) == examples/p/aws.toml ]] || fail "def_find: 上のディレクトリの定義を使わない"
  # 近い定義が優先
  : > examples/p/q/aws.toml
  [[ $(def_find examples/p/q/in.toml) == examples/p/q/aws.toml ]] || fail "def_find: 近い定義を使わない"
  if (def_find examples/other/in.toml) > /dev/null 2>&1; then fail "def_find: 定義がなくても止まらない"; fi
) || fail "def_find"

# 渡された定義を読み、道具のファイルがあるか（リポジトリの直下から走らせたときだけ）確かめる
for file in "$@"; do
  (
    def_parse "$(cat "$file")"
    if [[ -f flake.nix ]]; then
      for f in "${DEF_TOOLS[@]}" "$DEF_ENTRY"; do [[ -f $f ]] || fail "$file: $f がありません"; done
    fi
  ) || fail "$file を読めない"
done

if ((failures)); then
  echo "$failures 件の失敗" >&2
  exit 1
fi
echo "aws_definition: ok"
